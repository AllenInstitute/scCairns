import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")

from sccairns import summarize as sr
from sccairns.cli import COMMANDS
from sccairns.record_filter import infer_cluster_key, record_filter


def _with_embeddings(adata):
    adata.obsm["X_scVI"] = np.ones((adata.n_obs, 2))
    adata.obsm["X_umap"] = np.ones((adata.n_obs, 2))
    return adata


def _write_input_output(tmp_path, adata, kept_ids):
    input_path = tmp_path / "input.h5ad"
    output_path = tmp_path / "filtered.h5ad"
    adata.write_h5ad(input_path)
    adata[kept_ids, :].copy().write_h5ad(output_path)
    return input_path, output_path


def test_record_filter_registered_in_cli():
    assert COMMANDS["record-filter"] == "sccairns.record_filter"


def test_infer_cluster_key_from_embedding():
    assert infer_cluster_key("X_scVI") == "leiden"
    assert infer_cluster_key("X_scVI_small_gene_nb") == "leiden_small_gene_nb"
    assert infer_cluster_key("X_pca_harmony") == "leiden_harmony"
    assert infer_cluster_key("X_custom") is None


def test_record_filter_writes_manifest_and_sidecars(tmp_path, tiny_adata):
    adata = _with_embeddings(tiny_adata)
    input_path, output_path = _write_input_output(
        tmp_path, adata, ["cell0", "cell1"]
    )
    output_dir = tmp_path / "round_02"
    output_dir.mkdir()
    (output_dir / "round_manifest.json").write_text(json.dumps({
        "integration": {"n_cells": 2},
        "inspection": {"n_clusters": 1},
    }))

    manifest_path = record_filter(
        input_h5ad=str(input_path),
        output_h5ad=str(output_path),
        output_dir=str(output_dir),
        embedding_key="X_scVI",
        filter_used="manual neuronal subset",
        notes="interactive review",
        seed=0,
        batch_key="data_origin",
    )

    assert Path(manifest_path) == output_dir / "round_manifest.json"
    manifest = json.loads(Path(manifest_path).read_text())
    dec = manifest["decisions"]
    assert manifest["integration"] == {"n_cells": 2}
    assert manifest["inspection"] == {"n_clusters": 1}
    assert dec["input_cells"] == 4
    assert dec["output_cells"] == 2
    assert dec["removed"] == 2
    assert dec["removal_pct"] == 50.0
    assert dec["input_fingerprint"]["sha256"]
    assert dec["output_fingerprint"]["sha256"]
    assert dec["filtered_on"] == {
        "variant": "default",
        "cluster_key": "leiden",
        "latent_key": "X_scVI",
        "umap_key": "X_umap",
        "source": "interactive",
    }
    assert dec["actions_summary"] == ["interactive_filter: cluster 1 (-2)"]

    kept = pd.read_csv(output_dir / "cells_to_keep.csv")
    removed = pd.read_csv(output_dir / "cells_removed.csv")
    assert kept["obs_name"].tolist() == ["cell0", "cell1"]
    assert removed["obs_name"].tolist() == ["cell2", "cell3"]
    assert (output_dir / "decisions_applied.yaml").exists()


def test_record_filter_generic_action_when_cluster_key_missing(tmp_path, tiny_adata):
    adata = _with_embeddings(tiny_adata)
    del adata.obs["leiden"]
    input_path, output_path = _write_input_output(
        tmp_path, adata, ["cell0", "cell1"]
    )

    manifest_path = record_filter(
        input_h5ad=str(input_path),
        output_h5ad=str(output_path),
        output_dir=str(tmp_path / "round_02"),
        embedding_key="X_scVI",
        filter_used="manual lasso",
        notes="",
    )

    manifest = json.loads(Path(manifest_path).read_text())
    assert manifest["decisions"]["actions_summary"] == [
        "interactive_filter: manual lasso (-2)"
    ]


def test_record_filter_rejects_output_cells_not_in_input(tmp_path, tiny_adata):
    adata = _with_embeddings(tiny_adata)
    input_path = tmp_path / "input.h5ad"
    output_path = tmp_path / "filtered.h5ad"
    adata.write_h5ad(input_path)
    output = adata.copy()
    output.obs_names = ["cell0", "cell1", "cell2", "new_cell"]
    output.write_h5ad(output_path)

    with pytest.raises(ValueError, match="not present in input"):
        record_filter(
            input_h5ad=str(input_path),
            output_h5ad=str(output_path),
            output_dir=str(tmp_path / "round_02"),
            embedding_key="X_scVI",
            filter_used="manual lasso",
        )


def test_record_filter_manifest_is_visible_to_summarize(tmp_path, tiny_adata):
    adata = _with_embeddings(tiny_adata)
    round_dir = tmp_path / "round_02"
    input_path, output_path = _write_input_output(
        tmp_path, adata, ["cell0", "cell1"]
    )
    record_filter(
        input_h5ad=str(input_path),
        output_h5ad=str(output_path),
        output_dir=str(round_dir),
        embedding_key="X_scVI",
        filter_used="manual neuronal subset",
    )

    manifests = [sr.load_manifest(p) for p in sr.find_manifests(str(tmp_path))]
    ledger = sr.extract_decisions_ledger(manifests)
    assert ledger == [{
        "round": "round_02",
        "filtered_on": "default",
        "cluster_key": "leiden",
        "action": "interactive_filter: cluster 1 (-2)",
    }]

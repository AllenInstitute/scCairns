import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))


anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from inspect_integration import (  # noqa: E402
    apply_decisions,
    cluster_qc_summary,
    find_sweep_architecture_keys,
)


def _adata():
    adata = anndata.AnnData(np.ones((4, 3)))
    adata.obs_names = [f"cell{i}" for i in range(4)]
    adata.var_names = [f"gene{i}" for i in range(3)]
    adata.obs["leiden"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obs["data_origin"] = pd.Categorical(["a", "a", "b", "b"])
    adata.obs["n_genes_by_counts"] = [600, 700, 300, 800]
    return adata


def test_apply_decisions_filters_valid_query(tmp_path):
    filtered = apply_decisions(
        _adata(),
        {"keep_cells": [{"query": "n_genes_by_counts > 500"}]},
        batch_key="data_origin",
        output_dir=str(tmp_path),
    )

    assert filtered.n_obs == 3
    assert (tmp_path / "filtered.h5ad").exists()
    assert (tmp_path / "filtering_retention_summary.csv").exists()
    assert (tmp_path / "round_manifest.json").exists()

    summary = pd.read_csv(tmp_path / "filtering_retention_summary.csv")
    overall = summary.loc[summary["grouping"] == "overall"].iloc[0]
    assert overall["pre_filter_cells"] == 4
    assert overall["post_filter_cells"] == 3
    assert overall["removed_cells"] == 1

    batch_b = summary.loc[
        (summary["grouping"] == "batch") & (summary["batch"] == "b")
    ].iloc[0]
    assert batch_b["pre_filter_cells"] == 2
    assert batch_b["post_filter_cells"] == 1


def test_apply_decisions_raises_on_invalid_query(tmp_path):
    with pytest.raises(ValueError, match="Keep query failed"):
        apply_decisions(
            _adata(),
            {"keep_cells": [{"query": "missing_column > 0"}]},
            output_dir=str(tmp_path),
        )


def test_cluster_qc_summary_reports_latent_and_umap_silhouette():
    adata = _adata()
    adata.obsm["X_scVI"] = np.array(
        [[0.0, 0.0], [0.0, 0.2], [3.0, 3.0], [3.0, 3.2]]
    )
    adata.obsm["X_umap"] = np.array(
        [[0.0, 0.0], [0.2, 0.0], [8.0, 0.0], [8.2, 0.0]]
    )

    summary = cluster_qc_summary(
        adata,
        cluster_key="leiden",
        batch_key="data_origin",
        latent_key="X_scVI",
        umap_key="X_umap",
    )

    assert "silhouette_latent" in summary.columns
    assert "silhouette_umap" in summary.columns
    assert "silhouette" not in summary.columns
    assert np.isfinite(summary["silhouette_latent"]).all()
    assert np.isfinite(summary["silhouette_umap"]).all()


def test_find_sweep_architecture_keys_uses_config_order():
    adata = _adata()
    for name in ["small", "medium"]:
        adata.obs[f"leiden_{name}"] = pd.Categorical(["0", "0", "1", "1"])
        adata.obsm[f"X_scVI_{name}"] = np.ones((4, 2))
        adata.obsm[f"X_umap_{name}"] = np.ones((4, 2))

    config = {
        "integration": {
            "sweep": [
                {"name": "medium"},
                {"name": "small"},
            ]
        }
    }

    runs = find_sweep_architecture_keys(adata, config)

    assert [run["name"] for run in runs] == ["medium", "small"]
    assert runs[0]["cluster_key"] == "leiden_medium"
    assert runs[1]["latent_key"] == "X_scVI_small"

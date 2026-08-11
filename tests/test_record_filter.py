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

    with pytest.raises(ValueError, match="not present in the input"):
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


# ── Cell-ID lists and ingest-time filtering ────────────────────────────────────
#
# The filtered side is only ever read for the IDs that survived, so a tool that
# cannot write AnnData (Seurat, Loupe) can hand back a CSV instead.

def test_cell_id_list_stands_in_for_the_filtered_h5ad(tmp_path, tiny_adata):
    adata = _with_embeddings(tiny_adata)
    input_path = tmp_path / "input.h5ad"
    adata.write_h5ad(input_path)
    kept = tmp_path / "kept.csv"
    kept.write_text("obs_name\ncell0\ncell1\n")

    manifest_path = record_filter(
        input_h5ad=str(input_path),
        output_cell_ids=str(kept),
        output_dir=str(tmp_path / "round_02"),
        embedding_key="X_scVI",
        filter_used="Seurat QC",
    )
    dec = json.loads(Path(manifest_path).read_text())["decisions"]

    assert dec["input_cells"] == 4 and dec["output_cells"] == 2
    assert dec["removed"] == 2
    assert dec["compared"] == {"input": "h5ad", "output": "cell_ids"}
    assert dec["output_h5ad"] is None
    assert dec["output_cell_ids"] == str(kept)
    # The input h5ad still supplies genes and the per-cluster breakdown.
    assert dec["input_genes"] == 3
    assert dec["output_genes"] is None
    assert dec["actions_summary"] == ["interactive_filter: cluster 1 (-2)"]


@pytest.mark.parametrize("body,label", [
    ("obs_name\ncell0\ncell1\n", "our own cells_to_keep.csv header"),
    ("cell0\ncell1\n", "bare one-per-line list"),
    ('"","orig.ident"\n"cell0","A"\n"cell1","A"\n', "R write.csv of meta.data"),
    ("barcode,nCount\ncell0,10\ncell1,20\n", "named first column"),
    # write.csv(data.frame(x = Cells(obj))) prepends row numbers; reading column
    # 0 blindly yields "1","2" instead of barcodes.
    ('"","x"\n"1","cell0"\n"2","cell1"\n', "R write.csv of Cells()"),
    ("1,cell0\n2,cell1\n", "headerless row numbers then IDs"),
])
def test_cell_id_list_accepts_common_export_shapes(tmp_path, tiny_adata, body, label):
    input_path = tmp_path / "input.h5ad"
    tiny_adata.write_h5ad(input_path)
    kept = tmp_path / "kept.csv"
    kept.write_text(body)

    manifest_path = record_filter(
        input_h5ad=str(input_path),
        output_cell_ids=str(kept),
        output_dir=str(tmp_path / "out"),
        filter_used="QC",
    )
    dec = json.loads(Path(manifest_path).read_text())["decisions"]
    assert dec["output_cells"] == 2, f"failed for {label}"


def test_both_sides_as_cell_id_lists(tmp_path):
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    before.write_text("\n".join(f"cell{i}" for i in range(10)) + "\n")
    after.write_text("\n".join(f"cell{i}" for i in range(6)) + "\n")

    manifest_path = record_filter(
        input_cell_ids=str(before),
        output_cell_ids=str(after),
        output_dir=str(tmp_path / "round_00"),
        filter_used="Seurat ingest QC: nFeature_RNA > 500",
        notes="no h5ad exists yet",
    )
    dec = json.loads(Path(manifest_path).read_text())["decisions"]

    assert (dec["input_cells"], dec["output_cells"], dec["removed"]) == (10, 6, 4)
    assert dec["removal_pct"] == 40.0
    assert dec["compared"] == {"input": "cell_ids", "output": "cell_ids"}
    assert dec["input_genes"] is None and dec["output_genes"] is None
    assert dec["parent_round"] is None          # an ID list has no round beside it
    assert dec["actions_summary"] == [
        "interactive_filter: Seurat ingest QC: nFeature_RNA > 500 (-4)"]
    # The sidecars are still written, so the removed cells are enumerable.
    kept = pd.read_csv(tmp_path / "round_00" / "cells_to_keep.csv")
    removed = pd.read_csv(tmp_path / "round_00" / "cells_removed.csv")
    assert len(kept) == 6 and len(removed) == 4


def test_ingest_filter_needs_no_embedding_or_clustering(tmp_path, capsys):
    """At ingest there is no embedding yet; that is normal, not a warning."""
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    before.write_text("a\nb\nc\n")
    after.write_text("a\nb\n")

    manifest_path = record_filter(
        input_cell_ids=str(before),
        output_cell_ids=str(after),
        output_dir=str(tmp_path / "round_00"),
        filter_used="percent.mt < 10",
    )
    dec = json.loads(Path(manifest_path).read_text())["decisions"]

    assert dec["filtered_on"] == {
        "variant": None, "cluster_key": None, "latent_key": None,
        "umap_key": None, "source": "interactive",
    }
    assert dec["cluster_key"] is None
    assert "[WARN]" not in capsys.readouterr().err


def test_mismatched_cell_ids_explain_the_seurat_relabel_trap(tmp_path):
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    before.write_text("AAACCTG-1\nAAACGGT-1\n")
    after.write_text("AAACCTG-1_1\n")      # a Seurat merge suffix

    with pytest.raises(ValueError, match="survive the round trip unchanged"):
        record_filter(
            input_cell_ids=str(before),
            output_cell_ids=str(after),
            output_dir=str(tmp_path / "out"),
            filter_used="QC",
        )


def test_each_side_needs_exactly_one_source(tmp_path, tiny_adata):
    input_path = tmp_path / "input.h5ad"
    tiny_adata.write_h5ad(input_path)
    cells = tmp_path / "cells.csv"
    cells.write_text("cell0\n")

    with pytest.raises(ValueError, match="exactly one of output_h5ad"):
        record_filter(input_h5ad=str(input_path), output_dir=str(tmp_path / "o"),
                      filter_used="QC")
    with pytest.raises(ValueError, match="exactly one of input_h5ad"):
        record_filter(input_h5ad=str(input_path), input_cell_ids=str(cells),
                      output_cell_ids=str(cells),
                      output_dir=str(tmp_path / "o"), filter_used="QC")


def test_duplicate_cell_ids_are_rejected(tmp_path):
    dupes = tmp_path / "dupes.csv"
    dupes.write_text("cell0\ncell1\ncell0\n")
    with pytest.raises(ValueError, match="duplicate cell IDs"):
        record_filter(input_cell_ids=str(dupes), output_cell_ids=str(dupes),
                      output_dir=str(tmp_path / "o"), filter_used="QC")


def test_cli_accepts_cell_ids_and_omitted_embedding(tmp_path):
    from sccairns.record_filter import build_parser

    args = build_parser().parse_args([
        "--input-cell-ids", "before.csv", "--output-cell-ids", "after.csv",
        "--output-dir", str(tmp_path), "--filter-used", "QC",
    ])
    assert args.embedding_key is None
    assert args.input_h5ad is None and args.output_h5ad is None


def test_cli_rejects_both_forms_for_one_side(tmp_path):
    from sccairns.record_filter import build_parser

    with pytest.raises(SystemExit):
        build_parser().parse_args([
            "--input-h5ad", "a.h5ad", "--input-cell-ids", "a.csv",
            "--output-cell-ids", "b.csv",
            "--output-dir", str(tmp_path), "--filter-used", "QC",
        ])

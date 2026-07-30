import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from sccairns.config import DEFAULT_CONFIG  # noqa: E402
from sccairns.inspect import (  # noqa: E402
    apply_decisions,
    apply_inspection_cli_overrides,
    cluster_qc_summary,
    find_sweep_architecture_keys,
    variant_name_from_cluster_key,
    write_auto_flags_yaml,
)


def _adata():
    adata = anndata.AnnData(np.ones((4, 3)))
    adata.obs_names = [f"cell{i}" for i in range(4)]
    adata.var_names = [f"gene{i}" for i in range(3)]
    adata.obs["leiden"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obs["data_origin"] = pd.Categorical(["a", "a", "b", "b"])
    adata.obs["n_genes_by_counts"] = [600, 700, 300, 800]
    return adata


def _inspection_args(**overrides):
    """Namespace with every field apply_inspection_cli_overrides reads (all None,
    i.e. 'not provided', unless overridden)."""
    import argparse
    fields = dict(
        input=None, output_dir=None, batch_key=None, covariate_keys=None,
        cluster_key=None, latent_key=None, umap_key=None,
        all_sweep_architectures=False, markers=None, neuronal_markers=None,
        neuronal_cutoff=None, marker_threshold=None, pca_color_gene=None,
        mt_threshold=None, min_genes_threshold=None, min_cells=None,
        single_batch_threshold=None,
    )
    fields.update(overrides)
    return argparse.Namespace(**fields)


def test_all_architectures_config_field_default_and_cli_override():
    import copy
    # The config field exists and defaults to off.
    assert DEFAULT_CONFIG["inspection"]["all_architectures"] is False

    # The --all-sweep-architectures flag forces the config toggle on.
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    apply_inspection_cli_overrides(cfg, _inspection_args(all_sweep_architectures=True))
    assert cfg["inspection"]["all_architectures"] is True

    # Without the flag, a config value set to true is left untouched (config wins).
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["inspection"]["all_architectures"] = True
    apply_inspection_cli_overrides(cfg, _inspection_args(all_sweep_architectures=False))
    assert cfg["inspection"]["all_architectures"] is True


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


def test_variant_name_from_cluster_key():
    assert variant_name_from_cluster_key("leiden") == "default"
    assert variant_name_from_cluster_key("leiden_small_gene_nb") == "small_gene_nb"
    assert variant_name_from_cluster_key("leiden_harmony") == "harmony"
    assert variant_name_from_cluster_key(None) is None


def test_auto_flags_records_integration_variant(tmp_path):
    import yaml

    summary = pd.DataFrame({"n_cells": [120, 80]}, index=["3", "7"])
    out = tmp_path / "auto_flags.yaml"
    write_auto_flags_yaml(
        {"7": ["high mt"]},
        summary,
        str(out),
        {"mt": 15, "min_genes": 400, "min_cells": 20, "single_batch": 0.9},
        integration={
            "variant": "small_gene_nb",
            "cluster_key": "leiden_small_gene_nb",
            "latent_key": "X_scVI_small_gene_nb",
            "umap_key": "X_umap_small_gene_nb",
        },
    )
    parsed = yaml.safe_load(out.read_text())
    assert parsed["integration"]["variant"] == "small_gene_nb"
    assert parsed["integration"]["cluster_key"] == "leiden_small_gene_nb"
    assert parsed["integration"]["latent_key"] == "X_scVI_small_gene_nb"
    # The block must not disturb the rest of the editable template.
    assert parsed["remove_clusters"][0]["cluster"] == "7"


def test_apply_decisions_records_filtered_on(tmp_path):
    import json

    adata = _adata()
    adata.obs["leiden_small_gene_nb"] = pd.Categorical(["0", "0", "7", "7"])
    filtered_on = {
        "variant": "small_gene_nb",
        "cluster_key": "leiden_small_gene_nb",
        "latent_key": "X_scVI_small_gene_nb",
        "umap_key": "X_umap_small_gene_nb",
        "source": "decisions_file",
    }
    apply_decisions(
        adata,
        {"remove_clusters": [{"cluster": "7", "reasons": ["test"]}]},
        cluster_key="leiden_small_gene_nb",
        batch_key="data_origin",
        output_dir=str(tmp_path),
        filtered_on=filtered_on,
    )
    manifest = json.loads((tmp_path / "round_manifest.json").read_text())
    assert manifest["decisions"]["filtered_on"] == filtered_on
    assert manifest["decisions"]["cluster_key"] == "leiden_small_gene_nb"


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


def test_find_sweep_architecture_keys_includes_harmony():
    adata = _adata()
    for name in ["small", "medium"]:
        adata.obs[f"leiden_{name}"] = pd.Categorical(["0", "0", "1", "1"])
        adata.obsm[f"X_scVI_{name}"] = np.ones((4, 2))
        adata.obsm[f"X_umap_{name}"] = np.ones((4, 2))
    adata.obs["leiden_harmony"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obsm["X_pca_harmony"] = np.ones((4, 3))
    adata.obsm["X_umap_harmony"] = np.ones((4, 2))

    config = {
        "integration": {
            "sweep": [
                {"name": "medium"},
                {"name": "small"},
            ]
        }
    }

    runs = find_sweep_architecture_keys(adata, config)

    assert [run["name"] for run in runs] == ["medium", "small", "harmony"]
    assert runs[-1]["cluster_key"] == "leiden_harmony"
    assert runs[-1]["latent_key"] == "X_pca_harmony"
    assert runs[-1]["umap_key"] == "X_umap_harmony"


def test_find_sweep_architecture_keys_allows_harmony_only():
    adata = _adata()
    adata.obs["leiden_harmony"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obsm["X_pca_harmony"] = np.ones((4, 3))
    adata.obsm["X_umap_harmony"] = np.ones((4, 2))

    runs = find_sweep_architecture_keys(adata, {"integration": {}})

    assert runs == [
        {
            "name": "harmony",
            "cluster_key": "leiden_harmony",
            "latent_key": "X_pca_harmony",
            "umap_key": "X_umap_harmony",
        }
    ]

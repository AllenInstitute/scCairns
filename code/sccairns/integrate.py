#!/usr/bin/env python3
"""scVI integration with optional parameter sweep and scIB benchmarking.

Merges three workflows into a single CLI script:
  1. Default single-model integration (from 01_20260513 combined SNS notebook)
  2. Parameter sweep over multiple architectures (from 20260325 refinement notebook)
  3. Annotation-free scIB benchmarking (scib_metrics)

The input is a pre-concatenated h5ad with raw counts (in .X or .layers["counts"])
and at minimum a batch column (e.g. "data_origin") in .obs.

Usage
-----
    # Default single-model integration (large architecture)
    python integrate_scvi.py --input combined.h5ad

    # Override default architecture
    python integrate_scvi.py --input combined.h5ad \
        --n-hidden 128 --n-layers 2 --gene-likelihood zinb

    # Parameter sweep with built-in configs
    python integrate_scvi.py --input combined.h5ad --sweep

    # Parameter sweep with custom configs (JSON)
    python integrate_scvi.py --input combined.h5ad --sweep \
        --sweep-configs my_configs.json

    # Enable bio-conservation metrics using a proxy label
    python integrate_scvi.py --input combined.h5ad --sweep \
        --bench-label-key ganglion_group
"""

import argparse
import gc
import inspect
import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns

from .config import (
    ConfigError,
    collect_package_versions,
    fingerprint_file,
    load_pipeline_config,
    read_parent_provenance,
    set_global_seed,
    set_if_provided,
    update_round_manifest,
    validate_config,
    write_command_args,
    write_resolved_config,
)


# ═══════════════════════════════════════════════════════════════════════════════
#  Configuration
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ModelConfig:
    """Configuration for a single scVI model run."""
    name: str
    n_hidden: int = 256
    n_layers: int = 3
    n_latent: int = 32
    dispersion: str = "gene-cell"
    gene_likelihood: str = "nb"
    batch_key: str = "data_origin"
    categorical_covariate_keys: List[str] = field(default_factory=list)
    continuous_covariate_keys: List[str] = field(default_factory=list)
    counts_layer: str = "counts"
    hvg_batch_key: str = "tech"
    hvg_flavor: str = "seurat_v3"
    hvg_nbatches: Optional[int] = None  # min batches for HVG; None = no filter


# Built-in sweep configurations (from 20260325 sympathetic refinement notebook)
SWEEP_CONFIGS: Dict[str, dict] = {
    "custom2_tech": dict(
        n_hidden=16, n_layers=2, n_latent=32,
        dispersion="gene-cell", gene_likelihood="nb",
        batch_key="tech", categorical_covariate_keys=[],
        hvg_batch_key="tech", hvg_flavor="seurat", hvg_nbatches=3,
    ),
    "custom2_platform": dict(
        n_hidden=16, n_layers=2, n_latent=32,
        dispersion="gene-cell", gene_likelihood="nb",
        batch_key="platform_origin", categorical_covariate_keys=[],
        hvg_batch_key="platform_origin", hvg_flavor="seurat", hvg_nbatches=3,
    ),
    "medium_platform": dict(
        n_hidden=128, n_layers=2, n_latent=32,
        dispersion="gene-cell", gene_likelihood="nb",
        batch_key="platform_origin", categorical_covariate_keys=[],
        hvg_batch_key="platform_origin", hvg_flavor="seurat", hvg_nbatches=3,
    ),
    "medium_platform_cov": dict(
        n_hidden=128, n_layers=2, n_latent=32,
        dispersion="gene-cell", gene_likelihood="nb",
        batch_key="platform_origin", categorical_covariate_keys=["data_origin"],
        hvg_batch_key="platform_origin", hvg_flavor="seurat", hvg_nbatches=3,
    ),
}


# ═══════════════════════════════════════════════════════════════════════════════
#  QC helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _normalise_gene_names_for_matching(var_names, gene_symbol_case="preserve"):
    """Return gene names transformed for configurable QC pattern matching."""
    if gene_symbol_case == "lower":
        return var_names.str.lower()
    if gene_symbol_case == "upper":
        return var_names.str.upper()
    return var_names


def should_apply_qc_filter(qc_cfg, is_first_round):
    """Decide whether the QC gene/cell filter runs for this round.

    The filter runs when QC is enabled and not globally skipped. Additionally,
    ``skip_filter_after_round_1`` (default True) suppresses it on any round >1
    (``is_first_round`` False), because re-filtering genes against a later
    round's smaller cell subset can drop rare markers the first full-data round
    correctly kept.

    Returns (apply_filter, skipped_because_later_round).
    """
    globally_off = (not qc_cfg.get("enabled", True)) or qc_cfg.get("skip_filter", False)
    skip_later = (qc_cfg.get("skip_filter_after_round_1", True)
                  and not is_first_round)
    apply_filter = (not globally_off) and (not skip_later)
    # "skipped because later round" only when it would otherwise have run.
    skipped_later = skip_later and not globally_off
    return apply_filter, skipped_later


def compute_qc_metrics(
    adata,
    mt_gene_patterns=None,
    ribo_gene_patterns=None,
    hb_gene_pattern="^hb[^(p)]",
    gene_symbol_case="preserve",
):
    """Annotate mt/ribo/hb genes and compute QC metrics in-place."""
    mt_gene_patterns = mt_gene_patterns or ["mt-"]
    ribo_gene_patterns = ribo_gene_patterns or ["rps", "rpl"]
    names = _normalise_gene_names_for_matching(adata.var_names, gene_symbol_case)
    adata.var["mt"] = names.str.startswith(tuple(mt_gene_patterns))
    adata.var["ribo"] = names.str.startswith(tuple(ribo_gene_patterns))
    adata.var["hb"] = names.str.contains(hb_gene_pattern, regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True,
    )


def plot_qc_violins(adata, output_dir, groupby="data_origin", suffix=""):
    """QC violin plots grouped by a batch variable."""
    if groupby not in adata.obs.columns:
        print(f"  [WARN] '{groupby}' not in obs — skipping QC violin plots.")
        return
    qc_cols = [c for c in ["n_genes_by_counts", "total_counts",
                            "pct_counts_mt", "log1p_total_counts"]
               if c in adata.obs.columns]
    if not qc_cols:
        print("  No QC columns found — skipping violin plots.")
        return

    n = len(qc_cols)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, col in zip(axes, qc_cols):
        sc.pl.violin(adata, keys=col, groupby=groupby, ax=ax,
                     show=False, rotation=45)
        ax.set_title(col)

    tag = f"_{suffix}" if suffix else ""
    fig.suptitle(f"QC — {adata.n_obs:,} cells", fontsize=14, y=1.02)
    fig.tight_layout()
    fname = f"qc_violins{tag}.png"
    fig.savefig(os.path.join(output_dir, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")

    # Cell counts per study
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    counts = adata.obs[groupby].value_counts()
    counts.plot(kind="barh", ax=ax2, color="steelblue", edgecolor="black")
    ax2.set_xlabel("Number of cells")
    ax2.set_title(f"Cells per {groupby}")
    for i, (label, count) in enumerate(counts.items()):
        ax2.text(count + 50, i, f"{count:,}", va="center", fontsize=9)
    fig2.tight_layout()
    fig2.savefig(os.path.join(output_dir, f"cells_per_{groupby}{tag}.png"),
                 dpi=150, bbox_inches="tight")
    plt.close(fig2)


# ═══════════════════════════════════════════════════════════════════════════════
#  Data preparation
# ═══════════════════════════════════════════════════════════════════════════════

def ensure_counts_layer(adata, counts_layer="counts"):
    """Ensure raw integer counts live in the configured layer and .X."""
    from scipy.sparse import issparse

    def _looks_like_counts(mat, n_check=200):
        sample = mat[:n_check]
        if issparse(sample):
            sample = sample.toarray()
        vals = sample[sample > 0]
        return len(vals) > 0 and np.allclose(vals, np.round(vals))

    if counts_layer in adata.layers and _looks_like_counts(adata.layers[counts_layer]):
        print(f"  Using existing '{counts_layer}' layer (integer-valued).")
        adata.X = adata.layers[counts_layer].copy()
        return

    if _looks_like_counts(adata.X):
        print(f"  .X contains raw counts — copying to layers['{counts_layer}'].")
        adata.layers[counts_layer] = adata.X.copy()
        return

    # Fallback: assume X is usable even if not perfectly integer
    print("  [WARN] .X does not look like raw integer counts.")
    if counts_layer not in adata.layers:
        print(f"         No '{counts_layer}' layer found. Using .X as-is — verify input data.")
        adata.layers[counts_layer] = adata.X.copy()
    else:
        print(f"         Using existing '{counts_layer}' layer despite non-integer values.")
        adata.X = adata.layers[counts_layer].copy()


def select_hvgs(adata, n_top_genes=3000, batch_key="tech", flavor="seurat_v3",
                hvg_nbatches=None, counts_layer="counts"):
    """Select HVGs on a temporary normalised copy; return gene name list.

    The caller's adata is not modified.
    """
    tmp = adata.copy()
    tmp.X = tmp.layers[counts_layer].copy()

    if flavor == "seurat":
        # seurat flavor expects log-normalised .X (internally calls expm1)
        sc.pp.normalize_total(tmp)
        sc.pp.log1p(tmp)
    # seurat_v3 / cell_ranger expect raw counts — leave as-is

    sc.pp.highly_variable_genes(
        tmp, n_top_genes=n_top_genes,
        batch_key=batch_key, flavor=flavor, subset=False,
    )

    if hvg_nbatches is not None and "highly_variable_nbatches" in tmp.var.columns:
        mask = (tmp.var["highly_variable"]
                & (tmp.var["highly_variable_nbatches"] >= hvg_nbatches))
    else:
        mask = tmp.var["highly_variable"]

    hvgs = tmp.var_names[mask].tolist()
    del tmp; gc.collect()
    return hvgs


# ═══════════════════════════════════════════════════════════════════════════════
#  scVI training
# ═══════════════════════════════════════════════════════════════════════════════

def clean_batch_column(adata, col):
    """Fill NaN values in a categorical/batch column with 'other'."""
    if col not in adata.obs.columns:
        return
    n_nan = adata.obs[col].isna().sum()
    if n_nan == 0:
        return
    print(f"  [WARN] {n_nan} NaN in '{col}' — filling with 'other'.")
    if adata.obs[col].dtype.name != "category":
        adata.obs[col] = adata.obs[col].astype("category")
    if "other" not in adata.obs[col].cat.categories:
        adata.obs[col] = adata.obs[col].cat.add_categories("other")
    adata.obs[col] = adata.obs[col].fillna("other")


def train_with_num_workers(model, train_kwargs, num_workers=None):
    """Train with scvi-tools-version-compatible DataLoader worker kwargs."""
    if num_workers is None:
        model.train(**train_kwargs)
        return

    for kw_name in ("datasplitter_kwargs", "datamodule_kwargs"):
        try:
            model.train(**train_kwargs, **{kw_name: {"num_workers": num_workers}})
            return
        except TypeError:
            continue

    print("  [WARN] Could not set num_workers, training with default.")
    model.train(**train_kwargs)


def setup_and_train(adata_hvg, config, max_epochs=200,
                    early_stopping_patience=20, batch_size=256,
                    num_workers=None, output_dir=None):
    """Setup AnnData, build, and train one scVI model. Returns the model."""
    clean_batch_column(adata_hvg, config.batch_key)

    # Filter covariate keys to columns actually present
    cat_covs = [c for c in config.categorical_covariate_keys
                if c in adata_hvg.obs.columns]
    cont_covs = [c for c in config.continuous_covariate_keys
                 if c in adata_hvg.obs.columns]
    for c in cat_covs:
        clean_batch_column(adata_hvg, c)

    scvi.model.SCVI.setup_anndata(
        adata_hvg, layer=config.counts_layer,
        batch_key=config.batch_key,
        categorical_covariate_keys=cat_covs or None,
        continuous_covariate_keys=cont_covs or None,
    )

    model = scvi.model.SCVI(
        adata_hvg,
        n_hidden=config.n_hidden,
        n_layers=config.n_layers,
        n_latent=config.n_latent,
        dispersion=config.dispersion,
        gene_likelihood=config.gene_likelihood,
    )

    print(f"  Architecture: n_hidden={config.n_hidden}, n_layers={config.n_layers}, "
          f"n_latent={config.n_latent}")
    print(f"  Likelihood: {config.gene_likelihood}, dispersion: {config.dispersion}")
    print(f"  batch_key='{config.batch_key}', categorical_covariates={cat_covs or 'none'}")
    print(f"  continuous_covariates={cont_covs or 'none'}")

    train_kwargs = dict(
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=0.9,
        early_stopping=True,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor="elbo_validation",
        check_val_every_n_epoch=1,
    )
    train_with_num_workers(model, train_kwargs, num_workers=num_workers)

    n_epochs = len(model.history["elbo_train"])
    final_val = model.history["elbo_validation"].iloc[-1]
    if hasattr(final_val, "values"):
        final_val = final_val.values[0]
    print(f"  Trained {n_epochs} epochs — final val ELBO: {final_val:.1f}")

    if output_dir:
        model_dir = os.path.join(output_dir, f"scvi_model_{config.name}")
        os.makedirs(model_dir, exist_ok=True)
        if "cell_id" in model.adata.obs.columns:
            model.adata.obs["cell_id"] = model.adata.obs["cell_id"].astype(str)
        model.save(model_dir, save_anndata=True, overwrite=True)
        print(f"  Model saved → {model_dir}")

    return model


def resolve_harmony_hvg_spec(harmony_cfg, integration_cfg, configs, sweep_mode):
    """Resolve which HVG selection Harmony should use, as an explicit record.

    Returns a dict describing the HVG source and parameters. By default Harmony
    uses the top-level ``integration.hvg`` spec (identical to single-model scVI).
    In sweep mode, ``harmony.hvg_from`` may name a sweep entry whose HVG
    selection Harmony borrows instead. Raises ConfigError if ``hvg_from`` does
    not match an available sweep config. The returned ``source`` field makes the
    choice visible in the manifest rather than implicit in the code.
    """
    top_hvg = integration_cfg["hvg"]
    hvg_from = harmony_cfg.get("hvg_from")
    if hvg_from:
        sweep_cfg = configs.get(hvg_from) if sweep_mode else None
        if sweep_cfg is None:
            available = list(configs) if sweep_mode else []
            raise ConfigError(
                f"harmony.hvg_from='{hvg_from}' does not match a sweep config. "
                f"Available: {available or 'none (not a sweep run)'}."
            )
        return {
            "source": f"sweep:{hvg_from}",
            "n_top_genes": top_hvg["n_top_genes"],
            "batch_key": sweep_cfg.hvg_batch_key,
            "flavor": sweep_cfg.hvg_flavor,
            "min_batches": sweep_cfg.hvg_nbatches,
        }
    return {
        "source": "integration.hvg",
        "n_top_genes": top_hvg["n_top_genes"],
        "batch_key": top_hvg["batch_key"],
        "flavor": top_hvg["flavor"],
        "min_batches": top_hvg.get("min_batches"),
    }


def run_harmony(adata_hvg, batch_key, counts_layer, n_pcs=30, seed=None):
    """Compute a Harmony-integrated embedding from the HVG counts.

    Builds a standard log-normalized PCA on the HVG subset and runs
    ``scanpy.external.pp.harmony_integrate`` to batch-correct it, returning the
    corrected coordinates (Scanpy stores them in obsm["X_pca_harmony"]). Works
    on a copy so the caller's AnnData is untouched. The seed is threaded into
    PCA and Harmony's internal KMeans for best-effort reproducibility; harmonypy
    is only partially seed-controllable, so the result is deterministic on a
    given platform/library set but not guaranteed bit-identical across them.

    Raises ImportError (with an actionable message) if harmonypy is missing.
    """
    try:
        import scanpy.external as sce  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised via integration
        raise ImportError(
            "Harmony integration requires harmonypy. Install it with "
            "`pip install harmonypy` or add it to the environment."
        ) from exc

    tmp = adata_hvg.copy()
    # Source from raw counts so Harmony's PCA input is independent of whatever
    # .X currently holds.
    tmp.X = tmp.layers[counts_layer].copy()
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.pp.scale(tmp, max_value=10)
    n_comps = min(n_pcs, tmp.n_vars - 1, tmp.n_obs - 1)
    sc.pp.pca(tmp, n_comps=n_comps, random_state=seed if seed is not None else 0)

    clean_batch_column(tmp, batch_key)
    sce.pp.harmony_integrate(
        tmp,
        key=batch_key,
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        random_state=seed if seed is not None else 0,
    )
    return tmp.obsm["X_pca_harmony"]


# ═══════════════════════════════════════════════════════════════════════════════
#  Diagnostic plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_training_curves(models, output_dir):
    """ELBO convergence curves for one or more models, side by side."""
    n = len(models)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, (name, model) in zip(axes, models.items()):
        hist = model.history
        train = hist["elbo_train"].values.flatten()
        val = hist["elbo_validation"].values.flatten()
        epochs = np.arange(1, len(train) + 1)

        skip = min(5, len(train) - 1)
        ax.plot(epochs[skip:], train[skip:], label="Train", linewidth=1.5)
        ax.plot(epochs[skip:], val[skip:], label="Validation",
                linewidth=1.5, linestyle="--")

        ax.annotate(f"Train: {train[-1]:.1f}",
                    xy=(epochs[-1], train[-1]), xytext=(-60, 8),
                    textcoords="offset points", fontsize=8, color="C0")
        ax.annotate(f"Val: {val[-1]:.1f}",
                    xy=(epochs[-1], val[-1]), xytext=(-60, -14),
                    textcoords="offset points", fontsize=8, color="C1")

        max_ep = max(len(h["elbo_train"]) for h in
                     [m.history for m in models.values()])
        if len(train) < max_ep:
            ax.axvline(len(train), color="grey", linestyle=":", linewidth=1,
                       label=f"Stopped @ {len(train)}")

        ax.set_title(name, fontsize=9, fontweight="bold")
        ax.set_xlabel("Epoch")
        ax.set_ylabel("ELBO")
        ax.set_yscale("log")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        delta = val[-1] - val[max(0, len(val) - 20)]
        print(f"  {name:30s}  final_val={val[-1]:.1f}  "
              f"epochs={len(val)}  delta_val(last20)={delta:+.1f}")

    fig.suptitle("scVI training convergence", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "training_convergence.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved training_convergence.png")


def plot_umap_grid(adata, umap_keys, leiden_keys, color_vars, output_dir,
                   cell_size=8, dpi=150):
    """UMAP comparison grid: configs (rows) x colour variables (columns).

    Used in sweep mode to visually compare integration quality across configs.
    """
    n_rows = len(umap_keys)
    n_cols = len(color_vars)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5 * n_cols, 4 * n_rows),
                             squeeze=False)

    for row, (name, umap_key) in enumerate(umap_keys.items()):
        for col, color in enumerate(color_vars):
            ax = axes[row][col]
            resolved = (leiden_keys.get(name, color)
                        if color == "leiden" else color)
            if (resolved not in adata.obs.columns
                    and resolved not in adata.var_names):
                ax.set_visible(False)
                continue
            sc.pl.embedding(
                adata, basis=umap_key, color=resolved, size=cell_size,
                legend_loc="on data" if color == "leiden" else "right margin",
                legend_fontsize=6, frameon=False, show=False, ax=ax,
            )
            if row == 0:
                ax.set_title(resolved, fontsize=10, fontweight="bold")
            if col == 0:
                ax.set_ylabel(name, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "umap_config_comparison.png"),
                dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("  Saved umap_config_comparison.png")


def plot_single_umaps(adata, output_dir, basis="X_umap", output_name="umap_integration.png"):
    """Standard UMAP panels coloured by available metadata (single-model mode)."""
    candidates = ["data_origin", "tech", "platform_origin",
                   "ganglion_group", "ganglion", "leiden"]
    colors = [c for c in candidates if c in adata.obs.columns]
    if not colors:
        return
    n = len(colors)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, c in zip(axes, colors):
        sc.pl.embedding(
            adata, basis=basis, color=c, size=10,
            legend_loc="on data" if c == "leiden" else "right margin",
            legend_fontsize=7, frameon=False, show=False, ax=ax,
        )
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, output_name),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {output_name}")


_MISSING = object()


def _resolve_neighbor_graph_keys(adata, neighbors_key):
    """Return obsp keys backing a Scanpy neighbors graph."""
    neighbors = adata.uns.get(neighbors_key, {})
    if not isinstance(neighbors, dict):
        neighbors = {}

    connectivities_key = neighbors.get("connectivities_key")
    distances_key = neighbors.get("distances_key")

    if not connectivities_key:
        candidates = ["connectivities"] if neighbors_key == "neighbors" else []
        candidates.append(f"{neighbors_key}_connectivities")
        connectivities_key = next(
            (key for key in candidates if key in adata.obsp),
            candidates[0],
        )

    if not distances_key:
        candidates = ["distances"] if neighbors_key == "neighbors" else []
        candidates.append(f"{neighbors_key}_distances")
        distances_key = next(
            (key for key in candidates if key in adata.obsp),
            candidates[0],
        )

    if connectivities_key not in adata.obsp:
        raise KeyError(
            f'No "{connectivities_key}" in .obsp for neighbors_key='
            f'"{neighbors_key}"'
        )

    return connectivities_key, distances_key


def _restore_mapping_key(mapping, key, prior_value):
    if prior_value is _MISSING:
        if key in mapping:
            del mapping[key]
    else:
        mapping[key] = prior_value


def _normalise_neighbor_graph_metadata(adata, neighbors_key):
    """Ensure a named neighbors graph records its backing obsp keys."""
    connectivities_key, distances_key = _resolve_neighbor_graph_keys(
        adata, neighbors_key
    )
    neighbors = adata.uns.get(neighbors_key, {})
    neighbors = dict(neighbors) if isinstance(neighbors, dict) else {}
    neighbors["connectivities_key"] = connectivities_key
    if distances_key in adata.obsp:
        neighbors["distances_key"] = distances_key
    adata.uns[neighbors_key] = neighbors


def run_neighbors_with_key(adata, neighbors_key, seed=None, **kwargs):
    """Run neighbors and keep named graph metadata consistent across versions.

    ``seed`` is threaded into Scanpy's ``random_state`` so the configured
    ``reproducibility.seed`` governs the graph rather than Scanpy's implicit
    default. ``None`` normalizes to 0 — same convention as ``run_harmony`` — so
    an unseeded run is still deterministic here, and the manifest's recorded
    seed fully describes this step.
    """
    kwargs.setdefault("random_state", seed if seed is not None else 0)
    try:
        sc.pp.neighbors(adata, key_added=neighbors_key, **kwargs)
        _normalise_neighbor_graph_metadata(adata, neighbors_key)
        return
    except TypeError as e:
        if "key_added" not in str(e):
            raise
        print(
            "  [WARN] sc.pp.neighbors does not support key_added; "
            "using fallback."
        )

    prior_neighbors = adata.uns["neighbors"] if "neighbors" in adata.uns else _MISSING
    prior_connectivities = (
        adata.obsp["connectivities"] if "connectivities" in adata.obsp else _MISSING
    )
    prior_distances = (
        adata.obsp["distances"] if "distances" in adata.obsp else _MISSING
    )

    sc.pp.neighbors(adata, **kwargs)
    _normalise_neighbor_graph_metadata(adata, "neighbors")

    if neighbors_key == "neighbors":
        return

    connectivities_key, distances_key = _resolve_neighbor_graph_keys(
        adata, "neighbors"
    )
    target_connectivities_key = f"{neighbors_key}_connectivities"
    target_distances_key = f"{neighbors_key}_distances"
    source_neighbors = adata.uns.get("neighbors", {})
    target_neighbors = (
        dict(source_neighbors) if isinstance(source_neighbors, dict) else {}
    )
    target_neighbors["connectivities_key"] = target_connectivities_key

    adata.obsp[target_connectivities_key] = adata.obsp[connectivities_key]
    if distances_key in adata.obsp:
        adata.obsp[target_distances_key] = adata.obsp[distances_key]
        target_neighbors["distances_key"] = target_distances_key
    adata.uns[neighbors_key] = target_neighbors

    _restore_mapping_key(adata.uns, "neighbors", prior_neighbors)
    _restore_mapping_key(adata.obsp, "connectivities", prior_connectivities)
    _restore_mapping_key(adata.obsp, "distances", prior_distances)


def run_umap_with_key(adata, neighbors_key, umap_key, spread=1.0, min_dist=0.5,
                      seed=None):
    """Run UMAP with a requested obsm key across Scanpy versions.

    Newer Scanpy versions support sc.tl.umap(..., key_added=...). Some older
    versions always write obsm["X_umap"]. This wrapper tries key_added first,
    then falls back to copying/restoring X_umap so multiple embeddings can
    coexist safely.

    ``seed`` is threaded into ``random_state`` on both the primary and fallback
    paths (``None`` -> 0), so the layout does not depend on which path a given
    Scanpy version takes.
    """
    random_state = seed if seed is not None else 0
    try:
        sc.tl.umap(
            adata,
            neighbors_key=neighbors_key,
            key_added=umap_key,
            spread=spread,
            min_dist=min_dist,
            random_state=random_state,
        )
        return
    except TypeError as e:
        unsupported_args = [
            arg for arg in ("key_added", "neighbors_key") if arg in str(e)
        ]
        if not unsupported_args:
            raise
        print(
            "  [WARN] sc.tl.umap does not support "
            f"{', '.join(unsupported_args)}; using fallback."
        )
    except KeyError as e:
        if "connectivities" not in str(e) and "distances" not in str(e):
            raise
        print(
            "  [WARN] sc.tl.umap could not resolve the requested neighbors "
            "graph; using fallback."
        )

    prior_x_umap = adata.obsm["X_umap"].copy() if "X_umap" in adata.obsm else None
    prior_neighbors = adata.uns["neighbors"] if "neighbors" in adata.uns else _MISSING
    prior_connectivities = (
        adata.obsp["connectivities"] if "connectivities" in adata.obsp else _MISSING
    )
    prior_distances = (
        adata.obsp["distances"] if "distances" in adata.obsp else _MISSING
    )

    connectivities_key, distances_key = _resolve_neighbor_graph_keys(
        adata, neighbors_key
    )
    source_neighbors = adata.uns.get(neighbors_key, {})
    legacy_neighbors = (
        dict(source_neighbors) if isinstance(source_neighbors, dict) else {}
    )
    legacy_neighbors["connectivities_key"] = "connectivities"

    try:
        adata.obsp["connectivities"] = adata.obsp[connectivities_key]
        if distances_key in adata.obsp:
            adata.obsp["distances"] = adata.obsp[distances_key]
            legacy_neighbors["distances_key"] = "distances"
        adata.uns["neighbors"] = legacy_neighbors

        sc.tl.umap(
            adata,
            spread=spread,
            min_dist=min_dist,
            random_state=random_state,
        )
    finally:
        _restore_mapping_key(adata.uns, "neighbors", prior_neighbors)
        _restore_mapping_key(adata.obsp, "connectivities", prior_connectivities)
        _restore_mapping_key(adata.obsp, "distances", prior_distances)

    if umap_key != "X_umap":
        adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        if prior_x_umap is not None:
            adata.obsm["X_umap"] = prior_x_umap
        else:
            del adata.obsm["X_umap"]


def _call_leiden(adata, leiden_key, resolution=1.0, flavor="igraph", **graph_kwargs):
    """Call sc.tl.leiden while tolerating older optional kwargs."""
    kwargs = dict(graph_kwargs)
    kwargs["resolution"] = resolution
    kwargs["key_added"] = leiden_key
    if flavor is not None:
        kwargs["flavor"] = flavor

    copy_default_key = False
    prior_leiden = _MISSING
    while True:
        try:
            sc.tl.leiden(adata, **kwargs)
            break
        except TypeError as e:
            message = str(e)
            if "flavor" in message and "flavor" in kwargs:
                print("  [WARN] sc.tl.leiden does not support flavor; omitting it.")
                del kwargs["flavor"]
                continue
            if "key_added" in message and "key_added" in kwargs:
                print(
                    "  [WARN] sc.tl.leiden does not support key_added; "
                    "using fallback."
                )
                prior_leiden = (
                    adata.obs["leiden"].copy()
                    if "leiden" in adata.obs.columns
                    else _MISSING
                )
                del kwargs["key_added"]
                copy_default_key = leiden_key != "leiden"
                continue
            raise

    if copy_default_key:
        adata.obs[leiden_key] = adata.obs["leiden"].copy()
        _restore_mapping_key(adata.obs, "leiden", prior_leiden)


def run_leiden_with_key(
    adata,
    neighbors_key,
    leiden_key,
    resolution=1.0,
    flavor="igraph",
    seed=None,
):
    """Run Leiden using a named neighbors graph across Scanpy versions.

    ``seed`` is threaded into ``random_state`` (``None`` -> 0). This matters
    more than for the other steps: Leiden cluster IDs are what decisions.yaml
    files reference by number, so the recorded seed must describe them.
    """
    random_state = seed if seed is not None else 0
    try:
        _call_leiden(
            adata,
            neighbors_key=neighbors_key,
            resolution=resolution,
            flavor=flavor,
            leiden_key=leiden_key,
            random_state=random_state,
        )
        return
    except (KeyError, TypeError, ValueError) as e:
        message = str(e)
        graph_error = (
            "pp.neighbors" in message
            or "neighborhood graph" in message
            or "connectivities" in message
            or "neighbors_key" in message
        )
        if not graph_error:
            raise
        print(
            "  [WARN] sc.tl.leiden could not resolve the requested neighbors "
            "graph; using adjacency fallback."
        )

    connectivities_key, _ = _resolve_neighbor_graph_keys(adata, neighbors_key)
    _call_leiden(
        adata,
        adjacency=adata.obsp[connectivities_key],
        resolution=resolution,
        flavor=flavor,
        leiden_key=leiden_key,
        random_state=random_state,
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  scIB benchmarking
# ═══════════════════════════════════════════════════════════════════════════════

def _supported_metric_kwargs(cls, desired, aliases=None, name=None):
    """Return kwargs supported by a scib-metrics config class signature."""
    aliases = aliases or {}
    label = name or getattr(cls, "__name__", "metric config")
    sig = inspect.signature(cls)
    params = sig.parameters
    supports_kwargs = any(
        p.kind == inspect.Parameter.VAR_KEYWORD
        for p in params.values()
    )

    resolved = {}
    skipped = []
    for key, value in desired.items():
        candidates = aliases.get(key, [key])
        target = None
        for candidate in candidates:
            if supports_kwargs or candidate in params:
                target = candidate
                break
        if target is None:
            skipped.append(key)
            continue
        resolved[target] = value

    if skipped:
        print(f"  [WARN] {label} does not support: {', '.join(skipped)}")
    return resolved


def run_benchmarking(adata, embedding_keys, bench_batch_key,
                     bench_label_key=None, output_dir=None, n_pca_comps=50,
                     counts_layer="counts", seed=None):
    """Run scIB benchmarking across embedding keys.

    Parameters
    ----------
    bench_label_key : str or None
        If None, only batch-correction metrics are computed and a leiden proxy
        is generated internally so the Benchmarker constructor is satisfied.
        If provided, a subset of bio-conservation metrics is also enabled.
    seed : int or None
        Threaded into the PCA baseline and the leiden proxy (``None`` -> 0).
        Both feed the reported metrics, so the recorded seed has to cover them
        for scib_benchmark_results.csv to be reproducible from the manifest.
    """
    try:
        from scib_metrics.benchmark import (Benchmarker, BatchCorrection,
                                            BioConservation)
    except ImportError:
        print("  [SKIP] scib_metrics not installed — run:")
        print("         pip install scib-metrics")
        return None, None

    # PCA baseline (pre-integration reference)
    if "X_pca" not in adata.obsm:
        print("  Computing PCA baseline for benchmarking...")
        adata_tmp = adata.copy()
        adata_tmp.X = adata_tmp.layers[counts_layer].copy()
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)
        sc.pp.pca(adata_tmp, n_comps=n_pca_comps,
                  random_state=seed if seed is not None else 0)
        adata.obsm["X_pca"] = adata_tmp.obsm["X_pca"]
        del adata_tmp; gc.collect()

    # Label key handling
    if bench_label_key and bench_label_key in adata.obs.columns:
        label_key = bench_label_key
        use_bio = True
        print(f"  Label key: '{label_key}' — bio-conservation metrics enabled.")
    else:
        if bench_label_key:
            print(f"  [WARN] '{bench_label_key}' not in obs — falling back to "
                  f"leiden proxy.")
        # Generate a proxy so the Benchmarker constructor doesn't error
        if "leiden_proxy" not in adata.obs.columns:
            run_neighbors_with_key(adata, "_proxy_neighbors", use_rep="X_pca",
                                   seed=seed)
            run_leiden_with_key(
                adata,
                neighbors_key="_proxy_neighbors",
                leiden_key="leiden_proxy",
                resolution=0.5,
                flavor="igraph",
                seed=seed,
            )
        label_key = "leiden_proxy"
        use_bio = False
        print(f"  No annotation label — using '{label_key}' as proxy "
              f"(bio-conservation metrics disabled).")

    # Metric configuration
    if use_bio:
        bio = BioConservation(
            **_supported_metric_kwargs(
                BioConservation,
                {
                    "isolated_labels": False,
                    "silhouette_label": True,
                    "clisi_knn": True,
                    "nmi_ari_cluster_labels_kmeans": False,
                    "nmi_ari_cluster_labels_leiden": False,
                },
                name="BioConservation",
            )
        )
    else:
        bio = None

    batch = BatchCorrection(
        **_supported_metric_kwargs(
            BatchCorrection,
            {
                # scib-metrics <=0.5.1 used silhouette_batch; newer versions
                # use BRAS for the batch-silhouette-style metric.
                "batch_silhouette": True,
                "graph_connectivity": True,
                "pcr_comparison": True,
                "ilisi_knn": False,
                "kbet_per_label": False,
            },
            aliases={"batch_silhouette": ["silhouette_batch", "bras"]},
            name="BatchCorrection",
        )
    )

    benchmark_adata = adata
    if not adata.obs_names.is_unique:
        print("  [WARN] Observation names are not unique; using a benchmark copy with unique names.")
        benchmark_adata = adata.copy()
        benchmark_adata.obs_names_make_unique()

    bm = Benchmarker(
        benchmark_adata,
        batch_key=bench_batch_key,
        label_key=label_key,
        embedding_obsm_keys=embedding_keys,
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=bio,
        batch_correction_metrics=batch,
        progress_bar=True,
    )

    print("  Running benchmark...")
    bm.benchmark()

    # Save results
    if output_dir:
        try:
            result_fig = bm.plot_results_table(save_dir=output_dir)
            if result_fig is not None:
                plt.close(result_fig)
        except Exception as e:
            print(f"  [WARN] plot_results_table failed: {e}")
            # Fallback: manual save
            try:
                result_fig = bm.plot_results_table()
                result_fig.savefig(
                    os.path.join(output_dir, "scib_results.png"),
                    dpi=150, bbox_inches="tight")
                plt.close(result_fig)
            except Exception:
                pass
        print("  Saved scIB results table.")

    df = bm.get_results(min_max_scale=False)
    if output_dir:
        df.to_csv(os.path.join(output_dir, "scib_benchmark_results.csv"))
        print("  Saved scib_benchmark_results.csv")
    print("\n  Benchmark results:")
    print(df.to_string())

    return bm, df


def apply_integration_cli_overrides(config, args):
    """Apply backwards-compatible CLI overrides onto the resolved config."""
    set_if_provided(config, ["data", "input_h5ad"], args.input)
    set_if_provided(config, ["data", "output_dir"], args.output_dir)
    set_if_provided(config, ["data", "batch_key"], args.batch_key)
    set_if_provided(config, ["data", "categorical_covariate_keys"], args.covariate_keys)
    set_if_provided(config, ["qc", "min_genes"], args.min_genes)
    set_if_provided(config, ["qc", "min_cells"], args.min_cells)
    set_if_provided(config, ["integration", "n_hidden"], args.n_hidden)
    set_if_provided(config, ["integration", "n_layers"], args.n_layers)
    set_if_provided(config, ["integration", "n_latent"], args.n_latent)
    set_if_provided(config, ["integration", "dispersion"], args.dispersion)
    set_if_provided(config, ["integration", "gene_likelihood"], args.gene_likelihood)
    set_if_provided(config, ["integration", "max_epochs"], args.max_epochs)
    set_if_provided(
        config,
        ["integration", "early_stopping_patience"],
        args.early_stopping_patience,
    )
    set_if_provided(config, ["integration", "batch_size"], args.batch_size)
    set_if_provided(config, ["integration", "num_workers"], args.num_workers)
    set_if_provided(config, ["integration", "hvg", "n_top_genes"], args.n_hvgs)
    set_if_provided(config, ["integration", "hvg", "batch_key"], args.hvg_batch_key)
    set_if_provided(config, ["integration", "hvg", "flavor"], args.hvg_flavor)
    set_if_provided(config, ["embedding", "n_neighbors"], args.n_neighbors)
    set_if_provided(config, ["embedding", "umap_spread"], args.umap_spread)
    set_if_provided(config, ["embedding", "umap_min_dist"], args.umap_min_dist)
    set_if_provided(
        config, ["embedding", "leiden_resolution"], args.leiden_resolution
    )
    set_if_provided(config, ["benchmark", "batch_key"], args.bench_batch_key)
    set_if_provided(config, ["benchmark", "label_key"], args.bench_label_key)

    if args.skip_qc_filter:
        config["qc"]["skip_filter"] = True
    if args.skip_benchmark:
        config["benchmark"]["enabled"] = False


def config_to_model_config(name, cfg, config):
    """Create a ModelConfig from a dict plus shared data settings."""
    hvg_cfg = cfg.get("hvg", {}) or {}
    return ModelConfig(
        name=name,
        n_hidden=cfg.get("n_hidden", config["integration"]["n_hidden"]),
        n_layers=cfg.get("n_layers", config["integration"]["n_layers"]),
        n_latent=cfg.get("n_latent", config["integration"]["n_latent"]),
        dispersion=cfg.get("dispersion", config["integration"]["dispersion"]),
        gene_likelihood=cfg.get(
            "gene_likelihood", config["integration"]["gene_likelihood"]
        ),
        batch_key=cfg.get("batch_key", config["data"]["batch_key"]),
        categorical_covariate_keys=cfg.get(
            "categorical_covariate_keys",
            config["data"].get("categorical_covariate_keys", []),
        )
        or [],
        continuous_covariate_keys=cfg.get(
            "continuous_covariate_keys",
            config["data"].get("continuous_covariate_keys", []),
        )
        or [],
        counts_layer=config["data"]["counts_layer"],
        hvg_batch_key=cfg.get(
            "hvg_batch_key",
            hvg_cfg.get("batch_key", config["integration"]["hvg"]["batch_key"]),
        ),
        hvg_flavor=cfg.get(
            "hvg_flavor",
            hvg_cfg.get("flavor", config["integration"]["hvg"]["flavor"]),
        ),
        hvg_nbatches=cfg.get(
            "hvg_nbatches",
            hvg_cfg.get("min_batches", config["integration"]["hvg"].get("min_batches")),
        ),
    )


def load_sweep_configs_from_args_or_config(args, config):
    """Return sweep configs from JSON, YAML config, or built-in defaults."""
    if args.sweep_configs:
        print(f"  Loading configs from {args.sweep_configs}")
        with open(args.sweep_configs) as f:
            raw = json.load(f)
    elif config["integration"].get("sweep"):
        raw = config["integration"]["sweep"]
        if isinstance(raw, list):
            raw = {entry["name"]: {k: v for k, v in entry.items() if k != "name"}
                   for entry in raw}
    else:
        print("  Using built-in sweep configs.")
        raw = SWEEP_CONFIGS
    return {name: config_to_model_config(name, params, config)
            for name, params in raw.items()}


def train_scanvi_from_scvi(model, adata_hvg, adata_full, config, output_dir):
    """Train scANVI from a fitted SCVI model and store annotation outputs."""
    annotation = config["annotation"]
    labels_key = annotation["labels_key"]
    unlabeled = annotation["unlabeled_category"]
    prediction_key = annotation["prediction_key"]
    confidence_key = annotation["confidence_key"]

    if labels_key not in adata_hvg.obs.columns:
        raise ValueError(
            f"annotation.labels_key '{labels_key}' is not present in adata.obs."
        )

    if adata_hvg.obs[labels_key].isna().any():
        print(
            f"  [WARN] NaN values in '{labels_key}' — filling with "
            f"'{unlabeled}' for scANVI."
        )
        if adata_hvg.obs[labels_key].dtype.name != "category":
            adata_hvg.obs[labels_key] = adata_hvg.obs[labels_key].astype("category")
        if unlabeled not in adata_hvg.obs[labels_key].cat.categories:
            adata_hvg.obs[labels_key] = (
                adata_hvg.obs[labels_key].cat.add_categories(unlabeled)
            )
        adata_hvg.obs[labels_key] = adata_hvg.obs[labels_key].fillna(unlabeled)

    print("\n  Training scANVI annotation model")
    print(f"  labels_key='{labels_key}', unlabeled_category='{unlabeled}'")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key=labels_key,
        unlabeled_category=unlabeled,
    )
    train_kwargs = dict(
        max_epochs=config["integration"]["max_epochs"],
        batch_size=config["integration"]["batch_size"],
        train_size=0.9,
        early_stopping=True,
        early_stopping_patience=config["integration"]["early_stopping_patience"],
        early_stopping_monitor="elbo_validation",
        check_val_every_n_epoch=1,
    )
    train_with_num_workers(
        scanvi_model,
        train_kwargs,
        num_workers=config["integration"].get("num_workers"),
    )

    latent = scanvi_model.get_latent_representation()
    adata_full.obsm["X_scANVI"] = latent

    predictions = scanvi_model.predict()
    adata_full.obs[prediction_key] = pd.Series(
        predictions, index=adata_hvg.obs_names
    ).reindex(adata_full.obs_names).astype("category")

    try:
        probabilities = scanvi_model.predict(soft=True)
        if hasattr(probabilities, "max"):
            confidence = probabilities.max(axis=1)
            if hasattr(confidence, "values"):
                confidence = confidence.values
        else:
            confidence = np.max(probabilities, axis=1)
        adata_full.obs[confidence_key] = pd.Series(
            confidence, index=adata_hvg.obs_names
        ).reindex(adata_full.obs_names).astype(float)
    except Exception as e:
        print(f"  [WARN] Could not compute scANVI prediction confidence: {e}")

    model_dir = os.path.join(output_dir, "scanvi_model_default")
    os.makedirs(model_dir, exist_ok=True)
    scanvi_model.save(model_dir, save_anndata=True, overwrite=True)
    print(f"  scANVI latent stored → obsm['X_scANVI']")
    print(f"  scANVI predictions stored → obs['{prediction_key}']")
    if confidence_key in adata_full.obs:
        print(f"  scANVI confidence stored → obs['{confidence_key}']")
    print(f"  scANVI model saved → {model_dir}")
    return scanvi_model, model_dir


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="SNS scVI integration with optional parameter sweep "
                    "and scIB benchmarking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default single-model run
  python integrate_scvi.py --input combined.h5ad

  # Parameter sweep
  python integrate_scvi.py --input combined.h5ad --sweep

  # Custom sweep configs from JSON
  python integrate_scvi.py --input combined.h5ad --sweep \\
      --sweep-configs my_configs.json

  # Enable bio-conservation metrics with a proxy label
  python integrate_scvi.py --input combined.h5ad --sweep \\
      --bench-label-key ganglion_group

  # Override default architecture
  python integrate_scvi.py --input combined.h5ad \\
      --n-hidden 128 --n-layers 2 --gene-likelihood zinb

Sweep config JSON format:
  {
    "config_name": {
      "n_hidden": 128,
      "n_layers": 2,
      "n_latent": 32,
      "dispersion": "gene-cell",
      "gene_likelihood": "nb",
      "batch_key": "platform_origin",
      "categorical_covariate_keys": ["data_origin"],
      "hvg_batch_key": "platform_origin",
      "hvg_flavor": "seurat",
      "hvg_nbatches": 3
    }
  }
""",
    )

    # ── I/O ──
    io = parser.add_argument_group("input / output")
    io.add_argument("--config", default=None,
                    help="Versioned pipeline YAML config.")
    io.add_argument("--input", required=False, default=None,
                    help="Path to pre-concatenated h5ad.")
    io.add_argument("--output-dir", default=None,
                    help="Output directory (relative to the current directory; "
                         "default: ./results without --config).")

    # ── QC ──
    qc = parser.add_argument_group("quality control")
    qc.add_argument("--min-genes", type=int, default=None,
                    help="Min genes per cell (default: 500).")
    qc.add_argument("--min-cells", type=int, default=None,
                    help="Min cells per gene (default: 3).")
    qc.add_argument("--skip-qc-filter", action="store_true",
                    help="Skip cell/gene filtering (data already filtered).")

    # ── HVG ──
    hvg = parser.add_argument_group("HVG selection")
    hvg.add_argument("--n-hvgs", type=int, default=None,
                     help="Number of HVGs (default: 3000).")
    hvg.add_argument("--hvg-batch-key", default=None,
                     help="Batch key for HVG selection (default: tech).")
    hvg.add_argument("--hvg-flavor", default=None,
                     choices=["seurat_v3", "seurat", "cell_ranger"],
                     help="HVG flavor (default: seurat_v3).")

    # ── Architecture (single-model mode) ──
    arch = parser.add_argument_group("model architecture (single-model mode)")
    arch.add_argument("--batch-key", default=None,
                      help="Batch key for scVI (default: data_origin).")
    arch.add_argument("--covariate-keys", nargs="*", default=None,
                      help="Categorical covariate keys (default: tech). "
                           "Pass without values to disable.")
    arch.add_argument("--n-hidden", type=int, default=None)
    arch.add_argument("--n-layers", type=int, default=None)
    arch.add_argument("--n-latent", type=int, default=None)
    arch.add_argument("--dispersion", default=None,
                      choices=["gene", "gene-batch", "gene-label", "gene-cell"])
    arch.add_argument("--gene-likelihood", default=None,
                      choices=["zinb", "nb", "poisson"])

    # ── Training ──
    tr = parser.add_argument_group("training")
    tr.add_argument("--max-epochs", type=int, default=None)
    tr.add_argument("--early-stopping-patience", type=int, default=None)
    tr.add_argument("--batch-size", type=int, default=None)
    tr.add_argument("--num-workers", type=int, default=None,
                    help="DataLoader workers for scVI/scANVI training.")

    # ── Sweep ──
    sw = parser.add_argument_group("parameter sweep")
    sw.add_argument("--sweep", action="store_true",
                    help="Run parameter sweep over multiple configs.")
    sw.add_argument("--sweep-configs", type=str, default=None,
                    help="JSON file with custom sweep configs. "
                         "Omit to use built-in defaults.")

    # ── UMAP / clustering ──
    um = parser.add_argument_group("UMAP / clustering")
    um.add_argument("--n-neighbors", type=int, default=None)
    um.add_argument("--umap-spread", type=float, default=None)
    um.add_argument("--umap-min-dist", type=float, default=None)
    um.add_argument("--leiden-resolution", type=float, default=None)

    # ── Benchmarking ──
    bm = parser.add_argument_group("scIB benchmarking")
    bm.add_argument("--bench-batch-key", default=None,
                    help="Batch key for scIB (default: same as --batch-key).")
    bm.add_argument("--bench-label-key", default=None,
                    help="Label key for bio-conservation metrics (optional). "
                         "Omit to compute batch-correction metrics only.")
    bm.add_argument("--skip-benchmark", action="store_true",
                    help="Skip scIB benchmarking entirely.")

    args = parser.parse_args()
    try:
        config = load_pipeline_config(args.config)
        apply_integration_cli_overrides(config, args)
        if args.sweep and config["annotation"].get("enabled"):
            raise ConfigError("scANVI annotation is only supported for single-model runs in V1.")
        validate_config(config, mode="integration")
    except ConfigError as e:
        parser.error(str(e))

    input_path = config["data"]["input_h5ad"]
    output_dir = config["data"]["output_dir"]
    counts_layer = config["data"]["counts_layer"]
    data_cfg = config["data"]
    qc_cfg = config["qc"]

    # Round context: the input of a round >1 sits beside its parent round's
    # round_manifest.json, which read_parent_provenance detects. This is the
    # signal used to honor qc.skip_filter_after_round_1 (skip gene/cell filtering
    # on later rounds) without a manual per-round flag. Computed once and reused
    # for the manifest's parent_round pointer.
    parent_provenance = read_parent_provenance(input_path)
    is_first_round = parent_provenance is None
    integration_cfg = config["integration"]
    embedding_cfg = config["embedding"]
    benchmark_cfg = config["benchmark"]
    annotation_cfg = config["annotation"]

    # Seed all RNGs before any stochastic step (scVI/scANVI training, neighbors,
    # UMAP, Leiden) so the integrated object — and therefore the cluster IDs that
    # downstream decisions.yaml files reference — is reproducible.
    applied_seed = set_global_seed(config.get("reproducibility", {}).get("seed"))
    if applied_seed is not None:
        print(f"  Global seed: {applied_seed}")
    else:
        print("  Global seed: unset (non-reproducible run)")

    os.makedirs(output_dir, exist_ok=True)
    write_resolved_config(output_dir, config)
    write_command_args(output_dir, "integrate", vars(args))

    # ══════════════════════════════════════════════════════════════════════
    #  1. Load data
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("1. LOADING DATA")
    print(f"{'='*60}")
    adata = ad.read_h5ad(input_path)
    print(f"  Shape : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  obs   : {list(adata.obs.columns)}")

    # Cell counts at each attrition step, so round_manifest.json reconciles with
    # the QC plots (which are drawn post-obs_filter, pre-QC-filter).
    n_cells_input = adata.n_obs
    n_genes_input = adata.n_vars

    obs_filter = data_cfg.get("obs_filter")
    n_cells_after_obs_filter = None
    if obs_filter:
        n_before = adata.n_obs
        try:
            mask = adata.obs.eval(obs_filter)
        except Exception as e:
            raise ValueError(
                f"data.obs_filter query failed: {obs_filter!r}\n  {e}"
            ) from e
        adata = adata[mask].copy()
        n_cells_after_obs_filter = adata.n_obs
        print(f"  obs_filter: {obs_filter!r}")
        print(f"  Cells after filter: {adata.n_obs:,} (removed {n_before - adata.n_obs:,})")

    ensure_counts_layer(adata, counts_layer=counts_layer)

    # ══════════════════════════════════════════════════════════════════════
    #  2. QC metrics and optional filtering
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("2. QUALITY CONTROL")
    print(f"{'='*60}")
    compute_qc_metrics(
        adata,
        mt_gene_patterns=qc_cfg.get("mt_gene_patterns"),
        ribo_gene_patterns=qc_cfg.get("ribo_gene_patterns"),
        hb_gene_pattern=qc_cfg.get("hb_gene_pattern"),
        gene_symbol_case=data_cfg.get("gene_symbol_case", "preserve"),
    )
    plot_qc_violins(adata, output_dir, groupby=data_cfg["batch_key"],
                    suffix="prefilter")
    n_cells_prefilter_plot = adata.n_obs

    # Honor skip_filter_after_round_1: re-filtering genes against a later
    # round's smaller cell subset can drop rare markers the first (full-data)
    # round correctly kept.
    qc_filter_applied, skipped_later_round = should_apply_qc_filter(
        qc_cfg, is_first_round)
    if qc_filter_applied:
        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=qc_cfg["min_genes"])
        sc.pp.filter_genes(adata, min_cells=qc_cfg["min_cells"])
        print(f"  Filtered: {n_before:,} -> {adata.n_obs:,} cells")
    elif skipped_later_round:
        print("  Skipping QC filtering (round >1; skip_filter_after_round_1). "
              "Genes/cells are not re-filtered against this round's subset.")
        # Re-save counts after filtering
        adata.layers[counts_layer] = adata.X.copy()
        plot_qc_violins(adata, output_dir, groupby=data_cfg["batch_key"],
                        suffix="postfilter")
    else:
        print("  Skipping QC filtering.")

    # ══════════════════════════════════════════════════════════════════════
    #  3. Preserve full-gene object
    # ══════════════════════════════════════════════════════════════════════
    adata_full = adata.copy()
    print(f"\n  Full object: {adata_full.n_obs:,} cells x "
          f"{adata_full.n_vars:,} genes (preserved for output)")

    # ══════════════════════════════════════════════════════════════════════
    #  4. Integration
    # ══════════════════════════════════════════════════════════════════════
    models = {}
    embedding_keys = []
    embedding_key_by_name = {}
    model_dirs = {}
    configs = {}
    sweep_mode = args.sweep or bool(integration_cfg.get("sweep"))

    if sweep_mode:
        # ── SWEEP MODE ───────────────────────────────────────────────
        print(f"\n{'='*60}")
        print("4. INTEGRATION — PARAMETER SWEEP")
        print(f"{'='*60}")

        configs = load_sweep_configs_from_args_or_config(args, config)

        for name, cfg in configs.items():
            print(f"\n{'─'*50}")
            print(f"  Config: {name}")
            print(f"{'─'*50}")

            hvgs = select_hvgs(
                adata, n_top_genes=integration_cfg["hvg"]["n_top_genes"],
                batch_key=cfg.hvg_batch_key,
                flavor=cfg.hvg_flavor,
                hvg_nbatches=cfg.hvg_nbatches,
                counts_layer=counts_layer,
            )
            print(f"  HVGs: {len(hvgs)}")

            adata_hvg = adata[:, hvgs].copy()
            adata_hvg.layers[counts_layer] = adata_hvg.X.copy()

            model = setup_and_train(
                adata_hvg, cfg,
                max_epochs=integration_cfg["max_epochs"],
                early_stopping_patience=integration_cfg["early_stopping_patience"],
                batch_size=integration_cfg["batch_size"],
                num_workers=integration_cfg.get("num_workers"),
                output_dir=output_dir,
            )

            latent = model.get_latent_representation()
            latent_key = f"X_scVI_{name}"
            adata_full.obsm[latent_key] = latent
            embedding_keys.append(latent_key)
            embedding_key_by_name[name] = latent_key
            models[name] = model
            model_dirs[name] = os.path.join(output_dir, f"scvi_model_{name}")
            print(f"  Latent stored → obsm['{latent_key}']")

            del adata_hvg; gc.collect()

    else:
        # ── SINGLE MODEL MODE ────────────────────────────────────────
        print(f"\n{'='*60}")
        print("4. INTEGRATION — SINGLE MODEL")
        print(f"{'='*60}")

        cfg = ModelConfig(
            name="default",
            n_hidden=integration_cfg["n_hidden"],
            n_layers=integration_cfg["n_layers"],
            n_latent=integration_cfg["n_latent"],
            dispersion=integration_cfg["dispersion"],
            gene_likelihood=integration_cfg["gene_likelihood"],
            batch_key=data_cfg["batch_key"],
            categorical_covariate_keys=data_cfg.get("categorical_covariate_keys", []) or [],
            continuous_covariate_keys=data_cfg.get("continuous_covariate_keys", []) or [],
            counts_layer=counts_layer,
            hvg_batch_key=integration_cfg["hvg"]["batch_key"],
            hvg_flavor=integration_cfg["hvg"]["flavor"],
            hvg_nbatches=integration_cfg["hvg"].get("min_batches"),
        )

        hvgs = select_hvgs(
            adata, n_top_genes=integration_cfg["hvg"]["n_top_genes"],
            batch_key=cfg.hvg_batch_key, flavor=cfg.hvg_flavor,
            hvg_nbatches=cfg.hvg_nbatches,
            counts_layer=counts_layer,
        )
        print(f"  HVGs: {len(hvgs)}")

        pd.DataFrame({"gene": hvgs}).to_csv(
            os.path.join(output_dir, "hvg_genes.csv"), index=False)
        adata_full.var["highly_variable"] = adata_full.var_names.isin(hvgs)

        adata_hvg = adata[:, hvgs].copy()
        adata_hvg.layers[counts_layer] = adata_hvg.X.copy()

        model = setup_and_train(
            adata_hvg, cfg,
            max_epochs=integration_cfg["max_epochs"],
            early_stopping_patience=integration_cfg["early_stopping_patience"],
            batch_size=integration_cfg["batch_size"],
            num_workers=integration_cfg.get("num_workers"),
            output_dir=output_dir,
        )

        latent = model.get_latent_representation()
        adata_full.obsm["X_scVI"] = latent
        embedding_keys.append("X_scVI")
        embedding_key_by_name["default"] = "X_scVI"
        models["default"] = model
        model_dirs["default"] = os.path.join(output_dir, "scvi_model_default")

        scanvi_model_dir = None
        if annotation_cfg.get("enabled"):
            scanvi_model, scanvi_model_dir = train_scanvi_from_scvi(
                model, adata_hvg, adata_full, config, output_dir
            )
            embedding_keys.append("X_scANVI")
            embedding_key_by_name["scanvi"] = "X_scANVI"
            model_dirs["scanvi"] = scanvi_model_dir

        del adata_hvg; gc.collect()

    # ══════════════════════════════════════════════════════════════════════
    #  4b. Optional Harmony integration (both single-model and sweep modes)
    # ══════════════════════════════════════════════════════════════════════
    # Harmony runs once per round on its own HVG subset, independent of the
    # scVI model(s). In sweep mode there is no single canonical HVG set, so the
    # source is an EXPLICIT, recorded decision: by default Harmony uses the
    # top-level integration.hvg spec (identical to what single-model scVI uses);
    # set harmony.hvg_from to a sweep entry's name to borrow that architecture's
    # HVG selection instead. Either way the resolved spec is logged to the
    # manifest so the choice is visible, not implicit.
    harmony_cfg = integration_cfg.get("harmony", {}) or {}
    harmony_provenance = None
    if harmony_cfg.get("enabled"):
        harmony_batch_key = harmony_cfg.get("batch_key") or data_cfg["batch_key"]
        n_pcs = harmony_cfg.get("n_pcs", 30)
        hvg_spec = resolve_harmony_hvg_spec(
            harmony_cfg, integration_cfg,
            configs if sweep_mode else {}, sweep_mode,
        )

        print(f"\n  Running Harmony integration (batch_key='{harmony_batch_key}', "
              f"n_pcs={n_pcs}, HVG source={hvg_spec['source']})...")
        harmony_hvgs = select_hvgs(
            adata, n_top_genes=hvg_spec["n_top_genes"],
            batch_key=hvg_spec["batch_key"], flavor=hvg_spec["flavor"],
            hvg_nbatches=hvg_spec["min_batches"], counts_layer=counts_layer,
        )
        print(f"  Harmony HVGs: {len(harmony_hvgs)}")
        adata_harmony = adata[:, harmony_hvgs].copy()
        adata_harmony.layers[counts_layer] = adata_harmony.X.copy()
        harmony_emb = run_harmony(
            adata_harmony,
            batch_key=harmony_batch_key,
            counts_layer=counts_layer,
            n_pcs=n_pcs,
            seed=applied_seed,
        )
        adata_full.obsm["X_pca_harmony"] = harmony_emb
        embedding_keys.append("X_pca_harmony")
        embedding_key_by_name["harmony"] = "X_pca_harmony"
        harmony_provenance = {
            "batch_key": harmony_batch_key,
            "n_pcs": int(harmony_emb.shape[1]),
            "n_hvgs": len(harmony_hvgs),
            "hvg": hvg_spec,
        }
        print(f"  Harmony embedding stored → obsm['X_pca_harmony'] "
              f"({harmony_emb.shape[1]} dims)")
        del adata_harmony; gc.collect()

    del adata; gc.collect()

    # ══════════════════════════════════════════════════════════════════════
    #  5. Training diagnostics
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("5. TRAINING DIAGNOSTICS")
    print(f"{'='*60}")
    plot_training_curves(models, output_dir)

    # ══════════════════════════════════════════════════════════════════════
    #  6. UMAP + Leiden
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("6. UMAP + LEIDEN CLUSTERING")
    print(f"{'='*60}")
    umap_keys = {}
    leiden_keys = {}

    for name, latent_key in embedding_key_by_name.items():
        if sweep_mode:
            umap_key = f"X_umap_{name}"
            leiden_key = f"leiden_{name}"
            nbr_key = f"neighbors_{name}"
        elif name == "default":
            umap_key = "X_umap"
            leiden_key = "leiden"
            nbr_key = "neighbors"
        else:
            umap_key = f"X_umap_{name}"
            leiden_key = f"leiden_{name}"
            nbr_key = f"neighbors_{name}"

        n_latent = adata_full.obsm[latent_key].shape[1]
        print(f"  {name}: neighbors (n={embedding_cfg['n_neighbors']}) -> UMAP -> "
              f"leiden (res={embedding_cfg['leiden_resolution']})")

        run_neighbors_with_key(
            adata_full,
            nbr_key,
            use_rep=latent_key,
            n_neighbors=embedding_cfg["n_neighbors"],
            n_pcs=n_latent,
            seed=applied_seed,
        )
        run_umap_with_key(
            adata_full,
            neighbors_key=nbr_key,
            umap_key=umap_key,
            spread=embedding_cfg["umap_spread"],
            min_dist=embedding_cfg["umap_min_dist"],
            seed=applied_seed,
        )
        run_leiden_with_key(
            adata_full,
            neighbors_key=nbr_key,
            leiden_key=leiden_key,
            resolution=embedding_cfg["leiden_resolution"],
            flavor="igraph",
            seed=applied_seed,
        )

        umap_keys[name] = umap_key
        leiden_keys[name] = leiden_key

    # ══════════════════════════════════════════════════════════════════════
    #  7. UMAP plots
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("7. UMAP VISUALISATION")
    print(f"{'='*60}")
    if sweep_mode:
        color_candidates = ["data_origin", "tech", "platform_origin",
                            "ganglion_group", "leiden"]
        color_vars = [c for c in color_candidates
                      if c in adata_full.obs.columns or c == "leiden"]
        plot_umap_grid(adata_full, umap_keys, leiden_keys, color_vars,
                       output_dir)
    else:
        plot_single_umaps(adata_full, output_dir)
        if "X_umap_scanvi" in adata_full.obsm:
            plot_single_umaps(
                adata_full,
                output_dir,
                basis="X_umap_scanvi",
                output_name="umap_scanvi_annotation.png",
            )

    # ══════════════════════════════════════════════════════════════════════
    #  8. scIB benchmarking
    # ══════════════════════════════════════════════════════════════════════
    if benchmark_cfg.get("enabled", True):
        print(f"\n{'='*60}")
        print("8. scIB BENCHMARKING")
        print(f"{'='*60}")
        bench_batch_key = benchmark_cfg.get("batch_key")
        if bench_batch_key == "auto":
            bench_batch_key = data_cfg["batch_key"]
        run_benchmarking(
            adata_full,
            embedding_keys=embedding_keys,
            bench_batch_key=bench_batch_key,
            bench_label_key=benchmark_cfg.get("label_key"),
            output_dir=output_dir,
            counts_layer=counts_layer,
            seed=applied_seed,
        )
    else:
        print("\n  Skipping scIB benchmarking.")

    # ══════════════════════════════════════════════════════════════════════
    #  9. Save integrated object — ALL genes
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("9. SAVING INTEGRATED OBJECT")
    print(f"{'='*60}")

    if "cell_id" not in adata_full.obs.columns:
        adata_full.obs["cell_id"] = adata_full.obs_names

    # Normalise .X for downstream use; raw counts stay in layers
    adata_full.X = adata_full.layers[counts_layer].copy()
    sc.pp.normalize_total(adata_full)
    sc.pp.log1p(adata_full)

    out_path = os.path.join(output_dir, "integrated.h5ad")
    adata_full.write_h5ad(out_path)
    print(f"  Saved {out_path}")
    print(f"  Shape  : {adata_full.n_obs:,} cells x {adata_full.n_vars:,} genes")
    print(f"  .X     : normalised, log1p")
    print(f"  layers : {list(adata_full.layers.keys())}")
    print(f"  obsm   : {list(adata_full.obsm.keys())}")

    update_round_manifest(
        output_dir,
        "integration",
        {
            "input_h5ad": input_path,
            "input_fingerprint": fingerprint_file(input_path),
            "obs_filter": obs_filter or None,
            "parent_round": parent_provenance,
            "is_first_round": is_first_round,
            "output_h5ad": out_path,
            # Attrition chain. n_cells stays the final post-QC count so existing
            # readers (summarize_rounds) keep working; the rest explain the gap
            # between the qc_*_prefilter plots and that number.
            "n_cells_input": n_cells_input,
            "n_genes_input": n_genes_input,
            "n_cells_after_obs_filter": n_cells_after_obs_filter,
            "n_cells_prefilter_plot": n_cells_prefilter_plot,
            "qc_filter_applied": qc_filter_applied,
            "n_cells": adata_full.n_obs,
            "n_genes": adata_full.n_vars,
            "model_type": integration_cfg["model_type"],
            "annotation_enabled": annotation_cfg.get("enabled", False),
            "embedding_keys": embedding_keys,
            "umap_keys": umap_keys,
            "cluster_keys": leiden_keys,
            "model_dirs": model_dirs,
            "harmony": harmony_provenance,
            "counts_layer": counts_layer,
            "package_versions": collect_package_versions(
                [
                    "anndata",
                    "scanpy",
                    "scvi-tools",
                    "scib-metrics",
                    "torch",
                    "numpy",
                    "scikit-learn",
                    "leidenalg",
                    "umap-learn",
                    "harmonypy",
                ]
            ),
        },
        seed=applied_seed,
    )

    # Summary
    print(f"\n{'='*60}")
    print(f"  Done. All outputs in {output_dir}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

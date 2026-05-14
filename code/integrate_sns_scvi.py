#!/usr/bin/env python3
"""SNS scVI integration with optional parameter sweep and scIB benchmarking.

Merges three workflows into a single CLI script:
  1. Default single-model integration (from 01_20260513 combined SNS notebook)
  2. Parameter sweep over multiple architectures (from 20260325 refinement notebook)
  3. Annotation-free scIB benchmarking (scib_metrics)

The input is a pre-concatenated h5ad with raw counts (in .X or .layers["counts"])
and at minimum a batch column (e.g. "data_origin") in .obs.

Usage
-----
    # Default single-model integration (large architecture)
    python integrate_sns_scvi.py --input combined.h5ad

    # Override default architecture
    python integrate_sns_scvi.py --input combined.h5ad \
        --n-hidden 128 --n-layers 2 --gene-likelihood zinb

    # Parameter sweep with built-in configs
    python integrate_sns_scvi.py --input combined.h5ad --sweep

    # Parameter sweep with custom configs (JSON)
    python integrate_sns_scvi.py --input combined.h5ad --sweep \
        --sweep-configs my_configs.json

    # Enable bio-conservation metrics using a proxy label
    python integrate_sns_scvi.py --input combined.h5ad --sweep \
        --bench-label-key ganglion_group
"""

import argparse
import gc
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

def compute_qc_metrics(adata):
    """Annotate mt/ribo/hb genes and compute QC metrics in-place."""
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^hb[^(p)]", regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True,
    )


def plot_qc_violins(adata, output_dir, groupby="data_origin", suffix=""):
    """QC violin plots grouped by a batch variable."""
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

def ensure_counts_layer(adata):
    """Ensure raw integer counts live in adata.layers['counts'] and .X."""
    from scipy.sparse import issparse

    def _looks_like_counts(mat, n_check=200):
        sample = mat[:n_check]
        if issparse(sample):
            sample = sample.toarray()
        vals = sample[sample > 0]
        return len(vals) > 0 and np.allclose(vals, np.round(vals))

    if "counts" in adata.layers and _looks_like_counts(adata.layers["counts"]):
        print("  Using existing 'counts' layer (integer-valued).")
        adata.X = adata.layers["counts"].copy()
        return

    if _looks_like_counts(adata.X):
        print("  .X contains raw counts — copying to layers['counts'].")
        adata.layers["counts"] = adata.X.copy()
        return

    # Fallback: assume X is usable even if not perfectly integer
    print("  [WARN] .X does not look like raw integer counts.")
    if "counts" not in adata.layers:
        print("         No 'counts' layer found. Using .X as-is — verify input data.")
        adata.layers["counts"] = adata.X.copy()
    else:
        print("         Using existing 'counts' layer despite non-integer values.")
        adata.X = adata.layers["counts"].copy()


def select_hvgs(adata, n_top_genes=3000, batch_key="tech", flavor="seurat_v3",
                hvg_nbatches=None):
    """Select HVGs on a temporary normalised copy; return gene name list.

    The caller's adata is not modified.
    """
    tmp = adata.copy()
    tmp.X = tmp.layers["counts"].copy()

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


def setup_and_train(adata_hvg, config, max_epochs=200,
                    early_stopping_patience=20, batch_size=256,
                    output_dir=None):
    """Setup AnnData, build, and train one scVI model. Returns the model."""
    clean_batch_column(adata_hvg, config.batch_key)

    # Filter covariate keys to columns actually present
    cat_covs = [c for c in config.categorical_covariate_keys
                if c in adata_hvg.obs.columns]
    for c in cat_covs:
        clean_batch_column(adata_hvg, c)

    scvi.model.SCVI.setup_anndata(
        adata_hvg, layer="counts",
        batch_key=config.batch_key,
        categorical_covariate_keys=cat_covs or None,
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
    print(f"  batch_key='{config.batch_key}', covariates={cat_covs or 'none'}")

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=True,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor="elbo_validation",
        check_val_every_n_epoch=1,
        train_size=0.9,
    )

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


def plot_single_umaps(adata, output_dir, basis="X_umap"):
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
    fig.savefig(os.path.join(output_dir, "umap_integration.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved umap_integration.png")


# ═══════════════════════════════════════════════════════════════════════════════
#  scIB benchmarking
# ═══════════════════════════════════════════════════════════════════════════════

def run_benchmarking(adata, embedding_keys, bench_batch_key,
                     bench_label_key=None, output_dir=None, n_pca_comps=50):
    """Run scIB benchmarking across embedding keys.

    Parameters
    ----------
    bench_label_key : str or None
        If None, only batch-correction metrics are computed and a leiden proxy
        is generated internally so the Benchmarker constructor is satisfied.
        If provided, a subset of bio-conservation metrics is also enabled.
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
        adata_tmp.X = adata_tmp.layers["counts"].copy()
        sc.pp.normalize_total(adata_tmp)
        sc.pp.log1p(adata_tmp)
        sc.pp.pca(adata_tmp, n_comps=n_pca_comps)
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
            sc.pp.neighbors(adata, use_rep="X_pca", key_added="_proxy_neighbors")
            sc.tl.leiden(adata, neighbors_key="_proxy_neighbors", resolution=0.5,
                         flavor="igraph", key_added="leiden_proxy")
        label_key = "leiden_proxy"
        use_bio = False
        print(f"  No annotation label — using '{label_key}' as proxy "
              f"(bio-conservation metrics disabled).")

    # Metric configuration
    bio = BioConservation(
        isolated_labels=False,
        silhouette_label=use_bio,
        clisi_knn=use_bio,
        nmi_ari_cluster_labels_kmeans=False,
        nmi_ari_cluster_labels_leiden=False,
    )
    batch = BatchCorrection(
        silhouette_batch=True,
        graph_connectivity=True,
        pcr_comparison=True,
        ilisi_knn=False,
        kbet_per_label=False,
    )

    bm = Benchmarker(
        adata,
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
  python integrate_sns_scvi.py --input combined.h5ad

  # Parameter sweep
  python integrate_sns_scvi.py --input combined.h5ad --sweep

  # Custom sweep configs from JSON
  python integrate_sns_scvi.py --input combined.h5ad --sweep \\
      --sweep-configs my_configs.json

  # Enable bio-conservation metrics with a proxy label
  python integrate_sns_scvi.py --input combined.h5ad --sweep \\
      --bench-label-key ganglion_group

  # Override default architecture
  python integrate_sns_scvi.py --input combined.h5ad \\
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
    io.add_argument("--input", required=True,
                    help="Path to pre-concatenated h5ad.")
    io.add_argument("--output-dir", default="../results",
                    help="Output directory (default: ../results).")

    # ── QC ──
    qc = parser.add_argument_group("quality control")
    qc.add_argument("--min-genes", type=int, default=500,
                    help="Min genes per cell (default: 500).")
    qc.add_argument("--min-cells", type=int, default=3,
                    help="Min cells per gene (default: 3).")
    qc.add_argument("--skip-qc-filter", action="store_true",
                    help="Skip cell/gene filtering (data already filtered).")

    # ── HVG ──
    hvg = parser.add_argument_group("HVG selection")
    hvg.add_argument("--n-hvgs", type=int, default=3000,
                     help="Number of HVGs (default: 3000).")
    hvg.add_argument("--hvg-batch-key", default="tech",
                     help="Batch key for HVG selection (default: tech).")
    hvg.add_argument("--hvg-flavor", default="seurat_v3",
                     choices=["seurat_v3", "seurat", "cell_ranger"],
                     help="HVG flavor (default: seurat_v3).")

    # ── Architecture (single-model mode) ──
    arch = parser.add_argument_group("model architecture (single-model mode)")
    arch.add_argument("--batch-key", default="data_origin",
                      help="Batch key for scVI (default: data_origin).")
    arch.add_argument("--covariate-keys", nargs="*", default=["tech"],
                      help="Categorical covariate keys (default: tech). "
                           "Pass without values to disable.")
    arch.add_argument("--n-hidden", type=int, default=256)
    arch.add_argument("--n-layers", type=int, default=3)
    arch.add_argument("--n-latent", type=int, default=32)
    arch.add_argument("--dispersion", default="gene-cell",
                      choices=["gene", "gene-batch", "gene-label", "gene-cell"])
    arch.add_argument("--gene-likelihood", default="nb",
                      choices=["zinb", "nb", "poisson"])

    # ── Training ──
    tr = parser.add_argument_group("training")
    tr.add_argument("--max-epochs", type=int, default=200)
    tr.add_argument("--early-stopping-patience", type=int, default=20)
    tr.add_argument("--batch-size", type=int, default=256)

    # ── Sweep ──
    sw = parser.add_argument_group("parameter sweep")
    sw.add_argument("--sweep", action="store_true",
                    help="Run parameter sweep over multiple configs.")
    sw.add_argument("--sweep-configs", type=str, default=None,
                    help="JSON file with custom sweep configs. "
                         "Omit to use built-in defaults.")

    # ── UMAP / clustering ──
    um = parser.add_argument_group("UMAP / clustering")
    um.add_argument("--n-neighbors", type=int, default=30)
    um.add_argument("--umap-spread", type=float, default=3.0)
    um.add_argument("--umap-min-dist", type=float, default=0.4)
    um.add_argument("--leiden-resolution", type=float, default=0.3)

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
    os.makedirs(args.output_dir, exist_ok=True)
    if args.bench_batch_key is None:
        args.bench_batch_key = args.batch_key

    # ══════════════════════════════════════════════════════════════════════
    #  1. Load data
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("1. LOADING DATA")
    print(f"{'='*60}")
    adata = ad.read_h5ad(args.input)
    print(f"  Shape : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  obs   : {list(adata.obs.columns)}")

    ensure_counts_layer(adata)

    # ══════════════════════════════════════════════════════════════════════
    #  2. QC metrics and optional filtering
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("2. QUALITY CONTROL")
    print(f"{'='*60}")
    compute_qc_metrics(adata)
    plot_qc_violins(adata, args.output_dir, groupby=args.batch_key,
                    suffix="prefilter")

    if not args.skip_qc_filter:
        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=args.min_genes)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
        print(f"  Filtered: {n_before:,} -> {adata.n_obs:,} cells")
        # Re-save counts after filtering
        adata.layers["counts"] = adata.X.copy()
        plot_qc_violins(adata, args.output_dir, groupby=args.batch_key,
                        suffix="postfilter")
    else:
        print("  Skipping QC filtering (--skip-qc-filter).")

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

    if args.sweep:
        # ── SWEEP MODE ───────────────────────────────────────────────
        print(f"\n{'='*60}")
        print("4. INTEGRATION — PARAMETER SWEEP")
        print(f"{'='*60}")

        if args.sweep_configs:
            print(f"  Loading configs from {args.sweep_configs}")
            with open(args.sweep_configs) as f:
                raw = json.load(f)
            configs = {name: ModelConfig(name=name, **params)
                       for name, params in raw.items()}
        else:
            print("  Using built-in sweep configs.")
            configs = {name: ModelConfig(name=name, **params)
                       for name, params in SWEEP_CONFIGS.items()}

        for name, cfg in configs.items():
            print(f"\n{'─'*50}")
            print(f"  Config: {name}")
            print(f"{'─'*50}")

            hvgs = select_hvgs(
                adata, n_top_genes=args.n_hvgs,
                batch_key=cfg.hvg_batch_key,
                flavor=cfg.hvg_flavor,
                hvg_nbatches=cfg.hvg_nbatches,
            )
            print(f"  HVGs: {len(hvgs)}")

            adata_hvg = adata[:, hvgs].copy()
            adata_hvg.layers["counts"] = adata_hvg.X.copy()

            model = setup_and_train(
                adata_hvg, cfg,
                max_epochs=args.max_epochs,
                early_stopping_patience=args.early_stopping_patience,
                batch_size=args.batch_size,
                output_dir=args.output_dir,
            )

            latent = model.get_latent_representation()
            latent_key = f"X_scVI_{name}"
            adata_full.obsm[latent_key] = latent
            embedding_keys.append(latent_key)
            models[name] = model
            print(f"  Latent stored → obsm['{latent_key}']")

            del adata_hvg; gc.collect()

    else:
        # ── SINGLE MODEL MODE ────────────────────────────────────────
        print(f"\n{'='*60}")
        print("4. INTEGRATION — SINGLE MODEL")
        print(f"{'='*60}")

        cfg = ModelConfig(
            name="default",
            n_hidden=args.n_hidden, n_layers=args.n_layers,
            n_latent=args.n_latent, dispersion=args.dispersion,
            gene_likelihood=args.gene_likelihood,
            batch_key=args.batch_key,
            categorical_covariate_keys=args.covariate_keys or [],
            hvg_batch_key=args.hvg_batch_key,
            hvg_flavor=args.hvg_flavor,
        )

        hvgs = select_hvgs(
            adata, n_top_genes=args.n_hvgs,
            batch_key=cfg.hvg_batch_key, flavor=cfg.hvg_flavor,
        )
        print(f"  HVGs: {len(hvgs)}")

        pd.DataFrame({"gene": hvgs}).to_csv(
            os.path.join(args.output_dir, "hvg_genes.csv"), index=False)
        adata_full.var["highly_variable"] = adata_full.var_names.isin(hvgs)

        adata_hvg = adata[:, hvgs].copy()
        adata_hvg.layers["counts"] = adata_hvg.X.copy()

        model = setup_and_train(
            adata_hvg, cfg,
            max_epochs=args.max_epochs,
            early_stopping_patience=args.early_stopping_patience,
            batch_size=args.batch_size,
            output_dir=args.output_dir,
        )

        latent = model.get_latent_representation()
        adata_full.obsm["X_scVI"] = latent
        embedding_keys.append("X_scVI")
        models["default"] = model

        del adata_hvg; gc.collect()

    del adata; gc.collect()

    # ══════════════════════════════════════════════════════════════════════
    #  5. Training diagnostics
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("5. TRAINING DIAGNOSTICS")
    print(f"{'='*60}")
    plot_training_curves(models, args.output_dir)

    # ══════════════════════════════════════════════════════════════════════
    #  6. UMAP + Leiden
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("6. UMAP + LEIDEN CLUSTERING")
    print(f"{'='*60}")
    umap_keys = {}
    leiden_keys = {}

    for name in models:
        latent_key = f"X_scVI_{name}" if args.sweep else "X_scVI"
        umap_key = f"X_umap_{name}" if args.sweep else "X_umap"
        leiden_key = f"leiden_{name}" if args.sweep else "leiden"
        nbr_key = f"neighbors_{name}" if args.sweep else "neighbors"

        n_latent = adata_full.obsm[latent_key].shape[1]
        print(f"  {name}: neighbors (n={args.n_neighbors}) -> UMAP -> "
              f"leiden (res={args.leiden_resolution})")

        sc.pp.neighbors(adata_full, use_rep=latent_key,
                        n_neighbors=args.n_neighbors,
                        n_pcs=n_latent,
                        key_added=nbr_key)
        sc.tl.umap(adata_full,
                    neighbors_key=nbr_key,
                    key_added=umap_key,
                    spread=args.umap_spread,
                    min_dist=args.umap_min_dist)
        sc.tl.leiden(adata_full,
                     neighbors_key=nbr_key,
                     resolution=args.leiden_resolution,
                     flavor="igraph",
                     key_added=leiden_key)

        umap_keys[name] = umap_key
        leiden_keys[name] = leiden_key

    # ══════════════════════════════════════════════════════════════════════
    #  7. UMAP plots
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("7. UMAP VISUALISATION")
    print(f"{'='*60}")
    if args.sweep:
        color_candidates = ["data_origin", "tech", "platform_origin",
                            "ganglion_group", "leiden"]
        color_vars = [c for c in color_candidates
                      if c in adata_full.obs.columns or c == "leiden"]
        plot_umap_grid(adata_full, umap_keys, leiden_keys, color_vars,
                       args.output_dir)
    else:
        plot_single_umaps(adata_full, args.output_dir)

    # ══════════════════════════════════════════════════════════════════════
    #  8. scIB benchmarking
    # ══════════════════════════════════════════════════════════════════════
    if not args.skip_benchmark:
        print(f"\n{'='*60}")
        print("8. scIB BENCHMARKING")
        print(f"{'='*60}")
        run_benchmarking(
            adata_full,
            embedding_keys=embedding_keys,
            bench_batch_key=args.bench_batch_key,
            bench_label_key=args.bench_label_key,
            output_dir=args.output_dir,
        )
    else:
        print("\n  Skipping scIB benchmarking (--skip-benchmark).")

    # ══════════════════════════════════════════════════════════════════════
    #  9. Save integrated object — ALL genes
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("9. SAVING INTEGRATED OBJECT")
    print(f"{'='*60}")

    if "cell_id" not in adata_full.obs.columns:
        adata_full.obs["cell_id"] = adata_full.obs_names

    # Normalise .X for downstream use; raw counts stay in layers
    adata_full.X = adata_full.layers["counts"].copy()
    sc.pp.normalize_total(adata_full)
    sc.pp.log1p(adata_full)

    out_path = os.path.join(args.output_dir, "integrated.h5ad")
    adata_full.write_h5ad(out_path)
    print(f"  Saved {out_path}")
    print(f"  Shape  : {adata_full.n_obs:,} cells x {adata_full.n_vars:,} genes")
    print(f"  .X     : normalised, log1p")
    print(f"  layers : {list(adata_full.layers.keys())}")
    print(f"  obsm   : {list(adata_full.obsm.keys())}")

    # Summary
    print(f"\n{'='*60}")
    print(f"  Done. All outputs in {args.output_dir}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Post-integration inspection and iterative refinement.

Companion to integrate_sns_scvi.py.  Generates a diagnostic report from an
integrated h5ad, auto-flags suspicious clusters, and supports a
decisions.yaml workflow for reproducible iterative filtering.

Modes
-----
  # 1. Generate inspection report (default)
  python inspect_integration.py --input integrated.h5ad \
      --output-dir rounds/round_01/

  # 2. Apply decisions from a prior round and export filtered data
  python inspect_integration.py --input integrated.h5ad \
      --decisions rounds/round_01/decisions.yaml \
      --output-dir rounds/round_02/

  # decisions.yaml supports four actions (applied in order):
  #   keep_clusters  — whitelist: only these cluster IDs are retained
  #   keep_cells     — query whitelist: each query further restricts (AND)
  #   remove_clusters — blacklist: drop these clusters from the retained set
  #   remove_cells   — query blacklist: drop matching cells from the retained set

  # 3. Both at once (report on the current state, then filter)
  python inspect_integration.py --input integrated.h5ad \
      --decisions rounds/round_01/decisions.yaml \
      --output-dir rounds/round_02/ --report

Output
------
  Report mode:
    cluster_qc_summary.csv        Per-cluster QC table
    auto_flags.yaml               Suggested removals with reasons
    qc_summary_heatmap.png        Visual QC overview
    batch_composition.png         Batch breakdown per cluster
    marker_dotplot.png            Marker gene expression
    cluster_silhouettes.png       Per-cluster silhouette scores
    inspection_report.html        Self-contained HTML report

  Filter mode:
    cells_to_keep.csv             Obs names of retained cells
    filtered.h5ad                 Filtered AnnData ready for re-integration
    decisions_applied.yaml        Record of what was removed and why
    round_manifest.json           Provenance record for this round
"""

import argparse
import base64
import gc
import io
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


# ═══════════════════════════════════════════════════════════════════════════════
#  Model key extraction
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_setup_args(setup_args):
    """Build a keys dict from a scVI setup_args mapping."""
    result = {}
    bk = setup_args.get("batch_key")
    if bk:
        result["batch_key"] = bk
    cat = setup_args.get("categorical_covariate_keys") or []
    result["categorical_covariate_keys"] = cat if isinstance(cat, list) else [cat]
    cont = setup_args.get("continuous_covariate_keys") or []
    result["continuous_covariate_keys"] = cont if isinstance(cont, list) else [cont]
    return result


def _setup_args_from_pt(pt_path):
    """Try to extract setup_args from a raw model.pt torch save file.

    scvi-tools < 0.20 packed everything into a single model.pt instead of
    writing a separate attr_dict.json.  The nested structure varied across
    versions, so we try the most common paths.
    """
    try:
        import torch  # type: ignore[import-untyped]
    except ImportError:
        return None
    try:
        data = torch.load(pt_path, map_location="cpu", weights_only=False)
    except TypeError:
        data = torch.load(pt_path, map_location="cpu")

    # Common locations across scvi-tools 0.15–0.19
    candidate_paths = [
        ["attr_dict", "registry_", "setup_args"],          # 0.20 transition
        ["attr_dict", "adata_manager_", "registry", "setup_args"],  # 0.17-0.19
        ["registry_", "setup_args"],
        ["attr_dict", "setup_args"],
    ]
    for path in candidate_paths:
        d = data
        for key in path:
            if not isinstance(d, dict) or key not in d:
                d = None
                break
            d = d[key]
        if isinstance(d, dict) and "batch_key" in d:
            return d
    return None


def extract_keys_from_model(model_dir):
    """Read batch_key and covariate_keys from a saved scVI model directory.

    Tries attr_dict.json first (scvi-tools >= 0.20, the format produced by
    integrate_sns_scvi.py).  Falls back to loading model.pt directly via
    torch for models saved with older scvi-tools versions where attr_dict.json
    was not written.

    Returns a dict with keys batch_key, categorical_covariate_keys,
    continuous_covariate_keys, or None if neither source is readable.
    """
    model_dir = Path(model_dir)

    # ── attr_dict.json (scvi-tools >= 0.20) ──────────────────────────────────
    attr_path = model_dir / "attr_dict.json"
    if attr_path.exists():
        try:
            with open(attr_path) as f:
                attr = json.load(f)
            setup_args = attr.get("registry_", {}).get("setup_args", {})
            if setup_args:
                return _parse_setup_args(setup_args)
        except Exception:
            pass

    # ── model.pt fallback (older scvi-tools) ─────────────────────────────────
    pt_path = model_dir / "model.pt"
    if pt_path.exists():
        try:
            setup_args = _setup_args_from_pt(pt_path)
            if setup_args:
                print(f"  (Keys read from model.pt — no attr_dict.json found)")
                return _parse_setup_args(setup_args)
        except Exception as e:
            print(f"  [WARN] Could not parse model.pt: {e}")

    return None


# ═══════════════════════════════════════════════════════════════════════════════
#  Default marker genes (peripheral / sympathetic neurons)
# ═══════════════════════════════════════════════════════════════════════════════

DEFAULT_MARKERS = {
    "Pan-neuronal":     ["Snap25", "Tubb3", "Rbfox3", "Elavl4", "Isl1"],
    "Noradrenergic":    ["Th", "Dbh", "Ddc", "Slc6a2"],
    "Cholinergic":      ["Chat", "Slc18a3", "Slc5a7"],
    "Glutamatergic":    ["Slc17a6", "Slc17a7"],
    "GABAergic":        ["Slc32a1", "Gad1", "Gad2"],
    "Nitrergic":        ["Nos1"],
    "Neuropeptides":    ["Npy", "Sst", "Vip", "Pdyn", "Chga"],
    "Transcription":    ["Phox2b", "Phox2a", "Sox6", "Shox2"],
    "Satellite glia":   ["Sox10", "Fabp7", "S100b"],
    "Immune":           ["Ptprc", "Cd68"],
    "Mitochondrial":    ["mt-Co1", "mt-Co2", "mt-Cytb"],
}

# Markers used for neuronal cluster identification (fraction-expressed threshold)
DEFAULT_NEURONAL_MARKERS = [
    "Th", "Snap25", "Phox2b", "Dbh", "Chat", "Slc18a2", "Slc17a6",
]


# ═══════════════════════════════════════════════════════════════════════════════
#  Cluster QC summary
# ═══════════════════════════════════════════════════════════════════════════════

def cluster_qc_summary(adata, cluster_key="leiden", batch_key="data_origin",
                        latent_key=None, covariate_keys=None):
    """Compute per-cluster QC statistics.

    Returns a DataFrame with one row per cluster, columns:
        n_cells, median_genes, median_counts, median_pct_mt,
        n_batches, dominant_batch, dominant_batch_frac,
        dominant_{cov}_frac  (for each key in covariate_keys),
        silhouette (if latent_key provided)
    """
    obs = adata.obs.copy()
    clusters = sorted(obs[cluster_key].unique(),
                      key=lambda x: int(x) if x.isdigit() else x)

    rows = []
    for cl in clusters:
        mask = obs[cluster_key] == cl
        sub = obs[mask]
        n = int(mask.sum())

        row = {"cluster": cl, "n_cells": n}

        for col, label in [("n_genes_by_counts", "median_genes"),
                           ("total_counts", "median_counts"),
                           ("pct_counts_mt", "median_pct_mt"),
                           ("pct_counts_ribo", "median_pct_ribo")]:
            if col in sub.columns:
                row[label] = float(sub[col].median())

        if batch_key and batch_key in sub.columns:
            vc = sub[batch_key].value_counts()
            row["n_batches"] = int((vc > 0).sum())
            row["dominant_batch"] = vc.index[0]
            row["dominant_batch_frac"] = float(vc.iloc[0] / n)

        for cov_key in (covariate_keys or []):
            if cov_key in sub.columns:
                vc = sub[cov_key].value_counts()
                row[f"dominant_{cov_key}_frac"] = float(vc.iloc[0] / n)

        rows.append(row)

    df = pd.DataFrame(rows).set_index("cluster")

    # Silhouette per cluster (if latent space available)
    if latent_key and latent_key in adata.obsm:
        from sklearn.metrics import silhouette_samples
        X = adata.obsm[latent_key]
        labels = obs[cluster_key].values
        sil = silhouette_samples(X, labels, metric="euclidean")
        obs["_sil"] = sil
        sil_median = obs.groupby(cluster_key)["_sil"].median()
        df["silhouette"] = df.index.map(sil_median)

    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Automated cluster flagging
# ═══════════════════════════════════════════════════════════════════════════════

def auto_flag_clusters(summary_df, mt_threshold=15.0, min_genes_threshold=400,
                        min_cells=20, single_batch_threshold=0.90,
                        silhouette_threshold=-0.05, covariate_keys=None):
    """Flag clusters that meet common removal criteria.

    Returns a dict: {cluster_id: [reason1, reason2, ...]}
    """
    flags = {}

    for cl, row in summary_df.iterrows():
        reasons = []

        if "median_pct_mt" in row and row["median_pct_mt"] > mt_threshold:
            reasons.append(
                f"High MT fraction (median {row['median_pct_mt']:.1f}% "
                f"> {mt_threshold}%)")

        if "median_genes" in row and row["median_genes"] < min_genes_threshold:
            reasons.append(
                f"Low gene complexity (median {row['median_genes']:.0f} "
                f"< {min_genes_threshold})")

        if row["n_cells"] < min_cells:
            reasons.append(
                f"Very small cluster ({row['n_cells']} cells "
                f"< {min_cells})")

        if ("dominant_batch_frac" in row
                and row["dominant_batch_frac"] > single_batch_threshold):
            reasons.append(
                f"Single-batch dominated ({row['dominant_batch']}: "
                f"{row['dominant_batch_frac']:.0%} "
                f"> {single_batch_threshold:.0%})")

        for cov_key in (covariate_keys or []):
            col = f"dominant_{cov_key}_frac"
            if col in row and row[col] > single_batch_threshold:
                reasons.append(
                    f"Single-covariate dominated ({cov_key}: "
                    f"{row[col]:.0%} > {single_batch_threshold:.0%})")

        if ("silhouette" in row
                and row["silhouette"] < silhouette_threshold):
            reasons.append(
                f"Poor separation (silhouette {row['silhouette']:.3f} "
                f"< {silhouette_threshold})")

        if reasons:
            flags[str(cl)] = reasons

    return flags


# ═══════════════════════════════════════════════════════════════════════════════
#  Diagnostic plots
# ═══════════════════════════════════════════════════════════════════════════════

def _fig_to_base64(fig):
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _fmt_annot(value, col_name):
    """Format a heatmap annotation value with column-appropriate precision."""
    if pd.isna(value):
        return ""
    if col_name in ("n_cells", "n_batches"):
        return str(int(value))
    if col_name in ("median_genes", "median_counts"):
        return f"{int(value):,}"
    if col_name == "silhouette":
        return f"{value:.3f}"
    if col_name.endswith("_frac"):
        return f"{value:.2f}"
    return f"{value:.1f}"


def plot_qc_summary_heatmap(summary_df, flags, output_dir):
    """Heatmap of per-cluster QC metrics with flagged clusters highlighted."""
    cols = [c for c in ["n_cells", "median_genes", "median_counts",
                         "median_pct_mt", "median_pct_ribo",
                         "n_batches", "dominant_batch_frac", "silhouette"]
            if c in summary_df.columns]

    plot_df = summary_df[cols].copy()

    # Normalise each column to [0, 1] for the heatmap
    norm_df = plot_df.copy()
    for c in cols:
        cmin, cmax = norm_df[c].min(), norm_df[c].max()
        if cmax > cmin:
            norm_df[c] = (norm_df[c] - cmin) / (cmax - cmin)
        else:
            norm_df[c] = 0.5

    # Build per-column annotation strings to avoid overflow from large values
    annot = np.empty((len(plot_df), len(cols)), dtype=object)
    for j, col in enumerate(cols):
        for i in range(len(plot_df)):
            annot[i, j] = _fmt_annot(plot_df.iloc[i, j], col)

    n_clusters = len(plot_df)
    fig_h = max(4, n_clusters * 0.35 + 1)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    sns.heatmap(norm_df, annot=annot, fmt="",
                cmap="YlOrRd", linewidths=0.5, ax=ax,
                cbar_kws={"label": "Normalised (row-wise)", "shrink": 0.6})

    # Highlight flagged rows
    flagged_indices = [i for i, cl in enumerate(summary_df.index)
                       if str(cl) in flags]
    for idx in flagged_indices:
        ax.add_patch(plt.Rectangle((0, idx), len(cols), 1,
                                    fill=False, edgecolor="red",
                                    linewidth=2.5))

    ax.set_title("Per-cluster QC summary (red border = auto-flagged)")
    ax.set_ylabel("Cluster")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "qc_summary_heatmap.png"),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


def _plot_key_composition(adata, cluster_key, groupby_key, output_dir,
                           output_fname, title_label=None):
    """Stacked bar chart: composition of groupby_key per cluster."""
    ct = pd.crosstab(adata.obs[cluster_key], adata.obs[groupby_key],
                     normalize="index")
    clusters = sorted(ct.index, key=lambda x: int(x) if x.isdigit() else x)
    ct = ct.loc[clusters]

    fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.5), 5))
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.85, edgecolor="black",
            linewidth=0.3)
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("Cluster")
    label = title_label or groupby_key
    ax.set_title(f"{label} composition per cluster ({groupby_key})")
    ax.legend(title=groupby_key, bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=7)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, output_fname),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


def plot_batch_composition(adata, cluster_key, batch_key, output_dir):
    """Stacked bar chart: batch composition per cluster."""
    return _plot_key_composition(
        adata, cluster_key, batch_key, output_dir,
        "batch_composition.png", title_label="Batch")


def plot_covariate_compositions(adata, cluster_key, covariate_keys, output_dir):
    """One stacked bar chart per covariate key.  Returns {key: b64} dict."""
    results = {}
    for key in covariate_keys:
        if key not in adata.obs.columns:
            print(f"  [WARN] Covariate key '{key}' not in obs — skipping.")
            continue
        fname = f"covariate_{key}_composition.png"
        b64 = _plot_key_composition(
            adata, cluster_key, key, output_dir, fname,
            title_label=f"Covariate ({key})")
        results[key] = b64
    return results


def plot_marker_dotplot(adata, cluster_key, marker_dict, output_dir):
    """Dot plot of marker genes per cluster."""
    # Filter to genes actually in the dataset
    filtered = {}
    for group, genes in marker_dict.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            filtered[group] = present

    if not filtered:
        print("  [WARN] No marker genes found in dataset — skipping dot plot.")
        return None

    # Ensure .X is normalised for plotting
    if adata.X.max() > 50:  # likely raw counts
        tmp = adata.copy()
        sc.pp.normalize_total(tmp)
        sc.pp.log1p(tmp)
    else:
        tmp = adata

    fig = sc.pl.dotplot(tmp, var_names=filtered, groupby=cluster_key,
                        standard_scale="var", return_fig=True, show=False)
    fig.savefig(os.path.join(output_dir, "marker_dotplot.png"),
                dpi=150, bbox_inches="tight")
    # Read back the saved PNG for base64
    with open(os.path.join(output_dir, "marker_dotplot.png"), "rb") as f:
        b64 = base64.b64encode(f.read()).decode("utf-8")
    plt.close("all")
    return b64


def plot_cluster_silhouettes(summary_df, output_dir):
    """Bar chart of median silhouette score per cluster."""
    if "silhouette" not in summary_df.columns:
        return None

    clusters = summary_df.index.tolist()
    sil = summary_df["silhouette"].values

    colors = ["#CC3311" if s < 0 else "#EE8833" if s < 0.1 else "#4477AA"
              for s in sil]

    fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.5), 4))
    ax.bar(range(len(clusters)), sil, color=colors, edgecolor="black",
           linewidth=0.5)
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(0.1, color="grey", linestyle="--", linewidth=0.8,
               label="0.1 threshold")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Median silhouette score")
    ax.set_title("Cluster separation in latent space")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "cluster_silhouettes.png"),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


def plot_umap_overview(adata, cluster_key, batch_key, output_dir,
                        umap_key="X_umap"):
    """2-panel UMAP: clusters + batch."""
    if umap_key not in adata.obsm:
        print(f"  [WARN] '{umap_key}' not in obsm — skipping UMAP overview.")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.embedding(adata, basis=umap_key, color=cluster_key, ax=axes[0],
                    legend_loc="on data", legend_fontsize=7, show=False,
                    size=8, frameon=False)
    axes[0].set_title(f"Clusters ({cluster_key})")

    if batch_key and batch_key in adata.obs.columns:
        sc.pl.embedding(adata, basis=umap_key, color=batch_key, ax=axes[1],
                        legend_loc="right margin", legend_fontsize=7,
                        show=False, size=8, frameon=False)
        axes[1].set_title(f"Batch ({batch_key})")
    else:
        axes[1].axis("off")
        axes[1].set_title(f"Batch key '{batch_key}' not in obs")

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "umap_overview.png"),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


# ═══════════════════════════════════════════════════════════════════════════════
#  Neuronal marker analysis
# ═══════════════════════════════════════════════════════════════════════════════

def calculate_percent_expressed(adata, gene_list, groupby_key, threshold=0.0):
    """Fraction of cells per cluster expressing each gene above threshold.

    Uses layers["counts"] (raw integers) when available, otherwise .X.
    Returns a DataFrame: clusters (rows) × genes (columns), values in [0, 1].
    Clusters are sorted naturally by ID.
    """
    import scipy.sparse as sp

    available = [g for g in gene_list if g in adata.var_names]
    if not available:
        print("  [WARN] No target genes found in dataset — skipping fraction table.")
        return pd.DataFrame()

    X = (adata[:, available].layers["counts"]
         if "counts" in adata.layers
         else adata[:, available].X)
    if sp.issparse(X) and not isinstance(X, sp.csr_matrix):
        X = X.tocsr()

    expressed = (X > threshold).toarray() if sp.issparse(X) else (np.asarray(X) > threshold)

    df = pd.DataFrame(expressed, index=adata.obs_names, columns=available)
    df["_g"] = adata.obs[groupby_key].values

    num = df.groupby("_g")[available].sum()
    total = adata.obs[groupby_key].value_counts()
    frac = num.div(total, axis=0)

    frac.index = frac.index.astype(str)
    frac = frac.loc[sorted(frac.index, key=lambda x: int(x) if x.isdigit() else x)]
    return frac


def plot_marker_fraction_heatmap(fraction_table, output_dir):
    """Heatmap of fraction of cells expressing each marker gene per cluster."""
    if fraction_table.empty:
        return None

    n_r, n_c = len(fraction_table), len(fraction_table.columns)
    fig, ax = plt.subplots(
        figsize=(max(8, n_c * 0.65 + 2), max(5, n_r * 0.35 + 1)))
    sns.heatmap(fraction_table, annot=True, fmt=".2f", cmap="YlGnBu",
                vmin=0, vmax=1, linewidths=0.5, ax=ax,
                cbar_kws={"label": "Fraction expressing", "shrink": 0.6})
    ax.set_title("Fraction of cells expressing marker genes per cluster")
    ax.set_ylabel("Cluster")
    ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "marker_fraction_heatmap.png"),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


def plot_cluster_marker_pca(fraction_table, output_dir, color_gene="Snap25"):
    """PCA of clusters by marker expression profile, coloured by one gene's fraction."""
    if fraction_table.empty:
        return None

    mat = fraction_table.fillna(0).values
    n_comp = min(2, mat.shape[0] - 1, mat.shape[1])
    if n_comp < 1 or len(fraction_table) < 3:
        return None

    from sklearn.decomposition import PCA
    try:
        pca = PCA(n_components=n_comp).fit(mat)
        xy = pca.transform(mat)
        pct_var = pca.explained_variance_ratio_ * 100
    except Exception as e:
        print(f"  [WARN] Marker PCA failed: {e}")
        return None

    x = xy[:, 0]
    y = xy[:, 1] if xy.shape[1] > 1 else np.zeros(len(xy))

    if color_gene in fraction_table.columns:
        color_vals = fraction_table[color_gene].fillna(0).values
        color_label = color_gene
    else:
        color_vals = x
        color_label = "PC1"

    fig, ax = plt.subplots(figsize=(7, 6))
    sc_plot = ax.scatter(x, y, c=color_vals, cmap="viridis", s=80,
                         edgecolors="black", linewidths=0.5)
    for i, cl in enumerate(fraction_table.index):
        ax.annotate(str(cl), (x[i], y[i]), fontsize=8, ha="center",
                    va="bottom", xytext=(0, 5), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({pct_var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pct_var[1]:.1f}%)" if xy.shape[1] > 1 else "PC2")
    ax.set_title(f"Cluster PCA — marker expression (colour: {color_label})")
    plt.colorbar(sc_plot, ax=ax, label=f"Fraction {color_label}")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "cluster_marker_pca.png"),
                dpi=150, bbox_inches="tight")
    b64 = _fig_to_base64(fig)
    return b64


def identify_neuronal_clusters(fraction_table, markers, cutoff=0.50):
    """Clusters where ≥1 neuronal marker exceeds cutoff fraction of cells.

    Returns (cluster_id_list, fraction_subtable_for_those_clusters).
    The subtable contains all genes from fraction_table (not just the
    selection markers) so the full expression profile is available for
    the HTML report.
    """
    if fraction_table.empty:
        return [], pd.DataFrame()

    present = [m for m in markers if m in fraction_table.columns]
    if not present:
        print(f"  [WARN] None of the neuronal markers {markers} found in "
              f"fraction table — skipping neuronal cluster identification.")
        return [], pd.DataFrame()

    meets = (fraction_table[present] > cutoff).any(axis=1)
    neuronal_ids = fraction_table.index[meets].tolist()
    return neuronal_ids, fraction_table.loc[meets]


# ═══════════════════════════════════════════════════════════════════════════════
#  YAML I/O
# ═══════════════════════════════════════════════════════════════════════════════

def write_auto_flags_yaml(flags, summary_df, output_path, thresholds,
                           neuronal_clusters=None):
    """Write auto-flagged clusters to a YAML template the user can edit."""
    lines = [
        "# ──────────────────────────────────────────────────────────────────",
        "# Auto-generated cluster flags — edit this file to create your",
        "# decisions.yaml for the next integration round.",
        "#",
        "# Actions (applied in this order):",
        "#   keep_clusters:    only retain these cluster IDs (all others excluded)",
        "#   keep_cells:       restrict to cells matching each query (AND-chained)",
        "#   remove_clusters:  drop these cluster IDs from the retained set",
        "#   remove_cells:     drop cells matching each query from the retained set",
        "#   notes:            free-text notes for provenance",
        "#",
        "# To accept the auto-flags as-is, rename this file to decisions.yaml",
        "# and re-run with --decisions.",
        "# ──────────────────────────────────────────────────────────────────",
        "",
        f"generated: \"{datetime.now().isoformat()}\"",
        "",
        "thresholds:",
        f"  mt_pct: {thresholds['mt']}",
        f"  min_genes: {thresholds['min_genes']}",
        f"  min_cells: {thresholds['min_cells']}",
        f"  single_batch_frac: {thresholds['single_batch']}",
        "",
    ]

    if flags:
        lines.append("remove_clusters:")
        for cl, reasons in flags.items():
            lines.append(f"  - cluster: \"{cl}\"")
            lines.append(f"    n_cells: {int(summary_df.loc[cl, 'n_cells'])}")
            lines.append(f"    reasons:")
            for r in reasons:
                lines.append(f"      - \"{r}\"")
        lines.append("")
    else:
        lines.append("remove_clusters: []  # no clusters auto-flagged")
        lines.append("")

    if neuronal_clusters:
        lines.extend([
            "# ── Suggested neuronal clusters (marker fraction > threshold) ──────",
            "# These clusters showed strong expression of neuronal marker genes.",
            "# Uncomment keep_clusters below to restrict to neurons only:",
            "# keep_clusters:",
        ])
        for cl in neuronal_clusters:
            lines.append(f"#   - cluster: \"{cl}\"")
        lines.append("")
    else:
        lines.extend([
            "# Uncomment to keep only specific clusters (all others excluded):",
            "# keep_clusters:",
            "#   - cluster: \"1\"",
            "#     reason: \"Confirmed high-quality neurons\"",
            "",
        ])

    lines.extend([
        "# Uncomment to restrict to cells matching a query (AND-chained):",
        "# keep_cells:",
        "#   - query: \"n_genes_by_counts > 500\"",
        "#     reason: \"Minimum gene complexity\"",
        "",
        "# Uncomment and edit to remove specific cells by query:",
        "# remove_cells:",
        "#   - query: \"(tech == 'scale') & (neuron_type.isin(['GABAergic']))\"",
        "#     reason: \"Contaminating non-target neurons\"",
        "",
        "notes: \"\"",
    ])

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def load_decisions(path):
    """Load a decisions.yaml file. Supports both YAML and simplified format."""
    try:
        import yaml
        with open(path) as f:
            return yaml.safe_load(f)
    except ImportError:
        # Fallback: parse the simplified YAML subset we generate
        import re
        with open(path) as f:
            text = f.read()

        decisions = {
            "keep_clusters": [], "keep_cells": [],
            "remove_clusters": [], "remove_cells": [],
        }

        # Track which section each cluster/query belongs to by scanning
        # section headers and collecting entries that follow them.
        section_re = re.compile(
            r'^(keep_clusters|keep_cells|remove_clusters|remove_cells)\s*:',
            re.MULTILINE)
        cluster_re = re.compile(r'-\s+cluster:\s*["\']?(\w+)["\']?')
        query_re   = re.compile(r'-\s+query:\s*["\'](.+?)["\']')
        reason_re  = re.compile(r'reason:\s*["\'](.+?)["\']')

        # Split text into labelled sections
        sections = []
        for m in section_re.finditer(text):
            sections.append((m.group(1), m.end()))
        sections.append((None, len(text)))  # sentinel

        for i, (section, start) in enumerate(sections[:-1]):
            chunk = text[start:sections[i + 1][1]]
            if section in ("keep_clusters", "remove_clusters"):
                for cm in cluster_re.finditer(chunk):
                    decisions[section].append({"cluster": cm.group(1)})
            elif section in ("keep_cells", "remove_cells"):
                for qm in query_re.finditer(chunk):
                    q = qm.group(1)
                    rest = chunk[qm.end():]
                    rm = reason_re.search(rest[:200])
                    reason = rm.group(1) if rm else ""
                    decisions[section].append({"query": q, "reason": reason})

        return decisions


# ═══════════════════════════════════════════════════════════════════════════════
#  HTML report
# ═══════════════════════════════════════════════════════════════════════════════

def generate_html_report(summary_df, flags, images, adata, output_path,
                          cluster_key, batch_key, covariate_keys=None,
                          neuronal_clusters=None, neuronal_fraction_table=None):
    """Generate a self-contained HTML inspection report."""
    n_flagged = len(flags)
    n_flagged_cells = 0
    for cl in flags:
        if cl in summary_df.index:
            n_flagged_cells += int(summary_df.loc[cl, "n_cells"])

    flag_rows = ""
    for cl, reasons in flags.items():
        n = int(summary_df.loc[cl, "n_cells"]) if cl in summary_df.index else "?"
        reason_html = "<br>".join(f"&bull; {r}" for r in reasons)
        flag_rows += f"""
        <tr style="background-color: #fff3f3;">
          <td><strong>{cl}</strong></td>
          <td>{n}</td>
          <td>{reason_html}</td>
        </tr>"""

    summary_html = summary_df.round(2).to_html(
        classes="summary-table", border=0)

    img_sections = ""
    for title, b64 in images.items():
        if b64:
            img_sections += f"""
            <h2>{title}</h2>
            <img src="data:image/png;base64,{b64}"
                 style="max-width:100%; border:1px solid #ddd; border-radius:4px;">
            """

    # Extra stat boxes for each covariate key
    cov_stats_html = ""
    for cov_key in (covariate_keys or []):
        if cov_key in adata.obs.columns:
            n_vals = adata.obs[cov_key].nunique()
            cov_stats_html += f"""
  <div class="stat">
    <div class="stat-value">{n_vals}</div>
    <div class="stat-label">Covariate values ({cov_key})</div>
  </div>"""

    cov_footer = ""
    if covariate_keys:
        cov_footer = f" &middot; Covariate keys: {', '.join(covariate_keys)}"

    batch_stat = ""
    if batch_key and batch_key in adata.obs.columns:
        batch_stat = f"""
  <div class="stat">
    <div class="stat-value">{adata.obs[batch_key].nunique()}</div>
    <div class="stat-label">Batches ({batch_key})</div>
  </div>"""

    neuro_stat = ""
    if neuronal_clusters:
        neuro_stat = f"""
  <div class="stat">
    <div class="stat-value" style="color:#2a7a2a;">{len(neuronal_clusters)}</div>
    <div class="stat-label">Likely neuronal clusters</div>
  </div>"""

    neuronal_section = ""
    if neuronal_clusters:
        ft = neuronal_fraction_table
        if ft is not None and not ft.empty:
            header_cells = "".join(f"<th>{g}</th>" for g in ft.columns)
            body_rows = ""
            for cl in neuronal_clusters:
                if cl in ft.index:
                    cells = "".join(
                        f"<td>{ft.loc[cl, g]:.0%}</td>" for g in ft.columns)
                    body_rows += (
                        f"<tr><td style='text-align:left'><strong>{cl}"
                        f"</strong></td>{cells}</tr>")
            neuronal_section = (
                "<h2>Likely Neuronal Clusters</h2>"
                f"<p>{len(neuronal_clusters)} cluster(s) with ≥1 neuronal "
                f"marker above threshold. Suggested as <code>keep_clusters</code> "
                f"in auto_flags.yaml.</p>"
                f"<table><tr><th>Cluster</th>{header_cells}</tr>"
                f"{body_rows}</table>")
        else:
            neuronal_section = (
                "<h2>Likely Neuronal Clusters</h2>"
                f"<p>Clusters: {', '.join(neuronal_clusters)}</p>")

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Integration Inspection Report</title>
<style>
  body {{ font-family: -apple-system, 'Segoe UI', Roboto, sans-serif;
         max-width: 1200px; margin: 40px auto; padding: 0 20px;
         color: #333; background: #fafafa; }}
  h1 {{ border-bottom: 3px solid #4477AA; padding-bottom: 8px; }}
  h2 {{ color: #4477AA; margin-top: 40px; }}
  .summary-box {{ background: white; padding: 20px; border-radius: 8px;
                  box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin: 20px 0; }}
  .stat {{ display: inline-block; margin: 0 30px 10px 0; }}
  .stat-value {{ font-size: 28px; font-weight: 700; color: #1B2A4A; }}
  .stat-label {{ font-size: 13px; color: #888; }}
  .flag-warn {{ background: #FFF3CD; border: 1px solid #FFD93D;
                border-radius: 6px; padding: 15px; margin: 20px 0; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ padding: 6px 10px; border: 1px solid #ddd; text-align: right; }}
  th {{ background: #4477AA; color: white; }}
  tr:nth-child(even) {{ background: #f9f9f9; }}
  .summary-table td:first-child, .summary-table th:first-child {{
    text-align: left; }}
  img {{ margin: 10px 0; }}
</style>
</head>
<body>

<h1>Integration Inspection Report</h1>
<p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>

<div class="summary-box">
  <div class="stat">
    <div class="stat-value">{adata.n_obs:,}</div>
    <div class="stat-label">Total cells</div>
  </div>
  <div class="stat">
    <div class="stat-value">{adata.n_vars:,}</div>
    <div class="stat-label">Genes</div>
  </div>
  <div class="stat">
    <div class="stat-value">{len(summary_df)}</div>
    <div class="stat-label">Clusters</div>
  </div>
  {batch_stat}
  {cov_stats_html}
  {neuro_stat}
</div>

{"<div class='flag-warn'>" +
 f"<strong>{n_flagged} cluster(s) flagged</strong> for review " +
 f"({n_flagged_cells:,} cells). See auto_flags.yaml." +
 "</div>" if n_flagged else
 "<p style='color:green;'>No clusters auto-flagged.</p>"}

{img_sections}

<h2>Cluster QC Summary Table</h2>
{summary_html}

{"<h2>Auto-flagged Clusters</h2><table>" +
 "<tr><th>Cluster</th><th>N cells</th><th>Reasons</th></tr>" +
 flag_rows + "</table>" if flag_rows else ""}

{neuronal_section}

<hr>
<p style="font-size:12px; color:#999;">
  Cluster key: {cluster_key} &middot;
  Batch key: {batch_key or "none"} &middot;
  Cells: {adata.n_obs:,} &middot;
  Genes: {adata.n_vars:,}{cov_footer}
</p>
</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)


# ═══════════════════════════════════════════════════════════════════════════════
#  Apply decisions
# ═══════════════════════════════════════════════════════════════════════════════

def apply_decisions(adata, decisions, cluster_key="leiden", output_dir="."):
    """Apply a decisions dict to filter an AnnData.

    Order of operations:
      1. keep_clusters  — restrict to listed cluster IDs
      2. keep_cells     — restrict to cells matching each query (AND-chained)
      3. remove_clusters — drop cluster IDs from the retained set
      4. remove_cells   — drop cells matching each query from the retained set

    Writes:
      - cells_to_keep.csv
      - filtered.h5ad
      - decisions_applied.yaml
      - round_manifest.json
    """
    n_before = adata.n_obs
    keep = pd.Series(True, index=adata.obs_names)
    log_entries = []

    # ── Keep clusters (restrict to listed IDs) ──────────────────────────────
    keep_cluster_ids = [str(e["cluster"])
                        for e in decisions.get("keep_clusters", [])]
    if keep_cluster_ids:
        mask = adata.obs[cluster_key].astype(str).isin(keep_cluster_ids)
        n_excluded = int((~mask & keep).sum())
        keep &= mask
        reasons = "; ".join(
            e.get("reason", "") for e in decisions["keep_clusters"]
            if e.get("reason"))
        log_entries.append({
            "type": "keep_clusters",
            "clusters": keep_cluster_ids,
            "n_excluded": n_excluded,
            "reason": reasons or "keep listed clusters only",
        })
        print(f"  Keep clusters {keep_cluster_ids}: "
              f"{n_excluded:,} cells outside excluded")

    # ── Keep cells (restrict by query, AND-chained) ──────────────────────────
    for entry in decisions.get("keep_cells", []):
        query = entry["query"]
        reason = entry.get("reason", "")
        try:
            query_mask = adata.obs.eval(query)
            n_excluded = int((~query_mask & keep).sum())
            keep &= query_mask
            log_entries.append({
                "type": "keep_query",
                "query": query,
                "n_excluded": n_excluded,
                "reason": reason,
            })
            print(f"  Keep query '{query}': "
                  f"{n_excluded:,} cells excluded ({reason})")
        except Exception as e:
            print(f"  [WARN] Keep query failed: {query} — {e}")

    # ── Remove clusters ──────────────────────────────────────────────────────
    for entry in decisions.get("remove_clusters", []):
        cl = str(entry["cluster"])
        mask = adata.obs[cluster_key].astype(str) == cl
        n_removed = int((mask & keep).sum())
        keep[mask] = False
        reasons = entry.get("reasons", entry.get("reason", "unspecified"))
        if isinstance(reasons, list):
            reasons = "; ".join(reasons)
        log_entries.append({
            "type": "remove_cluster",
            "cluster": cl,
            "n_removed": n_removed,
            "reasons": reasons,
        })
        print(f"  Remove cluster {cl}: {n_removed:,} cells ({reasons})")

    # ── Remove cells by query ────────────────────────────────────────────────
    for entry in decisions.get("remove_cells", []):
        query = entry["query"]
        reason = entry.get("reason", "")
        try:
            query_mask = adata.obs.eval(query)
            n_removed = int((query_mask & keep).sum())
            keep[query_mask] = False
            log_entries.append({
                "type": "remove_query",
                "query": query,
                "n_removed": n_removed,
                "reason": reason,
            })
            print(f"  Remove query '{query}': {n_removed:,} cells ({reason})")
        except Exception as e:
            print(f"  [WARN] Query failed: {query} — {e}")

    n_after = int(keep.sum())
    n_total_removed = n_before - n_after
    print(f"\n  Summary: {n_before:,} → {n_after:,} cells "
          f"({n_total_removed:,} removed, {n_total_removed/n_before:.1%})")

    # Export cells to keep
    kept_names = adata.obs_names[keep].tolist()
    cells_path = os.path.join(output_dir, "cells_to_keep.csv")
    pd.DataFrame({"obs_name": kept_names}).to_csv(cells_path, index=False)
    print(f"  Saved {cells_path} ({len(kept_names):,} cells)")

    # Export filtered h5ad
    adata_filtered = adata[keep].copy()
    filtered_path = os.path.join(output_dir, "filtered.h5ad")
    adata_filtered.write_h5ad(filtered_path)
    print(f"  Saved {filtered_path}")

    # Write applied decisions log
    applied = {
        "timestamp": datetime.now().isoformat(),
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "n_removed": n_total_removed,
        "actions": log_entries,
        "notes": decisions.get("notes", ""),
    }

    applied_path = os.path.join(output_dir, "decisions_applied.yaml")
    with open(applied_path, "w") as f:
        f.write(f"timestamp: \"{applied['timestamp']}\"\n")
        f.write(f"n_cells_before: {n_before}\n")
        f.write(f"n_cells_after: {n_after}\n")
        f.write(f"n_removed: {n_total_removed}\n")
        f.write(f"notes: \"{applied.get('notes', '')}\"\n")
        f.write("actions:\n")
        for a in log_entries:
            f.write(f"  - type: \"{a['type']}\"\n")
            if "cluster" in a:
                f.write(f"    cluster: \"{a['cluster']}\"\n")
            if "clusters" in a:
                f.write(f"    clusters: {a['clusters']}\n")
            if "query" in a:
                f.write(f"    query: \"{a['query']}\"\n")
            # keep actions use n_excluded; remove actions use n_removed
            count = a.get("n_excluded", a.get("n_removed", 0))
            count_key = "n_excluded" if "n_excluded" in a else "n_removed"
            f.write(f"    {count_key}: {count}\n")
            r = a.get("reasons", a.get("reason", ""))
            f.write(f"    reason: \"{r}\"\n")
    print(f"  Saved {applied_path}")

    # Round manifest
    def _action_summary(a):
        if a["type"] == "keep_clusters":
            return f"keep_clusters {a['clusters']} (-{a.get('n_excluded', 0)})"
        if a["type"] == "keep_query":
            return f"keep_query: {a.get('query', '?')} (-{a.get('n_excluded', 0)})"
        return (f"{a['type']}: {a.get('cluster', a.get('query', '?'))} "
                f"(-{a.get('n_removed', 0)})")

    manifest_path = os.path.join(output_dir, "round_manifest.json")
    manifest = {
        "timestamp": applied["timestamp"],
        "input_cells": n_before,
        "output_cells": n_after,
        "removed": n_total_removed,
        "removal_pct": round(n_total_removed / n_before * 100, 1),
        "actions_summary": [_action_summary(a) for a in log_entries],
        "output_files": {
            "filtered_h5ad": "filtered.h5ad",
            "cells_to_keep": "cells_to_keep.csv",
        },
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Saved {manifest_path}")

    return adata_filtered


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Post-integration inspection and iterative refinement.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate inspection report
  python inspect_integration.py --input integrated.h5ad \\
      --output-dir rounds/round_01/

  # Apply decisions and export filtered data
  python inspect_integration.py --input integrated.h5ad \\
      --decisions rounds/round_01/decisions.yaml \\
      --output-dir rounds/round_02/

  # Custom marker genes (JSON: {"group": ["Gene1", "Gene2"]})
  python inspect_integration.py --input integrated.h5ad \\
      --markers my_markers.json

Typical iterative cycle:
  1. python integrate_sns_scvi.py --input data.h5ad --output-dir round_01/
  2. python inspect_integration.py --input round_01/integrated.h5ad \\
         --output-dir round_01/
  3. Edit round_01/auto_flags.yaml → save as round_01/decisions.yaml
  4. python inspect_integration.py --input round_01/integrated.h5ad \\
         --decisions round_01/decisions.yaml --output-dir round_02/
  5. python integrate_sns_scvi.py --input round_02/filtered.h5ad \\
         --output-dir round_02/ --skip-qc-filter
""",
    )

    # ── I/O ──
    parser.add_argument("--input", required=True,
                        help="Path to integrated h5ad.")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory (default: current dir).")

    # ── Cluster / batch / covariate keys ──
    parser.add_argument("--cluster-key", default=None,
                        help="Obs column for clusters (default: auto-detect "
                             "'leiden' or first 'leiden_*').")
    parser.add_argument("--batch-key", default=None,
                        help="Obs column for batch (default: auto-detect from "
                             "model or 'data_origin').")
    parser.add_argument("--covariate-keys", nargs="*", default=None,
                        help="Additional obs columns treated as covariates "
                             "(default: auto-detect from model).")
    parser.add_argument("--model-dir", default=None,
                        help="scVI model directory containing attr_dict.json. "
                             "Used to auto-detect batch_key and covariate_keys. "
                             "Auto-searched as scvi_model_* near the input "
                             "h5ad if omitted.")
    parser.add_argument("--latent-key", default=None,
                        help="Obsm key for latent space, used for silhouette "
                             "(default: auto-detect 'X_scVI' or first "
                             "'X_scVI_*').")
    parser.add_argument("--umap-key", default=None,
                        help="Obsm key for UMAP (default: auto-detect).")

    # ── Markers ──
    parser.add_argument("--markers", default=None,
                        help="JSON file with marker gene dict. "
                             "Omit to use built-in peripheral neuron markers.")
    parser.add_argument("--neuronal-markers", nargs="*", default=None,
                        help="Genes used to identify neuronal clusters by "
                             "fraction-expressed threshold. Default: "
                             "Th Snap25 Phox2b Dbh Chat Slc18a2 Slc17a6")
    parser.add_argument("--neuronal-cutoff", type=float, default=0.50,
                        help="Fraction threshold for neuronal cluster calling "
                             "(default: 0.50).")
    parser.add_argument("--marker-threshold", type=float, default=0.0,
                        help="Expression threshold used to call a cell "
                             "'expressing' in the percent-expressed calculation "
                             "(default: 0 = any positive value; use 1 for "
                             "raw count > 1).")
    parser.add_argument("--pca-color-gene", default="Snap25",
                        help="Gene to colour the cluster marker PCA by "
                             "(default: Snap25).")

    # ── Flagging thresholds ──
    thr = parser.add_argument_group("auto-flagging thresholds")
    thr.add_argument("--mt-threshold", type=float, default=15.0,
                     help="Max median %%MT per cluster (default: 15).")
    thr.add_argument("--min-genes-threshold", type=float, default=400,
                     help="Min median genes per cluster (default: 400).")
    thr.add_argument("--min-cells", type=int, default=20,
                     help="Min cells per cluster (default: 20).")
    thr.add_argument("--single-batch-threshold", type=float, default=0.90,
                     help="Fraction above which a cluster is flagged as "
                          "single-batch (default: 0.90).")

    # ── Decisions ──
    parser.add_argument("--decisions", default=None,
                        help="Path to decisions.yaml to apply.")
    parser.add_argument("--report", action="store_true",
                        help="Generate report even when --decisions is set.")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # ══════════════════════════════════════════════════════════════════════
    #  Load
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("LOADING DATA")
    print(f"{'='*60}")
    adata = ad.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Resolve batch_key / covariate_keys from model or CLI ──
    model_info = None
    if args.model_dir:
        model_info = extract_keys_from_model(args.model_dir)
        if model_info is None:
            print(f"  [WARN] Could not read attr_dict.json from {args.model_dir}")
        else:
            print(f"  Read model config from {args.model_dir}")
    else:
        for search_dir in [Path(args.input).parent, Path(args.output_dir)]:
            for candidate in sorted(search_dir.glob("scvi_model_*")):
                info = extract_keys_from_model(candidate)
                if info:
                    model_info = info
                    print(f"  Auto-detected model: {candidate}")
                    break
            if model_info:
                break

    if args.batch_key is None:
        if model_info and model_info.get("batch_key"):
            args.batch_key = model_info["batch_key"]
            print(f"  Batch key from model: {args.batch_key}")
        else:
            args.batch_key = "data_origin"

    if args.covariate_keys is None:
        cat_covs = (model_info or {}).get("categorical_covariate_keys", [])
        args.covariate_keys = [k for k in cat_covs
                                if k and k != args.batch_key]
    if args.covariate_keys:
        print(f"  Covariate keys: {args.covariate_keys}")

    # ── Auto-detect cluster key ──
    cluster_key = args.cluster_key
    if cluster_key is None:
        for candidate in ["leiden"] + [c for c in adata.obs.columns
                                        if c.startswith("leiden_")]:
            if candidate in adata.obs.columns:
                cluster_key = candidate
                break
    if cluster_key is None or cluster_key not in adata.obs.columns:
        print(f"  [ERROR] No cluster column found. Use --cluster-key.")
        return
    print(f"  Cluster key: {cluster_key}")

    latent_key = args.latent_key
    if latent_key is None:
        for candidate in ["X_scVI"] + [k for k in adata.obsm
                                        if k.startswith("X_scVI_")]:
            if candidate in adata.obsm:
                latent_key = candidate
                break
    print(f"  Latent key: {latent_key or 'none (silhouette disabled)'}")

    umap_key = args.umap_key
    if umap_key is None:
        for candidate in ["X_umap"] + [k for k in adata.obsm
                                        if k.startswith("X_umap_")]:
            if candidate in adata.obsm:
                umap_key = candidate
                break
    print(f"  UMAP key: {umap_key or 'none'}")

    # Ensure QC metrics exist
    if "n_genes_by_counts" not in adata.obs.columns:
        print("  Computing QC metrics...")
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"],
                                   inplace=True, log1p=True)

    # ══════════════════════════════════════════════════════════════════════
    #  Report mode
    # ══════════════════════════════════════════════════════════════════════
    generate_report = (args.decisions is None) or args.report

    if generate_report:
        print(f"\n{'='*60}")
        print("GENERATING INSPECTION REPORT")
        print(f"{'='*60}")

        # 1. Load marker dict (needed early for fraction analysis)
        if args.markers:
            with open(args.markers) as f:
                marker_dict = json.load(f)
        else:
            marker_dict = DEFAULT_MARKERS

        # 2. Cluster QC summary
        print("  Computing cluster QC summary...")
        summary = cluster_qc_summary(
            adata, cluster_key=cluster_key, batch_key=args.batch_key,
            latent_key=latent_key, covariate_keys=args.covariate_keys)
        summary.to_csv(os.path.join(args.output_dir,
                                     "cluster_qc_summary.csv"))
        print(f"  Saved cluster_qc_summary.csv ({len(summary)} clusters)")
        print(summary.to_string())

        # 3. Marker fraction table + neuronal cluster identification
        print("\n  Computing marker expression fractions...")
        flat_markers = [g for genes in marker_dict.values() for g in genes]
        fraction_table = calculate_percent_expressed(
            adata, flat_markers, cluster_key,
            threshold=args.marker_threshold)

        neuronal_markers = args.neuronal_markers or DEFAULT_NEURONAL_MARKERS
        neuronal_clusters, neuronal_frac = identify_neuronal_clusters(
            fraction_table, neuronal_markers, args.neuronal_cutoff)

        if neuronal_clusters:
            print(f"  Likely neuronal clusters "
                  f"({len(neuronal_clusters)}, "
                  f"cutoff={args.neuronal_cutoff:.0%}): "
                  f"{neuronal_clusters}")
        else:
            print("  No clusters met neuronal marker criteria.")

        # 4. Auto-flag
        print("\n  Auto-flagging clusters...")
        thresholds = {
            "mt": args.mt_threshold,
            "min_genes": args.min_genes_threshold,
            "min_cells": args.min_cells,
            "single_batch": args.single_batch_threshold,
        }
        flags = auto_flag_clusters(
            summary,
            mt_threshold=args.mt_threshold,
            min_genes_threshold=args.min_genes_threshold,
            min_cells=args.min_cells,
            single_batch_threshold=args.single_batch_threshold,
            covariate_keys=args.covariate_keys,
        )

        if flags:
            for cl, reasons in flags.items():
                print(f"    Cluster {cl}:")
                for r in reasons:
                    print(f"      - {r}")
        else:
            print("    No clusters flagged.")

        write_auto_flags_yaml(
            flags, summary,
            os.path.join(args.output_dir, "auto_flags.yaml"),
            thresholds,
            neuronal_clusters=neuronal_clusters)
        print("  Saved auto_flags.yaml")

        # 5. Plots
        print("\n  Generating diagnostic plots...")
        images = {}

        images["UMAP Overview"] = plot_umap_overview(
            adata, cluster_key, args.batch_key, args.output_dir,
            umap_key=umap_key or "X_umap")

        images["QC Summary Heatmap"] = plot_qc_summary_heatmap(
            summary, flags, args.output_dir)

        if args.batch_key and args.batch_key in adata.obs.columns:
            images["Batch Composition"] = plot_batch_composition(
                adata, cluster_key, args.batch_key, args.output_dir)

        if args.covariate_keys:
            print(f"  Generating covariate composition plots "
                  f"({args.covariate_keys})...")
            cov_plots = plot_covariate_compositions(
                adata, cluster_key, args.covariate_keys, args.output_dir)
            for key, b64 in cov_plots.items():
                images[f"Covariate Composition ({key})"] = b64

        images["Marker Gene Expression"] = plot_marker_dotplot(
            adata, cluster_key, marker_dict, args.output_dir)

        if not fraction_table.empty:
            images["Marker Fraction Heatmap"] = plot_marker_fraction_heatmap(
                fraction_table, args.output_dir)
            images["Cluster Marker PCA"] = plot_cluster_marker_pca(
                fraction_table, args.output_dir,
                color_gene=args.pca_color_gene)

        images["Cluster Silhouettes"] = plot_cluster_silhouettes(
            summary, args.output_dir)

        # 6. HTML report
        print("  Assembling HTML report...")
        generate_html_report(
            summary, flags, images, adata,
            os.path.join(args.output_dir, "inspection_report.html"),
            cluster_key, args.batch_key,
            covariate_keys=args.covariate_keys,
            neuronal_clusters=neuronal_clusters,
            neuronal_fraction_table=neuronal_frac if not neuronal_frac.empty else None)
        print("  Saved inspection_report.html")

    # ══════════════════════════════════════════════════════════════════════
    #  Filter mode
    # ══════════════════════════════════════════════════════════════════════
    if args.decisions:
        print(f"\n{'='*60}")
        print("APPLYING DECISIONS")
        print(f"{'='*60}")
        decisions = load_decisions(args.decisions)
        print(f"  Loaded {args.decisions}")
        apply_decisions(adata, decisions, cluster_key=cluster_key,
                        output_dir=args.output_dir)

    print(f"\n{'='*60}")
    print(f"  Done. Outputs in {args.output_dir}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

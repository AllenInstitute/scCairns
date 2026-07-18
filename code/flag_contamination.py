#!/usr/bin/env python3
"""Marker-based contamination flagging via per-cell Z-scores.

Reproducible, cluster-independent way to catch cells expressing off-target
lineage markers (e.g. endothelial or hepatic genes in a sympathetic-neuron
dataset).  For each named panel, every panel gene is z-scored across cells
on the log-normalized expression matrix; a cell is flagged when it clears a
per-gene z threshold on at least ``min_genes`` panel members.  The per-cell
panel score returned alongside is the mean of those per-gene z-scores.

Typical use
-----------
    import scanpy as sc
    from flag_contamination import flag_contamination, DEFAULT_CONTAM_PANELS

    adata = sc.read_h5ad("integrated.h5ad")
    flags = flag_contamination(
        adata,
        panels=DEFAULT_CONTAM_PANELS,     # or your own dict
        z_thresh=2.0,
        min_genes=2,
        layer=None,                       # use adata.X (log-normalized)
    )
    # flags is a DataFrame indexed by adata.obs_names with columns:
    #   score_endothelial, n_hits_endothelial, flag_endothelial,
    #   score_hepatic,     n_hits_hepatic,     flag_hepatic,
    #   score_mesenchymal, n_hits_mesenchymal, flag_mesenchymal,
    #   flag_any_contam

    # Write out the flagged cells (matches the Seurat WhichCells CSV format)
    flags.query("flag_any_contam").index.to_series().to_csv(
        "../results/260708_contamination_zscore.csv", header=False, index=False)

Comparison with cluster-based selection
---------------------------------------
    from flag_contamination import compare_selections
    cluster_cells = pd.read_csv("260708_hepatic_endothelial_cells.csv",
                                header=None).squeeze().tolist()
    zscore_cells  = flags.query("flag_any_contam").index.tolist()
    compare_selections(cluster_cells, zscore_cells,
                       labels=("cluster5+6", "zscore"))
"""

from __future__ import annotations

import warnings
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ═══════════════════════════════════════════════════════════════════════════════
#  Default contamination panels (mouse; sympathetic-ganglion context)
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Rationale for choices: keep panels tight (3-5 genes), lineage-specific,
#  and expressed at high enough level in the contaminant to give a clear
#  z-score signal. Drop broadly-expressed "endothelial-ish" genes (Sparc,
#  Osmr, Piezo2) that also fire in fibroblasts, glia, or sensory neurons.

DEFAULT_CONTAM_PANELS: Dict[str, list] = {
    # Pan-vascular endothelium; Flt1+Kdr+Emcn is the textbook triad,
    # Cdh5/Pecam1 add specificity if present in the object.
    "endothelial": ["Flt1", "Kdr", "Emcn", "Cdh5", "Pecam1"],
    # Hepatocyte-restricted: Slco1a1/1a4 sinusoidal transporters + G6pc
    # gluconeogenesis + Mug1 murine plasma protein + Ugt2b1 phase-II
    # metabolism. Together these are essentially not expressed outside liver.
    "hepatic":     ["Slco1a1", "Slco1a4", "G6pc", "Mug1", "Ugt2b1"],
    # Mesenchymal / fibroblast (endoneurial fibroblasts or Schwann-cell
    # precursors — both neural-crest-derived, easy to co-purify with
    # sympathetic neurons on standard dissociation). Panel choice:
    #   Pdgfra  — canonical tissue-resident fibroblast TF-target
    #   Col1a2  — fibrillar collagen, fibroblast-specific (Col4* shared
    #             with endothelium so intentionally not used here)
    #   Ddr2    — collagen-binding RTK, mesenchymal-restricted
    #   Cfh     — complement factor H, secreted by stromal fibroblasts
    #   Tgfbr2  — TGF-β signaling dominant in fibroblasts/stroma
    # None of these are expressed at appreciable level in adult sympathetic
    # neurons; coincident z-score elevation on ≥2 flags a mesenchymal cell.
    "mesenchymal": ["Pdgfra", "Col1a2", "Ddr2", "Cfh", "Tgfbr2"],
}


# ═══════════════════════════════════════════════════════════════════════════════
#  Core scoring
# ═══════════════════════════════════════════════════════════════════════════════

def _get_expression(adata, gene: str, layer: Optional[str]) -> Optional[np.ndarray]:
    """Return a dense 1-D expression vector for `gene`, or None if absent."""
    if gene not in adata.var_names:
        return None
    j = adata.var_names.get_loc(gene)
    X = adata.X if layer is None else adata.layers[layer]
    col = X[:, j]
    if sp.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col.astype(np.float64, copy=False)


def _zscore(x: np.ndarray) -> np.ndarray:
    """Robustly z-score a vector; returns zeros if variance is 0."""
    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if sd == 0.0 or not np.isfinite(sd):
        return np.zeros_like(x)
    return (x - mu) / sd


def score_panel(adata, genes: Sequence[str], *,
                layer: Optional[str] = None,
                z_thresh: float = 2.0
                ) -> Tuple[np.ndarray, np.ndarray, list, list]:
    """Score one contamination panel.

    Parameters
    ----------
    adata : AnnData
        Expression matrix; ``adata.X`` (or ``adata.layers[layer]``) must be
        log-normalized (e.g. ``sc.pp.normalize_total`` + ``sc.pp.log1p``).
    genes : sequence of str
        Panel gene symbols. Missing genes are skipped with a warning.
    layer : str or None
        Layer to score. None uses ``adata.X``.
    z_thresh : float
        Per-gene z-score above which a gene counts as a "hit" in a cell.

    Returns
    -------
    score : ndarray, shape (n_cells,)
        Mean per-gene z-score across the panel's *present* genes.
    n_hits : ndarray, shape (n_cells,)
        Number of panel genes for which the cell's z >= z_thresh.
    used : list of str
        Panel genes actually present in ``adata.var_names``.
    missing : list of str
        Panel genes absent from ``adata.var_names``.
    """
    n = adata.n_obs
    used, missing = [], []
    z_stack, hit_stack = [], []
    for g in genes:
        x = _get_expression(adata, g, layer)
        if x is None:
            missing.append(g)
            continue
        used.append(g)
        z = _zscore(x)
        z_stack.append(z)
        hit_stack.append(z >= z_thresh)

    if not used:
        warnings.warn(f"No panel genes found in adata.var_names (asked for {list(genes)}).")
        return np.zeros(n), np.zeros(n, dtype=int), used, missing

    score = np.mean(np.stack(z_stack, axis=1), axis=1)
    n_hits = np.sum(np.stack(hit_stack, axis=1), axis=1).astype(int)
    return score, n_hits, used, missing


def flag_contamination(adata,
                       panels: Dict[str, Sequence[str]] = DEFAULT_CONTAM_PANELS,
                       *,
                       z_thresh: float = 2.0,
                       min_genes: int = 2,
                       layer: Optional[str] = None,
                       verbose: bool = True
                       ) -> pd.DataFrame:
    """Flag likely contaminants using per-panel Z-score co-expression.

    A cell is flagged for panel ``P`` when it has ``>= min_genes`` panel
    genes whose per-cell z-score is ``>= z_thresh``.  Requiring co-expression
    of multiple panel genes is what separates a real off-lineage cell from a
    single-gene outlier (dropout, sequencing noise, ambient RNA).

    Parameters
    ----------
    adata : AnnData
        Log-normalized expression matrix.  Cell IDs come from ``adata.obs_names``.
    panels : dict of {label: [gene, ...]}
        One entry per contamination lineage.
    z_thresh : float, default 2.0
        Per-gene z threshold (~top 2.5% one-tailed if expression were normal;
        for scRNA-seq counts the tail is heavier, so 2.0 typically catches
        a percent or two of cells per panel).
    min_genes : int, default 2
        Minimum number of panel genes clearing ``z_thresh`` required to flag.
    layer : str or None
        Which matrix to score. None = ``adata.X``.
    verbose : bool
        Print a per-panel summary.

    Returns
    -------
    pd.DataFrame indexed by ``adata.obs_names`` with, for each panel P:
        score_P, n_hits_P, flag_P
    plus a combined ``flag_any_contam`` column.
    """
    cols = {}
    for label, genes in panels.items():
        score, n_hits, used, missing = score_panel(
            adata, genes, layer=layer, z_thresh=z_thresh)
        flag = n_hits >= min_genes
        cols[f"score_{label}"]  = score
        cols[f"n_hits_{label}"] = n_hits
        cols[f"flag_{label}"]   = flag
        if verbose:
            n_flag = int(flag.sum())
            pct = 100.0 * n_flag / adata.n_obs
            miss_txt = f", missing: {missing}" if missing else ""
            print(f"[{label}] used {used}{miss_txt} — "
                  f"{n_flag} cells flagged ({pct:.2f}%, "
                  f"z>={z_thresh} on >={min_genes} genes)")

    df = pd.DataFrame(cols, index=adata.obs_names)
    flag_cols = [c for c in df.columns if c.startswith("flag_")]
    df["flag_any_contam"] = df[flag_cols].any(axis=1)
    if verbose:
        n_any = int(df["flag_any_contam"].sum())
        print(f"[combined] flag_any_contam = {n_any} "
              f"({100.0 * n_any / adata.n_obs:.2f}%)")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Comparison with cluster-based selection (regression test hook)
# ═══════════════════════════════════════════════════════════════════════════════

def compare_selections(a: Iterable[str], b: Iterable[str],
                       labels: Tuple[str, str] = ("A", "B")
                       ) -> Dict[str, object]:
    """Set-compare two cell-ID lists and print a small confusion summary.

    Useful as a regression check: does the z-score panel recover the cells
    you picked by cluster ID, and does it pick up any additional cells that
    the cluster call missed?
    """
    A, B = set(a), set(b)
    inter = A & B
    only_a = A - B
    only_b = B - A
    union = A | B
    jaccard = len(inter) / len(union) if union else float("nan")
    la, lb = labels
    print(f"|{la}|             = {len(A)}")
    print(f"|{lb}|             = {len(B)}")
    print(f"|{la} ∩ {lb}|      = {len(inter)}")
    print(f"|{la} \\ {lb}|      = {len(only_a)}   (in {la}, missed by {lb})")
    print(f"|{lb} \\ {la}|      = {len(only_b)}   (in {lb}, missed by {la})")
    print(f"Jaccard         = {jaccard:.3f}")
    if A:
        print(f"Recall of {la} by {lb} = {len(inter)/len(A):.3f}")
    if B:
        print(f"Recall of {lb} by {la} = {len(inter)/len(B):.3f}")
    return {"intersection": inter, "only_a": only_a, "only_b": only_b,
            "jaccard": jaccard}


# ═══════════════════════════════════════════════════════════════════════════════
#  Convenience: cluster contamination profile (which clusters do flags land in?)
# ═══════════════════════════════════════════════════════════════════════════════

def flag_by_cluster(adata, flag_df: pd.DataFrame,
                    cluster_key: str = "cl.res.0.6") -> pd.DataFrame:
    """Cross-tab flag rate by cluster — sanity check that flags concentrate.

    Returns one row per cluster with n_cells, flagged counts per panel,
    and the flagged fraction.  If z-score flags concentrate in the same
    clusters you'd have removed manually, the two approaches agree.
    """
    obs = adata.obs[[cluster_key]].join(flag_df)
    flag_cols = [c for c in flag_df.columns if c.startswith("flag_")]
    g = obs.groupby(cluster_key, observed=True)
    out = g.size().to_frame("n_cells")
    for c in flag_cols:
        out[c + "_n"]   = g[c].sum().astype(int)
        out[c + "_pct"] = 100.0 * g[c].mean()
    return out.sort_index(key=lambda idx: [int(x) if str(x).isdigit() else str(x)
                                           for x in idx])


# ═══════════════════════════════════════════════════════════════════════════════
#  CLI (optional)
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    import scanpy as sc

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", required=True, help="path to integrated h5ad")
    ap.add_argument("--output", required=True, help="output CSV of per-cell flags")
    ap.add_argument("--cluster-key", default="cl.res.0.6",
                    help="obs column for cluster-vs-flag cross-tab (optional)")
    ap.add_argument("--z-thresh", type=float, default=2.0)
    ap.add_argument("--min-genes", type=int, default=2)
    ap.add_argument("--layer", default=None)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    flags = flag_contamination(adata, z_thresh=args.z_thresh,
                               min_genes=args.min_genes, layer=args.layer)
    flags.to_csv(args.output)
    if args.cluster_key in adata.obs.columns:
        print("\nFlag rate by cluster:")
        print(flag_by_cluster(adata, flags, cluster_key=args.cluster_key))

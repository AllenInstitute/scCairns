#!/usr/bin/env python3
"""Trace whether specific genes survive the integration pipeline's stages.

Given a directory of pipeline outputs, this checks each ``.h5ad`` (and any
``hvg_genes.csv``) for a list of gene symbols and reports *where* each gene is
present or absent — so you can tell which stage removed a gene, and whether a
``FALSE`` from ``gene %in% rownames(x)`` is a genuine loss or just a name/case
mismatch.

Why this matters
----------------
A gene can be "missing" for several distinct reasons, each with a different fix:

  1. Name/case mismatch      -> the gene IS there under a different string
                                (e.g. var_names are Ensembl IDs like
                                'ENSMUSG...', with the symbol in a var column;
                                or 'TAC1' vs 'Tac1'). No filtering happened.
  2. QC gene filter          -> sc.pp.filter_genes(min_cells=...) dropped it
                                because it was detected in < min_cells cells.
                                It will be ABSENT from integrated.h5ad.
                                Fix: relax qc.min_cells.
  3. HVG selection           -> the gene is PRESENT in integrated.h5ad but not
                                in hvg_genes.csv, so scVI never trained on it.
                                Relaxing min_cells does NOT help; you need a
                                gene allowlist for HVG selection.
  4. Never in the raw input  -> absent even from the earliest h5ad.

For every gene that IS found, the script reports the number of cells in which
it is detected (X > 0). Compare that to your qc.min_cells (default 3): a gene
detected in fewer than min_cells cells is exactly what the QC filter removes.

Stages (in pipeline order), if present in the directory:
  <input>.h5ad            raw / pre-integration input (if you point at one)
  integrated.h5ad         full-gene object AFTER QC filtering (adata_full)
  filtered.h5ad           after decision-based CELL filtering (genes unchanged)
  hvg_genes.csv           the ~3000 genes scVI actually trained on

Usage
-----
    python check_gene_survival.py [DIR] \
        --genes Tac1 Vip Nts Cck \
        [--min-cells 3] [--recursive] [--counts-layer counts] \
        [--out gene_survival_report.csv]

DIR defaults to the current directory. With no --genes, the four defaults above
are used. Requires the pipeline conda env (anndata / scanpy / pandas).
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

try:
    import anndata as ad
except ImportError:  # pragma: no cover
    sys.exit("This script needs anndata. Run it in the pipeline conda env "
             "(the one with scanpy/anndata installed).")


# Canonical pipeline stage order, for sorting the report sensibly.
STAGE_ORDER = ["input", "raw", "adata", "integrated", "filtered", "other"]


def classify_stage(filename: str) -> str:
    """Map a filename to a rough pipeline stage label for ordering."""
    name = filename.lower()
    if "integrated" in name:
        return "integrated"
    if "filtered" in name:
        return "filtered"
    if "raw" in name or "unfiltered" in name or "input" in name:
        return "raw"
    return "other"


def find_h5ad_files(directory: str, recursive: bool) -> list[str]:
    pattern = "**/*.h5ad" if recursive else "*.h5ad"
    files = sorted(glob.glob(os.path.join(directory, pattern), recursive=recursive))
    # Order by pipeline stage, then by path, so the report reads in flow order.
    files.sort(key=lambda p: (STAGE_ORDER.index(classify_stage(os.path.basename(p))),
                              p))
    return files


def find_hvg_csvs(directory: str, recursive: bool) -> list[str]:
    pattern = "**/hvg_genes.csv" if recursive else "hvg_genes.csv"
    return sorted(glob.glob(os.path.join(directory, pattern), recursive=recursive))


def _lower_map(names) -> dict[str, str]:
    """Map lower-cased name -> first original name (for case-insensitive hits)."""
    out: dict[str, str] = {}
    for n in names:
        key = str(n).lower()
        if key not in out:
            out[key] = str(n)
    return out


def _scan_var_columns(var: pd.DataFrame, gene: str) -> list[str]:
    """Return 'column=value' hits where a var column holds the gene symbol.

    Handles the common case of Ensembl-ID var_names with the symbol living in a
    column such as 'gene_symbol', 'feature_name', 'gene_ids', 'symbol', etc.
    Matching is case-insensitive.
    """
    hits = []
    gl = gene.lower()
    for col in var.columns:
        try:
            series = var[col].astype(str).str.lower()
        except Exception:
            continue
        if (series == gl).any():
            hits.append(col)
    return hits


def detection_counts(adata, present_names: list[str], counts_layer: str | None):
    """Cells-detected count per gene name. Reads only the requested columns.

    Prefers a counts layer when available (matches what QC filtered on); else
    uses X. Returns {gene_name: n_cells_detected}.
    """
    if not present_names:
        return {}
    try:
        sub = adata[:, present_names]
        sub = sub.to_memory() if adata.isbacked else sub.copy()
    except Exception as exc:
        print(f"    [warn] could not slice genes for detection counts: {exc}")
        return {name: None for name in present_names}

    if counts_layer and counts_layer in sub.layers:
        mat = sub.layers[counts_layer]
    else:
        mat = sub.X

    # (X > 0) column sums, sparse- or dense-safe.
    try:
        from scipy.sparse import issparse
        if issparse(mat):
            detected = np.asarray((mat > 0).sum(axis=0)).ravel()
        else:
            detected = np.asarray((np.asarray(mat) > 0).sum(axis=0)).ravel()
    except Exception as exc:
        print(f"    [warn] detection-count computation failed: {exc}")
        return {name: None for name in present_names}

    return {name: int(c) for name, c in zip(sub.var_names, detected)}


def check_one_h5ad(path: str, genes: list[str], counts_layer: str | None):
    """Return a list of per-gene record dicts for one h5ad file."""
    print(f"\n{'='*70}\n{path}\n{'='*70}")
    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as exc:
        print(f"  [error] could not open (trying non-backed): {exc}")
        try:
            adata = ad.read_h5ad(path)
        except Exception as exc2:
            print(f"  [error] failed to read: {exc2}")
            return []

    n_obs, n_vars = adata.n_obs, adata.n_vars
    print(f"  {n_obs:,} cells x {n_vars:,} genes")
    print(f"  var_names sample: {list(map(str, adata.var_names[:3]))}")
    if len(adata.var.columns):
        print(f"  var columns: {list(adata.var.columns)}")

    var_names = list(adata.var_names)
    exact_set = set(map(str, var_names))
    lower_map = _lower_map(var_names)

    raw_lower_map = {}
    if getattr(adata, "raw", None) is not None:
        raw_lower_map = _lower_map(list(adata.raw.var_names))

    # First pass: resolve each gene to a present var_name (or None).
    resolved = {}
    for gene in genes:
        if gene in exact_set:
            resolved[gene] = gene
        elif gene.lower() in lower_map:
            resolved[gene] = lower_map[gene.lower()]  # case-insensitive hit
        else:
            resolved[gene] = None

    present_names = [v for v in resolved.values() if v is not None]
    det = detection_counts(adata, present_names, counts_layer)

    records = []
    for gene in genes:
        actual = resolved[gene]
        var_col_hits = _scan_var_columns(adata.var, gene) if actual is None else []
        in_raw = (actual is None) and (gene.lower() in raw_lower_map)

        if actual == gene:
            status = "present (exact)"
        elif actual is not None:
            status = f"present (case: '{actual}')"
        elif var_col_hits:
            status = f"present via var[{var_col_hits[0]}] (var_names are IDs)"
        elif in_raw:
            status = "present in .raw only"
        else:
            status = "ABSENT"

        records.append({
            "file": os.path.basename(path),
            "path": path,
            "stage": classify_stage(os.path.basename(path)),
            "gene": gene,
            "status": status,
            "matched_name": actual or (var_col_hits[0] if var_col_hits else ""),
            "n_cells_detected": det.get(actual) if actual else None,
        })

    # Console table for this file.
    print(f"\n  {'gene':<10} {'status':<42} {'cells_detected':>14}")
    print(f"  {'-'*10} {'-'*42} {'-'*14}")
    for r in records:
        nc = r["n_cells_detected"]
        nc_str = f"{nc:,}" if isinstance(nc, int) else "-"
        print(f"  {r['gene']:<10} {r['status']:<42} {nc_str:>14}")

    return records


def check_hvg_csvs(hvg_files: list[str], genes: list[str]) -> list[dict]:
    records = []
    for path in hvg_files:
        try:
            df = pd.read_csv(path)
        except Exception as exc:
            print(f"  [warn] could not read {path}: {exc}")
            continue
        col = "gene" if "gene" in df.columns else df.columns[0]
        hvg_set = set(df[col].astype(str))
        hvg_lower = {g.lower() for g in hvg_set}
        print(f"\n{'='*70}\n{path}  ({len(hvg_set):,} HVGs)\n{'='*70}")
        print(f"  {'gene':<10} {'in HVG set?'}")
        print(f"  {'-'*10} {'-'*20}")
        for gene in genes:
            if gene in hvg_set:
                inset = "yes (exact)"
            elif gene.lower() in hvg_lower:
                inset = "yes (case-insensitive)"
            else:
                inset = "NO"
            print(f"  {gene:<10} {inset}")
            records.append({
                "file": os.path.basename(path), "path": path, "stage": "hvg",
                "gene": gene, "status": f"HVG: {inset}",
                "matched_name": "", "n_cells_detected": None,
            })
    return records


def print_diagnosis(records: list[dict], genes: list[str], min_cells: int):
    print(f"\n{'#'*70}\n# DIAGNOSIS (min_cells reference = {min_cells})\n{'#'*70}")
    for gene in genes:
        gene_recs = [r for r in records if r["gene"] == gene]
        h5ad_recs = [r for r in gene_recs if r["stage"] != "hvg"]
        hvg_recs = [r for r in gene_recs if r["stage"] == "hvg"]

        present_files = [r for r in h5ad_recs if not r["status"].startswith("ABSENT")]
        absent_files = [r for r in h5ad_recs if r["status"].startswith("ABSENT")]

        print(f"\n  {gene}:")
        if not h5ad_recs:
            print("    no h5ad files scanned.")
            continue

        # Detection counts where known.
        det_vals = [r["n_cells_detected"] for r in present_files
                    if isinstance(r["n_cells_detected"], int)]

        if present_files and not absent_files:
            msg = "    present in ALL h5ad files scanned"
            if det_vals:
                lo = min(det_vals)
                msg += f" (min cells-detected = {lo:,})"
                if lo < min_cells:
                    msg += f"  <-- BELOW min_cells={min_cells}: QC WOULD drop it"
            print(msg)
        elif present_files and absent_files:
            print("    present in some files, ABSENT in others -> lost between stages:")
            for r in h5ad_recs:
                print(f"       {r['file']:<22} {r['status']}")
        else:
            print("    ABSENT from every h5ad scanned -> either removed by the QC")
            print("    gene filter on the raw input, or never in the raw data.")
            print("    Check the ORIGINAL pre-pipeline input to distinguish.")

        # HVG verdict — the most common real cause.
        for r in hvg_recs:
            verdict = r["status"]
            if "NO" in verdict and present_files:
                print(f"       {r['file']:<22} {verdict}  <-- present in h5ad but "
                      "NOT trained on (HVG drop). Relaxing min_cells will NOT help; "
                      "use a gene allowlist.")
            else:
                print(f"       {r['file']:<22} {verdict}")

    print(f"\n{'#'*70}")
    print("# How to read this:")
    print("#  - 'present (case:...)' or 'via var[...]'  -> NOT filtered; a name")
    print("#      mismatch. Your R rownames() check used the wrong string.")
    print("#  - ABSENT from integrated.h5ad             -> QC min_cells drop.")
    print("#      Fix: relax qc.min_cells (integrate.py:1380).")
    print("#  - present in integrated.h5ad but HVG=NO   -> HVG selection drop.")
    print("#      Fix: gene allowlist in select_hvgs; min_cells is irrelevant.")
    print(f"{'#'*70}")


def main():
    ap = argparse.ArgumentParser(
        description="Trace gene survival across pipeline h5ad stages.")
    ap.add_argument("directory", nargs="?", default=".",
                    help="Directory of pipeline outputs (default: cwd).")
    ap.add_argument("--genes", nargs="+",
                    default=["Tac1", "Vip", "Nts", "Cck"],
                    help="Gene symbols to check (default: Tac1 Vip Nts Cck).")
    ap.add_argument("--min-cells", type=int, default=3,
                    help="Reference qc.min_cells for the diagnosis (default: 3).")
    ap.add_argument("--counts-layer", default="counts",
                    help="Layer to use for detection counts if present "
                         "(default: counts; falls back to X).")
    ap.add_argument("--recursive", action="store_true",
                    help="Recurse into subdirectories.")
    ap.add_argument("--out", default=None,
                    help="Optional CSV path for the full report table.")
    args = ap.parse_args()

    h5ad_files = find_h5ad_files(args.directory, args.recursive)
    hvg_files = find_hvg_csvs(args.directory, args.recursive)

    if not h5ad_files:
        sys.exit(f"No .h5ad files found in {args.directory!r} "
                 f"({'recursive' if args.recursive else 'non-recursive'}).")

    print(f"Genes to check: {args.genes}")
    print(f"Found {len(h5ad_files)} h5ad file(s), {len(hvg_files)} hvg_genes.csv.")

    records: list[dict] = []
    for path in h5ad_files:
        records += check_one_h5ad(path, args.genes, args.counts_layer)
    records += check_hvg_csvs(hvg_files, args.genes)

    print_diagnosis(records, args.genes, args.min_cells)

    if args.out:
        pd.DataFrame(records).to_csv(args.out, index=False)
        print(f"\nWrote report table -> {args.out}")


if __name__ == "__main__":
    main()

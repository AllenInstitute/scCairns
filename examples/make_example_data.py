#!/usr/bin/env python3
"""Generate a small, synthetic multi-batch dataset for the tutorial.

This writes a tiny (~2,000-cell) AnnData with **raw integer counts** in ``.X``
so the integration pipeline accepts it directly (see docs/getting-started.md for
the input contract). It is deliberately synthetic — biologically fake but
structurally realistic — so the tutorial runs end-to-end on a laptop CPU in a
few minutes with no downloads.

What it builds in, so the loop has something to do:
  * A clear batch effect across three ``data_origin`` batches (two ``tech``
    platforms) for scVI/Harmony to correct.
  * Distinct cell populations with mouse marker programs (noradrenergic,
    cholinergic, satellite glia, immune) so the marker dotplot shows signal.
  * A low-quality population (few genes, high mitochondrial fraction) so the
    inspection auto-flagging fires on a cluster.
  * Small endothelial and hepatic "contaminant" populations expressing the
    default mouse contamination panels, so the per-cell contamination flagging
    and a cells-based filtering decision have real targets.

Usage
-----
    python examples/make_example_data.py \
        --output examples/example_data/combined_example.h5ad

Then point ``data.input_h5ad`` at the written file (examples/pipeline_tutorial.yml
already does).
"""

from __future__ import annotations

import argparse
import os

import numpy as np


# Marker programs mirror code/pipeline_config.py:SNS_MARKER_SETS and
# code/flag_contamination.py:DEFAULT_CONTAM_PANELS so the tutorial exercises the
# real defaults without any config overrides.
PROGRAMS = {
    "noradrenergic": ["Snap25", "Tubb3", "Rbfox3", "Th", "Dbh", "Ddc", "Slc6a2"],
    "cholinergic": ["Snap25", "Tubb3", "Rbfox3", "Chat", "Slc18a3", "Slc5a7"],
    "satellite_glia": ["Sox10", "Fabp7", "S100b"],
    "immune": ["Ptprc", "Cd68"],
    "endothelial": ["Flt1", "Kdr", "Emcn", "Cdh5", "Pecam1"],
    "hepatic": ["Slco1a1", "Slco1a4", "G6pc", "Mug1", "Ugt2b1"],
}

# Housekeeping / QC genes present in every cell type.
MT_GENES = ["mt-Co1", "mt-Co2", "mt-Cytb", "mt-Nd1", "mt-Atp6"]
RIBO_GENES = [f"Rps{i}" for i in range(2, 15)] + [f"Rpl{i}" for i in range(2, 15)]

# Population plan: (label, n_cells, is_low_quality). Contaminant populations are
# small on purpose so filtering removes a visible-but-minority slice.
POPULATIONS = [
    ("noradrenergic", 620, False),
    ("cholinergic", 560, False),
    ("satellite_glia", 340, False),
    ("immune", 220, False),
    ("endothelial", 90, False),
    ("hepatic", 70, False),
    ("lowq", 160, True),
]

BATCHES = [("batch1", "tech1"), ("batch2", "tech1"), ("batch3", "tech2")]


def build(seed: int = 0):
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(seed)

    # Assemble the gene universe: markers + QC genes + generic background.
    marker_genes = sorted({g for genes in PROGRAMS.values() for g in genes})
    background = [f"Gene{i:04d}" for i in range(1, 1201)]
    var_names = marker_genes + MT_GENES + RIBO_GENES + background
    # de-duplicate while preserving order
    seen = set()
    var_names = [g for g in var_names if not (g in seen or seen.add(g))]
    gene_index = {g: j for j, g in enumerate(var_names)}
    n_genes = len(var_names)

    rows = []
    obs_records = []
    for label, n_cells, is_lowq in POPULATIONS:
        program = [] if label == "lowq" else PROGRAMS[label]
        for _ in range(n_cells):
            batch, tech = BATCHES[rng.integers(len(BATCHES))]
            # Base per-gene rate (log-space mean), low background everywhere.
            rate = np.full(n_genes, 0.08)
            # Housekeeping + ribosomal always modestly on.
            for g in RIBO_GENES:
                rate[gene_index[g]] = 0.9
            # Marker program strongly on for this population.
            for g in program:
                rate[gene_index[g]] = 6.0
            # Mitochondrial load: high in the low-quality population so its
            # cluster clears the inspection auto-flag mt_threshold (default 15%)
            # while its genes stay above the QC min_genes cutoff (so the cells
            # survive QC and reach clustering to be flagged there).
            mt_level = 5.5 if is_lowq else 0.4
            for g in MT_GENES:
                rate[gene_index[g]] = mt_level
            # Modest library size for low-quality cells.
            lib = 700.0 if is_lowq else float(rng.integers(1800, 3200))
            # Batch effect: multiplicative per-batch gene shift so integration
            # has a real nuisance signal to remove.
            batch_factor = {"batch1": 1.0, "batch2": 1.35, "batch3": 0.7}[batch]
            gene_shift = np.exp(
                0.35 * rng.standard_normal(n_genes)
                * (1.0 if batch == "batch1" else 1.0)
            )
            if batch == "batch2":
                gene_shift *= np.exp(0.25 * rng.standard_normal(n_genes))
            probs = rate * batch_factor * gene_shift
            probs = probs / probs.sum()
            counts = rng.multinomial(int(lib), probs)
            rows.append(counts)
            obs_records.append((label, batch, tech))

    X = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(
        obs_records, columns=["true_celltype", "data_origin", "tech"]
    )
    obs.index = [f"cell_{i:05d}" for i in range(X.shape[0])]
    # Shuffle so populations/batches are interleaved (not blocked).
    order = rng.permutation(X.shape[0])
    X = X[order]
    obs = obs.iloc[order].copy()

    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = var_names
    adata.var_names_make_unique()
    return adata


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--output",
        default="examples/example_data/combined_example.h5ad",
        help="where to write the synthetic h5ad",
    )
    ap.add_argument("--seed", type=int, default=0, help="RNG seed")
    args = ap.parse_args()

    adata = build(seed=args.seed)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)
    print(
        f"Wrote {adata.n_obs} cells x {adata.n_vars} genes to {args.output}\n"
        f"  batches (data_origin): {sorted(adata.obs['data_origin'].unique())}\n"
        f"  platforms (tech):      {sorted(adata.obs['tech'].unique())}\n"
        f"  populations:           {sorted(adata.obs['true_celltype'].unique())}\n"
        "  .X holds raw integer counts; feed it straight to integrate_scvi.py."
    )


if __name__ == "__main__":
    main()

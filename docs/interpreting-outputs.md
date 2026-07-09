# Interpreting the outputs

You've run integrate + inspect. Now: **is the integration any good, what should I
remove, and when do I stop?** This page reads the artifacts for you.

## The inspection report (`inspection_report.html`)

A self-contained HTML you open in a browser. Read it top to bottom:

1. **UMAPs** colored by batch, covariate, and cluster. What you *want*: clusters that
   mix batches (integration worked) but stay biologically coherent (real cell types
   aren't smeared together). What's *suspect*: a cluster that is one batch only, or a
   diffuse smear with no marker identity.
2. **Marker dotplot / fraction heatmap** — per-cluster expression of your marker
   sets. Use it to name clusters and to spot clusters with no coherent identity
   (candidates for removal).
3. **Cluster QC** — the table described below, rendered with highlights.
4. **Contamination** — UMAP and marker boxplots of the per-cell flags.
5. **Run Provenance footer** — seed, code commit (with a `+uncommitted changes`
   marker if the tree was dirty), input SHA-256, package versions. This is what makes
   a given report reproducible.

## Cluster QC (`cluster_qc_summary.csv`)

One row per cluster. The columns that drive decisions:

| Column | What it tells you | Watch for |
|---|---|---|
| `n_cells` | Cluster size | Tiny clusters (< `min_cells`, default 20) are unstable |
| `median_genes` | Complexity | Low values (< `min_genes_threshold`, default 400) → low-quality/empty droplets |
| `median_pct_mt` | Mitochondrial fraction | High values (> `mt_threshold`, default 15%) → dying/stressed cells |
| `dominant_batch_frac` | Batch purity | > `single_batch_threshold` (default 0.90) → a batch artifact, not biology |
| `silhouette_latent` | Separation in the scVI latent space | Very negative (< `silhouette_threshold`, default −0.05) → poorly separated |
| `silhouette_umap` | Separation in UMAP | **Diagnostic only** — not used for flagging (see below) |

**Why flagging uses `silhouette_latent`, not `silhouette_umap`:** UMAP is a 2-D
visualization that can exaggerate or invent gaps. The latent space is what scVI
actually models, so it's the honest basis for a keep/remove call. UMAP compactness
stays a separate visual cue.

## Auto-flags (`auto_flags.yaml`)

The inspection writes a *suggested* set of removals — any cluster tripping one of the
thresholds above, with the reason recorded. Treat it as a starting point, not a
verdict: rename it to `decisions.yaml`, keep the flags you agree with, and add your
own. It carries an `integration:` block naming which embedding/clustering the flags
came from, so filtering later applies to the same variant. See
[decisions.md](decisions.md).

## Contamination flags

`flag_contamination` scores each cell against marker panels (default: mouse
endothelial and hepatic) and flags cells where ≥ `min_genes` panel genes clear a
per-gene z-score of `z_thresh` (default 2.0 on ≥2 genes). Two files:

- `contamination_zscore.csv` — every cell, its per-panel score, hit counts, and flags.
- `contamination_flagged_cells.csv` — just the flagged cell IDs, one per line. This is
  the file you reference from `decisions.yaml` via `remove_cells.cells_file`.

Plus `contamination_umap.png` (where flagged cells land) and
`contamination_marker_boxplots.png` (flagged vs background expression). This is a
**cluster-independent** signal — it catches off-lineage cells even when they don't
form their own cluster. Adapt the panels to your biology in
[adapting-to-your-data.md](adapting-to-your-data.md).

## Benchmark (`scib_benchmark_results.csv`)

When more than one embedding is present (scVI, optional Harmony, and PCA baseline),
[scib-metrics](https://scib-metrics.readthedocs.io/) scores each on batch-correction
metrics (and bio-conservation when a `benchmark.label_key` is set). Read it as a
**relative** comparison within one run — "did scVI beat the PCA baseline; how does
Harmony compare?" — not as an absolute pass/fail. A higher aggregate batch-correction
score with preserved biology is the goal. On real data, use it to choose which
embedding to cluster and filter on.

## Filtering retention (`filtering_retention_summary.csv`)

After a filter step, this accounts for exactly what left: pre/post cell counts
**overall, by cluster, by batch, and by cluster × batch**, with retention and removal
fractions. Use it to confirm a decision did what you intended — e.g. that removing a
"low-quality" cluster didn't also gut one batch you care about.

## When to stop iterating

There is no built-in stopping rule — it's your call, guided by:

- **`auto_flags.yaml` comes back empty** — no cluster trips a QC threshold.
- **Contamination flags collapse** — the flagged fraction drops to near zero (in the
  [tutorial](tutorial.md), ~160 cells → ~4 between rounds 1 and 2).
- **Cell count and markers stabilize** — each round removes fewer cells and the
  cluster/marker structure stops changing.
- **`summarize_rounds.py` shows no new signal** — the round table flattens out.

When successive rounds stop changing the biology, stop.

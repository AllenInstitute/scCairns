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
| `neighbor_purity` | Fraction of a typical cell's nearest neighbors sharing its cluster label | Below `min_neighbor_purity` (default 0.50) the cluster has no territory of its own — it is a split or a passenger, not a population |
| `neighbor_entropy` | Normalized (0–1) diversity of cluster labels in the neighborhood | Read *with* purity: low purity + low entropy = absorbed by one cluster; low purity + high entropy = diffuse across many. On its own it cannot detect an absorbed cluster, which scores 0 |
| `dominant_neighbor` / `_frac` | Which other cluster the neighborhood belongs to, and how much of the foreign neighbors it accounts for | A high fraction names the merge candidate |
| `<key>_depleted` / `_depletion_ratio` / `_depletion_q` | The platform, batch, or donor most missing from this cluster, as a fraction of what chance predicts, with an FDR-corrected p-value | A ratio near 0 at low `q` means the group is absent beyond chance — e.g. a cluster with no SSv4 cells |
| `<key>_top` / `_enrichment` / `_enrich_q` | The group most over-represented, as a multiple of its expected share | Enrichment ≫ 1 at low `q` on the donor key means one donor carries the cluster |
| `<key>_n_groups_present` | How many levels of that key appear at all | Necessary but not sufficient — with large donors every cluster contains every donor |

### Compositional bias: why a test, not a fraction

A cluster's composition alone cannot tell you whether it is anomalous. With 15% SSv4
overall, a **20-cell** cluster containing no SSv4 has p ≈ 0.04 — that happens by chance —
while a **200-cell** cluster containing none has p ≈ 2e-15. Both are "100% one platform"
and any `dominant_*_frac > 0.9` rule scores them identically. The hypergeometric test
asks the right question: *given this cluster's size, how surprising is this makeup?*

Flagging requires **both** a significance gate (`composition_q`) and an effect-size gate
(`min_enrichment`, `max_depletion_ratio`), because once clusters are large almost any
deviation becomes statistically significant. In testing, 400-cell clusters differing from
background by 1.1–1.3× reached q ≈ 1e-3 while being entirely unremarkable biologically.

**A significant result is not a verdict.** The null is that cells land in clusters
without regard to platform or donor. If one platform was FACS-sorted and the other was
not, every cluster outside the sorted population is genuinely depleted of it — correct
biology, not junk. The flag says "explain this", not "remove this", which is why its
message ends with a prompt to check the sampling design.
| `median_genes` | Complexity | Low values (< `min_genes_threshold`, default 400) → low-quality/empty droplets |
| `median_pct_mt` | Mitochondrial fraction | High values (> `mt_threshold`, default 15%) → dying/stressed cells |
| `dominant_batch_frac` | Batch purity | > `single_batch_threshold` (default 0.90) → a batch artifact, not biology |
| `silhouette_latent` | Separation in the scVI latent space | Very negative (< `silhouette_threshold`, default −0.05) → poorly separated |
| `silhouette_umap` | Separation in UMAP | **Diagnostic only** — not used for flagging (see below) |

**Why flagging uses `silhouette_latent`, not `silhouette_umap`:** UMAP is a 2-D
visualization that can exaggerate or invent gaps. The latent space is what scVI
actually models, so it's the honest basis for a keep/remove call. UMAP compactness
stays a separate visual cue.

## Doublet score by cluster (`doublet_score_by_cluster.png`)

Rendered only when a per-cell `doublet_score` is present in `obs`. Nothing in scCairns
computes it — the panel plots whatever an upstream step produced (DoubletFinder in R,
`sc.pp.scrublet`, anything else), so it appears when the column exists and is silently
skipped when it doesn't. See
[the doublet-detection assessment](../archive/DOUBLET_DETECTION_ASSESSMENT.md) for why
native scoring is not implemented.

Doublet callers model co-encapsulation in a droplet, so a plate-based arm (SSv4) is
normally left unscored. **Unscored cells are excluded, not treated as zero** — the title
reports how many of the total were scored and names the platforms that carry a score, and
each cluster's tick label shows its scored `n`. A cluster whose median sits well above the
rest is a doublet candidate; if an optional call column (`predicted_doublet`, or
DoubletFinder's `DF.classifications`) is present, the fraction called is overlaid on the
right axis.

Read it alongside `neighbor_purity`: a doublet-rich cluster sitting between two real
populations usually shows low purity with the two parents as its `dominant_neighbor`.

## Stress score by cluster (`stress_score_by_cluster.png`)

The heat-shock / dissociation signature, scored per cell as the mean z-score across the
HSP genes present in the object and summarized per cluster as `median_stress` in
`cluster_qc_summary.csv`. A cluster sitting well above the rest is usually an artifact
of how long the tissue sat in protease rather than a cell type — worth resolving before
it gets annotated, because a "novel Hspa1a-high population" is a recurring way to
publish a handling artifact.

**The zero line is not "unstressed."** The score is a mean of z-scores, so 0 is the
average cell in *this* object. If the whole dissociation went badly, every cluster sits
near 0 and the panel shows nothing — it is a relative measure, and a uniformly stressed
prep is invisible to it. The absolute check is still the HSP genes' expression in the
marker dotplot.

**Read it beside `median_pct_mt`, not instead of it.** The two capture different failure
modes — handling stress versus mitochondrial burden — and a cluster can be high in one
and normal in the other. That is also why the panel is heat-shock-only: `pct_counts_mt`
already exists, is already auto-flagged by `mt_threshold`, and a second mitochondrial
number on a different scale would only disagree with the first.

Unlike the contamination panels, a high stress score does **not** put cells into
`contamination_flagged_cells.csv`. Stress is normally a whole-cluster verdict
(`remove_clusters`), whereas contamination is a per-cell, wrong-lineage call — so the
two stay separate and `flag_any_contam` is unaffected. Immediate-early genes (Fos, Jun,
Egr1) are excluded by default because in neurons they are genuine activity markers; add
them through `inspection.stress.gene_prefixes` where that ambiguity doesn't apply.

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

## Cell counts in `round_manifest.json`

The `integration` block records the whole attrition chain, because the QC plots and the
manifest are snapshots of different moments:

| Field | Moment |
| --- | --- |
| `n_cells_input` | Cells in `data.input_h5ad`, before anything is dropped |
| `n_cells_after_obs_filter` | After `data.obs_filter` (`null` when no filter is set) |
| `n_cells_prefilter_plot` | What `qc_violins_prefilter.png` and `cells_per_<batch>_prefilter.png` show |
| `n_cells` | Final count written to `integrated.h5ad`, after QC cell filtering |

So `n_cells` is smaller than the `*_prefilter` plots by however many cells
`qc.min_genes` removed, and smaller than `n_cells_input` by that plus any `obs_filter`
selection. `qc_filter_applied` tells you whether the QC filter ran at all
(`qc.enabled: false` or `qc.skip_filter: true` makes the pre- and post-filter counts
identical).

`summarize_rounds` reads this chain: the round table gains an **In cells** column
(`n_cells_input`) beside the final **Cells** count, retention percentages are measured
against the first round's input rather than its post-QC count, and `--verify` emits an
info-level `pre_integration_drop` finding naming what did the dropping. Rounds
integrated before these fields existed fall back to `n_cells` and report no drop.
A `warn`-level `input_cells_mismatch` fires when a round's integration read a different
number of cells than its own decision stage wrote — i.e. it was pointed at something
other than that round's `filtered.h5ad`.

## When to stop iterating

There is no built-in stopping rule — it's your call, guided by:

- **`auto_flags.yaml` comes back empty** — no cluster trips a QC threshold.
- **Contamination flags collapse** — the flagged fraction drops to near zero (in the
  [tutorial](tutorial.md), ~160 cells → ~4 between rounds 1 and 2).
- **Cell count and markers stabilize** — each round removes fewer cells and the
  cluster/marker structure stops changing.
- **`summarize_rounds.py` shows no new signal** — the round table flattens out.

When successive rounds stop changing the biology, stop.

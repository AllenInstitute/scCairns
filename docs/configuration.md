# Configuration reference

The YAML config is the source of truth for a run; CLI flags override individual
fields. Configs are **versioned** (`pipeline_version: 1`) and validated after
overrides — an invalid value fails fast with a clear error. The fully merged config
(defaults + your file + CLI) is written to `run_config_resolved.yml` in every output
directory, so you can always see exactly what ran.

This page lists each field with its default and, where it matters, **when to change
it**. For the raw templates see `examples/pipeline_scvi.yml` (mouse/SNS),
`examples/pipeline_generic.yml` (non-SNS), `examples/pipeline_scanvi.yml`, and
`examples/pipeline_sweep.yml`.

## `reproducibility`

| Field | Default | Notes |
|---|---|---|
| `seed` | `0` | Applied via `scvi.settings.seed` (seeds Python/NumPy/Torch). Makes scVI/scANVI, neighbors, UMAP, and **Leiden cluster IDs** deterministic across re-runs of the same config — which is what lets `decisions.yaml` reference cluster IDs safely. Set `null` for an intentionally unseeded run. |

## `data`

| Field | Default | When to change |
|---|---|---|
| `input_h5ad` | `null` | Required. Path to the combined counts h5ad (resolved relative to the config file). |
| `output_dir` | `results` | Where this round writes. Give each round its own directory. A relative path here is resolved against **the config file's location**, so a config works from any working directory; the bare default applies only when no config is given, and then it is relative to the current directory. |
| `obs_filter` | `null` | A pandas query applied to `adata.obs` right after load, e.g. `"condition in ['Control', 'Saline']"`. Integrates a subset without writing a pre-filtered h5ad. Recorded in `round_manifest.json`, and the cells it drops are accounted for by `n_cells_after_obs_filter` — see [interpreting outputs](interpreting-outputs.md#cell-counts-in-round_manifestjson). |
| `sample_key` | `null` | Donor / biological-replicate column (e.g. `donor_id`). **Not** used for integration — it is the grouping the per-cluster donor-bias test runs over, so a cluster carried by a single donor can be told from one present across all of them. Set it if you want that test. |
| `counts_layer` | `counts` | Only if your raw counts live in a differently named layer. |
| `batch_key` | `data_origin` | The `obs` column scVI integrates over. **Set this to your batch/sample variable.** |
| `categorical_covariate_keys` | `[tech]` | Extra nuisance factors to condition on (platform, chemistry). Set `[]` if none. |
| `continuous_covariate_keys` | `[]` | Continuous nuisance factors (e.g. percent-mito) if you want scVI to condition on them. |
| `gene_symbol_case` | `preserve` | `upper`/`lower` normalizes `var_names` before QC and marker matching. See [adapting to your data](adapting-to-your-data.md). |
| `species` | `mouse` | **Informational label only** — stored in the config but not read by the pipeline. QC is driven entirely by the `qc.*_gene_patterns` regex, so set *those* for your organism (setting `species` alone changes nothing). |

## `qc`

| Field | Default | When to change |
|---|---|---|
| `enabled` | `true` | Rarely off. |
| `min_genes` | `500` | Lower for shallow/nuclei data; raise to be stricter. Cells below this are dropped. |
| `min_cells` | `1` | Genes detected in fewer cells are dropped. Default `1` keeps every observed gene (only all-zero columns go), so rare markers (e.g. neuropeptides) survive; raise to be stricter. |
| `mt/ribo/hb_gene_patterns` | mouse regex | **Match your organism** (`^MT-` etc. for human). QC metrics are computed from these. |
| `skip_filter_after_round_1` | `true` | Runs the QC gene/cell filter on the **first round only**. Later rounds (input sits beside a parent `round_manifest.json`) skip it, so genes are never re-filtered against the smaller cell subsets later rounds operate on — which would otherwise drop rare genes the first round kept. In later rounds you filter by decision, not by QC. |

## `integration`

Only `model_type: scvi` is supported in V1. scANVI is an *annotation stage on top*
(see `annotation`), not a separate `model_type`. To tune the architecture, compare
architectures with a sweep, or swap in a different integration method (Scanorama,
etc.), see [Changing the model architecture or integration method](integration-methods.md).

| Field | Default | When to change |
|---|---|---|
| `n_hidden` | `256` | Width of the encoder/decoder. Larger = more capacity; raise for very large/complex datasets, lower (64–128) for small ones or speed. |
| `n_layers` | `3` | Depth. 1–2 for small data, 2–3 for typical, more rarely helps. |
| `n_latent` | `32` | Latent dimensionality. 10–15 for simple data, 20–32 for heterogeneous atlases. Too high can overfit batch. |
| `dispersion` | `gene-cell` | `gene` is the simplest/most stable; `gene-cell` and `gene-batch` add flexibility at some cost. |
| `gene_likelihood` | `nb` | `nb` (negative binomial) is the safe default; `zinb` if you have strong zero-inflation; `poisson` for UMI data with low overdispersion. |
| `max_epochs` | `200` | With early stopping, this is a ceiling. Lower for quick experiments; the default is fine for real data. |
| `early_stopping_patience` | `20` | Epochs without validation improvement before stopping. |
| `batch_size`, `num_workers` | `256`, `4` | Tune for your hardware; `num_workers: 0` avoids dataloader issues on some machines. |
| `hvg.n_top_genes` | `3000` | Highly-variable genes used for modeling. 2000–5000 is typical; more for atlases. |
| `hvg.batch_key` | `tech` | Stratify HVG selection so genes variable *within* each batch are chosen. Set `null` to pool. |
| `hvg.flavor` | `seurat_v3` | `seurat_v3` operates on raw counts (matches the counts input). |
| `hvg.min_batches` | `null` | Optionally require a gene be HV in ≥ N batches. |
| `sweep` | `[]` | A list of named architecture variants to train in one run — see [Sweeps](#sweeps). |

### Harmony (optional)

Set `integration.harmony.enabled: true` to also run
[Harmony](https://github.com/slowkow/harmonypy) on a log-normalized PCA of the HVGs,
as a same-run baseline that flows through the same UMAP/Leiden/benchmark stages.
Requires the `harmonypy` package. `batch_key` defaults to `data.batch_key`; `n_pcs`
(default 30) is the PCA dimensionality it corrects. In sweep mode, `hvg_from` names
the sweep entry whose HVGs Harmony should borrow (it has no HVG set of its own).

## `annotation` (scANVI)

Off by default. To transfer labels, set `enabled: true`, `method: scanvi`, and point
`labels_key` at an `obs` column of seed labels; unlabeled cells use
`unlabeled_category` (default `Unknown`). Outputs go to `obsm["X_scANVI"]`,
`obs["scanvi_label"]`, and `obs["scanvi_confidence"]`, and you can filter on
confidence in `decisions.yaml`. **scANVI is single-model only** — it can't be combined
with `sweep` in V1 (choosing which swept model seeds annotation needs an explicit
design).

## `embedding`

| Field | Default | When to change |
|---|---|---|
| `n_neighbors` | `30` | Smaller = more local structure; larger = smoother. |
| `umap_min_dist`, `umap_spread` | `0.4`, `3.0` | Visualization only; don't affect clustering or filtering. |
| `leiden_resolution` | `0.3` | **The main clustering knob.** Higher = more, smaller clusters (finer filtering granularity); lower = fewer, broader clusters. |

## `inspection`

`cluster_key`/`latent_key`/`umap_key` default to `auto` (resolved from what's in the
object). The `auto_flag` block sets the thresholds documented in
[interpreting-outputs.md](interpreting-outputs.md#cluster-qc-cluster_qc_summarycsv):
`mt_threshold` 15.0, `min_genes_threshold` 400, `min_cells` 20,
`single_batch_threshold` 0.90, `silhouette_threshold` −0.05,
`min_neighbor_purity` 0.50, `entropy_threshold` null, `entropy_neighbors` 15.

`min_neighbor_purity` flags a cluster with no territory of its own in the latent
space: fewer than this fraction of a typical cell's nearest neighbors carry its own
cluster label. The flag names the cluster the neighborhood *does* belong to when one
dominates, which is the difference between "merge this into cluster 7" and "this is
diffuse, drop it". Set `null` to disable.

The **compositional bias** keys — `composition_q` 0.01, `composition_min_cells` 20,
`min_enrichment` 2.0, `max_depletion_ratio` 0.20, `composition_keys` null — test whether
each cluster's platform/batch/donor makeup is what chance would give, by hypergeometric
test with Benjamini–Hochberg correction. `composition_keys: null` tests `data.batch_key`,
`data.sample_key`, and the categorical covariates; give an explicit list to override.
Flagging needs **both** gates to trip: `composition_q` for significance and
`min_enrichment` / `max_depletion_ratio` for effect size, because once clusters are large
almost any deviation is significant. Clusters below `composition_min_cells` are not
tested and are excluded from the FDR correction rather than diluting it. Set
`composition_q: null` to disable. Full rationale in
[interpreting outputs](interpreting-outputs.md#compositional-bias-why-a-test-not-a-fraction).

`entropy_threshold` is a second, opt-in test on neighborhood label entropy, off by
default because entropy measures how *diverse* a neighborhood is rather than how much
of it agrees. A cluster wholly absorbed into one other cluster has a homogeneous
neighborhood — of someone else's label — and so scores near-zero entropy, the same as a
perfectly isolated cluster. Purity separates those two cases; entropy only describes
which kind of mixing is happening once purity has flagged it. The per-cluster values
(`neighbor_purity`, `neighbor_entropy`, `dominant_neighbor`, `dominant_neighbor_frac`)
are written to `cluster_qc_summary.csv` either way. The `contamination` block
(`z_thresh` 2.0, `min_genes` 2, `layer` null, `panels` null→built-in) controls per-cell
flagging. Marker sets and panels are covered in
[adapting to your data](adapting-to-your-data.md).

The **`stress`** block (`enabled` true, `gene_prefixes` null→built-in, `genes` null,
`z_thresh` 2.0, `layer` null, `score_key` `stress_score`) scores the heat-shock /
dissociation signature per cell, writing `median_stress` to `cluster_qc_summary.csv`
and `stress_score_by_cluster.png` to the report. The panel is matched by **prefix**
against `var_names` at run time (default `Hsp`, `HSP`, `Dnaj`, `DNAJ`), so one default
covers mouse and human symbol casing; set `genes` to pin an explicit list instead.

Two deliberate omissions. **Mitochondrial genes are not in it** — `pct_counts_mt` /
`median_pct_mt` / `auto_flag.mt_threshold` already cover that, and a second
mitochondrial number on a different scale would just disagree with the first; read the
stress score *beside* `median_pct_mt`, not instead of it. **Immediate-early genes are
not in it** either: in neurons Fos/Jun/Egr1 are real activity markers, so add them via
`gene_prefixes` only where the tissue makes the interpretation unambiguous.

Unlike contamination, a high stress score does **not** feed
`contamination_flagged_cells.csv` — it is a "drop this cluster" signal, not a
wrong-lineage call, so it stays out of `flag_any_contam`.

## `benchmark`

`enabled: true` runs scib-metrics across the available embeddings. Set
`benchmark.label_key` to an `obs` column to add bio-conservation metrics; leave
`null` for batch-correction metrics only. `batch_key: auto` uses `data.batch_key`.

## `decisions`

`ignore_failed_queries` (default `false`): a `decisions.yaml` query that references a
missing column or errors is a **hard failure** by default, so a typo can't silently
skip a filter. Set `true` to restore permissive behavior.

## `record-filter` (no config block)

`cairns record-filter` is deliberately **config-free** — it takes no `--config` and
reads nothing from `pipeline.yml`. It records a filtering step that already happened
somewhere else, so there is no shared state to inherit: everything it writes comes from
the two objects you point it at plus the flags below. See
[filtering outside scCairns](decisions.md#filtering-outside-sccairns-seurat-loupe-a-notebook)
for the workflow.

| Flag | Required | Notes |
|---|---|---|
| `--input-h5ad` / `--input-cell-ids` | one of the two | The cells **before** filtering. An h5ad also supplies the gene count, the cluster column for a per-cluster breakdown, and a parent-round pointer; an ID list gives only the IDs. |
| `--output-h5ad` / `--output-cell-ids` | one of the two | The cells that **survived**. An ID list is enough — a CSV column of barcodes, a bare one-per-line file, or R's `write.csv` output. |
| `--output-dir` | yes | Where `round_manifest.json` and the sidecars go. Make this the directory holding the filtered h5ad, so the next round's `data.input_h5ad` sits beside the manifest and the lineage links up. |
| `--filter-used` | yes | Human-readable description of what was removed. This is the only record of *why*, so make it specific — and make it match what actually ran. |
| `--embedding-key` | no | The embedding the filtering was done on. Sets `filtered_on.variant` and infers the cluster/UMAP keys. **Omit at ingest**, before any embedding or clustering exists. |
| `--cluster-key` | no | Overrides the cluster key inferred from `--embedding-key`, for the per-cluster removal breakdown. |
| `--umap-key` | no | Overrides the UMAP key inferred from `--embedding-key`; recorded for provenance only. |
| `--notes` | no | Free text. A good place for the path to an external QC report. |
| `--seed` | no | Stamped into the manifest's `reproducibility` block, if the external step used one. |
| `--batch-key` | no | Recorded for provenance; not used to compute anything. |

Because it is config-free, the flags are the whole interface — which means they belong
in your run script next to the command that produced the filtered object, not in a
config file that a later round might silently reuse.

## Sweeps

`integration.sweep` is a list of architecture variants trained in one round, each
producing its own `X_scVI_<name>`, `X_umap_<name>`, and `leiden_<name>`:

```yaml
integration:
  sweep:
    - name: small
      n_hidden: 64
      n_layers: 1
      n_latent: 10
    - name: large
      n_hidden: 256
      n_layers: 3
      n_latent: 32
```

Because each variant has its own cluster IDs, a filtering decision must name the
variant it applies to (via the `decisions.yaml` `integration:` block or `--cluster-key`).
Inspect all variants at once with `--all-sweep-architectures`. See
[decisions.md](decisions.md#sweeps-pinning-the-variant).

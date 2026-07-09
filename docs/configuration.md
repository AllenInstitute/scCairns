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
| `output_dir` | `../results` | Where this round writes. Give each round its own directory. |
| `counts_layer` | `counts` | Only if your raw counts live in a differently named layer. |
| `batch_key` | `data_origin` | The `obs` column scVI integrates over. **Set this to your batch/sample variable.** |
| `categorical_covariate_keys` | `[tech]` | Extra nuisance factors to condition on (platform, chemistry). Set `[]` if none. |
| `continuous_covariate_keys` | `[]` | Continuous nuisance factors (e.g. percent-mito) if you want scVI to condition on them. |
| `species`, `gene_symbol_case` | `mouse`, `preserve` | See [adapting to your data](adapting-to-your-data.md). `gene_symbol_case: upper/lower` normalizes symbols. |

## `qc`

| Field | Default | When to change |
|---|---|---|
| `enabled` | `true` | Rarely off. |
| `min_genes` | `500` | Lower for shallow/nuclei data; raise to be stricter. Cells below this are dropped. |
| `min_cells` | `3` | Genes in fewer cells are dropped. |
| `mt/ribo/hb_gene_patterns` | mouse regex | **Match your organism** (`^MT-` etc. for human). QC metrics are computed from these. |
| `skip_filter_after_round_1` | `true` | Keeps round-1 QC from re-dropping cells in later rounds (you're filtering by decision then, not by QC). |

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
`single_batch_threshold` 0.90, `silhouette_threshold` −0.05. The `contamination` block
(`z_thresh` 2.0, `min_genes` 2, `layer` null, `panels` null→built-in) controls per-cell
flagging. Marker sets and panels are covered in
[adapting to your data](adapting-to-your-data.md).

## `benchmark`

`enabled: true` runs scib-metrics across the available embeddings. Set
`benchmark.label_key` to an `obs` column to add bio-conservation metrics; leave
`null` for batch-correction metrics only. `batch_key: auto` uses `data.batch_key`.

## `decisions`

`ignore_failed_queries` (default `false`): a `decisions.yaml` query that references a
missing column or errors is a **hard failure** by default, so a typo can't silently
skip a filter. Set `true` to restore permissive behavior.

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

# Troubleshooting

Common errors and confusing behaviors, with the fix.

### `[WARN] adata.X looks like raw counts`

Contamination scoring (and some inspection steps) expect **log-normalized**
expression in `.X`. Integration writes `.X` normalized/log1p'd, so this warning
usually means you pointed inspection at a raw-counts file instead of an
`integrated.h5ad`. Run inspection on the integration output, or set
`inspection.contamination.layer` to a log-normalized layer.

### `[WARN] Could not read ... model` / `scvi_model_*` not found

Inspection auto-discovers a trained `scvi_model_*` directory (beside the input h5ad
or in the output directory) **only to infer `batch_key`/covariate keys when you
haven't supplied them**. It is *not* required to generate a report or to filter — both
work from the embeddings, clusters, and `obs` already in the h5ad. So a missing model
is at most a warning, and it's irrelevant if you pass `--config` (or `--batch-key`),
or if you brought your own embedding (see
[integration-methods.md](integration-methods.md#level-b--bring-your-own-embedding-no-changes-to-this-pipeline)).
To silence it, pass `--model-dir`, or supply the keys via config/CLI.

### `Unsupported pipeline_version`

Your config is missing `pipeline_version: 1` at the top level, or has a different
value. Add it. The version gates schema validation.

### A `decisions.yaml` query did nothing / errored the run

By default a query referencing a missing obs column is a **hard error** (so typos
can't silently skip a filter). Check the column exists in `adata.obs` (open the object
or read `cluster_qc_summary.csv` for available fields). To make queries permissive,
set `decisions.ignore_failed_queries: true`. See [decisions.md](decisions.md).

### Sweep filter: "multiple variants present, nothing pinned"

In a sweep each architecture has its own `leiden_<variant>` clusters, so filtering
refuses to guess which one your cluster IDs mean. Pin it with `--cluster-key
leiden_<variant>` or an `integration:` block in the decisions file. See
[Sweeps](decisions.md#sweeps-pinning-the-variant).

### scANVI with a sweep is rejected

scANVI annotation is single-model only in V1 (`annotation.enabled: true` cannot be
combined with a non-empty `integration.sweep`). Run the sweep and annotation as
separate configs.

### `ModuleNotFoundError: No module named 'harmonypy'`

You enabled `integration.harmony.enabled: true` without the optional `harmonypy`
package. Either `pip install harmonypy` (it's in `requirements.txt`) or set
`harmony.enabled: false`.

### QC dropped more cells than expected

`qc.min_genes` (default 500) removes low-complexity cells before clustering, so a
population of shallow/nuclei cells can disappear in round 1. Lower `min_genes` if
that's not what you want. Note `skip_filter_after_round_1: true` means QC filtering
only applies in round 1 — later rounds are shaped by your `decisions.yaml`, not QC.

### A specific gene is missing from `integrated.h5ad`

A gene can be "missing" for several distinct reasons, and only one is a QC issue:

- **Name/case mismatch.** If `gene %in% rownames(x)` returns `FALSE` for several
  genes at once, first suspect the gene-name convention. `var_names` may be
  Ensembl IDs (symbol in a `var` column like `gene_symbol`) or a different case
  (`TAC1` vs `Tac1`). The gene is present — your check used the wrong string.
- **QC gene filter.** `sc.pp.filter_genes(min_cells=...)` drops genes detected in
  fewer than `qc.min_cells` cells. The default is now **1** (keep every observed
  gene); if you raised it, rare markers can be removed. This filter also used to
  re-run every round against a *shrinking* cell subset, so a gene kept at round 1
  (many cells) could fall below the threshold at round 3 (few cells) — that is
  now prevented by `skip_filter_after_round_1`.
- **HVG selection.** scVI trains only on the top HVGs. A gene can be *present* in
  `integrated.h5ad` but not among the HVGs, so it doesn't drive the embedding.
  Relaxing `min_cells` does **not** help here.

Diagnose which one with `code/check_gene_survival.py`, pointed at a round
directory — it traces each gene across `filtered.h5ad` / `integrated.h5ad` /
`hvg_genes.csv`, reports cells-detected counts, and names the responsible stage:

```bash
python code/check_gene_survival.py path/to/rounds --recursive \
    --genes Tac1 Vip Nts Cck --out gene_survival_report.csv
```

### Cluster IDs changed between runs

Cluster IDs are deterministic **only when the seed is fixed** (`reproducibility.seed`,
default 0) *and* the config and input are unchanged. If you set `seed: null`, changed
the architecture, or filtered the input, cluster numbering can differ — so re-check
which cluster your report flagged before copying an old `decisions.yaml`.

### Benchmark table plot warning (`plot_results_table failed` / `FigureCanvasAgg`)

A cosmetic incompatibility between some `scib-metrics`/`plottable`/matplotlib
versions. The numeric results (`scib_benchmark_results.csv`) are still written
correctly; only the rendered table image is affected.

### Running tests

```bash
python -m pytest -q
```

If tests fail on import, your environment is missing a pinned dependency — re-run
`pip install -r requirements.txt` inside the activated `sccairns` env.

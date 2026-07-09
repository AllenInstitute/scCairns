# Troubleshooting

Common errors and confusing behaviors, with the fix.

### `[WARN] adata.X looks like raw counts`

Contamination scoring (and some inspection steps) expect **log-normalized**
expression in `.X`. Integration writes `.X` normalized/log1p'd, so this warning
usually means you pointed inspection at a raw-counts file instead of an
`integrated.h5ad`. Run inspection on the integration output, or set
`inspection.contamination.layer` to a log-normalized layer.

### Filter mode can't find the scVI model / `scvi_model_*` missing

When filtering, inspection looks for the trained model directory beside the input
h5ad and in the output directory. If you moved or renamed files, the auto-discovery
misses it. Keep each round's `integrated.h5ad` next to its `scvi_model_*` directory,
or pass `--model-dir` explicitly.

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
`pip install -r requirements.txt` inside the activated `scvi-loops` env.

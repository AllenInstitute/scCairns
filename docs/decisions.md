# Filtering with `decisions.yaml`

Between rounds you decide which cells and clusters to carry forward. That decision is
a small YAML file — reproducible, self-documenting, and recorded in the round
manifest. The inspection step generates an `auto_flags.yaml` you rename and edit into
`decisions.yaml`.

## The four actions

```yaml
keep_clusters:      # optional: restrict to these clusters (drop all others)
  - cluster: "1"
    reason: "High-quality target population"

keep_cells:         # optional: restrict to cells matching each query (AND-chained)
  - query: "n_genes_by_counts > 500"
    reason: "Minimum gene complexity"

remove_clusters:    # optional: drop these clusters from what's retained
  - cluster: "7"
    reason: "Low-quality, high-MT cluster"

remove_cells:       # optional: drop cells matching a query or listed in a file
  - query: "scanvi_confidence < 0.5"
    reason: "Low-confidence annotation"
  - cells_file: "contamination_flagged_cells.csv"
    reason: "Marker-based contamination flags"

notes: "Round 2 review."   # free text, kept for provenance
```

**Order of application:** `keep_clusters` → `keep_cells` → `remove_clusters` →
`remove_cells`. Keeps narrow the set first; removes subtract from what remains. Every
action needs a `reason` — it's recorded in `decisions_applied.yaml` so the *why* of a
filter survives with the data.

## Two ways to name cells

- **`query`** — a pandas expression evaluated against `adata.obs`. Any obs column
  works: QC metrics (`n_genes_by_counts`, `pct_counts_mt`), annotation
  (`scanvi_confidence`), or your own columns. Combine with `&`, `|`, and parentheses.
- **`cells_file`** — a path to a CSV of `obs_names`, **one per line, no header** (the
  format `contamination_flagged_cells.csv` is written in). Paths are resolved relative
  to the `decisions.yaml` file.

## Failed queries fail loud

By default a query that references a missing column or otherwise errors **stops the
run** — a typo can't silently no-op a filter you intended. To opt into permissive
behavior (skip failed queries with a warning), set in your config:

```yaml
decisions:
  ignore_failed_queries: true
```

## Running a filter

```bash
python code/inspect_integration.py \
  --config pipeline.yml \
  --input  rounds/round_01/integrated.h5ad \
  --decisions rounds/round_01/decisions.yaml \
  --output-dir rounds/round_02
```

Writes `rounds/round_02/filtered.h5ad`, `decisions_applied.yaml` (what actually ran,
with counts removed per action), and `filtering_retention_summary.csv` (pre/post
counts overall, by cluster, by batch, and by cluster × batch). Feed `filtered.h5ad`
into the next integration round.

## Sweeps: pinning the variant

In a [sweep](configuration.md#sweeps), each architecture has its own Leiden cluster
IDs (`leiden_<variant>`), so a cluster-based decision is only meaningful against the
variant its IDs came from. The `auto_flags.yaml` carries an `integration:` block
naming that variant; it travels with the file when you rename it:

```yaml
integration:
  variant: "small"
  cluster_key: "leiden_small"
```

Selection precedence is `--cluster-key` (CLI) > the decisions file's `integration:`
block > `inspection.cluster_key` in the config > `auto`. In filter mode, if multiple
variants are present and nothing pins one, the run **errors** rather than guessing.
The variant actually used is recorded in `round_manifest.json` under
`decisions.filtered_on`.

# Tutorial: one full integration loop, end to end

This walks through a complete round of the workflow —
**integrate → inspect → decide → filter → re-integrate → summarize** — on a small
synthetic dataset you generate locally. It runs on a laptop CPU in a few minutes
and needs no downloads or real data.

By the end you'll have seen every moving part: the integrated object, the HTML
inspection report, auto-flagged clusters, per-cell contamination flags, a
`decisions.yaml`, a filtered object, and a cross-round summary.

> Prerequisite: a working install (see [getting-started.md](getting-started.md)).
> All commands below assume you're in the repo root with the `scvi-loops`
> environment activated.

## Step 0 — Generate the example data

```bash
python examples/make_example_data.py --output examples/example_data/combined_example.h5ad
```

This writes ~2,060 cells × ~1,256 genes with **raw counts in `.X`**, three batches
(`data_origin` = batch1/2/3) across two platforms (`tech` = tech1/tech2). It builds
in structure the loop can act on: a clear batch effect, distinct marker-defined
populations, a **low-quality population** (few genes, high mitochondrial fraction),
and small **endothelial** and **hepatic contaminant** populations.

The tutorial config `examples/pipeline_tutorial.yml` already points at this file and
uses a small, fast scVI architecture. (It resolves `input_h5ad` relative to itself,
so the paths just work.)

## Step 1 — Integrate (round 1)

```bash
cd code
python integrate_scvi.py --config ../examples/pipeline_tutorial.yml
```

Outputs land in `examples/tutorial_run/rounds/round_01/`. The key file is
`integrated.h5ad` (log-normalized `.X`, raw counts in `layers["counts"]`, the scVI
embedding in `obsm["X_scVI"]`, a UMAP, and Leiden clusters in `obs["leiden"]`).
Alongside it: `run_config_resolved.yml`, `round_manifest.json` (seed, git commit,
input fingerprint), `scib_benchmark_results.csv`, and QC plots.

## Step 2 — Inspect (round 1)

```bash
python inspect_integration.py --config ../examples/pipeline_tutorial.yml
```

Open `examples/tutorial_run/rounds/round_01/inspection_report.html`. On this dataset
you'll see the inspection:

- **Auto-flag one cluster** for high mitochondrial fraction — the low-quality
  population (median ~18–19% MT, above the 15% threshold). It's written to
  `auto_flags.yaml` as a suggested `remove_clusters` entry.
- **Flag ~160 contamination cells** — roughly 90 endothelial + 70 hepatic, written
  to `contamination_flagged_cells.csv` (and scored in `contamination_zscore.csv`).

Other useful artifacts in the same folder: `cluster_qc_summary.csv`,
`marker_dotplot.png`, `cluster_silhouettes.png`, `contamination_umap.png`. See
[interpreting-outputs.md](interpreting-outputs.md) for how to read them.

## Step 3 — Decide what to remove

`auto_flags.yaml` is the machine's *suggestion*; you turn it into a `decisions.yaml`
you endorse. Create `examples/tutorial_run/rounds/round_01/decisions.yaml`:

```yaml
integration:
  variant: "default"        # which embedding/clustering the flags came from
  cluster_key: "leiden"

remove_clusters:
  - cluster: "2"            # the high-MT low-quality cluster (check YOUR report:
    reason: "Low-quality: high mitochondrial fraction (auto-flagged)"

remove_cells:
  - cells_file: "contamination_flagged_cells.csv"   # resolved next to this file
    reason: "Endothelial/hepatic contamination (per-cell z-score flags)"

notes: "Tutorial round 1 review."
```

> Cluster IDs are seed-stable but the *number* of a given cluster depends on the
> run — confirm which cluster your report flagged and edit `cluster: "2"` to match.
> The `integration:` block records which embedding the flags came from so filtering
> uses the same one. See [decisions.md](decisions.md) for the full syntax.

## Step 4 — Filter into round 2

```bash
python inspect_integration.py \
  --config ../examples/pipeline_tutorial.yml \
  --input ../examples/tutorial_run/rounds/round_01/integrated.h5ad \
  --decisions ../examples/tutorial_run/rounds/round_01/decisions.yaml \
  --output-dir ../examples/tutorial_run/rounds/round_02
```

This writes `round_02/filtered.h5ad` plus `decisions_applied.yaml` (what actually
ran) and `filtering_retention_summary.csv` (pre/post counts by cluster and batch).
You'll see roughly **2,060 → 1,740 cells** (the ~160-cell cluster + ~160 contaminant
cells removed; because the contaminants formed their own clusters, those clusters go
to zero retention).

## Step 5 — Re-integrate and inspect (round 2)

```bash
python integrate_scvi.py \
  --config ../examples/pipeline_tutorial.yml \
  --input ../examples/tutorial_run/rounds/round_02/filtered.h5ad \
  --output-dir ../examples/tutorial_run/rounds/round_02

python inspect_integration.py \
  --config ../examples/pipeline_tutorial.yml \
  --input ../examples/tutorial_run/rounds/round_02/integrated.h5ad \
  --output-dir ../examples/tutorial_run/rounds/round_02
```

Round 2 is the payoff: **no clusters are auto-flagged**, and contamination drops from
~160 cells to a handful (~4, <0.3%). That's what convergence looks like — the reason
to *stop* iterating (more on stopping criteria in
[interpreting-outputs.md](interpreting-outputs.md#when-to-stop-iterating)).

## Step 6 — Summarize the rounds

```bash
python summarize_rounds.py \
  --rounds-dir ../examples/tutorial_run/rounds \
  --verify \
  --output-dir ../examples/tutorial_run/summary
```

Open `examples/tutorial_run/summary/pipeline_summary.html`. You get a round table
(cell/cluster counts, seed, commit per round), a `decisions_ledger.csv` of every
keep/remove action, a Mermaid lineage diagram, and `--verify` integrity checks
(input-hash continuity, seed/commit drift). The round table here shows
round_01 (2,060 cells, 7 clusters) → round_02 (1,740 cells, 4 clusters).

## What you just learned

- The loop is **integrate → inspect → decide → filter → re-integrate**, with a human
  gate at the `decisions.yaml` step.
- Two complementary filtering signals: **cluster-level** auto-flags (QC) and
  **cell-level** contamination flags (marker z-scores).
- Provenance is automatic — every round records its seed, code commit, and input hash.
- You stop when new rounds stop flagging things.

Next: point it at your own data with the [configuration reference](configuration.md)
and, for non-mouse/non-SNS datasets, [adapting to your data](adapting-to-your-data.md).

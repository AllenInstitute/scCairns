# Changelog

All notable changes to **scCairns** are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `cairns record-filter` (`sccairns.record_filter`): records provenance for a filtering
  step performed **outside** scCairns — e.g. an interactive lasso or cluster selection in
  a notebook. It applies no filtering of its own; it diffs an input/output `.h5ad` pair by
  `obs_names` and writes the same `decisions` block `cairns summarize` already consumes,
  so an interactive round appears in the round table, decisions ledger, and lineage
  instead of breaking the chain. Writes `round_manifest.json` (input **and** output
  fingerprints, parent-round pointer, per-cluster `actions_summary`, removal count and
  percentage, `filtered_on.source: "interactive"`) plus `cells_to_keep.csv`,
  `cells_removed.csv`, and `decisions_applied.yaml`. Cluster and UMAP keys are inferred
  from `--embedding-key` when not passed; missing keys warn rather than fail, but an
  output containing cell IDs absent from the input is an error.
- `round_manifest.json` now records `is_first_round` (derived from the parent-round
  pointer) so it is explicit whether QC filtering ran for the round.
- `data.obs_filter`: an optional pandas query string applied to `adata.obs` right after
  load (e.g. `"condition in ['Control', 'Saline']"`), for integrating a subset without
  writing a pre-filtered `.h5ad`. Recorded in `round_manifest.json`; `null` disables it.
- `round_manifest.json` records the full cell-attrition chain for a round —
  `n_cells_input`, `n_genes_input`, `n_cells_after_obs_filter`,
  `n_cells_prefilter_plot`, and `qc_filter_applied` — alongside the existing
  post-QC `n_cells`. Previously only the final count was written, so the QC plots
  (drawn pre-QC-filter) and the manifest disagreed with nothing to reconcile them.
- `cairns summarize` surfaces the input count: an **In cells** column in the round
  table and `rounds_table.csv`, input → final counts on the lineage nodes and round
  cards, and a `pre_integration_drop` finding under `--verify`. A new
  `input_cells_mismatch` warning fires when a round's integration read a different
  cell count than its own decision stage wrote.

### Changed
- Default `qc.min_cells` lowered from `3` to `1`. The gene filter now keeps every gene
  observed in at least one cell (only all-zero columns are dropped), so rare but real
  markers — e.g. low-abundance neuropeptides such as `Tac1`/`Cck` — are no longer
  removed. Updated in `DEFAULT_CONFIG` and all shipped example/`code` configs. Raise it
  to be stricter; `0` is discouraged (it keeps all-zero columns).
- Retention percentages in the summary report are measured against the first round's
  *input* cell count, so `obs_filter` and QC losses are included. Rounds whose
  manifests predate `n_cells_input` fall back to the old post-QC baseline.

### Fixed
- `qc.skip_filter_after_round_1` is now honored — it was previously declared in the
  config but read by no code, so the QC gene/cell filter re-ran on every round. Because
  each round operates on a progressively smaller cell subset, a gene kept in round 1
  (detected in many cells) could fall below `min_cells` and be dropped in a later round.
  The filter now runs on the first round only (detected via the input's parent-round
  pointer); later rounds skip it, so genes are never re-filtered against a subset. This
  was the root cause of neuropeptide genes disappearing between rounds.

### Documentation
- README: `record-filter` added to the module/command table, with a runnable example in
  the provenance section for the interactive-filtering case.

## [0.1.0] — unreleased

First packaged release: the project became an installable Python package with a unified
CLI, renamed from `scvi_integration_loops` to **scCairns**.

### Added
- Installable `sccairns` package with a single `cairns` command
  (`cairns integrate | inspect | summarize | flag-contamination`).
- `inspection.all_architectures` config field to drive per-architecture sweep inspection
  from config (still available as the `--all-sweep-architectures` flag).
- `cairns summarize` recognizes inspection-only sweep rounds (a re-inspection that writes
  no new `.h5ad`) and splices them inline in the lineage instead of leaving them dangling.

### Changed
- **Renamed** `scvi_integration_loops` → **scCairns**: import package `sccairns`, PyPI
  distribution `scCairns`, conda env `sccairns`.
- The package lives at `code/sccairns/` (not the repo root) so Code Ocean reproducible
  runs — which ship only the `code/` folder — can import it; `pyproject.toml`
  (`package-dir`) maps it back to the `sccairns` import name for local/PyPI installs.
- `harmonypy` is now an optional extra (`pip install 'scCairns[harmony]'`, capped `<1`;
  2.x breaks `scanpy.external.pp.harmony_integrate`). Harmony remains an opt-in path.
- The `code/*.py` entry points are thin backward-compat shims over the package, so the
  Code Ocean capsule keeps working unchanged.
- Renamed `integrate_sns_scvi.py` → `integrate_scvi.py` (the engine is dataset-agnostic;
  SNS lives in config defaults), with a forwarding shim at the old path.

### Documentation
- `docs/` guide set: getting-started, tutorial, configuration, integration-methods,
  decisions, interpreting-outputs, adapting-to-your-data, troubleshooting, Code Ocean.
- `examples/make_example_data.py` (synthetic multi-batch generator) + `pipeline_tutorial.yml`
  so the tutorial runs on a laptop CPU with no data; `pipeline_generic.yml` non-SNS template.
- README restructured into an on-ramp linking the guides.
- Archived the point-in-time `TRACEABILITY_EVALUATION.md` under [`archive/`](archive/).

---

## Pre-0.1.0 development history

Changes before the project was packaged and renamed (dates approximate, from git history).

### 2026-07
- **Contamination flagging** (`flag_contamination`): marker-based, per-cell,
  cluster-independent contamination detection via z-scores, with default mouse
  endothelial/hepatic panels; contamination plots added to the inspection report.
- **Multi-round summary** (`summarize_rounds`): combined round table, decisions ledger,
  Mermaid lineage, and `--verify` integrity checks; handles rounds archived under
  arbitrary paths and records the terminal filtering architecture.
- scIB metric fix; HTML output improvements for the summary.

### 2026-06
- **Reproducibility & provenance:** configurable `reproducibility.seed` applied to all
  stochastic steps (stable Leiden cluster IDs across re-runs); git commit + dirty flag,
  input SHA-256 fingerprint, and explicit `parent_round` lineage pointer recorded in
  `round_manifest.json`; Run Provenance footer added to the HTML report; safe YAML
  serialization of applied decisions; append-only `command_args.json` history.
- **Harmony option:** optional same-run Harmony integration benchmarked against scVI and
  a PCA baseline via scib-metrics.
- **Sweeps:** inspection accepts multiple swept architectures; explicit variant pinning
  for filtering decisions (`--cluster-key` / decisions `integration:` block).

### Earlier
- Config-driven V1 pipeline: versioned YAML config, deep-merge with CLI overrides, schema
  validation, resolved-config persistence, and the
  `integrate → inspect → decide → filter → re-integrate` loop with scVI/scANVI.

[Unreleased]: https://github.com/AllenInstitute/scCairns/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/AllenInstitute/scCairns/releases/tag/v0.1.0

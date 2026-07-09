# Changelog

Notable changes to the pipeline, newest first. Dates are approximate, drawn from git
history; this project does not (yet) publish tagged releases.

## Unreleased

### Changed
- Renamed the entry script `integrate_sns_scvi.py` → `integrate_scvi.py` (the engine is
  dataset-agnostic; SNS lives in config defaults). A deprecation shim at the old path
  forwards to the new one so existing commands keep working.

### Documentation & usability
- Added a `docs/` guide set: getting-started (install + input contract), a runnable
  end-to-end tutorial, output interpretation, adapting-to-your-data (non-mouse/SNS),
  configuration reference, decisions syntax, troubleshooting, and Code Ocean usage.
- Added `docs/integration-methods.md`: tuning the scVI architecture (config/CLI/sweeps)
  and using a different integration method (Harmony baseline, bring-your-own-embedding
  with Scanorama, or adding a first-class method via the Harmony template).
- Restructured `README.md` into an on-ramp that links to the guides; removed the
  duplicated Quick Start block and Code Ocean absolute paths from the local flow.
- Added `examples/make_example_data.py` (synthetic multi-batch generator) and
  `examples/pipeline_tutorial.yml` so the tutorial runs on a laptop CPU with no data.
- Added `examples/pipeline_generic.yml`, a non-SNS template documenting every
  domain-specific override point.
- Annotated `TRACEABILITY_EVALUATION.md` to reflect the reproducibility work that has
  since been implemented.

## 2026-07
- **Contamination flagging** (`code/flag_contamination.py`): marker-based, per-cell,
  cluster-independent contamination detection via z-scores, with default mouse
  endothelial/hepatic panels; contamination plots added to the inspection report.
- **Multi-round summary** (`code/summarize_rounds.py`): combined round table,
  decisions ledger, Mermaid lineage, and `--verify` integrity checks; handles rounds
  archived under arbitrary paths and records the terminal filtering architecture.
- scIB metric fix; HTML output improvements for the summary.

## 2026-06
- **Reproducibility & provenance:** configurable `reproducibility.seed` applied to all
  stochastic steps (stable Leiden cluster IDs across re-runs); git commit + dirty flag,
  input SHA-256 fingerprint, and explicit `parent_round` lineage pointer recorded in
  `round_manifest.json`; Run Provenance footer added to the HTML report; safe YAML
  serialization of applied decisions; append-only `command_args.json` history.
- **Harmony option:** optional same-run Harmony integration benchmarked against scVI
  and a PCA baseline via scib-metrics.
- **Sweeps:** inspection accepts multiple swept architectures; explicit variant pinning
  for filtering decisions (`--cluster-key` / decisions `integration:` block).

## Earlier
- Config-driven V1 pipeline: versioned YAML config, deep-merge with CLI overrides,
  schema validation, resolved-config persistence, and the
  `integrate → inspect → decide → filter → re-integrate` loop with scVI/scANVI.

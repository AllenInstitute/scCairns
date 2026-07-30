# Changelog

All notable changes to **scCairns** are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

_Nothing yet._

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

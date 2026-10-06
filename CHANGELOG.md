# Changelog

All notable changes to **scCairns** are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.1] — 2026-10-06

### Changed
- README: added `pip install scCairns` to Quick start; converted all relative doc links
  to absolute GitHub URLs so they resolve correctly on PyPI.
- `pyproject.toml`: added `[project.urls]` (Homepage, Source, Getting started, Changelog)
  for the PyPI sidebar.

## [0.3.0] — 2026-10-05 — 2026-10-05

### Added
- `sccairns.stress`: heat-shock / dissociation-stress module scoring, surfaced in the
  inspection report as `stress_score_by_cluster.png` and as a `median_stress` column in
  `cluster_qc_summary.csv`. A cluster defined by Hspa1a/Hspb1/Hsp90aa1 is usually an
  artifact of how long the tissue sat in protease rather than a cell type, and nothing
  in the report made that visible before. Configured under `inspection.stress`
  (`enabled`, `gene_prefixes`, `genes`, `z_thresh`, `layer`, `score_key`).

  **Heat-shock only, deliberately.** Mitochondrial burden — the other half of the usual
  "stress" bundle — is already covered end to end by `pct_counts_mt` → `median_pct_mt` →
  `auto_flag.mt_threshold`, so scoring it again would produce two differently-scaled
  numbers that can disagree about the same cluster. The panel is meant to be read
  *alongside* `median_pct_mt`, and the plot says so. Immediate-early genes (Fos, Jun,
  Egr1) are likewise excluded by default: in neurons they are genuine activity markers,
  so scoring them as stress would penalize real biology — add them via `gene_prefixes`
  where the tissue makes that unambiguous.

  The panel is resolved by **prefix** against `var_names` at run time rather than listed
  gene by gene, because HSP families are large and the informative member varies by
  dataset; the same default therefore covers mouse (`Hspa1a`) and human (`HSPA1A`)
  without restating it, which matters because `gene_symbol_case: preserve` means neither
  convention can be assumed. Resolution is `sorted()`, not `set()`-ordered — Python
  randomizes string hashing per process, so a set would make the recorded panel differ
  run to run.

  Scoring reuses `contamination.score_panel` (the same per-gene z-score path), so the
  two are read on one scale, but the stress panel is kept out of `DEFAULT_CONTAM_PANELS`
  on purpose: contamination panels feed `flag_any_contam` and therefore
  `contamination_flagged_cells.csv`, and a stressed cell is a "drop this cluster" call,
  not a wrong-lineage call. Scoring runs before `cluster_qc_summary()` so the column
  actually reaches the CSV — contamination runs after it, which is why no contamination
  column appears there.

## [0.2.0] — 2026-09-15

Diagnostics and reproducibility release. Adds three cluster-level diagnostic
families — doublet scoring, neighborhood entropy/purity, and compositional bias
testing — and closes the gap between what `reproducibility.seed` recorded and
what it actually governed. Upgrading changes UMAP/Leiden output for any config
with a seed other than 0; see the seed entry under **Fixed**.

### Added
- Per-cluster **doublet score** panel in the inspection report
  (`doublet_score_by_cluster.png`), rendered only when a per-cell `doublet_score` exists in
  `obs`. Nothing is computed here — it plots whatever an upstream step produced
  (DoubletFinder in R, `sc.pp.scrublet`, …). Because droplet-based callers leave a
  plate-based arm unscored, **unscored cells are excluded rather than counted as zero**:
  the title reports how many of the total were scored and names the scored platforms, and
  each cluster's tick shows its scored `n`. An optional call column is overlaid as the
  fraction called per cluster, accepting booleans, a 0/1 flag, or DoubletFinder's
  `"Doublet"`/`"Singlet"` strings. Configured under `inspection.doublets`
  (`score_key`, `call_key`, `platform_key`).
- `sccairns.composition`: per-cluster compositional bias testing. `compositional_bias()`
  asks whether a cluster's platform/batch/donor makeup is what random assignment would
  give — hypergeometric, both tails (a missing platform is *depletion*, a dominating donor
  is *enrichment*), Bonferroni-corrected within a cluster for the group being chosen by
  inspection and Benjamini–Hochberg across clusters. `cairns inspect` runs it over
  `data.batch_key`, `data.sample_key`, and the categorical covariates by default, adding
  `<key>_depleted`, `<key>_depletion_ratio`, `<key>_depletion_q`, `<key>_top`,
  `<key>_enrichment`, `<key>_enrich_q`, and `<key>_n_groups_present` to
  `cluster_qc_summary.csv`, with new `inspection.auto_flag` keys `composition_q` (0.01),
  `composition_min_cells` (20), `min_enrichment` (2.0), `max_depletion_ratio` (0.20), and
  `composition_keys` (null → auto).

  This exists because composition alone is not evidence. With 15% SSv4 overall, a 20-cell
  cluster with no SSv4 (p ≈ 0.04, i.e. chance) and a 200-cell cluster with no SSv4
  (p ≈ 2e-15, i.e. impossible) are both "100% one platform" and score identically under
  any `dominant_*_frac` threshold. Flagging requires a significance **and** an effect-size
  gate, since at 400 cells a biologically unremarkable 1.1–1.3× deviation already reaches
  q ≈ 1e-3. The flag's message ends by prompting a check of the sampling design, because
  a platform that was FACS-sorted differently will legitimately be absent from whole
  clusters.
- `data.sample_key` is now read. It was declared in `DEFAULT_CONFIG` but consumed by no
  code; it names the donor / biological-replicate column the donor-bias test groups over.
- `sccairns.entropy`: neighborhood label entropy and neighbor purity, computed on a k-NN
  graph over any `obsm` embedding. `neighborhood_entropy()` scores each cell (useful for
  batch mixing or annotation disagreement); `cluster_label_coherence()` aggregates per
  cluster. Ported from an R implementation, with the neighbor graph built on the full
  embedding rather than a single component, vectorized (100k cells × 30 dims, k=15, two
  annotations in ~4.5s), and missing labels excluded from a neighborhood rather than
  counted as their own category.
- `cairns inspect` flags clusters that are **not well supported in the embedding**:
  `cluster_qc_summary.csv` gains `neighbor_purity`, `neighbor_entropy`,
  `dominant_neighbor`, and `dominant_neighbor_frac`, and `inspection.auto_flag`
  gains `min_neighbor_purity` (default 0.50), `entropy_threshold` (default null, opt-in),
  and `entropy_neighbors` (default 15). The flag names the absorbing cluster when one
  dominates — "only 10% of neighbors share its label — 100% of its neighborhood is
  cluster 0" — which distinguishes a merge candidate from a diffuse cluster to drop.

  The flag is keyed on **purity rather than entropy** deliberately. Entropy measures how
  diverse a neighborhood is, not how much of it agrees, so a cluster wholly absorbed into
  another has a homogeneous neighborhood — of someone else's label — and scores near-zero
  entropy, tying with a perfectly isolated cluster. In testing, an entropy ranking placed
  the absorbed cluster at or below the healthy cluster hosting it, so no entropy threshold
  could catch it without also flagging good clusters. Entropy is still reported, and is
  the right statistic for characterizing *which kind* of mixing a flagged cluster has.

### Fixed
- The cluster QC summary in `inspection_report.html` was rendered with
  `summary_df.round(2)`, so every q-value below 0.005 displayed as `0.0` — a decisive
  q = 6.8e-20 and a marginal q = 0.003 both read as zero. The table now uses per-column
  formatters: `_q` columns keep scientific notation (`6.81e-20`), other floats get two
  decimals, and an underflowed q renders as `<1e-300` rather than implying an exact zero.
  Applied to the console table as well so the two agree. `cluster_qc_summary.csv` was
  never affected — it is written from the unrounded frame.
- `reproducibility.seed` now governs **every** stochastic step, not just model
  training. `set_global_seed()` seeded scvi-tools (and therefore scVI/scANVI), but the
  neighbor-graph, UMAP, and Leiden wrappers called Scanpy without `random_state`, so
  they silently ran on Scanpy's implicit default of 0 regardless of the config:
  `seed: 42` changed the latent space and left UMAP/Leiden on 0. The seed is now
  threaded into all three wrappers — on their older-Scanpy fallback paths as well as
  the primary ones, since which path executes depends on the installed Scanpy version
  and an unseeded fallback would reintroduce the gap invisibly — and into the
  benchmarking stage's PCA baseline and Leiden proxy label, both of which feed
  `scib_benchmark_results.csv`. `seed: null` normalizes to 0 throughout, matching the
  convention `run_harmony` already used, so an unseeded run stays deterministic here.

  Leiden is why this mattered: cluster IDs are what `decisions.yaml` files reference
  by number, so a manifest recording `seed: 42` did not describe the RNG state that
  produced the IDs its own decisions cite. **This changes output** for any config with
  a seed other than 0 — UMAP coordinates and Leiden labels differ from pre-0.2.0
  rounds, so a `decisions.yaml` written against the old cluster numbering should be
  re-checked against a fresh inspection report before it is reapplied.
- `cairns flag-contamination` raised `AttributeError` instead of running. The CLI
  dispatcher calls `module.main()`, but `contamination.py` defined its argparse only
  under `if __name__ == "__main__"`, so the subcommand had never worked since the
  package was introduced in 0.1.0. Refactored into a module-level `build_parser()` +
  `main(argv)` matching `record_filter`. No round record was affected: the
  functionality itself was never unavailable, because contamination flagging runs
  automatically inside `cairns inspect` report mode (`contamination_zscore.csv`,
  `contamination_flagged_cells.csv`). The existing CLI test only asserted that command
  names appeared in the usage string, which is why this went unnoticed — every command
  in `COMMANDS` is now asserted to export a callable `main()`.

### Documentation
- [`archive/DOUBLET_DETECTION_ASSESSMENT.md`](archive/DOUBLET_DETECTION_ASSESSMENT.md) —
  point-in-time feasibility review of adding Scrublet to the QC stage. **Not implemented:**
  `sc.pp.scrublet` ships inside scanpy 1.10.4 (no `scrublet` package needed) but hard-fails
  with `ModuleNotFoundError: No module named 'skimage'`, and `scikit-image` is pinned
  nowhere. That is a one-line fix, so the blocker is design rather than environment: the
  method assumes droplet co-encapsulation and this dataset mixes 10x Multiome with
  plate-based SSv4. Records the batch-key, first-round-only, and score-don't-remove
  decisions, and notes that `decisions.yaml` can already filter on `predicted_doublet`
  with no new code.
- [`archive/PIANO_INTEGRATION_ASSESSMENT.md`](archive/PIANO_INTEGRATION_ASSESSMENT.md) —
  point-in-time feasibility review of adding [PIANO](https://github.com/NingWang1729/piano)
  as an integration method. **Not implemented:** the code seam is ~1 day (Harmony is the
  precedent), but PIANO requires Python ≥3.11 and numpy ≥2 while the Code Ocean image is
  `python3.10.12` with `numpy==1.26.4` pinned, so the real cost is an environment
  migration. Also records the GPL-3.0 licensing question and a zero-code path for
  evaluating PIANO out-of-band today.
- The tutorial, decisions, and integration-methods guides now invoke the `cairns` CLI
  instead of `python code/<script>.py`. Those calls predate the packaging; they also
  required a `cd code`, so every path in the tutorial was relative to a directory the
  reader had no other reason to be in — and the tutorial's own preamble said to run from
  the repo root. The `code/*.py` shims remain for the Code Ocean capsule.
- README's "one full round, locally" ends with `cairns summarize --rounds-dir rounds
  --verify`, so the quick start covers the step that actually checks the round rather
  than stopping at re-integration.
- [Configuration reference](docs/configuration.md) documents `cairns record-filter`,
  which is deliberately config-free — every flag, which of them are required, and why the
  flags belong in a run script rather than a config a later round could silently reuse.
- Same page: added the missing `data.obs_filter` field, and corrected the `data.output_dir`
  default (`../results` → `results`) with a note that a relative path in a config resolves
  against the config file's own location.

## [0.1.0] — 2026-08-11

First tagged release: the project became an installable Python package with a unified
CLI, renamed from `scvi_integration_loops` to **scCairns**, and is now consumable as a
pinned dependency (`pip install git+https://github.com/AllenInstitute/scCairns.git@v0.1.0`)
rather than by cloning the repo into each Code Ocean capsule.

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
- `record-filter` accepts **cell-ID lists** in place of `.h5ad` files:
  `--input-cell-ids` / `--output-cell-ids`. The filtered side is only ever read for the
  IDs that survived, so a tool that cannot write AnnData — Seurat, Loupe — can hand back
  a CSV of barcodes instead. Reads a bare one-per-line list, its own `cells_to_keep.csv`,
  and both common R exports (`write.csv(data.frame(x = Cells(obj)))`, which prepends a
  row-number column, and `write.csv(obj@meta.data)`, where the barcodes are rownames).
  Duplicate IDs are rejected, and mismatched IDs name the likely cause — Seurat merges
  appending `_1`, sample-name prefixes, and `-1` suffix changes all break ID matching.
  `--embedding-key` is optional, for a filter applied at ingest before any embedding or
  clustering exists. The manifest records a `compared` block
  (`{"input": "h5ad"|"cell_ids", "output": …}`) so a null `input_genes`/`output_genes` or
  a missing per-cluster breakdown reads as a consequence rather than a gap.
- `round_manifest.json` now records `is_first_round` (derived from the parent-round
  pointer) so it is explicit whether QC filtering ran for the round.
- `data.obs_filter`: an optional pandas query string applied to `adata.obs` right after
  load (e.g. `"condition in ['Control', 'Saline']"`), for integrating a subset without
  writing a pre-filtered `.h5ad`. Recorded in `round_manifest.json`; `null` disables it.
- The manifest's `code` block identifies the pipeline by **version as well as commit**:
  new `version` (installed distribution), `source` (`checkout` / `installed` /
  `unknown`), and `repo_root` fields. Git metadata only exists when scCairns runs from a
  checkout, so a wheel installed into a Code Ocean image previously recorded an all-null
  `code` block — the round was untraceable to any code. `repo_root` names the repository
  a commit came from, so a commit picked up from a surrounding project checkout is
  visible rather than mistaken for the pipeline's own.
- `cairns --version` (also `-V`) reports the running scCairns using the same identity
  the manifest records — `cairns 0.1.0 (installed)` for a pinned release, or
  `cairns 0.1.0 (checkout <sha> in <repo>, dirty)` from a working tree. This is how you
  confirm which build a Code Ocean image actually got.
- `cairns summarize` falls back to the version when there is no commit: the **Commit**
  column shows `v0.1.0`, `code_version`/`code_source` are exported to `rounds_table.csv`,
  `--verify` reports version drift across rounds when no commits are recorded, and a new
  `code_unidentified` warning fires for a round that records neither identifier.
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
- The no-config output directory default is `./results` instead of `../results`, which
  only made sense when the CLI was invoked from `code/`. Paths *inside* a config file are
  unaffected — they already resolve against the config file's own directory, so a config
  is portable no matter where `cairns` is invoked from. `cairns integrate` no longer
  passes the legacy `../results`; every invocation in `code/run` supplies `--config`, so
  capsule behavior is unchanged.

### Fixed
- `qc.skip_filter_after_round_1` is now honored — it was previously declared in the
  config but read by no code, so the QC gene/cell filter re-ran on every round. Because
  each round operates on a progressively smaller cell subset, a gene kept in round 1
  (detected in many cells) could fall below `min_cells` and be dropped in a later round.
  The filter now runs on the first round only (detected via the input's parent-round
  pointer); later rounds skip it, so genes are never re-filtered against a subset. This
  was the root cause of neuropeptide genes disappearing between rounds.

### Packaging and rename
The project became an installable Python package with a unified CLI, renamed from
`scvi_integration_loops` to **scCairns**.

- Installable `sccairns` package with a single `cairns` command
  (`cairns integrate | inspect | summarize | record-filter | flag-contamination`).
- `inspection.all_architectures` config field to drive per-architecture sweep inspection
  from config (still available as the `--all-sweep-architectures` flag).
- `cairns summarize` recognizes inspection-only sweep rounds (a re-inspection that writes
  no new `.h5ad`) and splices them inline in the lineage instead of leaving them dangling.
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
- [Filtering outside scCairns](docs/decisions.md) — recording a filter applied in Seurat,
  Loupe, or a notebook: the accepted cell-ID export shapes, what an ID list costs you
  versus an `.h5ad`, and why cell IDs must survive the round trip unchanged.
- [Getting started](docs/getting-started.md) now leads with installing a **pinned
  release as a dependency** (`pip install "git+…@v0.1.0"`) for analysis projects, with
  cloning presented as the develop-scCairns path. Covers why the tag must be exact, that
  the repo is private so pip needs GitHub credentials, and that `pyproject.toml` deps are
  unpinned by design while `requirements.txt` holds the reproducible pins.
- [Code Ocean](docs/codeocean.md) gains a section on how `cairns` reaches the capsule:
  vendored under `code/` versus installed from a pinned tag in `postInstall`, why the
  install cannot live in `run` (reproducible runs are offline), and why it must not go in
  the Code-Ocean-generated `Dockerfile`.
- Fixed a stale verification command in getting-started: `compileall` referenced a
  root-level `sccairns/` that has not existed since the package moved to `code/sccairns/`.
- README: `record-filter` added to the module/command table, with a runnable example in
  the provenance section for the interactive-filtering case.
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

[Unreleased]: https://github.com/AllenInstitute/scCairns/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/AllenInstitute/scCairns/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/AllenInstitute/scCairns/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/AllenInstitute/scCairns/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/AllenInstitute/scCairns/releases/tag/v0.1.0

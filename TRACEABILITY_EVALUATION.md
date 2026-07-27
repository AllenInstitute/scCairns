# Traceability & Decision-Documentation Evaluation

**Repository:** `scCairns`
**Scope of evaluation:** How well the repository achieves its stated aim of
documenting decisions for complete traceability of the integration pipeline.
**Reviewed:** code, config system, provenance helpers, example round outputs,
environment pinning, tests, and git state.

---

## Update (2026-07): most high-priority gaps are now closed

This evaluation was written against an earlier state of the repo. The reproducibility
and run-identity gaps it flagged as highest-leverage have since been implemented — so
read the scorecard and "Gaps" section below as **historical context**, not the current
state. Mapping of the original gaps to what now exists:

| Original gap | Status | Where |
|---|---|---|
| 1. No seed / determinism | **Fixed** | `set_global_seed()` sets `scvi.settings.seed` before every stochastic step; recorded under `reproducibility.seed` in `round_manifest.json`. Leiden cluster IDs are now stable across re-runs of the same config. |
| 2. No code-version capture | **Fixed** | `collect_code_provenance()` records git commit/branch/`dirty` flag into `manifest["code"]`. |
| 3. Input identified by path only | **Fixed** | `fingerprint_file()` records path + size + mtime + SHA-256 per stage (`input_fingerprint`). |
| 4. Implicit cross-round lineage | **Fixed** | `read_parent_provenance()` writes a structured `parent_round` pointer; `summarize_rounds.py --verify` validates hash continuity. |
| 6. `command_args.json` overwrites | **Fixed** | Now keeps an append-only `history` list of every invocation. |
| 7. Report decoupled from provenance | **Fixed** | `inspection_report.html` ends with a Run Provenance footer (seed, commit, input SHA-256, package versions). |
| 8. `decisions_applied.yaml` hand-serialized | **Fixed** | `write_yaml_record()` uses `yaml.safe_dump`. |
| 5, 9, 10. Package-version breadth, enforced rationale, naive timestamps | **Partly open** | Lower-priority; see below. |

With items 1–4, 6–8 addressed, the pipeline has moved from "well-documented decisions"
to substantially reproducible traceability. The remaining open items are the
lower-priority ones. The original analysis is retained below for reference.

---

## Summary judgment

The repository is **deliberately and competently designed around decision
traceability**, and on the *decision-documentation* axis specifically it is
strong: filtering choices flow through an auditable
`auto_flags.yaml → decisions.yaml → decisions_applied.yaml` chain, each step
carries human-readable reasons and exact cell counts, and every run emits a
resolved config plus a cumulative manifest. The design intent is real, not
cosmetic.

Where it falls short of *complete* traceability is **reproducibility of the
artifacts the decisions are made against**. There is no random-seed control, no
code-version capture, and no input-data fingerprinting. As a result the pipeline
documents *what configuration and what rationale* produced a round, but it
cannot guarantee that re-running that configuration reproduces the same
embedding, the same Leiden cluster IDs, or therefore the same target of a
decision like `remove_clusters: cluster "7"`. That is the central gap between
"well-documented" and "fully traceable."

**Overall: B / B-plus.** Excellent decision-rationale capture and config
provenance; incomplete reproducibility and run-identity capture.

---

## Scorecard

| Dimension | Grade | One-line basis |
|---|---|---|
| Decision capture & rationale | **A-** | `decisions.yaml` + `decisions_applied.yaml` + retention CSV; reasons + counts logged |
| Config provenance | **A-** | Fully resolved config written per stage; defaults+overrides merged and persisted |
| Run/command provenance | **B** | `command_args.json` (argv + timestamp); but overwrites per stage, no history |
| Environment provenance | **B** | Dockerfile fully pinned (sha256); runtime manifest captures only 4 packages |
| Input-data provenance | **C+** | Code Ocean `datasets.json` pins dataset UUIDs; manifest records only file *path*, no hash |
| Code-version provenance | **D** | No git commit/revision captured in any output |
| Reproducibility / determinism | **D** | No seed set for scVI/scANVI/UMAP/Leiden; outputs not reproducible |
| Cross-round lineage | **C** | Round chaining is implicit in paths/run script, not recorded as a parent pointer |
| Human-facing report linkage | **C** | HTML report shows counts + timestamp only; no config/version/decision provenance embedded |
| Record integrity | **B-** | Fail-loud on bad queries (good); but `decisions_applied.yaml` is hand-serialized, fragile |
| Tests of the traceability path | **B** | Decisions-apply and config-merge are unit-tested |

---

## What the repository does well

**1. Decision rationale is first-class and auditable.**
The filtering workflow is the strongest part of the design. `inspect_integration.py`
auto-generates `auto_flags.yaml` — an editable template that records the
thresholds used (`mt_pct`, `min_genes`, `min_cells`, `single_batch_frac`), the
flagged clusters, their cell counts, and machine-suggested *reasons*. The user
promotes this to `decisions.yaml`, and `apply_decisions()` then writes
`decisions_applied.yaml` recording, per action, the type, the cluster/query, the
**exact number of cells excluded or removed**, and the reason string. This is a
genuine closed loop: a reviewer can read *what was removed, how many, and why*
without rerunning anything.

**2. Quantified, multi-axis retention accounting.**
`filtering_retention_summary.csv` reports pre/post counts overall, by cluster, by
batch, and by cluster×batch. This makes the *consequence* of a decision
inspectable, not just the decision itself — important for catching a filter that
silently nukes one batch.

**3. Config is the single source of truth, and it is persisted resolved.**
`run_config_resolved.yml` is written after defaults are merged and CLI overrides
applied, so the persisted file is the *effective* configuration, not the
user-edited fragment. Defaults live in one place (`DEFAULT_CONFIG`), are
deep-merged, validated (`validate_config`), and version-gated
(`pipeline_version: 1`). The example `round_01/run_config_resolved.yml` confirms
this captures even implicit defaults (marker sets, auto-flag thresholds).

**4. Cumulative per-round manifest.**
`round_manifest.json` accretes `integration`, `inspection`, and `decisions`
stage payloads into one file, with `created`/`updated` timestamps, Python
version + executable, and key output identifiers (embedding/UMAP/cluster keys,
model dirs, n_cells/n_genes).

**5. Fail-loud integrity default.**
Failed decision queries raise by default (`ignore_failed_queries: false`). A
mistyped `remove_cells` query fails the run rather than silently filtering
nothing — the right default for a record that is meant to be trusted.

**6. Environment is pinned at the image level.**
`environment/Dockerfile` is content-hashed (`sha256:e7f06…`) and pins every apt
and pip dependency to an exact version; `requirements.txt` mirrors it for local
runs. Input datasets are pinned by immutable UUID in `.codeocean/datasets.json`.

**7. The traceability code is tested.**
`test_apply_decisions_filters_valid_query`, `test_apply_decisions_raises_on_invalid_query`,
and the config-merge tests mean the decision-application and config-resolution
paths are guarded against regression.

---

## Gaps that prevent *complete* traceability

**1. No determinism control — the highest-impact gap.**
No seed is set anywhere (`scvi.settings.seed`, `torch.manual_seed`,
`np.random.seed`, or Scanpy `random_state` for `sc.tl.umap` / `sc.tl.leiden`).
scVI/scANVI training, UMAP, and Leiden are all stochastic. Two runs of the same
`run_config_resolved.yml` will produce different latent spaces, different UMAPs,
and — critically — **different Leiden cluster numberings**. Because
`decisions.yaml` references clusters by ID (`remove_clusters: cluster "7"`), a
decision is only meaningful against the *specific* integrated object it was
authored on. The manifest records the config but not a seed, so the chain
"config → embedding → cluster IDs → decision" is not reproducible. *Fix: add a
configurable `seed` to the config, set all RNGs from it, and record it in the
manifest.*

**2. No code-version capture.**
`collect_package_versions` records dependency versions but nothing records the
*pipeline's own* git commit. The git history shows active refactoring
("refactor inspect to accept sweeps", architecture changes). A round produced
last week and one produced today can be byte-identical in config yet run
different code, with no recorded way to tell. *Fix: capture `git rev-parse HEAD`
(+ dirty flag) into `round_manifest.json`.*

**3. Input data identified by path, not content.**
The manifest stores `input_h5ad` as a filesystem path. Paths get reused and
overwritten (the `run` script writes `round_03/filtered.h5ad` and then reads it
back). There is no checksum/size/n_obs fingerprint of the input, so the record
cannot prove the input was the file the path now points to. The Code Ocean
dataset UUIDs help for *attached* datasets but not for round-to-round
intermediates. *Fix: hash inputs (or record shape + a cheap content hash) in the
manifest.*

**4. Cross-round lineage is implicit.**
Each round directory has its own manifest, but nothing inside `round_03` records
that its input came from `round_02/filtered.h5ad` as a structured parent
pointer. The lineage lives only in the (commented, hand-edited) `run` script and
in directory naming. *Fix: record an explicit `parent_round` / `input_provenance`
field.*

**5. Runtime environment capture is narrower than the pinned environment.**
`collect_package_versions([...])` records only `anndata`, `scanpy`, `scvi-tools`,
`scib-metrics`. Numerically influential packages — `torch`, `numpy`,
`scikit-learn`, `leidenalg`, `umap-learn` — are pinned in the Dockerfile but not
confirmed in the per-run manifest, so a run outside the image (the supported
local conda path) leaves them unrecorded.

**6. `command_args.json` and resolved config overwrite within a round.**
Both are keyed/written per stage, so re-running `inspect` twice in the same
output dir overwrites the earlier record. Repeated invocations within a round
leave no history — only the last one is documented.

**7. The human-facing report is decoupled from provenance.**
`inspection_report.html` surfaces only a timestamp and cell/gene/cluster counts
plus the auto-flag summary. It does not embed the config, package versions,
input file, seed, or the applied decisions — so the artifact a human actually
reads carries no provenance trail back to the machine-readable record.

**8. `decisions_applied.yaml` is hand-serialized.**
It is written with manual f-string formatting rather than `yaml.safe_dump`.
Reason/query strings containing quotes, colons, or newlines can produce invalid
or misleading YAML, and the writer is not symmetric with `load_decisions`. For a
file whose entire purpose is to be a trustworthy record, this is an integrity
risk. *Fix: serialize with the same YAML library used to read.*

**9. Rationale is encouraged but not enforced.**
`keep_clusters` falls back to `"keep listed clusters only"` and `remove_clusters`
to `"unspecified"` when no reason is given. A decision can therefore be recorded
with no actual justification. *Fix: optionally require non-empty reasons.*

**10. Minor:** timestamps are naive local time (`datetime.now()`, no timezone);
`metadata.yml` description is a single terse line; the `run` reproducible
entrypoint is mostly commented-out hardcoded paths rather than a turnkey
sequence; stdout carries an informative per-step narrative that is not persisted
to a log file.

---

## Priority recommendations (highest leverage first)

1. **Set and record a seed** for all stochastic steps — without this, nothing
   downstream is reproducible, which undercuts the whole decision trail.
2. **Capture the git commit (+ dirty flag) in the manifest** — pairs the recorded
   config with the actual code that ran.
3. **Fingerprint inputs** (hash or shape+hash) and **record an explicit
   parent-round pointer** — turns the implicit round chain into a verifiable
   lineage graph.
4. **Serialize `decisions_applied.yaml` with `yaml.safe_dump`** and broaden
   `collect_package_versions` to the full numerics stack.
5. **Embed a provenance footer in the HTML report** (config digest, seed, git
   commit, package versions, input file, applied decisions) so the
   human-readable artifact is self-documenting.

With items 1–3 addressed, this pipeline would move from "well-documented
decisions" to genuinely "complete, reproducible traceability."

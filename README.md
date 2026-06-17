# Config-Driven scVI/scANVI Integration Pipeline

Iterative single-cell integration for QC, scVI/scANVI modeling, annotation review,
filtering decisions, and re-integration.


```text
integrate -> inspect -> edit decisions -> filter -> re-integrate
```

The primary interface is a versioned YAML config for integration/inspection and `decisions.yaml` files for filtering decisions. CLI arguments are supported as overrides for common fields.

## Workflow

```mermaid
flowchart TD
    cfg["pipeline.yml<br/>(versioned config)"]:::cfg

    subgraph INT["integrate_sns_scvi.py"]
        direction TB
        qc["QC filter<br/>(skippable after round 1)"] --> hvg["HVG selection<br/>(batch-aware)"]
        hvg --> scvi["scVI training<br/>obsm['X_scVI']"]
        scvi -.->|annotation.enabled| scanvi["scANVI annotation<br/>obsm['X_scANVI']"]
        hvg -.->|harmony.enabled| harmony["Harmony<br/>obsm['X_pca_harmony']"]
        scvi --> emb["neighbors → UMAP → Leiden<br/>per embedding"]
        scanvi --> emb
        harmony --> emb
        emb --> bench["scIB benchmarking<br/>scVI vs Harmony vs PCA"]
    end

    INT --> integ["integrated.h5ad<br/>+ round_manifest.json<br/>+ scib_benchmark_results.csv"]:::art

    subgraph INS["inspect_integration.py"]
        direction TB
        rpt["cluster QC + markers<br/>+ diagnostic plots"] --> flags["auto_flags.yaml<br/>(suggested removals + reasons)"]
        flags --> html["inspection_report.html<br/>(+ Run Provenance footer)"]:::art
    end

    integ --> INS
    html --> review{"human review"}:::dec
    review --> dec["decisions.yaml<br/>(keep/remove + reasons)"]:::cfg

    dec --> filter["apply_decisions()"]
    filter --> filt["filtered.h5ad<br/>+ decisions_applied.yaml<br/>+ filtering_retention_summary.csv"]:::art
    filt -->|next round| INT

    cfg --> INT
    cfg --> INS

    classDef cfg fill:#e8f0fe,stroke:#4477AA;
    classDef art fill:#eef7ee,stroke:#449944;
    classDef dec fill:#fff3e0,stroke:#cc7a00;
```

The loop is `integrate → inspect → edit decisions → filter → re-integrate`. Every
stage writes provenance into `round_manifest.json` (seed, git commit, input
fingerprint, parent-round pointer), and the round's filtered output feeds the
next round's integration.

## Main Scripts

| Script | Role |
|---|---|
| `code/integrate_sns_scvi.py` | QC, HVG selection, scVI training, optional scANVI annotation, UMAP/leiden, benchmarking, integrated h5ad output |
| `code/inspect_integration.py` | Integration report, cluster QC, marker/annotation diagnostics, auto flags, decisions-based filtering |
| `code/pipeline_config.py` | Shared V1 config defaults, validation, and provenance helpers |
| `code/summarize_rounds.py` | Cross-round summary: combines every `round_manifest.json` into a table, decisions ledger, Mermaid lineage diagram, and (optional) integrity report |
| `code/04_visualize.py` | Legacy plotting script; label-transfer plots are now handled by inspection when annotation columns exist |

## Quick Start

Copy one of the example configs and edit paths/metadata keys:

```bash
cp examples/pipeline_scvi.yml pipeline.yml
```

Run the first integration round:

```bash
python code/integrate_sns_scvi.py --config pipeline.yml
```

Inspect the integrated object:

```bash
python code/inspect_integration.py --config pipeline.yml
```

Review `inspection_report.html`, edit `auto_flags.yaml` into `decisions.yaml`,
then export a filtered object for the next round:

```bash
python code/inspect_integration.py \
  --config pipeline.yml \
  --input rounds/round_01/integrated.h5ad \
  --decisions rounds/round_01/decisions.yaml \
  --output-dir rounds/round_02
```

For the next round, update `data.input_h5ad` to `rounds/round_02/filtered.h5ad`
and set `data.output_dir` to `rounds/round_02`, then re-run integration.

Or, perform round_1 integration and inspection, then generate the `decisions.yaml` and re-run inspection to filter. Export the filtered object to round_2 of integration and inspection:

```
python -u ./inspect_integration.py \
  --config ./pipeline_scvi.yml \
  --input /root/capsule/results/rounds/round_01/integrated.h5ad \
  --decisions /root/capsule/results/rounds/round_01/decisions.yaml \
  --output-dir /root/capsule/results/rounds/round_02

python -u ./integrate_sns_scvi.py \
  --config ./pipeline_scvi.yml \
  --input /root/capsule/results/rounds/round_02/filtered.h5ad \
  --output-dir /root/capsule/results/rounds/round_02

python -u ./inspect_integration.py \
  --config ./pipeline_scvi.yml \
  --input /root/capsule/results/rounds/round_02/integrated.h5ad \
  --output-dir /root/capsule/results/rounds/round_02
  ```


## Local Test Environment

For local integrated testing, use a Python 3.10 conda/mamba environment and the
repo-level `requirements.txt`. The requirements file mirrors the package pins in
`environment/Dockerfile` where practical, but uses CPU/Mac-friendly
`scvi-tools==1.3.3` instead of the Docker image's CUDA extra.

```bash
mamba create -n scvi-loops -c conda-forge python=3.10 pip scikit-misc -y
mamba activate scvi-loops
python -m pip install -U pip
python -m pip install -r requirements.txt
```

Run the lightweight checks:

```bash
python -m compileall -q code tests
python -m pytest -q
```

Run an integration pass with your edited config:

```bash
python code/integrate_sns_scvi.py --config pipeline.yml
```

Use the Docker environment for exact CUDA/container parity; the conda
environment is intended for local regression testing and reproducing Scanpy API
compatibility issues.

## Config Schema

Minimal scVI-only config:

```yaml
pipeline_version: 1

reproducibility:
  seed: 0

data:
  input_h5ad: data/combined_sns_adata.h5ad
  output_dir: rounds/round_01
  counts_layer: counts
  batch_key: data_origin
  categorical_covariate_keys: [tech]
  continuous_covariate_keys: []
  species: mouse
  gene_symbol_case: preserve

qc:
  enabled: true
  min_genes: 500
  min_cells: 3
  mt_gene_patterns: ["mt-"]
  ribo_gene_patterns: ["rps", "rpl"]
  hb_gene_pattern: "^hb[^(p)]"

integration:
  model_type: scvi
  n_hidden: 256
  n_layers: 3
  n_latent: 32
  dispersion: gene-cell
  gene_likelihood: nb
  max_epochs: 200
  early_stopping_patience: 20
  batch_size: 256
  num_workers: 4
  hvg:
    n_top_genes: 3000
    batch_key: tech
    flavor: seurat_v3
    min_batches: null
  sweep: []
  harmony:
    enabled: false       # run Harmony alongside scVI (single-model or sweep)
    batch_key: null      # defaults to data.batch_key when null
    n_pcs: 30            # PCA dimensionality Harmony corrects (>= 2)
    hvg_from: null       # sweep mode: name a sweep entry to borrow its HVGs;
                         # null = use integration.hvg (the default)

annotation:
  enabled: false
  method: scanvi
  labels_key: null
  unlabeled_category: Unknown
  prediction_key: scanvi_label
  confidence_key: scanvi_confidence
  min_confidence: 0.5

embedding:
  n_neighbors: 30
  umap_min_dist: 0.4
  umap_spread: 3.0
  leiden_resolution: 0.3

inspection:
  cluster_key: auto
  latent_key: auto
  umap_key: auto
  markers_json: null
  marker_threshold: 0.0
  neuronal_cutoff: 0.5
  auto_flag:
    mt_threshold: 15.0
    min_genes_threshold: 400
    min_cells: 20
    single_batch_threshold: 0.90

benchmark:
  enabled: true
  batch_key: auto
  label_key: null

decisions:
  ignore_failed_queries: false
```

See:

- `examples/pipeline_scvi.yml` for scVI-only integration.
- `examples/pipeline_scanvi.yml` for scVI pretraining plus scANVI labels.
- `examples/pipeline_sweep.yml` for config-driven model sweeps.

## scANVI Annotation

Enable scANVI by setting:

```yaml
annotation:
  enabled: true
  method: scanvi
  labels_key: seed_cell_type
  unlabeled_category: Unknown
  prediction_key: scanvi_label
  confidence_key: scanvi_confidence
  min_confidence: 0.5
```

The integration script trains scVI first, initializes scANVI from that model,
then writes:

| Slot | Contents |
|---|---|
| `.obsm["X_scANVI"]` | scANVI latent representation |
| `.obs["scanvi_label"]` | Predicted label, or configured `prediction_key` |
| `.obs["scanvi_confidence"]` | Maximum predicted class probability, or configured `confidence_key` |
| `.obsm["X_umap_scanvi"]` | UMAP from the scANVI latent |
| `.obs["leiden_scanvi"]` | Leiden clusters from the scANVI latent |

Inspection adds annotation UMAPs, confidence distribution, label composition by
cluster, and low-confidence cluster summaries when those columns are present.

## Harmony Integration (optional)

Set `integration.harmony.enabled: true` to run
[Harmony](https://github.com/slowkow/harmonypy) alongside scVI, as an additional
integration to compare against. It is computed on a standard log-normalized PCA
of the HVGs and batch-corrected with `scanpy.external.pp.harmony_integrate`.

| Slot | Contents |
|---|---|
| `.obsm["X_pca_harmony"]` | Harmony-corrected embedding |
| `.obsm["X_umap_harmony"]` | UMAP from the Harmony embedding |
| `.obs["leiden_harmony"]` | Leiden clusters from the Harmony embedding |

The Harmony embedding flows through the same downstream stages as scVI —
neighbors, UMAP, Leiden, and scIB benchmarking — so `scib_benchmark_results.csv`
reports a **same-run scVI vs Harmony vs PCA** comparison with no extra steps.
Inspection auto-resolves the Harmony keys, or pass them explicitly with
`--latent-key X_pca_harmony --cluster-key leiden_harmony --umap-key X_umap_harmony`.

**Runs in both single-model and sweep modes.** Harmony runs **once per round**,
independent of the scVI model(s) — so in a sweep, the single Harmony embedding is
benchmarked against *every* swept scVI architecture, which is often the most
useful comparison. Because a sweep has no single canonical HVG set, the HVG
source is an **explicit, recorded decision**: by default Harmony uses the
top-level `integration.hvg` spec (identical to single-model scVI); set
`harmony.hvg_from: <sweep_entry_name>` to borrow a specific swept architecture's
HVG selection instead. The resolved choice (`source`, batch key, flavor, HVG
count) is written to `round_manifest.json` under `integration.harmony`, so the
decision is auditable rather than implicit. `hvg_from` is only valid in sweep
mode.

`batch_key` defaults to `data.batch_key` when null; `n_pcs` (default 30) sets the
PCA dimensionality Harmony corrects. The configured `reproducibility.seed` is
threaded into both the PCA and Harmony's internal KMeans; harmonypy is only
partially seed-controllable, so results are deterministic on a fixed
platform/library set but not guaranteed bit-identical across them. Requires
`harmonypy` (pinned in the Dockerfile and `requirements.txt`).

## CLI Overrides

YAML is the source of truth, but common CLI overrides still work:

```bash
python code/integrate_sns_scvi.py \
  --config pipeline.yml \
  --input rounds/round_02/filtered.h5ad \
  --output-dir rounds/round_02 \
  --skip-benchmark
```

No-config legacy usage remains supported:

```bash
python code/integrate_sns_scvi.py \
  --input data/combined_sns_adata.h5ad \
  --output-dir rounds/round_01

python code/inspect_integration.py \
  --input rounds/round_01/integrated.h5ad \
  --output-dir rounds/round_01
```

## Decisions

`decisions.yaml` supports:

```yaml
keep_clusters:
  - cluster: "1"
    reason: "High-quality target cells"

keep_cells:
  - query: "n_genes_by_counts > 500"
    reason: "Minimum complexity"

remove_clusters:
  - cluster: "7"
    reason: "Low-quality cluster"

remove_cells:
  - query: "scanvi_confidence < 0.5"
    reason: "Low-confidence annotation"

notes: "Manual review after round 1."
```

Failed decision queries are errors by default. To keep the previous permissive
behavior, set:

```yaml
decisions:
  ignore_failed_queries: true
```

## Outputs and Provenance

Each run writes a cumulative round record:

| File | Description |
|---|---|
| `integrated.h5ad` | All genes, normalized/log `.X`, raw counts in configured counts layer, latent embeddings, UMAPs, clusters |
| `run_config_resolved.yml` | Fully merged config after defaults and CLI overrides |
| `command_args.json` | CLI arguments recorded by stage |
| `round_manifest.json` | Cumulative integration, inspection, and decisions metadata |
| `inspection_report.html` | Self-contained inspection report |
| `auto_flags.yaml` | Editable starting point for filtering decisions |
| `filtered.h5ad` | Filtered object produced by applying decisions |
| `filtering_retention_summary.csv` | Pre/post decision-filtering counts overall, by cluster, by batch, and by cluster x batch |

`cluster_qc_summary.csv` reports both `silhouette_latent` and
`silhouette_umap`. Auto-flagging uses `silhouette_latent` so UMAP visual
compactness remains a separate diagnostic rather than the filtering criterion.

### Reproducibility and lineage

Every `round_manifest.json` now records, in addition to the per-stage payloads:

- `reproducibility.seed` — the global RNG seed applied to the run (see below).
- `code` — the pipeline's own git revision: `commit`, `commit_short`, `branch`,
  `commit_time`, and a `dirty` flag set when uncommitted tracked changes were
  present. Fields are null when git metadata is unavailable.
- `input_fingerprint` (per stage) — resolved path, `size_bytes`, `mtime`, and a
  `sha256` content digest of the input h5ad, so a recorded input can be verified
  rather than only named.
- `parent_round` (per stage) — when the input lives beside another round's
  `round_manifest.json` (e.g. `round_02/filtered.h5ad`), a structured pointer to
  that manifest with its `commit`, `seed`, and timestamps. This makes the
  round-to-round chain an explicit lineage link rather than an inference from
  directory names.

The decisions stage additionally fingerprints the applied `decisions.yaml`
(`decisions_fingerprint`) so the exact filtering record is tied to its content.

**Seeding.** `reproducibility.seed` is applied before any stochastic step via
`scvi.settings.seed`, which seeds Python `random`, NumPy, and Torch (through
Lightning's `seed_everything`). This makes scVI/scANVI training, the neighbor
graph, UMAP, and Leiden clustering deterministic — and therefore makes the
Leiden **cluster IDs** that `decisions.yaml` files reference stable across
re-runs of the same config. Set `seed: null` to leave RNG state untouched for an
intentionally non-reproducible run; the manifest records `null` in that case.

**Provenance footer.** `inspection_report.html` ends with a **Run Provenance**
section that surfaces the same record written to `round_manifest.json` —
pipeline version, RNG seed, code commit (with a `+uncommitted changes` marker
when the working tree was dirty), input file path and SHA-256, resolved-config
SHA-256, the parent-round pointer, and key package versions. The human-readable
artifact is therefore self-documenting; the JSON manifest remains the full
machine-readable record.

**Invocation history.** `command_args.json` keeps an append-only `history` list
of every stage invocation (in order), in addition to `commands[stage]` (the most
recent invocation per stage). Re-running a stage in the same output directory no
longer overwrites the earlier record.

**Safe serialization.** `decisions_applied.yaml` is written with
`yaml.safe_dump`, so reason and query strings containing quotes, colons, or
newlines are escaped correctly and the file round-trips through `yaml.safe_load`.

For sweep outputs, inspect every complete architecture key triplet in one call.
If Harmony keys are present, this also creates an `inspect_harmony/` report:

```bash
python code/inspect_integration.py \
  --config pipeline.yml \
  --input rounds/round_02_sweep/integrated.h5ad \
  --output-dir rounds/round_02_sweep \
  --all-sweep-architectures
```

This writes one report directory per architecture or embedding, for example
`rounds/round_02_sweep/inspect_small_gene_nb/` and
`rounds/round_02_sweep/inspect_harmony/`. Filtering with `--decisions` still
uses a single selected cluster key and cannot be combined with
`--all-sweep-architectures`.

### Selecting which sweep variant to filter on

In a sweep, Leiden clusters are computed **independently per scVI architecture**
(`leiden_<variant>` on `X_scVI_<variant>`), so a filtering decision is only
meaningful against the specific variant its cluster IDs came from. To keep that
choice explicit and recorded rather than implied by directory placement:

- Each `auto_flags.yaml` carries an `integration:` block naming the variant
  and the `cluster_key`/`latent_key`/`umap_key` the flags were generated against.
  When you rename the file to `decisions.yaml`, that block travels with it, so
  `apply_decisions` filters on the **same** variant.
- The selection precedence is `--cluster-key` (CLI) > the decisions file's
  `integration:` block > config `inspection.cluster_key` > `auto`.
- In filter mode, if multiple variants are present and nothing pins one, the run
  **errors out** instead of silently guessing — declare the variant in
  `decisions.yaml` or pass `--cluster-key leiden_<variant>`.

The variant actually used is recorded in `round_manifest.json` under
`decisions.filtered_on` (`variant`, `cluster_key`, `latent_key`, `umap_key`,
and `source` = how it was chosen), so every round documents which architecture
its filtering was performed on.

## Summarizing a completed set of rounds

After running several rounds, combine every `round_manifest.json` into one view.
Either point at a parent directory of `round_*` subdirectories:

```bash
python code/summarize_rounds.py --rounds-dir results/rounds --verify
```
...or pass an explicit, ordered list of round directories. Use this when rounds
are archived under arbitrary names/locations — e.g., CodeOcean overwrites
`/results` each run, so rounds are copied out to an S3 mount as
`data/260603_cmg_round1`, `data/260605_cmg_round3`, etc.:

```bash
python code/summarize_rounds.py \
  --rounds data/260603_cmg_round1 data/260603_cmg_round2 data/260605_cmg_round3 \
  --output-dir results/summary --verify
```

In discovery mode rounds are ordered by their recorded lineage (the
`parent_round` pointers). With `--rounds`, the **supplied order is taken as the
lineage** whenever those pointers can't be resolved — their run-time absolute
paths no longer exist after archiving. Lineage-continuity verification still
works either way: it hashes each parent's `integrated.h5ad` *as it sits on disk
now* and compares it to the child's recorded input fingerprint. Each `--rounds`
item may be a round directory or a `round_manifest.json` path. Outputs go to the
output directory (default `<rounds-dir>/summary`, or `./round_summary` with
`--rounds`):

| Output | Contents |
|---|---|
| `pipeline_summary.md` / `.html` | Round table + Mermaid lineage + integrity findings |
| `rounds_table.csv` | One row per round (cells, clusters, seed, commit, filtered-on variant, scIB selected vs best) |
| `decisions_ledger.csv` | Every keep/remove action across all rounds, with the variant it was applied on |
| `lineage.mmd` | Raw Mermaid lineage diagram (embeddable in a README or PR) |
| `pipeline_summary.json` | Full merged record + verification verdicts |

`--verify` adds integrity checks: round-to-round input-hash continuity, seed
drift, code-commit/dirty-tree drift, package-version drift, and any filtering
whose variant was chosen by `auto` (ambiguous). `--strict` exits non-zero when a
check raises an error, so the summary can gate CI. The scIB columns compare the
**selected** variant's score against the round's **best** scoring variant, so a
round filtered on a non-top architecture is visible at a glance.

## Notes

- V1 supports `integration.model_type: scvi`; scANVI is an annotation stage on top
  of the single-model scVI path.
- scANVI plus sweep mode is intentionally not supported in V1 because selecting
  which sweep model should seed annotation needs an explicit design.
- Existing SNS marker defaults are now part of the resolved config/template
  rather than implicit inspection logic.

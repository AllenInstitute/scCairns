# Config-Driven scVI/scANVI Integration Pipeline

Iterative single-cell integration for QC, scVI/scANVI modeling, annotation review,
filtering decisions, and re-integration.

The active workflow is still script-based:

```text
integrate -> inspect -> edit decisions -> filter -> re-integrate
```

The primary interface is now a versioned YAML config. Existing CLI arguments are
still supported as overrides for common fields.

## Main Scripts

| Script | Role |
|---|---|
| `code/integrate_sns_scvi.py` | QC, HVG selection, scVI training, optional scANVI annotation, UMAP/leiden, benchmarking, integrated h5ad output |
| `code/inspect_integration.py` | Integration report, cluster QC, marker/annotation diagnostics, auto flags, decisions-based filtering |
| `code/pipeline_config.py` | Shared V1 config defaults, validation, and provenance helpers |
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

## Config Schema

Minimal scVI-only config:

```yaml
pipeline_version: 1

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

## Notes

- V1 supports `integration.model_type: scvi`; scANVI is an annotation stage on top
  of the single-model scVI path.
- scANVI plus sweep mode is intentionally not supported in V1 because selecting
  which sweep model should seed annotation needs an explicit design.
- Existing SNS marker defaults are now part of the resolved config/template
  rather than implicit inspection logic.

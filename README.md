# SNS scVI Integration Pipeline

Iterative scVI integration and refinement for sympathetic nervous system (SNS) single-cell/single-nucleus RNA-seq data across multiple technologies (10x Chromium, ScaleBio, SMART-seq).

## Overview

Two scripts form an iterative cycle:

```
integrate  →  inspect  →  (edit decisions)  →  filter  →  re-integrate  →  ...
```

| Script | Role |
|---|---|
| `integrate_sns_scvi.py` | scVI model training, batch integration, UMAP/leiden, scIB benchmarking |
| `inspect_integration.py` | Post-integration QC report, automated cluster flagging, decision-based filtering |

The integration script produces an integrated h5ad. The inspection script generates a diagnostic report, auto-flags problematic clusters, and — once you've reviewed and edited the flags — exports a filtered dataset for the next round. Each round lives in its own directory with full provenance.

## Requirements

```
scanpy>=1.10
anndata>=0.10
scvi-tools>=1.1
scib-metrics         # optional, for benchmarking
scikit-learn
matplotlib
seaborn
leidenalg
igraph
pyyaml               # optional, for decisions.yaml parsing (fallback parser included)
```

## Quick Start

```bash
# Round 1: integrate
python integrate_sns_scvi.py \
    --input data/combined_sns_adata.h5ad \
    --output-dir rounds/round_01/

# Round 1: inspect
python inspect_integration.py \
    --input rounds/round_01/integrated.h5ad \
    --output-dir rounds/round_01/

# Review rounds/round_01/inspection_report.html
# Edit rounds/round_01/auto_flags.yaml → save as decisions.yaml

# Round 2: apply decisions + export filtered data
python inspect_integration.py \
    --input rounds/round_01/integrated.h5ad \
    --decisions rounds/round_01/decisions.yaml \
    --output-dir rounds/round_02/

# Round 2: re-integrate on filtered data
python integrate_sns_scvi.py \
    --input rounds/round_02/filtered.h5ad \
    --output-dir rounds/round_02/ \
    --skip-qc-filter

# Round 2: inspect again
python inspect_integration.py \
    --input rounds/round_02/integrated.h5ad \
    --output-dir rounds/round_02/
```

## Directory Structure

After two rounds the project looks like:

```
project/
├── data/
│   └── combined_sns_adata.h5ad          # Pre-concatenated input
├── rounds/
│   ├── round_01/
│   │   ├── integrated.h5ad              # All genes, scVI latent + UMAP + leiden
│   │   ├── scvi_model_default/          # Saved scVI model
│   │   ├── hvg_genes.csv
│   │   ├── training_convergence.png
│   │   ├── umap_integration.png
│   │   ├── scib_benchmark_results.csv
│   │   ├── qc_violins_prefilter.png
│   │   ├── inspection_report.html       # ← review this
│   │   ├── cluster_qc_summary.csv
│   │   ├── auto_flags.yaml              # ← auto-generated
│   │   ├── decisions.yaml               # ← you create this
│   │   ├── qc_summary_heatmap.png
│   │   ├── batch_composition.png
│   │   ├── marker_dotplot.png
│   │   └── cluster_silhouettes.png
│   └── round_02/
│       ├── cells_to_keep.csv            # From decisions applied
│       ├── filtered.h5ad                # Filtered input for this round
│       ├── decisions_applied.yaml       # Provenance log
│       ├── round_manifest.json
│       ├── integrated.h5ad
│       └── ...
└── scripts/
    ├── integrate_sns_scvi.py
    └── inspect_integration.py
```

---

## `integrate_sns_scvi.py`

### Modes

**Single model (default).** Trains one scVI model with the specified architecture. Good for initial exploration or when you've already identified the best config via a sweep.

```bash
python integrate_sns_scvi.py --input combined.h5ad
```

**Parameter sweep (`--sweep`).** Trains multiple scVI models with different architectures, batch keys, and covariate settings. Generates a UMAP comparison grid and scIB benchmarks across all configs.

```bash
python integrate_sns_scvi.py --input combined.h5ad --sweep
```

### Default Architecture

Matches the "large" config from the combined SNS notebook:

| Parameter | Default |
|---|---|
| `n_hidden` | 256 |
| `n_layers` | 3 |
| `n_latent` | 32 |
| `dispersion` | gene-cell |
| `gene_likelihood` | nb |
| `batch_key` | data_origin |
| `covariate_keys` | tech |

### Built-in Sweep Configs

From the 20260325 sympathetic refinement notebook:

| Config | n_hidden | batch_key | Covariates | HVG flavor |
|---|---|---|---|---|
| `custom2_tech` | 16 | tech | — | seurat, ≥3 batches |
| `custom2_platform` | 16 | platform_origin | — | seurat, ≥3 batches |
| `medium_platform` | 128 | platform_origin | — | seurat, ≥3 batches |
| `medium_platform_cov` | 128 | platform_origin | data_origin | seurat, ≥3 batches |

Override with `--sweep-configs my_configs.json` (see `--help` for format).

### scIB Benchmarking

Enabled by default. Since cell type annotations are not available, only batch-correction metrics are computed:

| Metric | What it measures |
|---|---|
| Silhouette (batch) | Whether batches overlap within clusters |
| Graph connectivity | Whether the kNN graph connects cells across batches |
| PCR comparison | Reduction in batch-driven variance vs. PCA baseline |

Pass `--bench-label-key ganglion_group` to additionally enable bio-conservation metrics (silhouette label, cLISI) using a proxy label.

### Key Output: `integrated.h5ad`

Contains **all genes** (not just HVGs):

| Slot | Contents |
|---|---|
| `.X` | Normalised, log1p |
| `.layers["counts"]` | Raw integer counts |
| `.obsm["X_scVI"]` | Latent representation (or `X_scVI_{config}` in sweep mode) |
| `.obsm["X_umap"]` | UMAP coordinates |
| `.obs["leiden"]` | Leiden cluster assignments |
| `.var["highly_variable"]` | Boolean HVG mask |

### Full CLI Reference

```
python integrate_sns_scvi.py --help
```

<details>
<summary>Key arguments</summary>

| Argument | Default | Description |
|---|---|---|
| `--input` | (required) | Path to pre-concatenated h5ad |
| `--output-dir` | `../results` | Output directory |
| `--batch-key` | `data_origin` | Batch key for scVI |
| `--covariate-keys` | `tech` | Categorical covariates |
| `--n-hidden` | 256 | Hidden layer size |
| `--n-layers` | 3 | Number of layers |
| `--n-latent` | 32 | Latent dimensions |
| `--gene-likelihood` | `nb` | `nb`, `zinb`, or `poisson` |
| `--n-hvgs` | 3000 | Number of HVGs |
| `--hvg-batch-key` | `tech` | Batch key for HVG selection |
| `--hvg-flavor` | `seurat_v3` | HVG method |
| `--max-epochs` | 200 | Training epochs |
| `--sweep` | off | Enable parameter sweep |
| `--sweep-configs` | built-in | JSON file with sweep configs |
| `--n-neighbors` | 30 | Neighbors for UMAP/leiden |
| `--leiden-resolution` | 0.3 | Leiden resolution |
| `--bench-batch-key` | (same as batch-key) | scIB batch key |
| `--bench-label-key` | None | scIB label key (optional proxy) |
| `--skip-qc-filter` | off | Skip cell/gene filtering |
| `--skip-benchmark` | off | Skip scIB benchmarking |

</details>

---

## `inspect_integration.py`

### Report Mode (default)

Generates a diagnostic report without modifying data.

```bash
python inspect_integration.py \
    --input integrated.h5ad \
    --output-dir rounds/round_01/
```

**Outputs:**

| File | Description |
|---|---|
| `inspection_report.html` | Self-contained HTML with all plots and tables |
| `cluster_qc_summary.csv` | Per-cluster: N cells, median genes/counts/%MT, batch composition, silhouette |
| `auto_flags.yaml` | Template of clusters that trip automated thresholds |
| `umap_overview.png` | 2-panel UMAP (clusters + batch) |
| `qc_summary_heatmap.png` | Heatmap of QC metrics, red borders on flagged clusters |
| `batch_composition.png` | Stacked bar of batch fractions per cluster |
| `marker_dotplot.png` | Marker gene expression per cluster |
| `cluster_silhouettes.png` | Per-cluster silhouette scores in latent space |

### Filter Mode (`--decisions`)

Applies a `decisions.yaml` file and exports filtered data.

```bash
python inspect_integration.py \
    --input integrated.h5ad \
    --decisions decisions.yaml \
    --output-dir rounds/round_02/
```

**Outputs:**

| File | Description |
|---|---|
| `cells_to_keep.csv` | Obs names of retained cells |
| `filtered.h5ad` | Filtered AnnData for re-integration |
| `decisions_applied.yaml` | Provenance: what was removed, why, how many |
| `round_manifest.json` | Machine-readable round summary |

### Auto-flagging Thresholds

| Criterion | Default | CLI flag |
|---|---|---|
| Median %MT | > 15% | `--mt-threshold` |
| Median genes | < 400 | `--min-genes-threshold` |
| Cluster size | < 20 cells | `--min-cells` |
| Dominant batch fraction | > 90% | `--single-batch-threshold` |
| Median silhouette | < −0.05 | (computed from latent) |

### Built-in Marker Genes

Default panel covers peripheral / sympathetic neurons:

- **Pan-neuronal:** Snap25, Tubb3, Rbfox3, Elavl4, Isl1
- **Noradrenergic:** Th, Dbh, Ddc, Slc6a2
- **Cholinergic:** Chat, Slc18a3, Slc5a7
- **Glutamatergic:** Slc17a6, Slc17a7
- **GABAergic:** Slc32a1, Gad1, Gad2
- **Nitrergic:** Nos1
- **Neuropeptides:** Npy, Sst, Vip, Pdyn
- **Transcription factors:** Phox2b, Phox2a, Sox6, Shox2
- **Satellite glia:** Sox10, Fabp7, S100b
- **Immune:** Ptprc, Cd68
- **Mitochondrial:** mt-Co1, mt-Co2, mt-Cytb

Override with `--markers my_markers.json`:

```json
{
  "My group": ["Gene1", "Gene2"],
  "Another group": ["Gene3", "Gene4"]
}
```

---

## `decisions.yaml` Format

The `auto_flags.yaml` output is a starting template. Copy/rename it to `decisions.yaml` and edit to reflect your actual decisions:

```yaml
remove_clusters:
  - cluster: "10"
    reasons:
      - "MT-enriched cluster (median 23.1%)"
  - cluster: "14"
    reasons:
      - "Very small, single-batch (ScaleBio 95%)"

remove_cells:
  - query: "(tech == 'scale') & (neuron_type.isin(['GABAergic', 'Glutamatergic']))"
    reason: "Contaminating non-noradrenergic neurons from ScaleBio"

notes: "Round 1 cleanup — removed dying cells and off-target neurons"
```

`remove_clusters` drops all cells in the named cluster. `remove_cells` uses pandas `.eval()` syntax on `.obs` for finer-grained filtering. Both are logged in `decisions_applied.yaml`.

---

## Tips

- **First round: skip the sweep.** The default large architecture is a good starting point for initial QC. Run a sweep after you've cleaned the data.
- **`--skip-qc-filter` on re-integration.** Filtering was already applied by the inspection script; don't double-filter.
- **Compare scIB scores across rounds.** The `scib_benchmark_results.csv` from each round lets you track whether filtering improved integration.
- **Custom sweep configs.** If the built-in configs don't match your data (different batch keys, platform columns), create a JSON file and pass `--sweep-configs`.
- **GPU.** scVI will automatically use a GPU if available. For large datasets (>50k cells), this significantly reduces sweep runtime.

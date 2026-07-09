# Getting started

This page takes you from a fresh clone to a verified install and explains the one
thing that trips up every new user: **what your input `.h5ad` must contain**.

New to the whole workflow? Read [the tutorial](tutorial.md) next — it runs the full
loop end-to-end on a synthetic dataset you generate locally, no real data required.

## 1. Prerequisites

- **Python 3.10** (the pins in `requirements.txt` and the Docker image target 3.10).
- **conda / mamba** for the environment (recommended), or any Python 3.10 venv.
- A C toolchain and **HDF5** headers are pulled in automatically by the pinned
  wheels on macOS/Linux; you rarely need to install them by hand.
- **GPU is optional.** Local runs use CPU-friendly `scvi-tools==1.3.3`. The
  [Docker image](codeocean.md) uses the CUDA extra for GPU parity — you only need
  that for large datasets or exact container reproducibility.

## 2. Install

```bash
git clone <this-repo> scvi_integration_loops
cd scvi_integration_loops

mamba create -n scvi-loops -c conda-forge python=3.10 pip scikit-misc -y
mamba activate scvi-loops
python -m pip install -U pip
python -m pip install -r requirements.txt
```

`harmonypy` is included in `requirements.txt` and is only needed if you enable the
optional [Harmony comparison](configuration.md#harmony-optional). Everything else
works without it.

## 3. Verify the install

```bash
python -m compileall -q code tests   # all modules import/compile
python -m pytest -q                   # the test suite passes
```

If both succeed, you're ready. If `pytest` reports import errors, your environment
is missing a dependency — re-run the `pip install` step inside the activated env.

## 4. The input data contract

Every run starts from **one pre-concatenated `.h5ad`** that combines all the samples
you want to integrate. The pipeline expects:

| Requirement | Detail |
|---|---|
| **Raw counts** | Integer counts in `.X` **or** in `.layers["counts"]`. The pipeline copies `.X → layers["counts"]` if no counts layer exists, then normalizes/log1p's internally. Do **not** pre-normalize `.X`. |
| **Batch column** | An `obs` column naming the batch to integrate over (default `data_origin`; set `data.batch_key`). This is the variable scVI corrects. |
| **Covariate columns** (optional) | Any `obs` columns listed in `data.categorical_covariate_keys` (default `[tech]`) / `continuous_covariate_keys` must exist. |
| **Label column** (only for scANVI) | If `annotation.enabled: true`, `annotation.labels_key` must point at an `obs` column of seed labels. |
| **Gene symbols** | `var_names` are gene symbols. QC mitochondrial/ribosomal/hemoglobin detection is **regex on `var_names`** (default patterns are mouse — see [adapting to your data](adapting-to-your-data.md)). |

### Building the combined input from per-sample files

There's no special format — just a standard AnnData concatenation. For example:

```python
import anndata as ad
import scanpy as sc

samples = {
    "sampleA": sc.read_10x_mtx("sampleA/"),
    "sampleB": sc.read_10x_mtx("sampleB/"),
}
for name, adata in samples.items():
    adata.obs["data_origin"] = name          # becomes batch_key
    adata.obs["tech"] = "10x_v3"              # a covariate, if useful

combined = ad.concat(samples.values(), join="outer", label=None)
combined.obs_names_make_unique()
combined.write_h5ad("data/my_combined.h5ad")  # raw counts in .X
```

Point `data.input_h5ad` at that file (or pass `--input`) and you're ready to run.

## 5. Two ways to run

**Config-driven (recommended).** A versioned YAML config is the source of truth;
CLI flags override individual fields. Copy an example and edit paths:

```bash
cp examples/pipeline_scvi.yml pipeline.yml   # mouse/SNS defaults
# or: cp examples/pipeline_generic.yml pipeline.yml   # non-SNS template
python code/integrate_scvi.py --config pipeline.yml
python code/inspect_integration.py --config pipeline.yml
```

**No-config legacy.** Every field has a default, so you can run straight off `--input`:

```bash
python code/integrate_scvi.py --input data/my_combined.h5ad --output-dir rounds/round_01
```

## Where to go next

- **[Tutorial](tutorial.md)** — do a complete, runnable integrate → filter → re-integrate loop.
- **[Configuration reference](configuration.md)** — every field, and *when to change it*.
- **[Integration methods](integration-methods.md)** — tune the scVI architecture or swap in another method.
- **[Adapting to your data](adapting-to-your-data.md)** — running on non-mouse / non-SNS datasets.
- **[Interpreting outputs](interpreting-outputs.md)** — is my integration any good?

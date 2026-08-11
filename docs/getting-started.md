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

Pick the mode that matches what you're doing. **Using** scCairns to analyze a dataset
does not require a clone — install a pinned release as a dependency and keep your
project's config and data in your own repo. Clone only if you intend to change
scCairns itself.

### A. As a pinned dependency (analysis projects)

```bash
mamba create -n myproject -c conda-forge python=3.10 pip scikit-misc -y
mamba activate myproject
python -m pip install -U pip
python -m pip install "git+https://github.com/AllenInstitute/scCairns.git@v0.1.0"
```

**Always pin an exact tag.** Installing from `@main` means a later rebuild of the same
environment silently gives you different code, and rounds produced by the old and new
builds are no longer comparable. The tag you pin is what every `round_manifest.json`
records as `code.version`, so it is the thing that makes a round reproducible.

The repository is private (Allen-internal), so pip needs GitHub credentials. Use a
credential helper or a PAT for the HTTPS form above, or SSH if that's how you
authenticate:

```bash
python -m pip install "git+ssh://git@github.com/AllenInstitute/scCairns.git@v0.1.0"
```

One caveat: the package's own dependencies are deliberately **unpinned** in
`pyproject.toml`, so that installing into an existing environment doesn't disturb a
working scientific stack. The reproducible pins live in `requirements.txt`. If you want
the exact stack this release was tested against, install those first — from a clone or
checkout of the same tag — and then the package:

```bash
python -m pip install -r requirements.txt   # pinned scientific stack (scanpy, scvi-tools, …)
python -m pip install "git+https://github.com/AllenInstitute/scCairns.git@v0.1.0"
```

For Code Ocean, the install belongs in the capsule's **`postInstall`**, not in `run` —
image builds have network access, reproducible runs do not. See
[Code Ocean](codeocean.md).

### B. From a clone (developing scCairns)

```bash
git clone https://github.com/AllenInstitute/scCairns.git
cd scCairns

mamba create -n sccairns -c conda-forge python=3.10 pip scikit-misc -y
mamba activate sccairns
python -m pip install -U pip
python -m pip install -r requirements.txt   # pinned scientific stack
python -m pip install -e .                  # the sccairns package + `cairns` CLI
```

The editable install puts the `cairns` command on your PATH and keeps your edits live.
The `code/*.py` scripts remain as backward-compat shims for the Code Ocean capsule.

### Harmony (optional, either mode)

Harmony is an **optional** integration path, so `harmonypy` is not a core dependency
of the package. It is pinned in `requirements.txt` (`harmonypy==0.0.10`), and you can
also install it via the extra: `python -m pip install -e '.[harmony]'` (from a clone) or
`python -m pip install "scCairns[harmony] @ git+https://github.com/AllenInstitute/scCairns.git@v0.1.0"`.
Only needed if you enable the optional
[Harmony comparison](configuration.md#harmony-optional); everything else works without
it. (Note: `harmonypy` 2.x is incompatible — the extra caps it `<1`.)

## 3. Verify the install

Either mode:

```bash
cairns --help       # the CLI is on your PATH
cairns --version    # which scCairns is actually running
```

`cairns --version` reports the same identity that gets stamped into every round:

```text
cairns 0.1.0 (installed)                                  # a pinned release
cairns 0.1.0 (checkout 76b1acde76a5 in scCairns, dirty)   # an editable clone
```

`(installed)` means the version is the only identifier available — there is no git
metadata in an installed distribution — which is exactly why the pin matters. A
`dirty` checkout means uncommitted changes are in play, and rounds produced from it
are flagged by `cairns summarize --verify`.

From a clone, also run the suite:

```bash
python -m compileall -q code tests   # all modules import/compile
python -m pytest -q                   # the test suite passes
```

If `pytest` reports import errors, your environment is missing a dependency — re-run
the `pip install` step inside the activated env.

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
cairns integrate --config pipeline.yml
cairns inspect --config pipeline.yml
```

**No-config legacy.** Every field has a default, so you can run straight off `--input`:

```bash
cairns integrate --input data/my_combined.h5ad --output-dir rounds/round_01
```

## Where to go next

- **[Tutorial](tutorial.md)** — do a complete, runnable integrate → filter → re-integrate loop.
- **[Configuration reference](configuration.md)** — every field, and *when to change it*.
- **[Integration methods](integration-methods.md)** — tune the scVI architecture or swap in another method.
- **[Adapting to your data](adapting-to-your-data.md)** — running on non-mouse / non-SNS datasets.
- **[Interpreting outputs](interpreting-outputs.md)** — is my integration any good?

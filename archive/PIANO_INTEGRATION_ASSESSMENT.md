> **📦 Archived — point-in-time assessment.** Feasibility review of adding
> [PIANO](https://github.com/NingWang1729/piano) as an integration method in scCairns.
> **Decision: not implemented at this time.** Facts below (PIANO's requirements, the
> pinned environment, the code seam) were true on the date given and will drift — re-check
> before acting. Assessed 2026-08-14.

# Adding PIANO as an integration option — feasibility

**Question asked:** how hard would it be to add a PIANO integration module/option?

**Answer:** the code is roughly a day. The environment is a multi-day migration and is a
hard blocker today. **Decision: deferred.** Evaluate PIANO out-of-band first (see
"Zero-cost path" below); revisit a native module only if it beats scVI on our data.

---

## What PIANO is

Probabilistic Inference Autoencoder Networks for multi-Omics — a PyTorch/Pyro
variational autoencoder for single-cell batch integration, AnnData-native, producing a
latent representation and optionally batch-corrected counts. Structurally the same shape
as scVI, which is why the code side is cheap.

| | |
|---|---|
| Source | https://github.com/NingWang1729/piano |
| PyPI | `piano-integration`, latest **0.1.9** (2026-07-23), 19 releases |
| License | **GPL-3.0-or-later** |
| Status | bioRxiv preprint; not peer-reviewed as of this assessment |
| Python | **`>=3.11`** |
| Core deps | `numpy>=2`, `torch>=2.2,<2.8`, `pyro-ppl`, `scanpy`, `anndata`, `scipy`, `pandas`, `fast-array-utils`, `scikit-misc` |
| Extras | `rapids`, `scvi-tools`, `torch`, `misc`, `all` |

API shape (maps almost 1:1 onto our existing config fields):

```python
from piano import Composer

pianist = Composer(adata=adata, batch_key=..., categorical_covariate_keys=[...],
                   n_top_genes=4096, latent_size=32, max_epochs=200, outdir=...)
pianist.run_pipeline()
latent = pianist.get_latent_representation()
corrected = pianist.get_counterfactual(covariates="marginal")   # optional
```

## The blocker: environment, not code

| Requirement | PIANO | scCairns today |
|---|---|---|
| Python | `>=3.11` | `requires-python = ">=3.10"`; **Code Ocean image is `python3.10.12`** |
| numpy | `>=2` | **`numpy==1.26.4`** pinned in `environment/Dockerfile` |

PIANO **cannot be installed into the current capsule image at all**. Installing it locally
forces a numpy major upgrade through the pinned, validated stack (scanpy 1.10.4,
scvi-tools 1.3.3, torch 2.1), which would invalidate the reproducibility of every round
produced on the current image. Adopting PIANO therefore means rebuilding the environment
on a 3.11+ base, re-pinning the whole scientific stack, and re-validating — days of work
and a new Code Ocean image, independent of any scCairns code.

## The code seam is small (~1 day)

Harmony is the working precedent: ~40 lines in `integrate.py` plus a `run_harmony()`
helper, a config block, and an optional dependency extra. The contract any integration
method must satisfy is only:

```python
adata_full.obsm["X_piano"] = latent          # (n_cells, d)
embedding_keys.append("X_piano")
embedding_key_by_name["piano"] = "X_piano"
```

Everything downstream iterates `embedding_key_by_name` generically, so a registered
embedding gets UMAP, Leiden, scIB benchmarking against scVI/Harmony/PCA, the inspection
report, decisions/filtering, `summarize`, and the neighbor-purity cluster flag for free.

Work items, if it were built:

1. `run_piano()` in `integrate.py`, mirroring `run_harmony()` (including the
   `ImportError` with an actionable install message).
2. `integration.piano` config block — `enabled`, `batch_key`, `latent_size`,
   `max_epochs`, `n_top_genes`, `adversarial`, `kld_weight` bounds.
3. Provenance dict alongside `harmony_provenance` in the round manifest.
4. Optional extra `piano = ["piano-integration==<pin>"]`, following the `harmony` extra.
5. Tests + a section in `docs/integration-methods.md`.

**`integration.model_type` is not involved.** `validate_config` still gates
`model_type != "scvi"`, but the Harmony-style *parallel arm* pattern
(`integration.piano.enabled: true`) sidesteps it — and is the better design anyway,
since the value of PIANO here is as a benchmarked comparison arm, not a replacement.

## Non-technical risks

- **License.** GPL-3.0-or-later against our Allen Institute license (BSD-2 + a
  non-commercial clause). A hard dependency on GPL code in something we distribute (pinned
  git installs today, possibly PyPI later) is a copyleft question that needs an
  institutional answer, not just a packaging decision. The `harmonypy` optional-extra
  pattern keeps it out of a base install, which mitigates but does not settle it.
- **Maturity.** 0.1.9, fast-moving, preprint-stage. Pin exactly; do not make it a default.

## Zero-cost path (available today, no code changes)

Run PIANO in its own Python 3.11 environment, write the latent and a UMAP/Leiden into the
h5ad, then use the bring-your-own-embedding route already documented for Scanorama in
[docs/integration-methods.md](../docs/integration-methods.md):

```bash
cairns inspect --config pipeline.yml \
  --input rounds/round_01/integrated.h5ad \
  --latent-key X_piano --umap-key X_umap_piano --cluster-key leiden_piano \
  --output-dir rounds/round_01
```

This keeps the whole inspection and filtering loop. The only thing given up is having
scIB benchmark PIANO head-to-head **in the same run** — which is exactly what a native
module would buy, and exactly what is worth paying the environment cost for *once PIANO
has demonstrated an advantage on our data*.

## Revisit when

- We migrate to Python 3.11+ / numpy 2 for other reasons (then the code cost is ~1 day).
- PIANO shows a real advantage over scVI on our data in the out-of-band evaluation.
- The license question gets an institutional answer.

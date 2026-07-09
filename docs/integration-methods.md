# Changing the model architecture or the integration method

Two different asks are covered here:

1. **[Tune the scVI architecture](#1-tune-the-scvi-architecture)** — deeper/wider network,
   different latent size, different likelihood. Pure config, no code.
2. **[Use a different integration method](#2-use-a-different-integration-method)** —
   e.g. Scanorama, fastMNN, scVI variants, or anything that produces a cell embedding.
   Three levels, from no-code to first-class.

---

## 1. Tune the scVI architecture

The scVI network is fully config-driven under `integration:` — no code change needed.
The knobs and *when to turn them* are in the
[configuration reference](configuration.md#integration); the essentials:

```yaml
integration:
  n_hidden: 256        # layer width  (64–128 small data, 256+ large/complex)
  n_layers: 3          # depth        (1–2 small, 2–3 typical)
  n_latent: 32         # latent dim   (10–15 simple, 20–32 heterogeneous atlases)
  dispersion: gene-cell        # gene | gene-batch | gene-label | gene-cell
  gene_likelihood: nb          # nb | zinb | poisson
  max_epochs: 200
  hvg:
    n_top_genes: 3000
```

Any of these can also be overridden per-run on the CLI (they win over the YAML):

```bash
python code/integrate_scvi.py --config pipeline.yml \
  --n-hidden 128 --n-layers 2 --n-latent 24 --gene-likelihood zinb
```

### Comparing architectures in one run (sweeps)

To evaluate several architectures side by side without hand-running each, list them
under `integration.sweep`. Each entry trains independently and produces its own
embedding, UMAP, and clustering (`X_scVI_<name>`, `X_umap_<name>`, `leiden_<name>`),
and all are scored together by scib-metrics:

```yaml
integration:
  sweep:
    - name: shallow
      n_hidden: 64
      n_layers: 1
      n_latent: 10
    - name: deep
      n_hidden: 256
      n_layers: 3
      n_latent: 32
```

See [`examples/pipeline_sweep.yml`](../examples/pipeline_sweep.yml) and
[Sweeps](configuration.md#sweeps). When you then filter, pin which variant your cluster
IDs came from — see [decisions.md](decisions.md#sweeps-pinning-the-variant).

---

## 2. Use a different integration method

**What's supported out of the box:** scVI (optionally with scANVI annotation) as the
integration model, plus **Harmony** as an optional same-run baseline you can toggle on
(`integration.harmony.enabled: true`) and compare against scVI and a PCA baseline via
scib-metrics. The config validator currently accepts only `integration.model_type:
scvi`, so a *different* method (Scanorama, fastMNN, scVI derivatives, …) is added one of
three ways, in increasing effort:

### Level A — Compare against Harmony today (no code)

If all you want is "how does scVI compare to a linear batch-correction baseline," it's
already there:

```yaml
integration:
  harmony:
    enabled: true      # requires the harmonypy package
    n_pcs: 30
```

Harmony runs once per round on a log-normalized PCA of the HVGs and flows through the
same UMAP/Leiden/benchmark stages, producing `X_pca_harmony`, `X_umap_harmony`,
`leiden_harmony`. `scib_benchmark_results.csv` then reports **scVI vs Harmony vs PCA**
in one table. See [Harmony](configuration.md#harmony-optional).

### Level B — Bring your own embedding (no changes to this pipeline)

The **inspect → decide → filter** half of the loop is method-agnostic: it works on any
embedding you put in the object, because you can point it at arbitrary keys with
`--latent-key`, `--umap-key`, and `--cluster-key`. So you can integrate with *any* tool
externally, then use this pipeline for the QC/inspection/filtering loop.

The object just needs, by convention:

- the embedding in `obsm["X_<method>"]`,
- a neighbors graph + UMAP in `obsm["X_umap_<method>"]`,
- Leiden clusters in `obs["leiden_<method>"]` (filtering by cluster needs cluster IDs).

Example with **Scanorama** (install `scanorama`; `scanpy.external` wraps it):

```python
import scanpy as sc
import scanpy.external as sce

adata = sc.read_h5ad("data/my_combined.h5ad")     # raw counts in .X

# standard preprocessing + PCA that Scanorama corrects
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=30)

# Scanorama expects cells grouped by batch; sort first
adata = adata[adata.obs.sort_values("data_origin").index].copy()
sce.pp.scanorama_integrate(adata, key="data_origin",
                           basis="X_pca", adjusted_basis="X_scanorama")

# make the embedding usable by the inspection loop (same convention as scVI)
sc.pp.neighbors(adata, use_rep="X_scanorama", key_added="neighbors_scanorama")
sc.tl.umap(adata, neighbors_key="neighbors_scanorama")
adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"]
sc.tl.leiden(adata, neighbors_key="neighbors_scanorama",
             key_added="leiden_scanorama", resolution=0.3)

adata.write_h5ad("rounds/round_01/integrated.h5ad")
```

Then inspect and filter against those keys — the rest of the workflow is unchanged:

```bash
python code/inspect_integration.py --config pipeline.yml \
  --input rounds/round_01/integrated.h5ad \
  --latent-key X_scanorama \
  --umap-key   X_umap_scanorama \
  --cluster-key leiden_scanorama \
  --output-dir rounds/round_01
```

Trade-off: you skip this pipeline's scVI training and its built-in scib benchmarking of
the new method, but you keep the full reproducible inspection/filtering loop. Contamination
flagging still works (it reads log-normalized `.X`, independent of the embedding).

### Level C — Add it as a first-class method (small code change)

To make a method a config toggle like Harmony — trained inside the pipeline, benchmarked
automatically, recorded in the manifest — follow the **Harmony template**, which is the
model to copy. Three edits in `code/integrate_scvi.py`:

1. **Write a `run_<method>()` function** that returns the corrected embedding, mirroring
   [`run_harmony()`](../code/integrate_scvi.py) (around line 394): copy the AnnData, build
   its input from `layers[counts_layer]`, run the method, return the coordinates. Raise a
   clear `ImportError` if the method's package is missing.

2. **Add a config block** (e.g. `integration.scanorama: {enabled: false, ...}`) in
   `DEFAULT_CONFIG` (`code/pipeline_config.py`) and, if it needs validation, extend
   `validate_config()`.

3. **Register the embedding** in `main()`, right after the Harmony block (around line
   1529), using the same two structures the downstream loop reads:

   ```python
   if scanorama_cfg.get("enabled"):
       emb = run_scanorama(adata_scan, batch_key=..., counts_layer=counts_layer, seed=seed)
       adata_full.obsm["X_scanorama"] = emb
       embedding_keys.append("X_scanorama")           # -> included in scib benchmark
       embedding_key_by_name["scanorama"] = "X_scanorama"  # -> gets neighbors/UMAP/Leiden
   ```

Because the round-6 loop iterates `embedding_key_by_name`, registering the embedding is
all it takes — it automatically gets `X_umap_scanorama`, `leiden_scanorama`, and a row in
`scib_benchmark_results.csv`, exactly like scVI and Harmony. Record any method parameters
into `round_manifest.json` (as the Harmony block does with `harmony_provenance`) to keep
the run traceable.

> Note: swapping the *primary* `integration.model_type` (rather than adding an alongside
> embedding) additionally means relaxing the `model_type == "scvi"` check in
> `validate_config()` and giving the new model its own train path. Adding it as an
> alongside embedding (above) is the lower-risk pattern and matches how Harmony is wired.

"""scCairns — reproducible single-cell integration with a documented decision trail.

Iterative scVI/scANVI (and Harmony/Scanorama) integration wrapped in a
reproducibility and decision-documentation layer: per-round manifests, a
decisions ledger, and cross-round lineage summaries.

Stage modules (also exposed as ``cairns`` subcommands):
  - :mod:`sccairns.integrate`     — QC, HVG, training, UMAP/Leiden, benchmarking
  - :mod:`sccairns.inspect`       — inspection report, auto-flags, decision filtering
  - :mod:`sccairns.summarize`     — cross-round summary, ledger, lineage, integrity
  - :mod:`sccairns.contamination` — marker-based per-cell contamination flagging
  - :mod:`sccairns.config`        — config load/merge/validate + provenance helpers
"""

__version__ = "0.2.0"

> **📦 Archived — point-in-time assessment.** Feasibility review of adding Scrublet-based
> doublet detection to the QC stage. **Decision: not implemented at this time.** The
> environment facts (what scanpy ships, what is pinned) were true on the date given and
> will drift — re-check before acting. Assessed 2026-08-17.

# Adding doublet detection (Scrublet) to the QC pipeline — feasibility

**Question asked:** how hard would it be to add Scrublet to the QC pipeline?

**Answer:** about a day of code and a one-line dependency addition. Unlike the
[PIANO assessment](PIANO_INTEGRATION_ASSESSMENT.md), there is no environment migration
blocking it. **Decision: deferred** — the plumbing is easy, but the design questions
below (which platforms, which batch key, which round) need answers first, and doublet
scores are currently produced upstream in R by DoubletFinder.

---

## Dependency situation: better than expected, but not free

`sc.pp.scrublet` **ships inside scanpy 1.10.4** — scanpy vendored the algorithm, so the
standalone `scrublet` package is *not* needed.

However it fails immediately in the current environment:

```
FAILED: ModuleNotFoundError: No module named 'skimage'
```

scanpy's implementation imports `skimage.filters.threshold_minimum` for automatic
doublet-score thresholding, and **`scikit-image` is not pinned in `requirements.txt`,
`environment/Dockerfile`, or `pyproject.toml`**.

| | |
|---|---|
| New algorithm dependency | none — `sc.pp.scrublet` is in scanpy 1.10.4 |
| Actually required | `scikit-image` (BSD-3), currently absent |
| Pulls in | networkx, imageio, tifffile, lazy_loader, pillow |
| Python / numpy constraints | none that conflict with the pinned 3.10 / numpy 1.26.4 stack |
| Cost | one line in `requirements.txt` + the Dockerfile pip block, then a Code Ocean image rebuild |

No compilation, no CUDA, no version wall. Contrast with PIANO, which needs Python ≥3.11
and numpy ≥2 and therefore a full re-pin.

## The filtering seam already exists — removal needs no new code

Scrublet writes `doublet_score` (float) and `predicted_doublet` (bool) into `obs`.
`decisions.yaml` `remove_cells` entries accept any `obs.eval` expression, and both forms
were verified to evaluate correctly:

```yaml
remove_cells:
  - query: "predicted_doublet"
    reason: "Scrublet-predicted doublet"
  # or a manual cut:
  - query: "doublet_score > 0.35"
    reason: "High doublet score"
```

So the work is **scoring and reporting**, not filtering.

## Work items (~1 day)

1. `run_scrublet()` in the QC section of `integrate.py` (near
   [line 1391](../code/sccairns/integrate.py#L1391)), mirroring `run_harmony()`: config
   gated, with an actionable `ImportError` naming `scikit-image`.
2. A `qc.doublets` config block — `enabled`, `batch_key`, `restrict_to`,
   `expected_doublet_rate`, `threshold`, `first_round_only`.
3. Doublet provenance in `round_manifest.json` (scores are an input to later decisions,
   so which parameters produced them must be recorded).
4. Per-cluster aggregation in `cluster_qc_summary` (median score, fraction predicted
   doublet) plus an `auto_flag` threshold — the same pattern used for neighbor purity and
   compositional bias.
5. Tests + a section in `docs/configuration.md` and `docs/interpreting-outputs.md`.

## Design decisions that must be settled first

**Scrublet is a droplet method and this object is not all droplets.** It simulates
doublets by summing random pairs of transcriptomes, modelling co-encapsulation. SSv4 is
plate-based with physically sorted cells, so the generative assumption does not hold and
scores are not comparable across arms. Scoring must be restricted to the droplet arm
(`tech == "10x_Multiome"`) rather than silently applied to everything. A `restrict_to`
config field should make that explicit.

**`batch_key` must be the sequencing library, not `data.batch_key`.** `data.batch_key` is
`tech`, and pooling every Multiome library into one simulation is wrong. `sample_key`
(`donor_id`) is closer; a per-lane key is better still.

**Score on the first round only.** Re-running on an already-filtered subset changes the
simulated background, so the same cell gets a different score in round 2. This should
mirror the existing `qc.skip_filter_after_round_1` mechanism (detected via the input's
parent-round pointer).

**Score and flag; do not auto-remove.** Doublet callers have real false-positive rates on
transitional states. Predicted doublets belong in `auto_flags.yaml` for endorsement in
`decisions.yaml`, consistent with the rest of the loop's human gate.

**For 10x Multiome specifically, ATAC-based detection (AMULET) generally outperforms
RNA-only Scrublet.** Scrublet on the RNA half is standard but is not the strongest option
available for this modality.

## Current state (2026-08)

Doublet scores are being produced **upstream in R by DoubletFinder** and arrive in `obs`
as `doublet_score`, present only for the droplet platform. A per-cluster doublet-score
plot that keys off whatever `doublet_score` is present — regardless of what produced it —
has been added to the inspection report, so the diagnostic is available without adopting
Scrublet at all. If DoubletFinder upstream stays, native Scrublet may never be needed;
the value of implementing it would be keeping scoring inside the reproducible round rather
than in a separate R step.

## Revisit when

- Doublet scoring needs to move out of R and into the reproducible round.
- `scikit-image` lands in the image for another reason (then the cost is code only).
- The droplet/plate split is resolved — e.g. the dataset becomes single-platform.

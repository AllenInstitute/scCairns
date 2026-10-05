#!/usr/bin/env python3
"""Heat-shock / stress module scoring, per cell and per cluster.

Dissociation and handling induce a heat-shock response that is *technical*, not
biological: a cluster defined by Hspa1a/Hspa1b/Hsp90aa1 is usually an artifact
of how long cells sat in the protease, not a cell type. This module scores that
signature so such clusters are visible in the inspection report.

Deliberately scoped to heat-shock genes only. Mitochondrial burden — the other
half of the usual "stress" bundle — is already covered end to end: QC computes
``pct_counts_mt`` from ``qc.mt_gene_patterns``, the inspection report surfaces
it as ``median_pct_mt``, and ``inspection.auto_flag.mt_threshold`` flags on it.
Adding a second mitochondrial number would give two metrics that can disagree
about the same cluster, so read the stress score *alongside* ``median_pct_mt``
rather than folding them together.

Scoring reuses :func:`sccairns.contamination.score_panel` — the same per-gene
z-score path the contamination panels use, so the two are read on one scale.
Keeping it out of ``DEFAULT_CONTAM_PANELS`` is intentional: contamination
panels feed ``flag_any_contam`` and therefore
``contamination_flagged_cells.csv``, and a stressed cell is a different
decision from an ambient/contaminating one — usually "drop this cluster", not
"drop these cells as the wrong lineage".

Panels are resolved by *prefix* against ``adata.var_names`` at run time rather
than listed gene by gene, because the HSP families are large and the useful
member varies by dataset. Prefix matching is also how the same default covers
mouse (``Hspa1a``) and human (``HSPA1A``) without restating the panel.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .contamination import score_panel

# Heat-shock families. `Hsp` also covers Hspa/Hspb/Hsph and their co-chaperones
# (Hspbp1); that breadth is intentional — the score is the mean z across
# whatever is present, so a few extra family members move it very little, while
# a missing Hspa1a would matter. Both cases are listed so the default works on
# mouse and human objects; `gene_symbol_case: preserve` means we cannot assume
# either has been normalized.
#
# Immediate-early genes (Fos, Jun, Egr1, Atf3) are the other classic
# dissociation signature and are deliberately NOT included by default: in
# neurons they are also genuine activity markers, so scoring them as "stress"
# would penalize real biology. Add them via `inspection.stress.gene_prefixes`
# when the tissue makes that unambiguous.
DEFAULT_STRESS_PREFIXES: List[str] = ["Hsp", "HSP", "Dnaj", "DNAJ"]


def resolve_stress_genes(adata, prefixes: Optional[Sequence[str]] = None) -> List[str]:
    """Panel genes present in ``adata``, matched by prefix and sorted.

    Sorted rather than set-ordered on purpose: ``list(set(...))`` varies run to
    run because Python randomizes string hashing per process, which would make
    the recorded panel — and any downstream diff of it — nondeterministic.
    """
    prefixes = tuple(DEFAULT_STRESS_PREFIXES if prefixes is None else prefixes)
    if not prefixes:
        return []
    return sorted(g for g in adata.var_names if str(g).startswith(prefixes))


def score_stress(adata,
                 *,
                 prefixes: Optional[Sequence[str]] = None,
                 genes: Optional[Sequence[str]] = None,
                 layer: Optional[str] = None,
                 z_thresh: float = 2.0) -> Dict[str, object]:
    """Score the stress panel for every cell.

    Parameters
    ----------
    genes : sequence of str or None
        Explicit panel. Overrides ``prefixes`` when given — an escape hatch for
        a curated list; prefix matching is the default path.
    layer : str or None
        Layer to score; None uses ``adata.X``, which must be log-normalized
        (the same assumption ``score_panel`` documents).
    z_thresh : float
        Per-gene z above which a gene counts as a hit in a cell.

    Returns
    -------
    dict with ``score`` (ndarray), ``n_hits`` (ndarray), ``used`` and
    ``missing`` (lists). ``used`` is empty when nothing matched — the caller
    decides whether that is worth reporting, since an object with no HSP genes
    is a legitimate state (e.g. a heavily subset gene space), not an error.
    """
    panel = list(genes) if genes else resolve_stress_genes(adata, prefixes)
    if not panel:
        return {"score": np.zeros(adata.n_obs), "n_hits": np.zeros(adata.n_obs, dtype=int),
                "used": [], "missing": []}

    score, n_hits, used, missing = score_panel(
        adata, panel, layer=layer, z_thresh=z_thresh)
    return {"score": score, "n_hits": n_hits, "used": used, "missing": missing}


def stress_by_cluster(adata, cluster_key: str,
                      score_key: str = "stress_score") -> pd.DataFrame:
    """Per-cluster median and spread of the stress score.

    Median rather than mean: the score is a mean of z-scores and a handful of
    extreme cells would otherwise drag a cluster's value, which is the opposite
    of what the flag wants to detect (a cluster that is *uniformly* stressed).
    """
    if score_key not in adata.obs.columns:
        return pd.DataFrame()
    obs = adata.obs[[cluster_key, score_key]].copy()
    g = obs.groupby(cluster_key, observed=True)[score_key]
    return pd.DataFrame({
        "median_stress": g.median(),
        "stress_q90": g.quantile(0.90),
    })

#!/usr/bin/env python3
"""Per-cluster compositional bias: is a cluster's makeup what chance would give?

Two questions with one test. A cluster missing a platform entirely (no SSv4
cells) is *depletion*; a cluster carried by a single donor is *enrichment*. Both
ask whether a cluster's composition departs from the dataset-wide composition by
more than random assignment would produce.

Why a test rather than a fraction threshold: composition alone cannot tell a
20-cell cluster with no SSv4 (unremarkable — that happens by chance) from a
200-cell one (impossible by chance). With 15% SSv4 overall those score p = 0.04
and p = 2e-15 respectively while looking identical to any "fraction > 0.9" rule.

The null is exchangeability — that cells land in clusters without regard to which
platform or donor produced them. **A significant result means the deviation is
likely real, not that the cluster may be an artifact.** If one platform was sorted for
neurons and the other was not, every non-neuronal cluster is genuinely depleted
of the sorted platform and that is correct biology, not junk to remove.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd


def _fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, ignoring NaN entries."""
    out = np.full(len(pvalues), np.nan)
    finite = np.flatnonzero(np.isfinite(pvalues))
    if finite.size == 0:
        return out
    try:
        from scipy.stats import false_discovery_control

        out[finite] = false_discovery_control(pvalues[finite])
    except ImportError:  # pragma: no cover - scipy < 1.11
        ordered = finite[np.argsort(pvalues[finite])]
        m = ordered.size
        adjusted = pvalues[ordered] * m / np.arange(1, m + 1)
        out[ordered] = np.minimum.accumulate(adjusted[::-1])[::-1]
    return np.clip(out, 0, 1)


def compositional_bias(
    obs: pd.DataFrame,
    cluster_key: str,
    group_key: str,
    *,
    min_cells: int = 20,
    prefix: Optional[str] = None,
) -> pd.DataFrame:
    """Test each cluster's composition over ``group_key`` against the whole dataset.

    For every cluster, the group most over-represented relative to its
    dataset-wide share is tested for enrichment (hypergeometric upper tail) and
    the most under-represented is tested for depletion (lower tail). Because the
    tested group is *chosen* by looking at the data, each p-value is Bonferroni-
    corrected by the number of groups before Benjamini-Hochberg is applied across
    clusters.

    Clusters smaller than ``min_cells`` are reported as NaN and excluded from the
    FDR correction — they carry no power and would only dilute it.

    Returns one row per cluster with ``{prefix}_`` columns: ``top``,
    ``top_frac``, ``enrichment`` (observed/expected for that group),
    ``enrich_q``, ``depleted``, ``depletion_ratio``, ``depletion_q``, and
    ``n_groups_present``.
    """
    from scipy.stats import hypergeom

    if cluster_key not in obs.columns:
        raise KeyError(f"cluster column not in obs: {cluster_key!r}")
    if group_key not in obs.columns:
        raise KeyError(f"group column not in obs: {group_key!r}")

    prefix = prefix or group_key
    table = pd.crosstab(obs[cluster_key].astype(str), obs[group_key].astype(str))
    # A category with no cells anywhere carries no information and would make
    # every cluster trivially "depleted" of it.
    table = table.loc[:, table.sum(axis=0) > 0]

    total = int(table.values.sum())
    group_totals = table.sum(axis=0)
    n_groups = table.shape[1]

    records = []
    for cluster in table.index:
        counts = table.loc[cluster]
        n_cells = int(counts.sum())
        record = {
            "cluster": cluster,
            f"{prefix}_n_groups_present": int((counts > 0).sum()),
            f"{prefix}_top": None,
            f"{prefix}_top_frac": np.nan,
            f"{prefix}_enrichment": np.nan,
            f"{prefix}_enrich_p": np.nan,
            f"{prefix}_depleted": None,
            f"{prefix}_depletion_ratio": np.nan,
            f"{prefix}_depletion_p": np.nan,
        }
        if n_groups < 2 or n_cells < min_cells or n_cells == 0:
            records.append(record)
            continue

        expected = group_totals * n_cells / total
        ratio = counts / expected

        top = ratio.idxmax()
        record[f"{prefix}_top"] = str(top)
        record[f"{prefix}_top_frac"] = float(counts[top] / n_cells)
        record[f"{prefix}_enrichment"] = float(ratio[top])
        # sf(k-1) is P(X >= k): the chance of seeing this many or more.
        record[f"{prefix}_enrich_p"] = min(1.0, float(
            hypergeom.sf(counts[top] - 1, total, group_totals[top], n_cells)
        ) * n_groups)

        low = ratio.idxmin()
        record[f"{prefix}_depleted"] = str(low)
        record[f"{prefix}_depletion_ratio"] = float(ratio[low])
        record[f"{prefix}_depletion_p"] = min(1.0, float(
            hypergeom.cdf(counts[low], total, group_totals[low], n_cells)
        ) * n_groups)

        records.append(record)

    frame = pd.DataFrame(records).set_index("cluster")
    frame[f"{prefix}_enrich_q"] = _fdr(frame[f"{prefix}_enrich_p"].to_numpy())
    frame[f"{prefix}_depletion_q"] = _fdr(frame[f"{prefix}_depletion_p"].to_numpy())
    return frame.drop(columns=[f"{prefix}_enrich_p", f"{prefix}_depletion_p"])


def resolve_composition_keys(
    obs: pd.DataFrame,
    batch_key: Optional[str] = None,
    sample_key: Optional[str] = None,
    covariate_keys: Optional[Sequence[str]] = None,
    explicit: Optional[Sequence[str]] = None,
) -> list:
    """Which obs columns to test, de-duplicated and filtered to what exists.

    ``explicit`` wins outright; otherwise the batch key, the sample/donor key,
    and any categorical covariates are tested, since those are the variables a
    cluster is not supposed to be a function of.
    """
    if explicit is not None:
        candidates = list(explicit)
    else:
        candidates = [batch_key, sample_key, *(covariate_keys or [])]

    keys, seen = [], set()
    for key in candidates:
        if not key or key in seen or key not in obs.columns:
            continue
        # A constant column has nothing to test against.
        if obs[key].nunique(dropna=True) < 2:
            continue
        seen.add(key)
        keys.append(key)
    return keys

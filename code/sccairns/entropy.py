#!/usr/bin/env python3
"""Neighborhood label entropy — how mixed each cell's local neighborhood is.

For every cell, look at its k nearest neighbors in an embedding and measure the
Shannon entropy of some categorical label across them. Entropy near 0 means the
neighborhood agrees (the cell sits inside a homogeneous island); high entropy
means the label is mixed there, which flags doublets, transitional states, cells
whose annotation disagrees with their neighbors, or — when the label is the batch
— residual batch mixing.

Ported from an R implementation (FNN + entropy). Three deliberate departures:

1. **Neighbors come from a multi-dimensional embedding.** The R version indexed a
   single meta.data column, so neighborhoods were computed along one principal
   component rather than in the embedding space. Pass ``n_comps`` to bound the
   dimensionality explicitly.
2. **Vectorized.** The R version assigned into a data.table row by row inside a
   per-cell loop, which does not scale past a few thousand cells.
3. **Missing labels are excluded** from a neighborhood rather than counted as a
   category, so an unannotated neighbor lowers the denominator instead of
   inventing a level.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd


def _neighbor_indices(
    coords: np.ndarray, n_neighbors: int, include_self: bool
) -> np.ndarray:
    """Indices of each row's ``n_neighbors`` nearest rows in ``coords``."""
    from sklearn.neighbors import NearestNeighbors

    n_cells = coords.shape[0]
    # One extra so that, after dropping each point's own index, every row still
    # has n_neighbors. Duplicate coordinates mean self is not reliably first.
    k = min(n_neighbors + (0 if include_self else 1), n_cells)
    finder = NearestNeighbors(n_neighbors=k).fit(coords)
    indices = finder.kneighbors(coords, return_distance=False)

    if include_self:
        return indices[:, :n_neighbors]

    is_self = indices == np.arange(n_cells)[:, None]
    # Stable sort pushes the self column to the end, keeping distance order.
    order = np.argsort(is_self, kind="stable", axis=1)
    return np.take_along_axis(indices, order, axis=1)[:, :n_neighbors]


def _entropy_per_row(
    codes: np.ndarray,
    indices: np.ndarray,
    n_categories: int,
    base: Optional[float],
    chunk: int = 8192,
) -> np.ndarray:
    """Shannon entropy of ``codes`` over each row of ``indices``.

    ``codes`` is a categorical encoding where -1 marks a missing label; those
    neighbors are dropped from both the counts and the denominator. Runs in
    row-blocks so the one-hot count matrix stays bounded regardless of how many
    categories the annotation has.
    """
    n_rows = indices.shape[0]
    out = np.empty(n_rows, dtype=float)

    for start in range(0, n_rows, chunk):
        stop = min(start + chunk, n_rows)
        block = codes[indices[start:stop]]
        n_block = block.shape[0]
        observed = block >= 0

        # Offset each row into its own slice of a flat bincount.
        flat = np.where(observed, block, 0) + np.arange(n_block)[:, None] * n_categories
        counts = np.bincount(
            flat.ravel(),
            weights=observed.ravel().astype(float),
            minlength=n_block * n_categories,
        ).reshape(n_block, n_categories)

        totals = counts.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            proportions = counts / totals[:, None]
            terms = np.where(proportions > 0, proportions * np.log(proportions), 0.0)
            entropy = -terms.sum(axis=1)
        # Negating a zero sum yields -0.0, which is equal to 0 but reads as a
        # negative entropy in tables and CSVs.
        entropy = np.where(entropy == 0, 0.0, entropy)
        # A neighborhood with no annotated neighbors has no entropy to report.
        entropy[totals == 0] = np.nan
        out[start:stop] = entropy

    if base is not None:
        out /= np.log(base)
    return out


def neighborhood_entropy(
    adata,
    annotation_keys: Union[str, Sequence[str]],
    *,
    use_rep: str = "X_pca",
    n_comps: Optional[int] = None,
    n_neighbors: int = 15,
    include_self: bool = False,
    normalize: bool = False,
    base: Optional[float] = None,
    key_suffix: str = "_entropy",
) -> pd.DataFrame:
    """Add per-cell neighborhood entropy for one or more ``obs`` annotations.

    Parameters
    ----------
    adata
        AnnData with an embedding in ``obsm[use_rep]``.
    annotation_keys
        One or more categorical ``obs`` columns to measure mixing for.
    use_rep, n_comps
        Embedding to find neighbors in, and how many of its dimensions to use
        (all by default). Use the *integrated* embedding — entropy computed on an
        uncorrected PCA mostly measures batch structure.
    n_neighbors
        Neighborhood size. Entropy is bounded by ``log(n_neighbors)``, so small
        values quantize the result coarsely.
    include_self
        Whether a cell counts as one of its own neighbors. Default False, which
        makes the score a statement about the *surroundings* rather than partly
        about the cell itself.
    normalize
        Divide by the maximum attainable entropy,
        ``log(min(n_neighbors, n_categories))``, putting the result on 0–1 so
        annotations with different category counts are comparable.
    base
        Log base. ``None`` gives nats (matching R's ``entropy::entropy``); pass 2
        for bits. Ignored when ``normalize`` is set, which is already unitless.

    Returns
    -------
    The entropy columns as a DataFrame indexed like ``adata.obs``. The same
    columns are written into ``adata.obs`` as ``f"{key}{key_suffix}"``.
    """
    if isinstance(annotation_keys, str):
        annotation_keys = [annotation_keys]
    annotation_keys = list(annotation_keys)

    target = adata

    if use_rep not in target.obsm:
        raise KeyError(
            f"embedding not found in obsm: {use_rep!r} "
            f"(available: {sorted(target.obsm.keys())})"
        )
    missing = [k for k in annotation_keys if k not in target.obs.columns]
    if missing:
        raise KeyError(f"annotation column(s) not in obs: {', '.join(missing)}")
    if n_neighbors < 1:
        raise ValueError(f"n_neighbors must be >= 1, got {n_neighbors}")
    if target.n_obs <= 1:
        raise ValueError("neighborhood entropy needs at least 2 cells")

    coords = np.asarray(target.obsm[use_rep])
    if coords.ndim == 1:
        coords = coords[:, None]
    if n_comps is not None:
        coords = coords[:, :n_comps]
    if coords.shape[1] == 1:
        # Legal, but almost never intended — the R original did this implicitly.
        print(f"  [WARN] {use_rep} has one dimension; neighborhoods are being "
              "computed along a single axis, not in the embedding space.")

    effective_k = min(n_neighbors, target.n_obs - (0 if include_self else 1))
    if effective_k < n_neighbors:
        print(f"  [WARN] only {target.n_obs} cells; using {effective_k} neighbors")

    indices = _neighbor_indices(coords, effective_k, include_self)

    results = {}
    for key in annotation_keys:
        categorical = pd.Categorical(target.obs[key])
        n_categories = len(categorical.categories)
        if n_categories == 0:
            raise ValueError(f"annotation column has no categories: {key!r}")

        values = _entropy_per_row(
            np.asarray(categorical.codes, dtype=np.intp),
            indices,
            n_categories,
            base=None if normalize else base,
        )
        if normalize:
            ceiling = np.log(min(effective_k, n_categories))
            values = values / ceiling if ceiling > 0 else np.zeros_like(values)

        column = f"{key}{key_suffix}"
        target.obs[column] = values
        results[column] = values

    return pd.DataFrame(results, index=target.obs_names)


def cluster_label_coherence(
    adata,
    cluster_key: str,
    *,
    use_rep: str = "X_pca",
    n_comps: Optional[int] = None,
    n_neighbors: int = 15,
) -> pd.DataFrame:
    """Per-cluster measures of how well each cluster holds together locally.

    Two complementary statistics, both computed from one k-NN graph built on
    ``use_rep`` and both reported per cluster:

    ``entropy_median`` / ``entropy_mean``
        Normalized (0–1) Shannon entropy of the cluster label across each cell's
        neighborhood. Answers *how many different things is this mixed with*.
    ``purity_median`` / ``purity_mean``
        Fraction of a cell's neighbors carrying its own cluster label. Answers
        *how much of the neighborhood is mine*, and is the more sensitive of the
        two to the commonest instability — a cluster split roughly in half with
        one neighbor. That case is only moderately *entropic* (two categories out
        of many) while being unambiguously impure.

    A well-supported cluster sits near entropy 0 / purity 1. Read them together:
    low purity with low entropy means one specific cluster is bleeding into this
    one; low purity with high entropy means it is diffuse across many.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"cluster column not in obs: {cluster_key!r}")
    if use_rep not in adata.obsm:
        raise KeyError(
            f"embedding not found in obsm: {use_rep!r} "
            f"(available: {sorted(adata.obsm.keys())})"
        )

    coords = np.asarray(adata.obsm[use_rep])
    if coords.ndim == 1:
        coords = coords[:, None]
    if n_comps is not None:
        coords = coords[:, :n_comps]

    effective_k = min(n_neighbors, adata.n_obs - 1)
    indices = _neighbor_indices(coords, effective_k, include_self=False)

    categorical = pd.Categorical(adata.obs[cluster_key])
    codes = np.asarray(categorical.codes, dtype=np.intp)
    n_categories = len(categorical.categories)

    entropy = _entropy_per_row(codes, indices, n_categories, base=None)
    ceiling = np.log(min(effective_k, n_categories))
    if ceiling > 0:
        entropy = entropy / ceiling

    # Fraction of each cell's neighbors that share its own label.
    own = codes[:, None]
    neighbor_codes = codes[indices]
    observed = neighbor_codes >= 0
    purity = np.where(
        observed.sum(axis=1) > 0,
        ((neighbor_codes == own) & observed).sum(axis=1) / np.maximum(
            observed.sum(axis=1), 1),
        np.nan,
    )

    frame = pd.DataFrame(
        {"cluster": adata.obs[cluster_key].astype(str).to_numpy(),
         "entropy": entropy, "purity": purity}
    )
    grouped = frame.groupby("cluster", observed=True)
    summary = pd.DataFrame({
        "entropy_median": grouped["entropy"].median(),
        "entropy_mean": grouped["entropy"].mean(),
        "purity_median": grouped["purity"].median(),
        "purity_mean": grouped["purity"].mean(),
    })
    summary.index.name = "cluster"

    # Which *other* cluster the neighborhood belongs to, pooled over the whole
    # cluster. Names the merge candidate for an absorbed cluster, and reads as a
    # low fraction when the cluster is spread across many rather than one.
    dominant, dominant_frac = [], []
    categories = [str(c) for c in categorical.categories]
    for cluster in summary.index:
        own_code = categories.index(cluster)
        rows = np.flatnonzero(codes == own_code)
        foreign = neighbor_codes[rows]
        foreign = foreign[(foreign >= 0) & (foreign != own_code)]
        if foreign.size == 0:
            dominant.append(None)
            dominant_frac.append(0.0)
            continue
        counts = np.bincount(foreign, minlength=n_categories)
        top = int(counts.argmax())
        dominant.append(categories[top])
        dominant_frac.append(float(counts[top] / foreign.size))

    summary["dominant_neighbor"] = dominant
    summary["dominant_neighbor_frac"] = dominant_frac
    return summary

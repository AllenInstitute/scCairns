"""Neighborhood label entropy, checked against hand-computable cases."""

import numpy as np
import pytest

pytest.importorskip("anndata")
pytest.importorskip("sklearn")

import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402

from sccairns.entropy import neighborhood_entropy  # noqa: E402


def _adata(coords, labels, key="label"):
    obj = ad.AnnData(np.ones((len(labels), 2)))
    obj.obs_names = [f"cell{i}" for i in range(len(labels))]
    obj.obs[key] = pd.Categorical(labels)
    obj.obsm["X_emb"] = np.asarray(coords, dtype=float)
    return obj


def test_uniform_neighborhood_has_zero_entropy():
    # Every cell is surrounded by its own label; nothing is mixed.
    coords = [[0, 0], [0, 1], [0, 2], [10, 0], [10, 1], [10, 2]]
    labels = ["a", "a", "a", "b", "b", "b"]
    obj = _adata(coords, labels)

    out = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=2)

    assert np.allclose(out["label_entropy"], 0.0)
    assert "label_entropy" in obj.obs        # written back to obs


def test_evenly_split_neighborhood_is_log_two():
    # Two neighbors, one of each label -> H = -2 * 0.5*ln(0.5) = ln 2.
    coords = [[0, 0], [1, 0], [2, 0], [3, 0]]
    labels = ["a", "b", "a", "b"]
    obj = _adata(coords, labels)

    out = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=2)

    # cell1's neighbors are cell0 ("a") and cell2 ("a") -> homogeneous.
    # cell0's are cell1 ("b") and cell2 ("a") -> evenly split.
    assert out["label_entropy"].iloc[0] == pytest.approx(np.log(2))
    assert out["label_entropy"].iloc[1] == pytest.approx(0.0)


def test_base_two_reports_bits_and_normalize_reports_a_fraction():
    coords = [[0, 0], [1, 0], [2, 0], [3, 0]]
    labels = ["a", "b", "a", "b"]
    obj = _adata(coords, labels)

    bits = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=2,
                                base=2, key_suffix="_bits")
    assert bits["label_bits"].iloc[0] == pytest.approx(1.0)

    fraction = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=2,
                                    normalize=True, key_suffix="_norm")
    # ln2 over the ceiling ln(min(k=2, categories=2)) = ln2 -> exactly 1.
    assert fraction["label_norm"].iloc[0] == pytest.approx(1.0)
    assert fraction["label_norm"].max() <= 1.0


def test_self_is_excluded_by_default():
    # cell0 is the lone "a" in a sea of "b". Excluding self, its neighborhood is
    # pure "b" (entropy 0); including self it is 1-of-3 "a" (entropy > 0).
    coords = [[0, 0], [1, 0], [2, 0], [3, 0]]
    labels = ["a", "b", "b", "b"]
    obj = _adata(coords, labels)

    without = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=3,
                                   key_suffix="_without")
    with_self = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=3,
                                     include_self=True, key_suffix="_with")

    assert without["label_without"].iloc[0] == pytest.approx(0.0)
    assert with_self["label_with"].iloc[0] > 0.0


def test_duplicate_coordinates_still_exclude_self():
    # With identical points, sklearn need not return self first; the row must
    # still contain exactly n_neighbors entries and never the cell itself.
    coords = [[0, 0]] * 6
    labels = ["a", "a", "a", "b", "b", "b"]
    obj = _adata(coords, labels)

    out = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=3)

    assert out["label_entropy"].notna().all()
    assert (out["label_entropy"] >= 0).all()


def test_missing_labels_shrink_the_denominator():
    # cell0's three neighbors are "b", NaN, NaN -> one observed label -> H = 0,
    # rather than NaN counting as its own category.
    coords = [[0, 0], [1, 0], [2, 0], [3, 0]]
    labels = ["a", "b", None, None]
    obj = _adata(coords, labels)

    out = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=3)

    assert out["label_entropy"].iloc[0] == pytest.approx(0.0)


def test_neighborhood_with_no_annotated_neighbors_is_nan():
    coords = [[0, 0], [1, 0], [2, 0]]
    labels = ["a", None, None]
    obj = _adata(coords, labels)

    out = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=2)

    assert np.isnan(out["label_entropy"].iloc[0])


def test_several_annotations_share_one_neighbor_graph():
    coords = [[0, 0], [1, 0], [2, 0], [3, 0]]
    obj = _adata(coords, ["a", "a", "b", "b"], key="cluster")
    obj.obs["batch"] = pd.Categorical(["x", "y", "x", "y"])

    out = neighborhood_entropy(obj, ["cluster", "batch"], use_rep="X_emb",
                               n_neighbors=2)

    assert list(out.columns) == ["cluster_entropy", "batch_entropy"]
    assert list(out.index) == list(obj.obs_names)


def test_chunking_does_not_change_the_result():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(300, 5))
    labels = rng.choice(list("abcde"), size=300)
    obj = _adata(coords, labels)

    from sccairns import entropy as mod

    full = neighborhood_entropy(obj, "label", use_rep="X_emb", n_neighbors=10)
    indices = mod._neighbor_indices(coords, 10, include_self=False)
    codes = np.asarray(pd.Categorical(obj.obs["label"]).codes, dtype=np.intp)
    chunked = mod._entropy_per_row(codes, indices, 5, base=None, chunk=7)

    assert np.allclose(full["label_entropy"].to_numpy(), chunked)


def test_missing_embedding_and_column_name_the_problem():
    obj = _adata([[0, 0], [1, 0]], ["a", "b"])

    with pytest.raises(KeyError, match="X_nope"):
        neighborhood_entropy(obj, "label", use_rep="X_nope")
    with pytest.raises(KeyError, match="nosuch"):
        neighborhood_entropy(obj, "nosuch", use_rep="X_emb")


# ── Per-cluster coherence and auto-flagging ───────────────────────────────────

def _structured_adata(seed=0):
    """Four clean clusters, one blob split arbitrarily in two, one diffuse."""
    rng = np.random.default_rng(seed)
    d = 10
    centers = rng.normal(scale=10, size=(6, d))
    blocks, labels = [], []
    for i in range(4):
        blocks.append(rng.normal(centers[i], 1.0, size=(120, d)))
        labels += [f"good{i}"] * 120
    blocks.append(rng.normal(centers[4], 1.0, size=(200, d)))
    labels += ["split_a"] * 100 + ["split_b"] * 100
    blocks.append(rng.normal(centers.mean(axis=0), 14.0, size=(80, d)))
    labels += ["diffuse"] * 80

    X = np.vstack(blocks)
    obj = ad.AnnData(np.ones((X.shape[0], 2)))
    obj.obsm["X_scVI"] = X
    obj.obs["leiden"] = pd.Categorical(labels)
    return obj


def test_well_separated_clusters_are_pure_and_ordered_clusters_are_not():
    from sccairns.entropy import cluster_label_coherence

    out = cluster_label_coherence(_structured_adata(), "leiden",
                                  use_rep="X_scVI", n_neighbors=15)

    for i in range(4):
        assert out.loc[f"good{i}", "purity_median"] == pytest.approx(1.0)
        assert out.loc[f"good{i}", "entropy_median"] == pytest.approx(0.0)
    # An arbitrarily split blob is about half its own label.
    assert 0.3 < out.loc["split_a", "purity_median"] < 0.7
    assert out.loc["diffuse", "purity_median"] < 0.5


def test_entropy_alone_misses_a_cluster_absorbed_into_another():
    """The case that decides purity-vs-entropy for flagging.

    A cluster with no territory of its own sits entirely inside a larger one —
    the clearest "not well supported" there is. Its neighborhoods are
    homogeneous, just homogeneously *someone else's* label, so entropy scores it
    no worse than a perfectly isolated cluster. Purity reads 0 against 1.

    This holds at every size tried: the absorbed cluster's entropy is never
    above the healthy cluster's, so no entropy threshold can catch it without
    also flagging the good ones.
    """
    from sccairns.entropy import cluster_label_coherence

    rng = np.random.default_rng(1)
    d = 10
    big = rng.normal(0, 1.0, size=(600, d))
    ghost = rng.normal(0, 1.0, size=(10, d))     # same distribution as big
    far = rng.normal(25, 1.0, size=(200, d))
    X = np.vstack([big, ghost, far])
    obj = ad.AnnData(np.ones((X.shape[0], 2)))
    obj.obsm["X_scVI"] = X
    obj.obs["leiden"] = pd.Categorical(["big"] * 600 + ["ghost"] * 10
                                       + ["far"] * 200)

    out = cluster_label_coherence(obj, "leiden", use_rep="X_scVI",
                                  n_neighbors=15)

    # Entropy is blind here: every cluster scores identically.
    assert out.loc["ghost", "entropy_median"] == pytest.approx(0.0)
    assert out.loc["far", "entropy_median"] == pytest.approx(0.0)
    assert out.loc["ghost", "entropy_median"] <= out.loc["big", "entropy_median"]
    # Purity separates them decisively.
    assert out.loc["ghost", "purity_median"] < 0.2
    assert out.loc["far", "purity_median"] == pytest.approx(1.0)
    assert out.loc["big", "purity_median"] > 0.9


def test_auto_flag_uses_purity_and_names_the_absorbing_cluster():
    pytest.importorskip("scanpy")
    from sccairns.inspect import auto_flag_clusters

    summary = pd.DataFrame(
        {"n_cells": [500, 500, 500],
         "neighbor_purity": [1.0, 0.10, 0.45],
         "neighbor_entropy": [0.0, 0.02, 0.60],
         "dominant_neighbor": [None, "7", "7"],
         "dominant_neighbor_frac": [0.0, 0.95, 0.30]},
        index=["clean", "absorbed", "diffuse"],
    )

    flags = auto_flag_clusters(summary, min_cells=1, min_neighbor_purity=0.5)

    assert "clean" not in flags
    # Naming the absorbing cluster is what tells you to merge rather than drop.
    assert "95% of its neighborhood is cluster 7" in flags["absorbed"][0]
    assert "10%" in flags["absorbed"][0]
    # No single cluster dominates the diffuse one, so none is named.
    assert "spread across several clusters" in flags["diffuse"][0]


def test_purity_flag_is_disablable_and_entropy_flag_is_opt_in():
    pytest.importorskip("scanpy")
    from sccairns.inspect import auto_flag_clusters

    summary = pd.DataFrame(
        {"n_cells": [500], "neighbor_purity": [0.10], "neighbor_entropy": [0.9]},
        index=["absorbed"],
    )

    assert auto_flag_clusters(summary, min_cells=1,
                              min_neighbor_purity=None) == {}
    reasons = auto_flag_clusters(summary, min_cells=1, min_neighbor_purity=None,
                                 entropy_threshold=0.5)["absorbed"]
    assert "High neighborhood entropy" in reasons[0]


def test_cluster_qc_summary_gains_the_coherence_columns():
    pytest.importorskip("scanpy")
    from sccairns.inspect import cluster_qc_summary

    obj = _structured_adata()
    obj.obs["data_origin"] = pd.Categorical(["a"] * obj.n_obs)

    summary = cluster_qc_summary(obj, cluster_key="leiden",
                                 batch_key="data_origin", latent_key="X_scVI")

    assert {"neighbor_purity", "neighbor_entropy"} <= set(summary.columns)
    assert summary.loc["good0", "neighbor_purity"] == pytest.approx(1.0)
    assert summary["neighbor_purity"].notna().all()


def test_coherence_columns_absent_without_a_latent_embedding():
    pytest.importorskip("scanpy")
    from sccairns.inspect import cluster_qc_summary

    obj = _structured_adata()
    obj.obs["data_origin"] = pd.Categorical(["a"] * obj.n_obs)

    summary = cluster_qc_summary(obj, cluster_key="leiden",
                                 batch_key="data_origin", latent_key=None)

    assert "neighbor_purity" not in summary.columns


def test_dominant_neighbor_names_the_absorbing_cluster():
    from sccairns.entropy import cluster_label_coherence

    rng = np.random.default_rng(2)
    d = 10
    X = np.vstack([
        rng.normal(0, 1.0, size=(400, d)),     # host
        rng.normal(0, 1.0, size=(30, d)),      # absorbed into the host
        rng.normal(30, 1.0, size=(200, d)),    # unrelated, far away
    ])
    obj = ad.AnnData(np.ones((X.shape[0], 2)))
    obj.obsm["X_scVI"] = X
    obj.obs["leiden"] = pd.Categorical(["host"] * 400 + ["ghost"] * 30
                                       + ["far"] * 200)

    out = cluster_label_coherence(obj, "leiden", use_rep="X_scVI",
                                  n_neighbors=15)

    assert out.loc["ghost", "dominant_neighbor"] == "host"
    assert out.loc["ghost", "dominant_neighbor_frac"] == pytest.approx(1.0)
    # An isolated cluster has no foreign neighbors at all to name.
    assert out.loc["far", "dominant_neighbor"] is None
    assert out.loc["far", "dominant_neighbor_frac"] == 0.0

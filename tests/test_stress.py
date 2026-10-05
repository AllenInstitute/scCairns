"""Heat-shock / dissociation-stress module scoring.

Covers the panel resolution rules (prefix matching, determinism, case), the
score itself, and the per-cluster aggregation. The panel is resolved against
`var_names` at run time rather than hard-coded, so most of the risk is in
resolution, not arithmetic.
"""

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")

from sccairns.stress import (  # noqa: E402
    DEFAULT_STRESS_PREFIXES,
    resolve_stress_genes,
    score_stress,
    stress_by_cluster,
)


def _adata(genes, n=60, seed=0):
    rng = np.random.RandomState(seed)
    X = rng.rand(n, len(genes)).astype("float32")
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=list(genes)),
    )


# ── panel resolution ────────────────────────────────────────────────────────
def test_resolves_hsp_genes_by_prefix():
    adata = _adata(["Hspa1a", "Hspb1", "Dnaja1", "Snap25", "Actb"])
    assert resolve_stress_genes(adata) == ["Dnaja1", "Hspa1a", "Hspb1"]


def test_resolution_covers_mouse_and_human_casing():
    """`gene_symbol_case: preserve` means we cannot assume either convention."""
    mouse = _adata(["Hspa1a", "Snap25"])
    human = _adata(["HSPA1A", "SNAP25"])
    assert resolve_stress_genes(mouse) == ["Hspa1a"]
    assert resolve_stress_genes(human) == ["HSPA1A"]


def test_resolution_is_sorted_not_set_ordered():
    """`list(set(...))` varies per process — string hashing is randomized."""
    adata = _adata(["Hspb1", "Hspa1a", "Hsp90aa1", "Actb"])
    out = resolve_stress_genes(adata)
    assert out == sorted(out)
    assert out == resolve_stress_genes(adata)


def test_custom_prefixes_override_the_default():
    adata = _adata(["Fos", "Jun", "Hspa1a"])
    assert resolve_stress_genes(adata, ["Fos", "Jun"]) == ["Fos", "Jun"]


def test_empty_prefix_list_resolves_to_nothing():
    adata = _adata(["Hspa1a"])
    assert resolve_stress_genes(adata, []) == []


def test_mitochondrial_genes_are_not_in_the_default_panel():
    """mt burden is pct_counts_mt's job; two metrics would disagree."""
    adata = _adata(["mt-Co1", "mt-Nd1", "Hspa1a"])
    assert resolve_stress_genes(adata) == ["Hspa1a"]


def test_immediate_early_genes_are_not_default():
    """In neurons IEGs are also activity markers — opt in, never automatic."""
    adata = _adata(["Fos", "Jun", "Egr1", "Atf3", "Hspa1a"])
    assert resolve_stress_genes(adata) == ["Hspa1a"]
    assert not any(p in ("Fos", "Jun", "Egr1") for p in DEFAULT_STRESS_PREFIXES)


# ── scoring ─────────────────────────────────────────────────────────────────
def test_score_tracks_stressed_cells():
    genes = ["Hspa1a", "Hspb1", "Hsp90aa1", "Actb"]
    adata = _adata(genes, n=100)
    # 10 of 100 cells stressed across the panel. The stressed fraction has to
    # stay under ~20%: a z-score caps at sqrt((1-p)/p) for a p-sized boosted
    # group, so at p=0.25 nothing can clear z_thresh=2.0 however large the
    # boost — the panel is only "hit" by a minority signal, by construction.
    for g in genes[:3]:
        adata.X[:10, adata.var_names.get_loc(g)] += 10.0

    out = score_stress(adata)

    assert out["used"] == ["Hsp90aa1", "Hspa1a", "Hspb1"]
    assert out["score"][:10].mean() > out["score"][10:].mean()
    assert out["n_hits"][:10].mean() > out["n_hits"][10:].mean()


def test_explicit_genes_override_prefixes():
    adata = _adata(["Hspa1a", "Hspb1", "Actb"])
    out = score_stress(adata, genes=["Actb"])
    assert out["used"] == ["Actb"]


def test_absent_panel_is_not_an_error():
    """A heavily subset gene space legitimately has no HSP genes."""
    adata = _adata(["Snap25", "Actb"])
    out = score_stress(adata)
    assert out["used"] == []
    assert out["score"].shape == (adata.n_obs,)
    assert not out["score"].any()


def test_score_reads_a_named_layer():
    adata = _adata(["Hspa1a", "Actb"])
    adata.layers["other"] = np.zeros_like(adata.X)
    out = score_stress(adata, layer="other")
    # Zero-variance input z-scores to zeros rather than blowing up.
    assert not out["score"].any()


# ── per-cluster aggregation ─────────────────────────────────────────────────
def test_stress_by_cluster_medians():
    adata = _adata(["Hspa1a"], n=40)
    adata.obs["leiden"] = pd.Categorical(["0"] * 20 + ["1"] * 20)
    adata.obs["stress_score"] = np.r_[np.ones(20) * 3.0, np.zeros(20)]

    out = stress_by_cluster(adata, "leiden")

    assert list(out.columns) == ["median_stress", "stress_q90"]
    assert out.loc["0", "median_stress"] == pytest.approx(3.0)
    assert out.loc["1", "median_stress"] == pytest.approx(0.0)


def test_stress_by_cluster_without_a_score_is_empty():
    adata = _adata(["Hspa1a"], n=10)
    adata.obs["leiden"] = pd.Categorical(["0"] * 10)
    assert stress_by_cluster(adata, "leiden").empty

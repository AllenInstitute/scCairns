"""Per-cluster compositional bias tests (platform depletion, donor dominance)."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scipy")

from sccairns.composition import compositional_bias, resolve_composition_keys


def _obs(spec, group_key="tech"):
    """spec: {cluster: {group: count}} -> a long obs frame."""
    rows = []
    for cluster, groups in spec.items():
        for group, count in groups.items():
            rows.extend([(cluster, group)] * count)
    return pd.DataFrame(rows, columns=["leiden", group_key])


def test_cluster_size_decides_significance_not_composition():
    """The whole reason this is a test and not a fraction threshold.

    Both clusters are 100% Multiome. Only the large one is beyond chance.
    """
    obs = _obs({
        "small": {"10x_Multiome": 20, "SSv4": 0},
        "big": {"10x_Multiome": 200, "SSv4": 0},
        "bulk1": {"10x_Multiome": 850, "SSv4": 150},
        "bulk2": {"10x_Multiome": 855, "SSv4": 150},
    })

    out = compositional_bias(obs, "leiden", "tech", min_cells=10)

    assert out.loc["small", "tech_depletion_ratio"] == 0.0
    assert out.loc["big", "tech_depletion_ratio"] == 0.0
    assert out.loc["small", "tech_depletion_q"] > 0.05      # unremarkable
    assert out.loc["big", "tech_depletion_q"] < 1e-9        # impossible by chance


def test_partial_depletion_is_caught_and_mild_deviation_is_not():
    obs = _obs({
        "partial": {"10x_Multiome": 190, "SSv4": 10},    # 10 seen, ~30 expected
        "mild": {"10x_Multiome": 175, "SSv4": 25},       # 25 seen, ~30 expected
        "bulk1": {"10x_Multiome": 850, "SSv4": 150},
        "bulk2": {"10x_Multiome": 855, "SSv4": 150},
    })

    out = compositional_bias(obs, "leiden", "tech", min_cells=10)

    assert out.loc["partial", "tech_depletion_q"] < 0.01
    assert out.loc["mild", "tech_depletion_q"] > 0.05


def test_donor_dominance_shows_up_as_enrichment():
    obs = _obs({
        "one_donor": {"D1": 235, "D2": 5, "D3": 5, "D4": 5},
        "even1": {"D1": 250, "D2": 250, "D3": 250, "D4": 250},
        "even2": {"D1": 250, "D2": 250, "D3": 250, "D4": 250},
    }, group_key="donor_id")

    out = compositional_bias(obs, "leiden", "donor_id", min_cells=20)

    assert out.loc["one_donor", "donor_id_top"] == "D1"
    assert out.loc["one_donor", "donor_id_enrichment"] > 2.0
    assert out.loc["one_donor", "donor_id_enrich_q"] < 1e-10

    # The even clusters can still come out *significant*, because the one-donor
    # cluster skews the dataset-wide donor shares that everything is compared
    # against, and n=1000 resolves a 1.1x deviation. Their effect size stays
    # near 1, which is what the flag's second gate keys on — significance alone
    # is not a usable criterion once clusters are large.
    assert out.loc["even1", "donor_id_enrichment"] < 1.2
    assert out.loc["even2", "donor_id_enrichment"] < 1.2


def test_small_clusters_are_not_tested_and_do_not_dilute_the_fdr():
    obs = _obs({
        "tiny": {"10x_Multiome": 5, "SSv4": 0},
        "big": {"10x_Multiome": 200, "SSv4": 0},
        "bulk": {"10x_Multiome": 1700, "SSv4": 300},
    })

    out = compositional_bias(obs, "leiden", "tech", min_cells=20)

    assert np.isnan(out.loc["tiny", "tech_depletion_q"])
    assert out.loc["tiny", "tech_top"] is None
    assert out.loc["big", "tech_depletion_q"] < 1e-9


def test_a_single_group_cannot_be_tested():
    obs = _obs({"a": {"10x_Multiome": 100}, "b": {"10x_Multiome": 100}})

    out = compositional_bias(obs, "leiden", "tech", min_cells=10)

    assert out["tech_depletion_q"].isna().all()
    assert out["tech_enrich_q"].isna().all()


def test_categories_with_no_cells_anywhere_are_ignored():
    # An unused categorical level would otherwise make every cluster
    # "depleted" of a group that does not exist.
    obs = _obs({"a": {"10x_Multiome": 100, "SSv4": 20},
                "b": {"10x_Multiome": 100, "SSv4": 20}})
    obs["tech"] = pd.Categorical(obs["tech"],
                                 categories=["10x_Multiome", "SSv4", "Ghost"])

    out = compositional_bias(obs, "leiden", "tech", min_cells=10)

    assert "Ghost" not in set(out["tech_depleted"].dropna())


def test_group_counts_are_reported():
    obs = _obs({"a": {"D1": 50, "D2": 50}, "b": {"D1": 100, "D2": 0}},
               group_key="donor_id")

    out = compositional_bias(obs, "leiden", "donor_id", min_cells=10)

    assert out.loc["a", "donor_id_n_groups_present"] == 2
    assert out.loc["b", "donor_id_n_groups_present"] == 1


def test_missing_columns_name_the_problem():
    obs = _obs({"a": {"D1": 30, "D2": 30}}, group_key="donor_id")

    with pytest.raises(KeyError, match="nosuch"):
        compositional_bias(obs, "nosuch", "donor_id")
    with pytest.raises(KeyError, match="nope"):
        compositional_bias(obs, "leiden", "nope")


def test_resolve_composition_keys_skips_absent_and_constant_columns():
    obs = pd.DataFrame({
        "tech": ["a", "b", "a", "b"],
        "donor_id": ["D1", "D2", "D1", "D2"],
        "constant": ["x", "x", "x", "x"],
    })

    keys = resolve_composition_keys(obs, batch_key="tech", sample_key="donor_id",
                                    covariate_keys=["constant", "absent", "tech"])

    assert keys == ["tech", "donor_id"]     # deduped, constant and absent dropped
    assert resolve_composition_keys(obs, explicit=["donor_id"]) == ["donor_id"]
    assert resolve_composition_keys(obs, explicit=[]) == []


def test_auto_flag_needs_both_significance_and_effect_size():
    pytest.importorskip("scanpy")
    from sccairns.inspect import auto_flag_clusters

    summary = pd.DataFrame({
        "n_cells": [400, 400, 400],
        "tech_depleted": ["SSv4", "SSv4", "SSv4"],
        "tech_depletion_ratio": [0.0, 0.95, 0.0],
        "tech_depletion_q": [1e-20, 1e-20, 0.5],
        "tech_top": ["10x", "10x", "10x"],
        "tech_enrichment": [1.1, 1.1, 1.1],
        "tech_enrich_q": [1.0, 1.0, 1.0],
    }, index=["real", "significant_but_tiny_effect", "big_effect_but_not_significant"])

    flags = auto_flag_clusters(summary, min_cells=1, composition_keys=["tech"],
                               min_neighbor_purity=None)

    assert "Depleted of a tech level" in flags["real"][0]
    assert "significant_but_tiny_effect" not in flags
    assert "big_effect_but_not_significant" not in flags


def test_composition_flagging_is_disablable():
    pytest.importorskip("scanpy")
    from sccairns.inspect import auto_flag_clusters

    summary = pd.DataFrame({
        "n_cells": [400],
        "tech_depleted": ["SSv4"], "tech_depletion_ratio": [0.0],
        "tech_depletion_q": [1e-20],
        "tech_top": ["10x"], "tech_enrichment": [1.1], "tech_enrich_q": [1.0],
    }, index=["c0"])

    assert auto_flag_clusters(summary, min_cells=1, composition_keys=["tech"],
                              composition_q=None,
                              min_neighbor_purity=None) == {}


# ── Display formatting ────────────────────────────────────────────────────────

def test_q_values_keep_their_exponent_in_the_report():
    """Regression: the summary table used to be `.round(2)`-ed for HTML, which
    rendered every significant q-value as 0.0 — wrong, and indistinguishable
    from a marginal one."""
    pytest.importorskip("scanpy")
    from sccairns.inspect import summary_table_formatters

    df = pd.DataFrame({
        "n_cells": [300, 250, 400],
        "median_pct_mt": [3.14159, 2.71828, 4.5],
        "tech_depleted": ["SSv4", "10x", "SSv4"],
        "tech_depletion_q": [6.81e-20, 1.0, 0.00229],
        "donor_id_enrich_q": [0.00299, 1.62e-97, 0.0],
    }, index=["2", "3", "0"])

    formatters = summary_table_formatters(df)
    rendered = df.to_string(formatters=formatters)

    assert "6.81e-20" in rendered
    assert "1.62e-97" in rendered
    assert "<1e-300" in rendered          # underflow is not reported as zero
    assert "0.002" in rendered            # above 1e-3 stays decimal
    assert "3.14" in rendered             # other floats still get 2 decimals
    # Non-float columns are left alone.
    assert "SSv4" in rendered
    assert "n_cells" not in formatters


def test_formatters_apply_to_the_html_report_too():
    pytest.importorskip("scanpy")
    from sccairns.inspect import summary_table_formatters

    df = pd.DataFrame({"tech_depletion_q": [6.81e-20]}, index=["2"])
    html = df.to_html(formatters=summary_table_formatters(df))

    assert "6.81e-20" in html
    assert ">0.0<" not in html

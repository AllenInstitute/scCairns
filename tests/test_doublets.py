"""Per-cluster doublet-score reporting.

Nothing here computes doublet scores — they arrive from upstream (DoubletFinder
in R, or a future in-pipeline caller). These cover the reporting path, whose main
hazard is that scores exist for only some platforms.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")
pytest.importorskip("scanpy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import anndata as ad  # noqa: E402

from sccairns.inspect import _as_doublet_calls, plot_doublet_score_by_cluster  # noqa: E402


def _adata(n=400, seed=0, scored_platform_only=True):
    """Droplet + plate object where only the droplet arm carries a score."""
    rng = np.random.default_rng(seed)
    obj = ad.AnnData(np.ones((n, 5)))
    obj.obs_names = [f"cell{i}" for i in range(n)]
    obj.obs["leiden"] = pd.Categorical(rng.choice(["0", "1", "2", "10"], n))
    obj.obs["tech"] = pd.Categorical(
        rng.choice(["10x_Multiome", "SSv4"], n, p=[0.75, 0.25]))

    score = np.full(n, np.nan)
    droplet = (obj.obs["tech"] == "10x_Multiome").to_numpy()
    target = droplet if scored_platform_only else np.ones(n, bool)
    score[target] = rng.beta(2, 20, target.sum())
    rich = target & (obj.obs["leiden"] == "10").to_numpy()
    score[rich] = rng.beta(6, 6, rich.sum())
    obj.obs["doublet_score"] = score
    return obj


@pytest.mark.parametrize("column,expected", [
    (pd.Series([True, False, True]), [True, False, True]),
    (pd.Series([1, 0, 1]), [True, False, True]),
    (pd.Series([0.0, 1.0]), [False, True]),
    # DoubletFinder writes strings; a numeric coercion would make these all False.
    (pd.Series(["Doublet", "Singlet", "Doublet"]), [True, False, True]),
    (pd.Series(["doublet", "singlet"]), [True, False]),
])
def test_doublet_call_columns_are_coerced_from_every_common_form(column, expected):
    assert list(_as_doublet_calls(column)) == expected


def test_absent_call_column_is_all_false():
    assert _as_doublet_calls(None) is False


def test_plot_is_skipped_when_no_score_exists(tmp_path, capsys):
    obj = _adata()
    del obj.obs["doublet_score"]

    assert plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path)) is None
    assert "No 'doublet_score' in obs" in capsys.readouterr().out


def test_plot_is_skipped_when_the_score_is_entirely_missing(tmp_path, capsys):
    obj = _adata()
    obj.obs["doublet_score"] = np.nan

    assert plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path)) is None
    assert "no finite values" in capsys.readouterr().out


def test_unscored_platforms_are_excluded_not_counted_as_zero(tmp_path):
    """The point of the panel: a plate-based arm has no score and must not
    dilute a cluster's distribution."""
    obj = _adata()
    scored = obj.obs["doublet_score"].notna()
    # Only droplet cells carry a score.
    assert set(obj.obs.loc[scored, "tech"].unique()) == {"10x_Multiome"}
    assert not scored.all()

    b64 = plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path),
                                        platform_key="tech")

    assert b64 is not None
    assert (tmp_path / "doublet_score_by_cluster.png").exists()


def test_call_key_path_runs_with_a_doubletfinder_string_column(tmp_path):
    """Regression: the call fraction was computed by indexing a cell-name-indexed
    Series with a positionally-indexed mask, which raised IndexingError on every
    real run (the report passes call_key by default)."""
    obj = _adata()
    score = obj.obs["doublet_score"].to_numpy()
    obj.obs["DF.classifications"] = np.where(
        np.isnan(score), "", np.where(score > 0.35, "Doublet", "Singlet"))

    b64 = plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path),
                                        platform_key="tech",
                                        call_key="DF.classifications")

    assert b64 is not None


def test_missing_call_column_does_not_prevent_the_plot(tmp_path):
    obj = _adata()

    b64 = plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path),
                                        call_key="not_a_column")

    assert b64 is not None


def test_non_numeric_score_column_is_coerced(tmp_path):
    obj = _adata()
    obj.obs["doublet_score"] = obj.obs["doublet_score"].astype(str)

    assert plot_doublet_score_by_cluster(obj, "leiden", str(tmp_path)) is not None

"""Tests for the QC gene/cell filter default and round-aware gating.

Covers two changes:
  1. Default qc.min_cells is 1 (keep every observed gene), so rare markers are
     not dropped by the gene filter.
  2. skip_filter_after_round_1 is honored: the QC filter runs on the first
     round only, so genes are not re-filtered against later rounds' smaller
     cell subsets.

should_apply_qc_filter lives in sccairns.integrate, which imports scvi at
module top; the function itself does not use scvi, so we stub it for import.
"""

import sys
import types

import pytest

from sccairns.config import DEFAULT_CONFIG


@pytest.fixture
def stubbed_scvi():
    had = "scvi" in sys.modules
    saved = sys.modules.get("scvi")
    fake = types.ModuleType("scvi")
    fake.settings = types.SimpleNamespace(seed=None)
    sys.modules["scvi"] = fake
    try:
        yield
    finally:
        sys.modules.pop("sccairns.integrate", None)
        if had:
            sys.modules["scvi"] = saved
        else:
            sys.modules.pop("scvi", None)


def test_default_min_cells_is_one():
    assert DEFAULT_CONFIG["qc"]["min_cells"] == 1
    assert DEFAULT_CONFIG["qc"]["skip_filter_after_round_1"] is True


def test_first_round_filter_applies(stubbed_scvi):
    from sccairns.integrate import should_apply_qc_filter

    qc = {"enabled": True, "skip_filter": False,
          "skip_filter_after_round_1": True}
    apply_filter, skipped_later = should_apply_qc_filter(qc, is_first_round=True)
    assert apply_filter is True
    assert skipped_later is False


def test_later_round_skips_when_flag_set(stubbed_scvi):
    from sccairns.integrate import should_apply_qc_filter

    qc = {"enabled": True, "skip_filter": False,
          "skip_filter_after_round_1": True}
    apply_filter, skipped_later = should_apply_qc_filter(qc, is_first_round=False)
    assert apply_filter is False
    assert skipped_later is True  # reports it was skipped for being a later round


def test_later_round_still_filters_when_flag_off(stubbed_scvi):
    from sccairns.integrate import should_apply_qc_filter

    qc = {"enabled": True, "skip_filter": False,
          "skip_filter_after_round_1": False}
    apply_filter, skipped_later = should_apply_qc_filter(qc, is_first_round=False)
    assert apply_filter is True
    assert skipped_later is False


def test_global_skip_filter_wins_on_first_round(stubbed_scvi):
    from sccairns.integrate import should_apply_qc_filter

    qc = {"enabled": True, "skip_filter": True,
          "skip_filter_after_round_1": True}
    apply_filter, skipped_later = should_apply_qc_filter(qc, is_first_round=True)
    assert apply_filter is False
    # Not attributed to "later round" — it was globally skipped.
    assert skipped_later is False


def test_disabled_qc_never_filters(stubbed_scvi):
    from sccairns.integrate import should_apply_qc_filter

    qc = {"enabled": False, "skip_filter": False,
          "skip_filter_after_round_1": True}
    for is_first in (True, False):
        apply_filter, skipped_later = should_apply_qc_filter(qc, is_first)
        assert apply_filter is False
        assert skipped_later is False

"""The configured seed must reach the embedding/clustering RNG.

`set_global_seed()` seeds scvi-tools (and therefore scVI/scANVI training), but
the neighbors/UMAP/Leiden wrappers used to call Scanpy without `random_state`,
so they silently ran on Scanpy's implicit default of 0 no matter what
`reproducibility.seed` said. A manifest recording `seed: 42` therefore did not
describe the RNG state of the steps that produce Leiden cluster IDs — the very
IDs a decisions.yaml references by number.

These tests capture what actually reaches Scanpy. Each fake raises a sentinel
so the wrapper short-circuits before touching graph state we have not built;
the wrappers catch only TypeError/KeyError/ValueError, so the sentinel
propagates cleanly.

sccairns.integrate imports scvi at module top (heavy, CUDA); it is stubbed for
the duration of these tests, mirroring tests/test_harmony.py.
"""

import sys
import types

import pytest


class _Captured(Exception):
    """Sentinel: stop the wrapper once we've seen the Scanpy kwargs."""


@pytest.fixture
def isi():
    """sccairns.integrate with the heavy scvi import stubbed out."""
    had = "scvi" in sys.modules
    saved = sys.modules.get("scvi")
    fake = types.ModuleType("scvi")
    fake.settings = types.SimpleNamespace(seed=None)
    sys.modules["scvi"] = fake
    try:
        from sccairns import integrate as module

        yield module
    finally:
        sys.modules.pop("sccairns.integrate", None)
        if had:
            sys.modules["scvi"] = saved
        else:
            sys.modules.pop("scvi", None)


@pytest.fixture
def adata():
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    ad = pytest.importorskip("anndata")

    n = 20
    obj = ad.AnnData(
        X=np.random.RandomState(0).rand(n, 5).astype("float32"),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
    )
    obj.obsm["X_pca"] = np.random.RandomState(1).rand(n, 4)
    return obj


def _capture(monkeypatch, target, attr):
    """Replace target.attr with a recorder that aborts via _Captured."""
    seen = {}

    def fake(*args, **kwargs):
        seen.update(kwargs)
        raise _Captured

    monkeypatch.setattr(target, attr, fake)
    return seen


# ── the configured seed is what reaches Scanpy ──────────────────────────────
@pytest.mark.parametrize("seed", [42, 7, 0])
def test_neighbors_threads_the_configured_seed(isi, adata, monkeypatch, seed):
    seen = _capture(monkeypatch, isi.sc.pp, "neighbors")
    with pytest.raises(_Captured):
        isi.run_neighbors_with_key(adata, "nbr", seed=seed, use_rep="X_pca")
    assert seen["random_state"] == seed


@pytest.mark.parametrize("seed", [42, 7, 0])
def test_umap_threads_the_configured_seed(isi, adata, monkeypatch, seed):
    seen = _capture(monkeypatch, isi.sc.tl, "umap")
    with pytest.raises(_Captured):
        isi.run_umap_with_key(adata, neighbors_key="nbr", umap_key="X_umap_x",
                              seed=seed)
    assert seen["random_state"] == seed


@pytest.mark.parametrize("seed", [42, 7, 0])
def test_leiden_threads_the_configured_seed(isi, adata, monkeypatch, seed):
    """Leiden matters most: its labels are what decisions.yaml cites."""
    seen = _capture(monkeypatch, isi.sc.tl, "leiden")
    with pytest.raises(_Captured):
        isi.run_leiden_with_key(adata, neighbors_key="nbr",
                                leiden_key="leiden_x", seed=seed)
    assert seen["random_state"] == seed


# ── seed: null normalizes to 0, matching run_harmony ────────────────────────
@pytest.mark.parametrize("wrapper,module_attr,attr,kwargs", [
    ("run_neighbors_with_key", "pp", "neighbors", {"use_rep": "X_pca"}),
    ("run_umap_with_key", "tl", "umap", {"umap_key": "X_umap_x"}),
    ("run_leiden_with_key", "tl", "leiden", {"leiden_key": "leiden_x"}),
])
def test_unset_seed_normalizes_to_zero(isi, adata, monkeypatch, wrapper,
                                       module_attr, attr, kwargs):
    """An unseeded run stays deterministic here rather than drifting.

    Same convention as run_harmony's `seed if seed is not None else 0`.
    """
    seen = _capture(monkeypatch, getattr(isi.sc, module_attr), attr)
    call = getattr(isi, wrapper)
    args = (adata, "nbr") if wrapper == "run_neighbors_with_key" else (adata,)
    if wrapper != "run_neighbors_with_key":
        kwargs = {"neighbors_key": "nbr", **kwargs}
    with pytest.raises(_Captured):
        call(*args, seed=None, **kwargs)
    assert seen["random_state"] == 0


# ── an explicit caller-supplied random_state still wins ─────────────────────
def test_explicit_random_state_is_not_overridden(isi, adata, monkeypatch):
    """run_neighbors_with_key forwards **kwargs; setdefault must not clobber."""
    seen = _capture(monkeypatch, isi.sc.pp, "neighbors")
    with pytest.raises(_Captured):
        isi.run_neighbors_with_key(adata, "nbr", seed=42, use_rep="X_pca",
                                   random_state=99)
    assert seen["random_state"] == 99

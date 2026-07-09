"""Regression tests for the optional Harmony integration path.

Covers:
  1. Config defaults + validation (single-model only; n_pcs sanity).
  2. run_harmony — produces a batch-corrected embedding, is deterministic
     under a fixed seed, and does not mutate its input AnnData.
  3. Inspection auto-resolves X_pca_harmony as a latent candidate.

run_harmony lives in integrate_scvi, which imports scvi at module top.
scvi is heavy/CUDA and not needed by run_harmony or the validators, so we stub
it in sys.modules before import. The test is skipped if harmonypy is absent.
"""

import sys
import types
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from pipeline_config import (  # noqa: E402
    ConfigError,
    DEFAULT_CONFIG,
    load_pipeline_config,
    validate_config,
)


@pytest.fixture
def stubbed_scvi():
    """Temporarily stub the heavy scvi import for run_harmony tests.

    integrate_scvi imports scvi at module top; run_harmony itself does not
    use it. We install a stub only for the duration of the test and remove it
    afterward so the fake module cannot leak into other test modules (e.g.
    set_global_seed, which imports scvi and would otherwise treat the stub as a
    real seeding backend and skip its numpy/random fallback).
    """
    had = "scvi" in sys.modules
    saved = sys.modules.get("scvi")
    fake = types.ModuleType("scvi")
    fake.settings = types.SimpleNamespace(seed=None)
    sys.modules["scvi"] = fake
    try:
        yield
    finally:
        # Drop the cached integrate module so a later real import re-runs.
        sys.modules.pop("integrate_scvi", None)
        if had:
            sys.modules["scvi"] = saved
        else:
            sys.modules.pop("scvi", None)


def _base_config():
    cfg = load_pipeline_config(None)
    cfg["data"]["input_h5ad"] = "x.h5ad"
    cfg["data"]["output_dir"] = "out"
    return cfg


# ── Config ──────────────────────────────────────────────────────────────────
def test_default_config_has_harmony_block():
    h = DEFAULT_CONFIG["integration"]["harmony"]
    assert h["enabled"] is False
    assert h["batch_key"] is None
    assert h["n_pcs"] == 30
    assert h["hvg_from"] is None


# ── HVG-source resolution (the explicit decision that lifts the sweep ban) ───
def test_resolve_harmony_hvg_spec_default_uses_top_level(stubbed_scvi):
    import integrate_scvi as isi  # noqa: E402

    integration_cfg = {
        "hvg": {"n_top_genes": 3000, "batch_key": "tech",
                "flavor": "seurat_v3", "min_batches": None},
    }
    spec = isi.resolve_harmony_hvg_spec(
        {"enabled": True}, integration_cfg, configs={}, sweep_mode=False)
    assert spec["source"] == "integration.hvg"
    assert spec["batch_key"] == "tech"
    assert spec["flavor"] == "seurat_v3"


def test_resolve_harmony_hvg_spec_borrows_sweep_entry(stubbed_scvi):
    import types as _types

    import integrate_scvi as isi  # noqa: E402

    integration_cfg = {
        "hvg": {"n_top_genes": 3000, "batch_key": "tech",
                "flavor": "seurat_v3", "min_batches": None},
    }
    # A stand-in for the sweep ModelConfig with its own HVG selection.
    sweep_entry = _types.SimpleNamespace(
        hvg_batch_key="platform_origin", hvg_flavor="seurat", hvg_nbatches=3)
    spec = isi.resolve_harmony_hvg_spec(
        {"enabled": True, "hvg_from": "medium"},
        integration_cfg, configs={"medium": sweep_entry}, sweep_mode=True)
    assert spec["source"] == "sweep:medium"
    assert spec["batch_key"] == "platform_origin"
    assert spec["flavor"] == "seurat"
    assert spec["min_batches"] == 3
    # n_top_genes still comes from the top-level spec.
    assert spec["n_top_genes"] == 3000


def test_resolve_harmony_hvg_spec_unknown_entry_raises(stubbed_scvi):
    import integrate_scvi as isi  # noqa: E402

    integration_cfg = {"hvg": {"n_top_genes": 3000, "batch_key": "tech",
                               "flavor": "seurat_v3", "min_batches": None}}
    with pytest.raises(ConfigError, match="does not match a sweep"):
        isi.resolve_harmony_hvg_spec(
            {"enabled": True, "hvg_from": "nope"},
            integration_cfg, configs={"medium": object()}, sweep_mode=True)


def test_harmony_validation_single_model_ok():
    cfg = _base_config()
    cfg["integration"]["harmony"]["enabled"] = True
    validate_config(cfg, mode="integration")  # must not raise


def test_harmony_allowed_with_sweep():
    # Harmony now runs once per round alongside a sweep (no restriction).
    cfg = _base_config()
    cfg["integration"]["harmony"]["enabled"] = True
    cfg["integration"]["sweep"] = [{"name": "small"}]
    validate_config(cfg, mode="integration")  # must not raise


def test_harmony_hvg_from_requires_sweep():
    # hvg_from names a sweep entry; it is invalid in single-model mode.
    cfg = _base_config()
    cfg["integration"]["harmony"]["enabled"] = True
    cfg["integration"]["harmony"]["hvg_from"] = "small"
    with pytest.raises(ConfigError, match="hvg_from"):
        validate_config(cfg, mode="integration")

    # With a matching sweep entry present it validates.
    cfg["integration"]["sweep"] = [{"name": "small"}]
    validate_config(cfg, mode="integration")  # must not raise


def test_harmony_rejects_bad_n_pcs():
    cfg = _base_config()
    cfg["integration"]["harmony"]["enabled"] = True
    cfg["integration"]["harmony"]["n_pcs"] = 1
    with pytest.raises(ConfigError, match="n_pcs"):
        validate_config(cfg, mode="integration")


def test_harmony_disabled_by_default_validates():
    # The default config (harmony disabled) must validate even with a sweep.
    cfg = _base_config()
    cfg["integration"]["sweep"] = [{"name": "small"}]
    validate_config(cfg, mode="integration")  # must not raise


# ── run_harmony ───────────────────────────────────────────────────────────────
def test_run_harmony_embedding_deterministic_and_pure(stubbed_scvi):
    pytest.importorskip("harmonypy")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    ad = pytest.importorskip("anndata")
    import integrate_scvi as isi  # noqa: E402

    rng = np.random.RandomState(0)
    n, g = 100, 40
    X = rng.poisson(2.0, size=(n, g)).astype("float32")
    batch = np.array((["a"] * 50) + (["b"] * 50))
    X[batch == "b"] += rng.poisson(1.0, size=(50, g)).astype("float32")
    adata = ad.AnnData(
        X=X.copy(),
        obs=pd.DataFrame({"tech": batch}, index=[f"c{i}" for i in range(n)]),
    )
    adata.layers["counts"] = X.copy()

    emb1 = isi.run_harmony(adata, batch_key="tech", counts_layer="counts",
                           n_pcs=10, seed=7)
    emb2 = isi.run_harmony(adata, batch_key="tech", counts_layer="counts",
                           n_pcs=10, seed=7)

    assert emb1.shape == (n, 10)
    assert np.allclose(emb1, emb2)  # deterministic under fixed seed
    # run_harmony works on a copy — the caller's AnnData is not mutated.
    assert "X_pca_harmony" not in adata.obsm


# ── Inspection auto-resolution ───────────────────────────────────────────────
def test_inspection_resolves_harmony_latent():
    pytest.importorskip("scanpy")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    ad = pytest.importorskip("anndata")
    from inspect_integration import resolve_inspection_keys  # noqa: E402

    n = 30
    adata = ad.AnnData(
        X=np.random.RandomState(0).rand(n, 5).astype("float32"),
        obs=pd.DataFrame(
            {"leiden_harmony": pd.Categorical(["0", "1", "2"] * 10)},
            index=[f"c{i}" for i in range(n)],
        ),
    )
    adata.obsm["X_pca_harmony"] = np.random.RandomState(1).rand(n, 10)
    adata.obsm["X_umap_harmony"] = np.random.RandomState(2).rand(n, 2)

    cluster_key, latent_key, umap_key = resolve_inspection_keys(
        adata, {"cluster_key": "auto", "latent_key": "auto", "umap_key": "auto"},
        {"enabled": False})

    assert cluster_key == "leiden_harmony"
    assert latent_key == "X_pca_harmony"
    assert umap_key == "X_umap_harmony"

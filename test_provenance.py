"""Regression tests for the reproducibility/traceability provenance helpers.

Covers the three traceability fixes:
  1. set_global_seed       — deterministic RNG seeding
  2. collect_code_provenance — git revision capture stamped into the manifest
  3. fingerprint_file / read_parent_provenance — input fingerprints and
     explicit cross-round lineage pointers
"""

import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from pipeline_config import (  # noqa: E402
    DEFAULT_CONFIG,
    collect_code_provenance,
    fingerprint_file,
    read_parent_provenance,
    set_global_seed,
    update_round_manifest,
)


def test_default_config_has_seed():
    assert "reproducibility" in DEFAULT_CONFIG
    assert DEFAULT_CONFIG["reproducibility"]["seed"] == 0


def test_set_global_seed_is_deterministic_and_nullable():
    np = pytest.importorskip("numpy")

    assert set_global_seed(None) is None  # no-op, returns None

    assert set_global_seed(123) == 123
    first = np.random.rand(5)
    assert set_global_seed(123) == 123
    second = np.random.rand(5)
    assert np.allclose(first, second)


def test_collect_code_provenance_shape():
    prov = collect_code_provenance()
    # The repo is a git checkout, so a commit should be captured here. Even
    # outside a checkout the keys must always be present (values may be None).
    for key in ("commit", "commit_short", "branch", "commit_time", "dirty"):
        assert key in prov
    if prov["commit"] is not None:
        assert prov["commit_short"] == prov["commit"][:12]


def test_fingerprint_file_roundtrip(tmp_path):
    f = tmp_path / "input.h5ad"
    f.write_bytes(b"some bytes")
    fp = fingerprint_file(str(f))
    assert fp is not None
    assert fp["size_bytes"] == len(b"some bytes")
    assert fp["partial"] is False
    # SHA-256 of the exact content is stable and verifiable.
    import hashlib

    assert fp["sha256"] == hashlib.sha256(b"some bytes").hexdigest()


def test_fingerprint_file_missing_returns_none():
    assert fingerprint_file(None) is None
    assert fingerprint_file("/no/such/path.h5ad") is None


def test_read_parent_provenance_links_round(tmp_path):
    parent_dir = tmp_path / "round_01"
    parent_dir.mkdir()
    (parent_dir / "round_manifest.json").write_text(
        json.dumps(
            {
                "created": "t0",
                "updated": "t1",
                "code": {"commit": "deadbeef"},
                "reproducibility": {"seed": 7},
            }
        ),
        encoding="utf-8",
    )
    child_input = parent_dir / "filtered.h5ad"
    child_input.write_bytes(b"")

    ptr = read_parent_provenance(str(child_input))
    assert ptr is not None
    assert ptr["input_file"] == "filtered.h5ad"
    assert ptr["commit"] == "deadbeef"
    assert ptr["seed"] == 7

    # No sibling manifest -> no pointer.
    orphan = tmp_path / "elsewhere.h5ad"
    orphan.write_bytes(b"")
    assert read_parent_provenance(str(orphan)) is None
    assert read_parent_provenance(None) is None


def test_update_round_manifest_stamps_code_and_seed(tmp_path):
    out = tmp_path / "round_02"
    mp = update_round_manifest(str(out), "integration", {"n_cells": 10}, seed=42)
    manifest = json.loads(Path(mp).read_text(encoding="utf-8"))

    assert "code" in manifest
    assert manifest["reproducibility"]["seed"] == 42
    assert manifest["integration"]["n_cells"] == 10

    # A later seedless stage update must not erase the recorded seed.
    update_round_manifest(str(out), "inspection", {"n_clusters": 3})
    manifest2 = json.loads(Path(mp).read_text(encoding="utf-8"))
    assert manifest2["reproducibility"]["seed"] == 42
    assert manifest2["inspection"]["n_clusters"] == 3

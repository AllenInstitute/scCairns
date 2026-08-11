"""Regression tests for the reproducibility/traceability provenance helpers.

Covers the three traceability fixes:
  1. set_global_seed       — deterministic RNG seeding
  2. collect_code_provenance — git revision capture stamped into the manifest
  3. fingerprint_file / read_parent_provenance — input fingerprints and
     explicit cross-round lineage pointers
"""

import json
from pathlib import Path

import pytest

from sccairns.config import (
    DEFAULT_CONFIG,
    collect_code_provenance,
    fingerprint_file,
    read_parent_provenance,
    set_global_seed,
    update_round_manifest,
    write_command_args,
    write_yaml_record,
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
    for key in ("version", "source", "commit", "commit_short", "branch",
                "commit_time", "repo_root", "dirty"):
        assert key in prov
    if prov["commit"] is not None:
        assert prov["commit_short"] == prov["commit"][:12]
        assert prov["source"] == "checkout"


def test_collect_code_provenance_outside_a_checkout(tmp_path, monkeypatch):
    """A wheel in site-packages has no .git — the version must still identify it.

    This is the deployment where scCairns is pip-installed into a Code Ocean
    image rather than vendored in the capsule, so it is the case that must not
    silently produce an unidentifiable round.
    """
    from sccairns import config as cfg

    monkeypatch.setattr(cfg, "package_version", lambda: "0.1.0")
    prov = cfg.collect_code_provenance(repo_dir=str(tmp_path))

    assert prov["commit"] is None
    assert prov["repo_root"] is None
    assert prov["version"] == "0.1.0"
    assert prov["source"] == "installed"


def test_collect_code_provenance_unidentified(tmp_path, monkeypatch):
    """Neither git nor an installed distribution: say so rather than imply a match."""
    from sccairns import config as cfg

    monkeypatch.setattr(cfg, "package_version", lambda: None)
    prov = cfg.collect_code_provenance(repo_dir=str(tmp_path))

    assert prov["version"] is None
    assert prov["source"] == "unknown"


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


# ── Fix: append-only invocation history in command_args.json ────────────────
def test_command_args_keeps_invocation_history(tmp_path):
    out = tmp_path / "round_01"
    write_command_args(str(out), "inspect", {"input": "a.h5ad"})
    write_command_args(str(out), "inspect", {"input": "b.h5ad"})
    write_command_args(str(out), "integrate", {"input": "b.h5ad"})

    payload = json.loads(
        (out / "command_args.json").read_text(encoding="utf-8"))

    # history records every invocation in order (re-runs are not lost).
    assert len(payload["history"]) == 3
    assert [h["stage"] for h in payload["history"]] == [
        "inspect", "inspect", "integrate"]
    assert payload["history"][0]["args"]["input"] == "a.h5ad"
    assert payload["history"][1]["args"]["input"] == "b.h5ad"
    # commands[stage] still holds the most recent invocation per stage.
    assert payload["commands"]["inspect"]["args"]["input"] == "b.h5ad"


# ── Fix: decisions records serialized with safe YAML ────────────────────────
def test_write_yaml_record_escapes_special_chars(tmp_path):
    yaml = pytest.importorskip("yaml")
    path = tmp_path / "decisions_applied.yaml"
    record = {
        "timestamp": "2026-06-04T12:00:00",
        "n_removed": 13,
        "notes": 'review: contains a "quote", a colon: and\na newline',
        "actions": [
            {
                "type": "remove_query",
                "query": "neuron_type.isin(['GABAergic']) & (tech == 'scale')",
                "reason": 'dropped: it had "issues"',
                "n_removed": 13,
            }
        ],
    }
    write_yaml_record(str(path), record)

    # The whole point: it round-trips back to the exact structure, which the
    # old hand-built f-string serializer could not guarantee.
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded == record


# ── Fix: provenance footer rendered into the HTML report ────────────────────
def test_render_provenance_section(tmp_path):
    pytest.importorskip("scanpy")  # sccairns.inspect imports scanpy
    from sccairns.inspect import _render_provenance_section

    assert _render_provenance_section(None) == ""
    assert _render_provenance_section({}) == ""

    html = _render_provenance_section({
        "pipeline_version": 1,
        "seed": 7,
        "code": {"commit_short": "abc123def456", "branch": "main", "dirty": True},
        "input": {"path": "/data/integrated.h5ad",
                  "sha256": "a" * 64, "size_bytes": 1234},
        "config_sha256": "b" * 64,
        "parent_round": {"input_file": "filtered.h5ad",
                         "commit": "deadbeef0000", "seed": 99},
        "packages": {"scanpy": "1.10.4", "scvi-tools": "1.3.3"},
        "generated": "2026-06-04T12:00:00",
    })

    assert "Run Provenance" in html
    assert "abc123def456" in html
    assert "+uncommitted changes" in html  # dirty flag surfaced
    assert "filtered.h5ad" in html         # parent lineage surfaced
    assert "scanpy 1.10.4" in html
    # seed None renders as a clear "unset" label rather than a blank.
    assert "unset" in _render_provenance_section({"seed": None})

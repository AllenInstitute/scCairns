import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import summarize_rounds as sr  # noqa: E402


def _write_manifest(d: Path, manifest: dict):
    d.mkdir(parents=True, exist_ok=True)
    (d / "round_manifest.json").write_text(json.dumps(manifest))


def _round1(rounds: Path):
    """Root round: sweep with two scVI variants, no filtering."""
    (rounds / "round_01").mkdir(parents=True, exist_ok=True)
    (rounds / "round_01" / "scib_benchmark_results.csv").write_text(
        "Embedding,Batch correction,Bio conservation,Total\n"
        "X_scVI_small_gene_nb,0.72,0.45,0.88\n"
        "X_scVI_big_zinb,0.72,0.23,0.43\n"
        "X_pca,0.10,0.05,0.07\n"
        "Metric Type,Aggregate score,Aggregate score,Aggregate score\n"
    )
    integ_path = rounds / "round_01" / "integrated.h5ad"
    integ_path.write_text("integrated-bytes")
    _write_manifest(rounds / "round_01", {
        "pipeline_version": 1,
        "reproducibility": {"seed": 0},
        "code": {"commit": "aaaa", "commit_short": "aaaa", "dirty": False},
        "integration": {
            "input_h5ad": "/external/combined.h5ad",
            "parent_round": None,
            "output_h5ad": str(integ_path),
            "n_cells": 48210, "n_genes": 21000, "model_type": "scvi",
            "embedding_keys": ["X_scVI_small_gene_nb", "X_scVI_big_zinb",
                               "X_pca_harmony"],
            "package_versions": {"scvi-tools": "1.1.0", "torch": "2.2.0"},
        },
        "inspection": {"n_clusters": 24, "n_auto_flagged": 3},
    })
    return integ_path


def _round2(rounds: Path, integ_path: Path, *, variant="big_zinb",
            source="decisions_file", input_sha=None):
    sha = input_sha or sr._sha256(str(integ_path))
    _write_manifest(rounds / "round_02", {
        "pipeline_version": 1,
        "reproducibility": {"seed": 0},
        "code": {"commit": "aaaa", "commit_short": "aaaa", "dirty": False},
        "decisions": {
            "input_h5ad": str(integ_path),
            "input_fingerprint": {"path": str(integ_path), "sha256": sha},
            "parent_round": {"output_dir": str(integ_path.parent)},
            "cluster_key": f"leiden_{variant}",
            "filtered_on": {"variant": variant, "cluster_key": f"leiden_{variant}",
                            "source": source},
            "input_cells": 48210, "output_cells": 44933,
            "removed": 3277, "removal_pct": 6.8,
            "actions_summary": ["remove_cluster: 7 (-2000)"],
        },
        "integration": {
            "input_h5ad": str(rounds / "round_02" / "filtered.h5ad"),
            "parent_round": None, "n_cells": 44933, "n_genes": 21000,
            "model_type": "scvi",
            "embedding_keys": ["X_scVI_small_gene_nb", "X_scVI_big_zinb"],
            "package_versions": {"scvi-tools": "1.1.0", "torch": "2.2.0"},
        },
        "inspection": {"n_clusters": 19, "n_auto_flagged": 0},
    })


def _load_ordered(tmp_path):
    manifests = [sr.load_manifest(str(p)) for p in
                 sorted(tmp_path.glob("*/round_manifest.json"))]
    parent_map = sr.build_parent_map(manifests)
    return sr.order_by_lineage(manifests, parent_map), parent_map


def test_lineage_ordering_and_edges(tmp_path):
    integ = _round1(tmp_path)
    _round2(tmp_path, integ)
    manifests, parent_map = _load_ordered(tmp_path)
    assert [m["_round_id"] for m in manifests] == ["round_01", "round_02"]

    by_id = {m["_round_id"]: m for m in manifests}
    child = by_id["round_02"]
    assert parent_map[child["_dir"]] == by_id["round_01"]["_dir"]


def test_explicit_order_sequential_fallback(tmp_path):
    # Archived rounds: parent pointer references a path/name not in the set,
    # so lineage must fall back to the supplied order.
    _round1(tmp_path)
    integ_missing = tmp_path / "round_01" / "integrated.h5ad"
    _round2(tmp_path, integ_missing)
    # Rename round dirs to arbitrary archive names and rewrite parent pointer
    # to a now-nonexistent location to defeat realpath + basename matching.
    a = tmp_path / "260603_cmg_round1"
    b = tmp_path / "260605_cmg_round3"
    (tmp_path / "round_01").rename(a)
    (tmp_path / "round_02").rename(b)
    mb = json.loads((b / "round_manifest.json").read_text())
    mb["decisions"]["parent_round"] = {"output_dir": "/gone/round_xyz"}
    mb["decisions"]["input_h5ad"] = "/gone/round_xyz/integrated.h5ad"
    (b / "round_manifest.json").write_text(json.dumps(mb))

    manifests = [sr.load_manifest(str(a / "round_manifest.json")),
                 sr.load_manifest(str(b / "round_manifest.json"))]
    order = [m["_dir"] for m in manifests]
    parent_map = sr.build_parent_map(manifests, explicit_order=order)
    assert parent_map[manifests[1]["_dir"]] == manifests[0]["_dir"]
    # Without the explicit order, the broken pointer resolves to nothing.
    assert sr.build_parent_map(manifests)[manifests[1]["_dir"]] is None


def test_scib_selected_vs_best_from_parent(tmp_path):
    integ = _round1(tmp_path)
    _round2(tmp_path, integ, variant="big_zinb")  # non-top variant
    manifests, parent_map = _load_ordered(tmp_path)
    by_id = {m["_round_id"]: m for m in manifests}

    parent = by_id["round_01"]
    child = by_id["round_02"]
    row1 = sr.extract_row(parent, None)
    row2 = sr.extract_row(child, parent)

    # Round 1 own-integration quality.
    assert row1["scib_best"] == 0.88
    assert row1["scib_best_variant"] == "small_gene_nb"
    # Round 2 was filtered on big_zinb — selected vs best sourced from parent.
    assert row2["filtered_on"] == "big_zinb"
    assert row2["scib_selected"] == 0.43
    assert row2["scib_selected_best"] == 0.88
    assert row2["scib_selected_is_best"] is False
    assert "✗" in row2["scib_sel_vs_best"]


def test_verify_flags_drift_and_changed_input(tmp_path):
    integ = _round1(tmp_path)
    _round2(tmp_path, integ, source="auto", input_sha="0" * 64)  # wrong hash
    manifests, parent_map = _load_ordered(tmp_path)
    findings = sr.verify(manifests, parent_map, scan_dir=str(tmp_path))
    checks = {f["check"] for f in findings}
    assert "input_changed" in checks            # parent's file ≠ recorded sha
    assert "auto_variant" in checks             # ambiguous selection
    severities = {f["check"]: f["severity"] for f in findings}
    assert severities["input_changed"] == "error"


def test_clean_run_has_no_findings(tmp_path):
    integ = _round1(tmp_path)
    _round2(tmp_path, integ)
    manifests, parent_map = _load_ordered(tmp_path)
    findings = sr.verify(manifests, parent_map, scan_dir=str(tmp_path))
    assert findings == []

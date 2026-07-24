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
    # Round 1's CSV carries a full Total → no metric tag.
    assert row1["scib_metric"] == "Total"


def test_scib_batch_only_aggregate_fallback(tmp_path):
    # Annotation-free benchmark: no Total/Bio conservation, only Batch correction
    # (matches the real CSV shape). scib_scores must fall back to that aggregate.
    d = tmp_path / "round_01"
    d.mkdir(parents=True, exist_ok=True)
    (d / "scib_benchmark_results.csv").write_text(
        "Embedding,BRAS,Graph connectivity,PCR comparison,Batch correction\n"
        "X_scVI_large_geneXcell_cell_nb,0.7,0.8,0.6,0.71\n"
        "X_pca_harmony,0.9,0.85,0.7,0.83\n"
    )
    scores, metric = sr.scib_scores(str(d))
    assert metric == "Batch correction"
    assert scores == {"large_geneXcell_cell_nb": 0.71, "harmony": 0.83}
    best_var, best_score = sr._best_integrated(scores)
    assert (best_var, best_score) == ("harmony", 0.83)
    # The batch-only aggregate is tagged so it isn't read as a Total.
    assert sr._scib_metric_tag(metric) == " (batch)"
    assert sr._scib_metric_tag("Total") == ""


def _sweep_inspection(base: Path, variant: str, n_clusters: int, n_flagged=0):
    """Write a per-architecture sweep inspection manifest (inspection stage)."""
    _write_manifest(base / f"inspect_{variant}", {
        "pipeline_version": 1,
        "code": {"commit": "aaaa", "commit_short": "aaaa", "dirty": False},
        "inspection": {"cluster_key": f"leiden_{variant}",
                       "n_clusters": n_clusters, "n_auto_flagged": n_flagged},
    })


def test_clusters_from_sweep_inspection_manifests(tmp_path):
    # Rounds whose own manifest carries integration+decisions but NO inspection
    # stage; clusters/flagged live in per-architecture inspect_* sweep manifests.
    integ = _round1(tmp_path)
    # Strip the inline inspection stage so resolution must fall to the sweep.
    m1 = json.loads((tmp_path / "round_01" / "round_manifest.json").read_text())
    del m1["inspection"]
    (tmp_path / "round_01" / "round_manifest.json").write_text(json.dumps(m1))
    # Round 1's sweep wrote inspect_* subdirs INSIDE the round dir.
    _sweep_inspection(tmp_path / "round_01", "small_gene_nb", 24, 3)
    _sweep_inspection(tmp_path / "round_01", "big_zinb", 30, 1)

    _round2(tmp_path, integ, variant="small_gene_nb")
    m2 = json.loads((tmp_path / "round_02" / "round_manifest.json").read_text())
    del m2["inspection"]
    (tmp_path / "round_02" / "round_manifest.json").write_text(json.dumps(m2))
    # Round 2's sweep wrote to a SIBLING <round>_sweep dir, one architecture.
    _sweep_inspection(tmp_path / "round_02_sweep", "harmony", 12, 0)

    manifests, parent_map = _load_ordered(tmp_path)
    by_dir = {m["_dir"]: m for m in manifests}
    acted_on = {}
    for m in manifests:
        pdir = parent_map.get(m["_dir"])
        v = sr._g(m, "decisions", "filtered_on", "variant")
        if pdir and v and pdir not in acted_on:
            acted_on[pdir] = v

    by_id = {m["_round_id"]: m for m in manifests}
    p, c = by_id["round_01"], by_id["round_02"]
    row1 = sr.extract_row(p, None, acted_on_variant=acted_on.get(p["_dir"]))
    row2 = sr.extract_row(c, p, acted_on_variant=acted_on.get(c["_dir"]))

    # Round 1: child filtered on small_gene_nb → that architecture's clusters.
    assert row1["n_clusters"] == 24
    assert row1["n_auto_flagged"] == 3
    assert row1["clusters_variant"] == "small_gene_nb"
    # Round 2 (terminal): no child; sole sweep architecture is the fallback.
    assert row2["n_clusters"] == 12
    assert row2["clusters_variant"] == "harmony"


def test_terminal_run_architecture_flag(tmp_path):
    # Terminal round with multiple sweep architectures and no scIB tiebreak:
    # clusters are ambiguous until the user pins one.
    integ = _round1(tmp_path)
    m1 = json.loads((tmp_path / "round_01" / "round_manifest.json").read_text())
    del m1["inspection"]
    (tmp_path / "round_01" / "round_manifest.json").write_text(json.dumps(m1))
    _sweep_inspection(tmp_path / "round_01", "small_gene_nb", 24, 3)

    # Filter this round on a variant NOT present in its own sweep, so neither
    # the own-filtered-variant nor the sole-architecture fallback can resolve it.
    _round2(tmp_path, integ, variant="big_zinb")
    m2 = json.loads((tmp_path / "round_02" / "round_manifest.json").read_text())
    del m2["inspection"]
    (tmp_path / "round_02" / "round_manifest.json").write_text(json.dumps(m2))
    _sweep_inspection(tmp_path / "round_02_sweep", "harmony", 12, 0)
    _sweep_inspection(tmp_path / "round_02_sweep", "other", 15, 1)

    manifests, parent_map = _load_ordered(tmp_path)
    by_id = {m["_round_id"]: m for m in manifests}
    p, c = by_id["round_01"], by_id["round_02"]

    # Without a pin, the terminal round is ambiguous (two archs, no tiebreak).
    assert sr.extract_row(c, p, acted_on_variant=None)["n_clusters"] is None

    # Pinned by variant name or cluster key — both resolve.
    for key in ("harmony", "leiden_harmony"):
        row = sr.extract_row(c, p, acted_on_variant=sr.variant_from_cluster_key(key))
        assert row["n_clusters"] == 12
        assert row["clusters_variant"] == "harmony"


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


def test_html_report_includes_full_featured_sections(tmp_path):
    integ = _round1(tmp_path)
    _round2(tmp_path, integ)
    manifests, parent_map = _load_ordered(tmp_path)
    by_dir = {m["_dir"]: m for m in manifests}

    rows = []
    for m in manifests:
        pdir = parent_map.get(m["_dir"])
        rows.append(sr.extract_row(m, by_dir.get(pdir) if pdir else None))
    rows_by_dir = {m["_dir"]: r for m, r in zip(manifests, rows)}
    mermaid_src = sr.render_mermaid(manifests, rows_by_dir, parent_map)
    ledger = sr.extract_decisions_ledger(manifests)

    html = sr.render_html(rows, mermaid_src, [], True, ledger)

    assert "Cells retained" in html
    assert "Retention" in html
    assert "Round details" in html
    assert "Decisions ledger" in html
    assert "remove_cluster: 7 (-2000)" in html
    assert "48,210 start cells to 44,933 final cells" in html
    assert "No integrity issues detected" in html


def test_synthesize_inspection_only_round(tmp_path):
    # A round that only re-inspected a prior round's output across a sweep: it
    # has inspect_<arch>/round_manifest.json sidecars but NO top-level
    # round_manifest.json. It should be recognized as a round, labeled by its
    # directory, with clusters sourced from the acted-on architecture.
    insp_dir = tmp_path / "round_03_v3"
    _sweep_inspection(insp_dir, "harmony", 14, 2)
    _sweep_inspection(insp_dir, "small_gene_nb", 20, 5)
    assert not (insp_dir / "round_manifest.json").exists()

    syn = sr.synthesize_inspection_round(str(insp_dir))
    assert syn is not None
    assert syn["_round_id"] == "round_03_v3"
    assert syn["_synthetic_inspection_only"] is True
    assert "inspection" not in syn  # so extract_row falls to the sweep sidecars

    # A child round filtered on harmony makes harmony the acted-on architecture.
    acted = sr.variant_from_cluster_key("harmony")
    row = sr.extract_row(syn, None, acted_on_variant=acted)
    assert row["round"] == "round_03_v3"
    assert row["n_clusters"] == 14
    assert row["n_auto_flagged"] == 2
    assert row["clusters_variant"] == "harmony"

    # A directory with no inspect_* sidecars is not a synthetic round.
    (tmp_path / "output").mkdir()
    assert sr.synthesize_inspection_round(str(tmp_path / "output")) is None


def test_discover_inspection_only_rounds(tmp_path):
    # Discovery must surface an inspection-only round alongside ordinary rounds,
    # without mistaking an ordinary round's own inspect_* sidecars (or a _sweep
    # sibling) for standalone rounds.
    _round1(tmp_path)                                   # ordinary round_01
    _sweep_inspection(tmp_path / "round_01", "harmony", 24, 3)   # its own sidecar
    _sweep_inspection(tmp_path / "round_02_sweep", "harmony", 19, 0)  # _sweep sib
    _sweep_inspection(tmp_path / "round_01b", "harmony", 10, 1)  # inspection-only

    real = [sr.load_manifest(str(p)) for p in
            sorted(tmp_path.glob("round_*/round_manifest.json"))]
    assert [m["_round_id"] for m in real] == ["round_01"]
    exclude = {m["_dir"] for m in real} | {m["_dir"] + "_sweep" for m in real}

    syn = sr.discover_inspection_only_rounds(str(tmp_path), exclude)
    # Only round_01b is a standalone inspection-only round; round_01's sidecar
    # and the round_02_sweep sibling must be excluded.
    assert [m["_round_id"] for m in syn] == ["round_01b"]


def test_explicit_rounds_accepts_inspection_only_dir(tmp_path):
    # --rounds should accept an inspection-only round directory (no top-level
    # manifest) the same as any other round path, via synthesis.
    insp_dir = tmp_path / "round_03_v3"
    _sweep_inspection(insp_dir, "harmony", 14, 2)
    syn = sr.synthesize_inspection_round(str(insp_dir))
    assert syn is not None and syn["_round_id"] == "round_03_v3"


def test_splice_inspection_only_round_into_lineage(tmp_path):
    # P (round_03) has BOTH a filtering data-child (round_04_filter_cells_v3,
    # whose fingerprint points back to round_03's own bytes) and an
    # inspection-only child (round_03_v3). The inspection round informed that
    # filtering, so it should be spliced onto the edge:
    #   round_03 → round_03_v3 → round_04_filter_cells_v3 → round_04
    p = str(tmp_path / "round_03")
    insp = str(tmp_path / "round_03_v3")
    filt = str(tmp_path / "round_04_filter_cells_v3")
    integ4 = str(tmp_path / "round_04")
    manifests = [
        {"_dir": p, "_round_id": "round_03"},
        {"_dir": insp, "_round_id": "round_03_v3",
         "_synthetic_inspection_only": True},
        {"_dir": filt, "_round_id": "round_04_filter_cells_v3"},
        {"_dir": integ4, "_round_id": "round_04"},
    ]
    parent_map = {p: None, insp: p, filt: p, integ4: filt}
    order = [p, insp, filt, integ4]

    sr.splice_inspection_only_rounds(manifests, parent_map, order)

    assert parent_map[insp] == p       # inspection still hangs off round_03
    assert parent_map[filt] == insp    # filtering now descends from it
    assert parent_map[integ4] == filt  # downstream integration unchanged


def test_splice_respects_supplied_order(tmp_path):
    # An inspection round listed AFTER the filtering it would attach to must not
    # be spliced onto that earlier edge.
    p = str(tmp_path / "round_03")
    filt = str(tmp_path / "round_04_filter")
    insp = str(tmp_path / "round_03_v3")
    manifests = [
        {"_dir": p, "_round_id": "round_03"},
        {"_dir": filt, "_round_id": "round_04_filter"},
        {"_dir": insp, "_round_id": "round_03_v3",
         "_synthetic_inspection_only": True},
    ]
    parent_map = {p: None, filt: p, insp: p}
    order = [p, filt, insp]  # inspection comes last

    sr.splice_inspection_only_rounds(manifests, parent_map, order)

    assert parent_map[filt] == p  # unchanged: inspection did not precede it

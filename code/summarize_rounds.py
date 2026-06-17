#!/usr/bin/env python3
"""Summarize a completed set of integration rounds.

Harvests every ``round_manifest.json`` under a rounds directory, orders the
rounds by their recorded lineage (the ``parent_round`` pointers, not directory
names), and emits a combined view of the whole experiment:

  - a per-round summary table (CSV + Markdown + HTML),
  - a flattened decisions ledger (every keep/remove action across all rounds),
  - a Mermaid lineage diagram (round-to-round, with seed/commit/cell deltas),
  - an integrity report (``--verify``): lineage hash continuity, seed/code/env
    drift, and ambiguous ("auto") variant selections.

The script only *reads* manifests and the sidecar artifacts each round already
writes (``scib_benchmark_results.csv``); it never recomputes anything, so the
summary is a faithful record rather than a new analysis.

Usage
-----
    python code/summarize_rounds.py --rounds-dir results/rounds
    python code/summarize_rounds.py --rounds-dir results/rounds --verify --strict
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime
from html import escape
from typing import Any, Dict, List, Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover - pandas is a pipeline dependency
    pd = None


# ═══════════════════════════════════════════════════════════════════════════════
#  Discovery & loading
# ═══════════════════════════════════════════════════════════════════════════════

def find_manifests(rounds_dir: str) -> List[str]:
    """Return every round_manifest.json at or one level below rounds_dir."""
    found = []
    direct = os.path.join(rounds_dir, "round_manifest.json")
    if os.path.exists(direct):
        found.append(direct)
    for entry in sorted(os.listdir(rounds_dir)):
        sub = os.path.join(rounds_dir, entry)
        if os.path.isdir(sub):
            candidate = os.path.join(sub, "round_manifest.json")
            if os.path.exists(candidate):
                found.append(candidate)
    return found


def load_manifest(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    # realpath so symlinked roots (e.g. macOS /tmp → /private/tmp) and the
    # abspaths recorded inside the manifest at run time compare equal.
    real = os.path.realpath(path)
    manifest["_path"] = real
    manifest["_dir"] = os.path.dirname(real)
    manifest["_round_id"] = os.path.basename(os.path.dirname(real))
    return manifest


# ═══════════════════════════════════════════════════════════════════════════════
#  Lineage ordering
# ═══════════════════════════════════════════════════════════════════════════════

def _parent_dir_of(manifest: Dict[str, Any]) -> Optional[str]:
    """Best-effort parent round directory for a manifest.

    The cross-round edge is the input the *decisions* (filtering) stage consumed
    — i.e. the previous round's integrated object. Fall back to the integration
    stage's parent pointer when a round has no filtering step (e.g. round 1 with
    an external input has neither and is a root).
    """
    for stage in ("decisions", "integration"):
        block = manifest.get(stage)
        if not isinstance(block, dict):
            continue
        parent = block.get("parent_round")
        if isinstance(parent, dict) and parent.get("output_dir"):
            return os.path.realpath(parent["output_dir"])
        in_path = block.get("input_h5ad")
        if in_path:
            cand = os.path.realpath(os.path.dirname(in_path))
            # A round's own re-integration reads its own filtered.h5ad; that is
            # a self-edge, not a parent. Ignore it.
            if cand != manifest["_dir"]:
                return cand
    return None


def build_parent_map(manifests: List[Dict[str, Any]],
                     explicit_order: Optional[List[str]] = None
                     ) -> Dict[str, Optional[str]]:
    """Map each round's _dir to its parent round's _dir.

    Resolution order:
      1. the recorded ``parent_round`` pointer, matched by realpath;
      2. failing that, by round-id basename (set was moved but names still chain);
      3. failing that — only when rounds were supplied as an explicit ordered
         list (``--rounds``) — the previous item in that list.

    Step 3 is what makes archived round sets work: after CodeOcean overwrites
    ``/results`` and rounds are copied to S3/data under new names, the abspaths
    and round-ids captured at run time no longer match, so the user-supplied
    order is the only reliable lineage signal.
    """
    by_dir = {m["_dir"]: m for m in manifests}
    by_id = {m["_round_id"]: m for m in manifests}
    order = explicit_order or []
    pos = {d: i for i, d in enumerate(order)}
    parent_map: Dict[str, Optional[str]] = {}
    for m in manifests:
        pdir = _parent_dir_of(m)
        resolved: Optional[str] = None
        if pdir is not None:
            if pdir in by_dir:
                resolved = pdir
            else:
                sib = by_id.get(os.path.basename(pdir))
                if sib is not None and sib["_dir"] != m["_dir"]:
                    resolved = sib["_dir"]
        if resolved is None and pos.get(m["_dir"], 0) > 0:
            resolved = order[pos[m["_dir"]] - 1]
        parent_map[m["_dir"]] = resolved
    return parent_map


def order_by_lineage(manifests: List[Dict[str, Any]],
                     parent_map: Dict[str, Optional[str]]
                     ) -> List[Dict[str, Any]]:
    """Topologically order rounds parent-before-child; ties broken by name."""
    by_dir = {m["_dir"]: m for m in manifests}
    depth_cache: Dict[str, int] = {}

    def depth(m: Dict[str, Any], seen=None) -> int:
        seen = seen or set()
        if m["_dir"] in depth_cache:
            return depth_cache[m["_dir"]]
        parent_dir = parent_map.get(m["_dir"])
        if parent_dir is None or parent_dir in seen or parent_dir not in by_dir:
            depth_cache[m["_dir"]] = 0
            return 0
        d = 1 + depth(by_dir[parent_dir], seen | {m["_dir"]})
        depth_cache[m["_dir"]] = d
        return d

    return sorted(manifests, key=lambda m: (depth(m), m["_round_id"]))


# ═══════════════════════════════════════════════════════════════════════════════
#  scIB benchmark harvesting
# ═══════════════════════════════════════════════════════════════════════════════

def variant_from_embedding_key(key: str) -> str:
    """Map an obsm embedding key to its integration-variant label."""
    if key == "X_scVI":
        return "default"
    if key == "X_scANVI":
        return "scanvi"
    if key == "X_pca_harmony":
        return "harmony"
    if key == "X_pca":
        return "baseline (PCA)"
    if key.startswith("X_scVI_"):
        return key[len("X_scVI_"):]
    return key


def read_scib(round_dir: str) -> Optional[pd.DataFrame]:
    """Read scib_benchmark_results.csv → {embedding: total_score}, if present."""
    if pd is None:
        return None
    path = os.path.join(round_dir, "scib_benchmark_results.csv")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, index_col=0)
    except Exception:
        return None
    # scib-metrics appends a "Metric Type" descriptor row; drop it.
    df = df[df.index.astype(str) != "Metric Type"]
    return df


def scib_scores(round_dir: str) -> Dict[str, float]:
    """Return {variant: Total score} for each benchmarked embedding."""
    df = read_scib(round_dir)
    if df is None:
        return {}
    total_col = None
    for cand in ("Total", "total"):
        if cand in df.columns:
            total_col = cand
            break
    if total_col is None:
        return {}
    scores = {}
    for emb, val in df[total_col].items():
        try:
            scores[variant_from_embedding_key(str(emb))] = round(float(val), 4)
        except (TypeError, ValueError):
            continue
    return scores


# ═══════════════════════════════════════════════════════════════════════════════
#  Per-round row extraction
# ═══════════════════════════════════════════════════════════════════════════════

def _g(d: Any, *keys, default=None):
    """Nested .get with dict-safety."""
    cur = d
    for k in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
    return cur if cur is not None else default


def _best_integrated(scores: Dict[str, float]):
    """(variant, score) of the top non-baseline embedding, or (None, None)."""
    integrated = {v: s for v, s in scores.items() if v != "baseline (PCA)"}
    if not integrated:
        return None, None
    best = max(integrated, key=integrated.get)
    return best, integrated[best]


def extract_row(m: Dict[str, Any], parent: Optional[Dict[str, Any]] = None
                ) -> Dict[str, Any]:
    """Build one table row.

    ``parent`` is the manifest of the round this one was filtered from. The
    filtering variant was clustered/scored in that parent's integration, so the
    "selected vs best" scIB comparison is sourced from the PARENT benchmark,
    while "scIB best" reflects THIS round's own integration quality (the trend).
    """
    integ = m.get("integration") or {}
    insp = m.get("inspection") or {}
    dec = m.get("decisions") or {}
    filtered_on = dec.get("filtered_on") or {}

    embedding_keys = integ.get("embedding_keys") or []
    scvi_variants = [variant_from_embedding_key(k) for k in embedding_keys
                     if k.startswith("X_scVI_") or k == "X_scVI"]

    own_best_variant, own_best_score = _best_integrated(scib_scores(m["_dir"]))

    # Selected-variant comparison, against the parent benchmark where the
    # filtering choice was actually made.
    selected_variant = filtered_on.get("variant")
    sel_score = sel_best_variant = sel_best_score = sel_is_best = None
    if selected_variant and parent is not None:
        parent_scores = scib_scores(parent["_dir"])
        sel_score = parent_scores.get(selected_variant)
        sel_best_variant, sel_best_score = _best_integrated(parent_scores)
        if sel_best_variant is not None:
            sel_is_best = (selected_variant == sel_best_variant)

    # Compact "selected / parent-best" display; ✗ flags a non-top choice.
    if sel_score is not None and sel_best_score is not None:
        sel_vs_best = f"{sel_score} / {sel_best_score}" + ("" if sel_is_best else " ✗")
    elif sel_score is not None:
        sel_vs_best = str(sel_score)
    else:
        sel_vs_best = None

    return {
        "round": m["_round_id"],
        "updated": m.get("updated"),
        "n_cells": integ.get("n_cells"),
        "n_genes": integ.get("n_genes"),
        "model_type": integ.get("model_type"),
        "annotation": integ.get("annotation_enabled"),
        "n_variants": len(scvi_variants) or (1 if embedding_keys else None),
        "variants": ", ".join(scvi_variants) if scvi_variants else None,
        "seed": _g(m, "reproducibility", "seed"),
        "commit": _g(m, "code", "commit_short"),
        "dirty": _g(m, "code", "dirty"),
        "n_clusters": insp.get("n_clusters"),
        "n_auto_flagged": insp.get("n_auto_flagged"),
        "input_cells": dec.get("input_cells"),
        "output_cells": dec.get("output_cells"),
        "removed": dec.get("removed"),
        "removal_pct": dec.get("removal_pct"),
        "filtered_on": selected_variant,
        "filter_source": filtered_on.get("source"),
        "filter_cluster_key": filtered_on.get("cluster_key") or dec.get("cluster_key"),
        "scib_best": own_best_score,
        "scib_best_variant": own_best_variant,
        "scib_selected": sel_score,
        "scib_selected_best": sel_best_score,
        "scib_selected_best_variant": sel_best_variant,
        "scib_selected_is_best": sel_is_best,
        "scib_sel_vs_best": sel_vs_best,
    }


def extract_decisions_ledger(manifests: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for m in manifests:
        dec = m.get("decisions") or {}
        actions = dec.get("actions_summary") or []
        filtered_on = dec.get("filtered_on") or {}
        for action in actions:
            rows.append({
                "round": m["_round_id"],
                "filtered_on": filtered_on.get("variant"),
                "cluster_key": filtered_on.get("cluster_key") or dec.get("cluster_key"),
                "action": action,
            })
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
#  Verification
# ═══════════════════════════════════════════════════════════════════════════════

def _sha256(path: str, max_bytes: Optional[int] = None) -> Optional[str]:
    if not os.path.exists(path):
        return None
    digest = hashlib.sha256()
    read = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            if max_bytes is not None and read + len(chunk) > max_bytes:
                digest.update(chunk[: max_bytes - read])
                break
            digest.update(chunk)
            read += len(chunk)
    return digest.hexdigest()


def _parent_integrated_file(parent: Dict[str, Any]) -> Optional[str]:
    """The parent round's integrated.h5ad as it sits on disk now.

    Uses the recorded output filename but resolves it inside the parent's
    *current* directory, so the check survives the round set being archived to a
    new location (the recorded abspath would be stale).
    """
    out = _g(parent, "integration", "output_h5ad")
    base = os.path.basename(out) if out else "integrated.h5ad"
    for cand in (os.path.join(parent["_dir"], base),
                 os.path.join(parent["_dir"], "integrated.h5ad")):
        if os.path.exists(cand):
            return cand
    return None


def verify(manifests: List[Dict[str, Any]],
           parent_map: Dict[str, Optional[str]],
           scan_dir: Optional[str] = None) -> List[Dict[str, Any]]:
    """Return a list of findings: {round, severity, check, message}."""
    findings: List[Dict[str, Any]] = []
    by_dir = {m["_dir"]: m for m in manifests}

    def add(rnd, severity, check, message):
        findings.append({"round": rnd, "severity": severity,
                         "check": check, "message": message})

    # Round dirs that exist but carry no manifest (discovery mode only).
    if scan_dir:
        manifest_dirs = {m["_dir"] for m in manifests}
        for entry in sorted(os.listdir(scan_dir)):
            sub = os.path.realpath(os.path.join(scan_dir, entry))
            if os.path.isdir(sub) and sub not in manifest_dirs:
                if any(f.endswith(".h5ad") for f in os.listdir(sub)):
                    add(entry, "warn", "missing_manifest",
                        "Directory has h5ad output but no round_manifest.json.")

    seeds, commits = set(), set()
    pkg_baseline: Optional[Dict[str, Any]] = None

    for m in manifests:
        rid = m["_round_id"]
        seed = _g(m, "reproducibility", "seed")
        if seed is not None:
            seeds.add(seed)
        commit = _g(m, "code", "commit")
        if commit:
            commits.add(commit)
        if _g(m, "code", "dirty") is True:
            add(rid, "warn", "dirty_tree",
                "Round was produced from a dirty working tree (uncommitted changes).")

        # Ambiguous variant selection.
        dec = m.get("decisions") or {}
        filtered_on = dec.get("filtered_on") or {}
        if filtered_on.get("source") == "auto":
            add(rid, "warn", "auto_variant",
                "Filtering variant was resolved by 'auto' — selection is implicit, "
                "not declared. Pin it in decisions.yaml or via --cluster-key.")
        if dec and not filtered_on:
            add(rid, "info", "no_filtered_on",
                "Decisions stage predates filtered_on capture; variant is only "
                "implied by cluster_key.")

        # Environment drift across rounds.
        pkgs = _g(m, "integration", "package_versions")
        if isinstance(pkgs, dict):
            if pkg_baseline is None:
                pkg_baseline = pkgs
            else:
                changed = [f"{k}: {pkg_baseline.get(k)}→{v}"
                           for k, v in pkgs.items()
                           if k in pkg_baseline and pkg_baseline[k] != v]
                if changed:
                    add(rid, "warn", "env_drift",
                        "Package versions differ from the first round: "
                        + "; ".join(changed))

        # Lineage hash continuity: the content the decision consumed should match
        # the parent round's integrated output as it sits on disk now. Hashing
        # the parent's current file (rather than the stale recorded path) makes
        # the check robust to the round set being archived/relocated.
        fp = _g(m, "decisions", "input_fingerprint")
        if isinstance(fp, dict) and fp.get("sha256"):
            recorded = fp["sha256"]
            max_bytes = fp.get("hashed_bytes") if fp.get("partial") else None
            parent = by_dir.get(parent_map.get(m["_dir"]) or "")
            target = _parent_integrated_file(parent) if parent else None
            if target is None and fp.get("path") and os.path.exists(fp["path"]):
                target = fp["path"]  # fall back to recorded path if still present
            if target is None:
                add(rid, "info", "input_unverifiable",
                    "Parent integrated file not found on disk; cannot verify "
                    "lineage hash continuity (set may be archived).")
            else:
                current = _sha256(target, max_bytes=max_bytes)
                if current != recorded:
                    add(rid, "error", "input_changed",
                        f"Decision input content does not match the parent's "
                        f"integrated file ({target}): recorded sha256 differs.")

    if len(seeds) > 1:
        add("(global)", "warn", "seed_drift",
            f"Multiple seeds across rounds: {sorted(seeds)} — embeddings/cluster "
            "IDs are not comparable round-to-round.")
    if len(commits) > 1:
        add("(global)", "warn", "code_drift",
            f"Rounds were produced at {len(commits)} different commits — config "
            "alone does not pin the code that ran.")

    return findings


# ═══════════════════════════════════════════════════════════════════════════════
#  Rendering
# ═══════════════════════════════════════════════════════════════════════════════

TABLE_COLUMNS = [
    ("round", "Round"),
    ("n_cells", "Cells"),
    ("n_genes", "Genes"),
    ("model_type", "Model"),
    ("n_variants", "Variants"),
    ("seed", "Seed"),
    ("commit", "Commit"),
    ("n_clusters", "Clusters"),
    ("n_auto_flagged", "Flagged"),
    ("removed", "Removed"),
    ("removal_pct", "Rem%"),
    ("filtered_on", "Filtered on"),
    ("filter_source", "Sel. via"),
    ("scib_best", "scIB best"),
    ("scib_best_variant", "best var"),
    ("scib_sel_vs_best", "Sel/best (parent)"),
]


def _fmt(val: Any) -> str:
    if val is None:
        return "—"
    if isinstance(val, bool):
        return "yes" if val else "no"
    if isinstance(val, int):
        return f"{val:,}"
    return str(val)


def render_markdown_table(rows: List[Dict[str, Any]]) -> str:
    headers = [label for _, label in TABLE_COLUMNS]
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        cells = []
        for key, _ in TABLE_COLUMNS:
            v = _fmt(r.get(key))
            if key == "commit" and r.get("dirty") is True and r.get("commit"):
                v += " ⚠"
            if key == "filter_source" and r.get(key) == "auto":
                v += " ⚠"
            cells.append(v)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def render_mermaid(manifests: List[Dict[str, Any]],
                   rows_by_dir: Dict[str, Dict[str, Any]],
                   parent_map: Dict[str, Optional[str]]) -> str:
    lines = ["flowchart TD"]
    node_id = {m["_dir"]: f"r{i}" for i, m in enumerate(manifests)}
    for m in manifests:
        row = rows_by_dir[m["_dir"]]
        cls = "clean"
        if _g(m, "code", "dirty") is True or row.get("filter_source") == "auto":
            cls = "drift"
        label_bits = [
            f"<b>{row['round']}</b>",
            f"{_fmt(row.get('n_cells'))} cells · {_fmt(row.get('n_clusters'))} clusters",
            f"seed {_fmt(row.get('seed'))} · {_fmt(row.get('commit'))}",
        ]
        if row.get("scib_best") is not None:
            label_bits.append(f"scIB {row['scib_best']} ({row.get('scib_best_variant')})")
        lines.append(f'    {node_id[m["_dir"]]}["{"<br/>".join(label_bits)}"]:::{cls}')
    for m in manifests:
        parent_dir = parent_map.get(m["_dir"])
        if parent_dir in node_id:
            row = rows_by_dir[m["_dir"]]
            edge_bits = []
            if row.get("removed") is not None:
                edge_bits.append(f"−{_fmt(row.get('removed'))} ({_fmt(row.get('removal_pct'))}%)")
            if row.get("filtered_on"):
                edge_bits.append(f"on {row['filtered_on']}")
            label = "<br/>".join(edge_bits) or "next round"
            lines.append(f'    {node_id[parent_dir]} -->|"{label}"| {node_id[m["_dir"]]}')
    lines += [
        "    classDef clean fill:#eef7ee,stroke:#449944;",
        "    classDef drift fill:#fff3e0,stroke:#cc7a00;",
    ]
    return "\n".join(lines)


def render_findings_md(findings: List[Dict[str, Any]]) -> str:
    if not findings:
        return "_No integrity issues detected._"
    order = {"error": 0, "warn": 1, "info": 2}
    icon = {"error": "❌", "warn": "⚠️", "info": "ℹ️"}
    rows = ["| Severity | Round | Check | Detail |", "|---|---|---|---|"]
    for f in sorted(findings, key=lambda f: order.get(f["severity"], 9)):
        rows.append(f"| {icon.get(f['severity'], '')} {f['severity']} | "
                    f"{f['round']} | {f['check']} | {f['message']} |")
    return "\n".join(rows)


def render_markdown(rows, manifests, rows_by_dir, parent_map, findings,
                    verified) -> str:
    parts = [
        "# Pipeline summary",
        "",
        f"_Generated {datetime.now().isoformat(timespec='seconds')} · "
        f"{len(rows)} round(s)_",
        "",
        "## Rounds",
        "",
        render_markdown_table(rows),
        "",
        "## Lineage",
        "",
        "```mermaid",
        render_mermaid(manifests, rows_by_dir, parent_map),
        "```",
        "",
        "## Integrity" + ("" if verified else " (run with --verify)"),
        "",
        render_findings_md(findings) if verified else "_Not checked._",
        "",
    ]
    return "\n".join(parts)


def render_html(md_table_rows, mermaid_src, findings, verified) -> str:
    head = (
        "<!doctype html><meta charset='utf-8'>"
        "<title>Pipeline summary</title>"
        "<style>body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;"
        "margin:2rem;color:#222}table{border-collapse:collapse;margin:1rem 0}"
        "th,td{border:1px solid #ccc;padding:4px 8px;font-size:13px;text-align:right}"
        "th{background:#f3f6fb}td:first-child,th:first-child{text-align:left}"
        "h1,h2{color:#234}.warn{color:#cc7a00}.error{color:#c0392b}</style>"
        "<script src='https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js'>"
        "</script><script>mermaid.initialize({startOnLoad:true});</script>"
    )
    # Table
    thead = "".join(f"<th>{escape(lbl)}</th>" for _, lbl in TABLE_COLUMNS)
    body_rows = []
    for r in md_table_rows:
        tds = []
        for key, _ in TABLE_COLUMNS:
            v = _fmt(r.get(key))
            cls = ""
            if key == "commit" and r.get("dirty") is True:
                v += " ⚠"; cls = " class='warn'"
            if key == "filter_source" and r.get(key) == "auto":
                v += " ⚠"; cls = " class='warn'"
            tds.append(f"<td{cls}>{escape(v)}</td>")
        body_rows.append("<tr>" + "".join(tds) + "</tr>")
    table = f"<table><tr>{thead}</tr>{''.join(body_rows)}</table>"

    # Findings
    if not verified:
        findings_html = "<p><em>Not checked (run with --verify).</em></p>"
    elif not findings:
        findings_html = "<p><em>No integrity issues detected.</em></p>"
    else:
        rows = ["<tr><th>Severity</th><th>Round</th><th>Check</th><th>Detail</th></tr>"]
        for f in findings:
            rows.append(
                f"<tr><td class='{escape(f['severity'])}'>{escape(f['severity'])}</td>"
                f"<td>{escape(str(f['round']))}</td><td>{escape(f['check'])}</td>"
                f"<td>{escape(f['message'])}</td></tr>")
        findings_html = "<table>" + "".join(rows) + "</table>"

    return (f"{head}<h1>Pipeline summary</h1>"
            f"<p><em>Generated {escape(datetime.now().isoformat(timespec='seconds'))}</em></p>"
            f"<h2>Rounds</h2>{table}"
            f"<h2>Lineage</h2><pre class='mermaid'>{escape(mermaid_src)}</pre>"
            f"<h2>Integrity</h2>{findings_html}")


# ═══════════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Summarize a completed set of integration rounds.")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--rounds-dir",
                     help="Directory containing round_* subdirectories "
                          "(auto-discovers and lineage-orders them).")
    src.add_argument("--rounds", nargs="+", metavar="PATH",
                     help="Explicit, ordered list of round directories (or "
                          "round_manifest.json paths). Use this when rounds are "
                          "archived under arbitrary names/locations (e.g. an S3 "
                          "mount). The supplied order is taken as the lineage "
                          "when recorded parent pointers can't be resolved.")
    parser.add_argument("--output-dir", default=None,
                        help="Where to write outputs (default: <rounds-dir>/summary, "
                             "or ./round_summary with --rounds).")
    parser.add_argument("--verify", action="store_true",
                        help="Run integrity checks (lineage, seed/code/env drift).")
    parser.add_argument("--strict", action="store_true",
                        help="Exit non-zero if any verification finding is raised.")
    parser.add_argument("--format", default="md,html,csv,json,mmd",
                        help="Comma-separated subset of: md,html,csv,json,mmd.")
    args = parser.parse_args()
    formats = {f.strip() for f in args.format.split(",") if f.strip()}

    explicit_order = None
    scan_dir = None
    if args.rounds_dir:
        if not os.path.isdir(args.rounds_dir):
            parser.error(f"rounds-dir not found: {args.rounds_dir}")
        manifest_paths = find_manifests(args.rounds_dir)
        if not manifest_paths:
            parser.error(f"No round_manifest.json found under {args.rounds_dir}")
        manifests = [load_manifest(p) for p in manifest_paths]
        scan_dir = args.rounds_dir
        default_out = os.path.join(args.rounds_dir, "summary")
    else:
        manifest_paths = []
        for item in args.rounds:
            mp = os.path.join(item, "round_manifest.json") if os.path.isdir(item) else item
            if not os.path.exists(mp):
                parser.error(f"No round_manifest.json at: {item}")
            manifest_paths.append(mp)
        manifests = [load_manifest(p) for p in manifest_paths]
        explicit_order = [m["_dir"] for m in manifests]  # trust supplied order
        default_out = "round_summary"

    output_dir = args.output_dir or default_out

    parent_map = build_parent_map(manifests, explicit_order)
    # Discovery mode: topologically reorder. Explicit list: keep supplied order.
    if explicit_order is None:
        manifests = order_by_lineage(manifests, parent_map)
    print(f"  Found {len(manifests)} round(s): "
          f"{', '.join(m['_round_id'] for m in manifests)}")

    by_dir = {m["_dir"]: m for m in manifests}
    rows = []
    for m in manifests:
        pdir = parent_map.get(m["_dir"])
        rows.append(extract_row(m, by_dir.get(pdir) if pdir else None))
    rows_by_dir = {m["_dir"]: r for m, r in zip(manifests, rows)}
    ledger = extract_decisions_ledger(manifests)
    findings = verify(manifests, parent_map, scan_dir) if args.verify else []

    os.makedirs(output_dir, exist_ok=True)
    mermaid_src = render_mermaid(manifests, rows_by_dir, parent_map)

    written = []
    if "csv" in formats and pd is not None:
        p = os.path.join(output_dir, "rounds_table.csv")
        pd.DataFrame(rows).to_csv(p, index=False)
        written.append(p)
        p = os.path.join(output_dir, "decisions_ledger.csv")
        pd.DataFrame(ledger).to_csv(p, index=False)
        written.append(p)
    if "mmd" in formats:
        p = os.path.join(output_dir, "lineage.mmd")
        with open(p, "w", encoding="utf-8") as h:
            h.write(mermaid_src + "\n")
        written.append(p)
    if "md" in formats:
        p = os.path.join(output_dir, "pipeline_summary.md")
        with open(p, "w", encoding="utf-8") as h:
            h.write(render_markdown(rows, manifests, rows_by_dir, parent_map,
                                    findings, args.verify))
        written.append(p)
    if "html" in formats:
        p = os.path.join(output_dir, "pipeline_summary.html")
        with open(p, "w", encoding="utf-8") as h:
            h.write(render_html(rows, mermaid_src, findings, args.verify))
        written.append(p)
    if "json" in formats:
        p = os.path.join(output_dir, "pipeline_summary.json")
        payload = {
            "generated": datetime.now().isoformat(),
            "rounds_dir": os.path.abspath(args.rounds_dir) if args.rounds_dir else None,
            "round_dirs": [m["_dir"] for m in manifests],
            "rounds": rows,
            "decisions_ledger": ledger,
            "findings": findings,
            "verified": args.verify,
        }
        with open(p, "w", encoding="utf-8") as h:
            json.dump(payload, h, indent=2, default=str)
            h.write("\n")
        written.append(p)

    print("\n" + render_markdown_table(rows))
    if args.verify:
        n_err = sum(1 for f in findings if f["severity"] == "error")
        n_warn = sum(1 for f in findings if f["severity"] == "warn")
        print(f"\n  Integrity: {n_err} error(s), {n_warn} warning(s)")
        for f in findings:
            print(f"    [{f['severity']}] {f['round']}: {f['check']} — {f['message']}")
    print("\n  Wrote:")
    for p in written:
        print(f"    {p}")

    if args.strict and any(f["severity"] == "error" for f in findings):
        sys.exit(1)


if __name__ == "__main__":
    main()

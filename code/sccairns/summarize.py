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

Clusters/flagged counts come from the per-architecture sweep inspection
manifests (``inspect_<arch>/round_manifest.json``); each round shows the
architecture the human acted on (the variant its child round was filtered on).
The terminal round has no child, so pin its architecture explicitly:

    python code/summarize_rounds.py --rounds r1 r2 r3 \\
        --terminal-run-architecture leiden_harmony

A round that only *re-inspected* a prior round's output across a sweep has no
top-level ``round_manifest.json`` — just ``inspect_<arch>/`` sidecars. Such a
directory is still recognized as a round (labeled by its directory name), with
clusters/flagged sourced from the acted-on architecture; point ``--rounds`` or
``--rounds-dir`` at it exactly as for any other round.
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


def splice_inspection_only_rounds(
        manifests: List[Dict[str, Any]],
        parent_map: Dict[str, Optional[str]],
        order: Optional[List[str]] = None) -> Dict[str, Optional[str]]:
    """Insert inspection-only rounds onto the lineage edge they informed.

    An inspection-only round re-inspects its parent P's integrated output but
    writes no new h5ad, so nothing's input fingerprint points to it and it lands
    as a leaf dangling off P. Its purpose, though, is to update the
    flags/decisions that P's downstream filtering acts on. Re-point P's
    data-children (rounds whose resolved parent is P) through the inspection
    round so it sits inline — ``P → inspect → filtered-child`` — instead of
    dangling.

    Only fires when P has exactly one inspection-only child (otherwise the
    target edge is ambiguous). When an explicit ``order`` is supplied, a
    data-child is re-parented only if the inspection round precedes it in that
    order, so a later inspection is never spliced onto an earlier filtering.
    """
    pos = {d: i for i, d in enumerate(order)} if order else {}
    insp_children: Dict[str, List[str]] = {}
    for m in manifests:
        if m.get("_synthetic_inspection_only"):
            p = parent_map.get(m["_dir"])
            if p:
                insp_children.setdefault(p, []).append(m["_dir"])
    for p, kids in insp_children.items():
        if len(kids) != 1:
            continue
        insp_dir = kids[0]
        for m in manifests:
            d = m["_dir"]
            if d == insp_dir or m.get("_synthetic_inspection_only"):
                continue
            if parent_map.get(d) != p:
                continue
            if pos and pos.get(insp_dir, 0) >= pos.get(d, 0):
                continue
            parent_map[d] = insp_dir
    return parent_map


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


def variant_from_cluster_key(key: Optional[str]) -> Optional[str]:
    """Map a leiden cluster key to its integration-variant label."""
    if not key:
        return None
    if key == "leiden":
        return "default"
    if key.startswith("leiden_"):
        return key[len("leiden_"):]
    return key


def find_inspection_manifests(round_dir: str) -> Dict[str, Dict[str, Any]]:
    """Locate a round's per-architecture inspection manifests.

    In sweep mode ``inspect_integration.py`` writes one ``round_manifest.json``
    per architecture into ``inspect_<arch>/`` subdirectories (the only manifests
    that carry the ``inspection`` stage — n_clusters/n_auto_flagged). Depending
    on how the run was invoked these sit either inside the round dir itself or in
    a sibling ``<round>_sweep`` directory. Return ``{variant: manifest}`` keyed
    by the integration variant the inspection clustered on.
    """
    result: Dict[str, Dict[str, Any]] = {}
    for root in (round_dir, round_dir + "_sweep"):
        if not os.path.isdir(root):
            continue
        for entry in sorted(os.listdir(root)):
            if not entry.startswith("inspect_"):
                continue
            sub = os.path.join(root, entry)
            cand = os.path.join(sub, "round_manifest.json")
            if not (os.path.isdir(sub) and os.path.exists(cand)):
                continue
            try:
                man = load_manifest(cand)
            except Exception:
                continue
            insp = man.get("inspection")
            if not isinstance(insp, dict):
                continue
            variant = (variant_from_cluster_key(insp.get("cluster_key"))
                       or entry[len("inspect_"):])
            result.setdefault(variant, man)
    return result


def synthesize_inspection_round(round_dir: str) -> Optional[Dict[str, Any]]:
    """Build a round record for an inspection-only sweep directory.

    Some rounds only *re-inspect* a prior round's integrated output across a
    sweep of architectures — they produce no integration of their own, so the
    directory has no top-level ``round_manifest.json``, just
    ``inspect_<arch>/round_manifest.json`` sidecars. Treat such a directory as a
    round in its own right: label it by the directory name and leave the
    ``inspection`` stage absent so :func:`extract_row` sources clusters/flagged
    from the acted-on architecture's sidecar (the normal sweep path). Run-level
    provenance (seed, commit) is donated from a sidecar, since every stage of a
    run records the same values.

    Returns ``None`` when the directory has no inspection sidecars, so callers
    can fall back to erroring as before.
    """
    insp_by_variant = find_inspection_manifests(round_dir)
    if not insp_by_variant:
        return None
    donor = next(iter(insp_by_variant.values()))
    real = os.path.realpath(round_dir)
    manifest: Dict[str, Any] = {
        "pipeline_version": donor.get("pipeline_version"),
        "reproducibility": donor.get("reproducibility"),
        "code": donor.get("code"),
        "updated": donor.get("updated"),
        "_synthetic_inspection_only": True,
    }
    # Best-effort parent link: if a sidecar recorded the object it inspected,
    # expose it as the decision input so discovery-mode lineage can chain this
    # round to its parent. (In --rounds mode the supplied order is authoritative,
    # so this is only a convenience for --rounds-dir.)
    for man in insp_by_variant.values():
        src = _g(man, "inspection", "input_h5ad") or _g(man, "integration",
                                                         "input_h5ad")
        if src:
            manifest["decisions"] = {"input_h5ad": src}
            break
    manifest["_path"] = os.path.join(real, "round_manifest.json")
    manifest["_dir"] = real
    manifest["_round_id"] = os.path.basename(real)
    return manifest


def discover_inspection_only_rounds(rounds_dir: str,
                                    exclude_dirs: set) -> List[Dict[str, Any]]:
    """Discover inspection-only sweep rounds under rounds_dir.

    Looks at rounds_dir and each immediate subdirectory for a dir that has
    ``inspect_*`` sidecars but no top-level ``round_manifest.json`` (see
    :func:`synthesize_inspection_round`). ``exclude_dirs`` holds the realpaths of
    directories already accounted for as ordinary rounds (and their ``_sweep``
    siblings), so their sidecars are not mistaken for standalone rounds. Dirs
    whose names start with ``inspect_`` or end with ``_sweep`` are skipped for
    the same reason.
    """
    found: List[Dict[str, Any]] = []
    seen = set(exclude_dirs)
    candidates = [rounds_dir]
    for entry in sorted(os.listdir(rounds_dir)):
        sub = os.path.join(rounds_dir, entry)
        if os.path.isdir(sub):
            candidates.append(sub)
    for cand in candidates:
        real = os.path.realpath(cand)
        if real in seen:
            continue
        base = os.path.basename(real)
        if base.startswith("inspect_") or base.endswith("_sweep"):
            continue
        syn = synthesize_inspection_round(cand)
        if syn is not None:
            found.append(syn)
            seen.add(real)
    return found


def _scib_csv_path(round_dir: str) -> Optional[str]:
    """First existing scib_benchmark_results.csv for a round.

    Integration writes the CSV into the round dir, but relocated/archived sets
    may carry it under a sibling ``<round>_sweep`` or an ``inspect_*`` subdir;
    check those as fallbacks.
    """
    candidates = [os.path.join(round_dir, "scib_benchmark_results.csv")]
    for root in (round_dir, round_dir + "_sweep"):
        if not os.path.isdir(root):
            continue
        candidates.append(os.path.join(root, "scib_benchmark_results.csv"))
        for entry in sorted(os.listdir(root)):
            if entry.startswith("inspect_"):
                candidates.append(
                    os.path.join(root, entry, "scib_benchmark_results.csv"))
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def read_scib(round_dir: str) -> Optional[pd.DataFrame]:
    """Read scib_benchmark_results.csv → {embedding: total_score}, if present."""
    if pd is None:
        return None
    path = _scib_csv_path(round_dir)
    if path is None:
        return None
    try:
        df = pd.read_csv(path, index_col=0)
    except Exception:
        return None
    # scib-metrics appends a "Metric Type" descriptor row; drop it.
    df = df[df.index.astype(str) != "Metric Type"]
    return df


# scib-metrics aggregate columns, best-to-worst. "Total" only exists when both
# metric families ran; annotation-free runs (no label_key) compute batch metrics
# only, so their best available aggregate is "Batch correction".
_SCIB_AGGREGATE_COLUMNS = ("Total", "total", "Batch correction", "Bio conservation")


def scib_scores(round_dir: str) -> "tuple[Dict[str, float], Optional[str]]":
    """Return ({variant: aggregate score}, aggregate_column_name).

    Prefers the overall ``Total`` score; falls back to whichever aggregate the
    benchmark actually produced (annotation-free runs only have ``Batch
    correction``). Returns ``({}, None)`` when no aggregate column is present.
    """
    df = read_scib(round_dir)
    if df is None:
        return {}, None
    agg_col = next((c for c in _SCIB_AGGREGATE_COLUMNS if c in df.columns), None)
    if agg_col is None:
        return {}, None
    scores = {}
    for emb, val in df[agg_col].items():
        try:
            scores[variant_from_embedding_key(str(emb))] = round(float(val), 4)
        except (TypeError, ValueError):
            continue
    return scores, agg_col


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


def extract_row(m: Dict[str, Any], parent: Optional[Dict[str, Any]] = None,
                acted_on_variant: Optional[str] = None) -> Dict[str, Any]:
    """Build one table row.

    ``parent`` is the manifest of the round this one was filtered from. The
    filtering variant was clustered/scored in that parent's integration, so the
    "selected vs best" scIB comparison is sourced from the PARENT benchmark,
    while "scIB best" reflects THIS round's own integration quality (the trend).

    ``acted_on_variant`` is the variant the CHILD round was filtered on — i.e.
    the architecture whose clustering the human reviewed and acted upon in THIS
    round. It drives which per-architecture inspection manifest supplies the
    clusters/flagged counts (see below).
    """
    integ = m.get("integration") or {}
    dec = m.get("decisions") or {}
    filtered_on = dec.get("filtered_on") or {}
    selected_variant = filtered_on.get("variant")

    embedding_keys = integ.get("embedding_keys") or []
    scvi_variants = [variant_from_embedding_key(k) for k in embedding_keys
                     if k.startswith("X_scVI_") or k == "X_scVI"]

    own_scores, own_scib_metric = scib_scores(m["_dir"])
    own_best_variant, own_best_score = _best_integrated(own_scores)

    # Inspection (clusters/flagged) is written per-architecture into
    # inspect_<arch>/ sweep manifests, not the round's own manifest. Choose the
    # architecture the human acted on: the variant the CHILD round was filtered
    # on, then this round's own filtered variant, then its best-scoring variant,
    # then the sole architecture if there is only one.
    insp = m.get("inspection") if isinstance(m.get("inspection"), dict) else None
    clusters_variant = None
    if insp is None:
        insp_by_variant = find_inspection_manifests(m["_dir"])
        for cand in (acted_on_variant, selected_variant, own_best_variant):
            if cand and cand in insp_by_variant:
                insp = insp_by_variant[cand].get("inspection")
                clusters_variant = cand
                break
        if insp is None and len(insp_by_variant) == 1:
            clusters_variant, only_man = next(iter(insp_by_variant.items()))
            insp = only_man.get("inspection")
    insp = insp or {}

    # Selected-variant comparison, against the parent benchmark where the
    # filtering choice was actually made.
    sel_score = sel_best_variant = sel_best_score = sel_is_best = None
    if selected_variant and parent is not None:
        parent_scores, _ = scib_scores(parent["_dir"])
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
        "clusters_variant": clusters_variant,
        "input_cells": dec.get("input_cells"),
        "output_cells": dec.get("output_cells"),
        "removed": dec.get("removed"),
        "removal_pct": dec.get("removal_pct"),
        "filtered_on": selected_variant,
        "filter_source": filtered_on.get("source"),
        "filter_cluster_key": filtered_on.get("cluster_key") or dec.get("cluster_key"),
        "scib_best": own_best_score,
        "scib_best_variant": own_best_variant,
        "scib_metric": own_scib_metric,
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


def _fmt_pct(val: Any) -> str:
    if val is None:
        return "—"
    try:
        return f"{float(val):.1f}%"
    except (TypeError, ValueError):
        return _fmt(val)


def _fmt_date(val: Any) -> str:
    if not val:
        return "—"
    text = str(val)
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00")).strftime(
            "%Y-%m-%d %H:%M")
    except ValueError:
        return text


def _as_int(val: Any) -> Optional[int]:
    if val is None:
        return None
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return None


def _attr(val: Any) -> str:
    return escape(str(val), quote=True)


def _scib_metric_tag(metric: Optional[str]) -> str:
    """Suffix flagging a non-Total scIB aggregate (e.g. batch-only runs).

    The "scIB best" score is a full Total only when both metric families ran;
    annotation-free runs report "Batch correction" instead. Tag those so the
    number isn't mistaken for a Total.
    """
    if not metric or metric.lower() == "total":
        return ""
    short = {"Batch correction": "batch", "Bio conservation": "bio"}
    return f" ({short.get(metric, metric)})"


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
            if key == "scib_best" and r.get("scib_best") is not None:
                v += _scib_metric_tag(r.get("scib_metric"))
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
            label_bits.append(f"scIB {row['scib_best']}"
                              f"{_scib_metric_tag(row.get('scib_metric'))} "
                              f"({row.get('scib_best_variant')})")
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
            default = ("re-inspect" if m.get("_synthetic_inspection_only")
                       else "next round")
            label = "<br/>".join(edge_bits) or default
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


def _html_filter_text(*values: Any) -> str:
    return _attr(" ".join(str(v) for v in values if v is not None).lower())


def _badge(label: Any, tone: str = "neutral") -> str:
    return f"<span class='badge {tone}'>{escape(_fmt(label))}</span>"


def _stat_card(label: str, value: str, detail: str = "",
               tone: str = "neutral") -> str:
    detail_html = f"<p>{escape(detail)}</p>" if detail else ""
    return (
        f"<article class='stat-card {tone}'>"
        f"<span>{escape(label)}</span>"
        f"<strong>{escape(value)}</strong>"
        f"{detail_html}</article>"
    )


def _report_metrics(rows: List[Dict[str, Any]], findings: List[Dict[str, Any]],
                    verified: bool) -> Dict[str, Any]:
    start_cells = next((_as_int(r.get("n_cells")) for r in rows
                        if _as_int(r.get("n_cells")) is not None), None)
    final_cells = next((_as_int(r.get("n_cells")) for r in reversed(rows)
                        if _as_int(r.get("n_cells")) is not None), None)
    removed = sum(_as_int(r.get("removed")) or 0 for r in rows)
    flagged = sum(_as_int(r.get("n_auto_flagged")) or 0 for r in rows)
    max_clusters = max((_as_int(r.get("n_clusters")) or 0 for r in rows),
                       default=0)
    variant_counts = sorted({_as_int(r.get("n_variants")) for r in rows
                             if _as_int(r.get("n_variants")) is not None})
    retention_pct = None
    if start_cells and final_cells is not None:
        retention_pct = 100 * final_cells / start_cells
    severity_counts = {"error": 0, "warn": 0, "info": 0}
    for f in findings:
        sev = str(f.get("severity", "info"))
        severity_counts[sev] = severity_counts.get(sev, 0) + 1
    if not verified:
        integrity_label = "Not checked"
        integrity_detail = "Run with --verify to add integrity findings"
        integrity_tone = "pending"
    elif not findings:
        integrity_label = "Clean"
        integrity_detail = "No integrity issues detected"
        integrity_tone = "good"
    else:
        integrity_label = (
            f"{severity_counts.get('error', 0)} error / "
            f"{severity_counts.get('warn', 0)} warn"
        )
        integrity_detail = f"{len(findings)} total finding(s)"
        integrity_tone = "risk" if severity_counts.get("error") else "warn"
    return {
        "start_cells": start_cells,
        "final_cells": final_cells,
        "removed": removed,
        "flagged": flagged,
        "max_clusters": max_clusters,
        "variant_counts": variant_counts,
        "retention_pct": retention_pct,
        "integrity_label": integrity_label,
        "integrity_detail": integrity_detail,
        "integrity_tone": integrity_tone,
        "severity_counts": severity_counts,
    }


def _render_overview_cards(rows: List[Dict[str, Any]],
                           findings: List[Dict[str, Any]],
                           verified: bool) -> str:
    metrics = _report_metrics(rows, findings, verified)
    start = _fmt(metrics["start_cells"])
    final = _fmt(metrics["final_cells"])
    retention = _fmt_pct(metrics["retention_pct"])
    variants = metrics["variant_counts"]
    variant_value = "—"
    variant_detail = "No embedding variant count recorded"
    if variants:
        variant_value = _fmt(variants[-1])
        variant_detail = (
            "consistent across rounds" if len(variants) == 1
            else f"range {_fmt(variants[0])}-{_fmt(variants[-1])}")
    cards = [
        _stat_card("Rounds", _fmt(len(rows)),
                   f"{start} start cells to {final} final cells", "neutral"),
        _stat_card("Cells retained", retention,
                   f"{_fmt(metrics['removed'])} total removed", "good"),
        _stat_card("Clusters / flagged",
                   f"{_fmt(metrics['max_clusters'])} / {_fmt(metrics['flagged'])}",
                   "maximum clusters and total auto-flagged clusters", "warn"),
        _stat_card("Variants", variant_value, variant_detail, "neutral"),
        _stat_card("Integrity", metrics["integrity_label"],
                   metrics["integrity_detail"], metrics["integrity_tone"]),
    ]
    return "<section class='stat-grid' aria-label='Overview'>" + "".join(cards) + "</section>"


def _render_rounds_table(rows: List[Dict[str, Any]]) -> str:
    header_cells = "".join(
        f"<th scope='col'>{escape(lbl)}</th>" for _, lbl in TABLE_COLUMNS)
    body_rows = []
    text_columns = {"round", "model_type", "commit", "filtered_on",
                    "filter_source", "scib_best_variant", "scib_sel_vs_best"}
    for r in rows:
        filter_text = _html_filter_text(*[r.get(k) for k, _ in TABLE_COLUMNS],
                                        r.get("variants"),
                                        r.get("filter_cluster_key"))
        tds = []
        for key, _ in TABLE_COLUMNS:
            v = _fmt(r.get(key))
            cls = "text" if key in text_columns else "num"
            if key == "commit" and r.get("dirty") is True:
                v += " !"
                cls += " warn"
            if key == "filter_source" and r.get(key) == "auto":
                v += " !"
                cls += " warn"
            if key == "scib_best" and r.get("scib_best") is not None:
                v += _scib_metric_tag(r.get("scib_metric"))
            title = ""
            if key == "n_variants" and r.get("variants"):
                title = f" title='{_attr(r.get('variants'))}'"
            tds.append(f"<td class='{cls}'{title}>{escape(v)}</td>")
        body_rows.append(
            f"<tr data-filter='{filter_text}'>" + "".join(tds) + "</tr>")
    return (
        "<div class='table-wrap'><table id='rounds-table'>"
        f"<thead><tr>{header_cells}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table></div>"
    )


def _render_retention(rows: List[Dict[str, Any]]) -> str:
    cell_counts = [_as_int(r.get("n_cells")) for r in rows]
    cell_counts = [v for v in cell_counts if v is not None]
    if not cell_counts:
        return "<p class='empty-state'>No cell counts recorded.</p>"
    start = cell_counts[0] or 0
    max_cells = max(cell_counts) or 1
    items = []
    for r in rows:
        cells = _as_int(r.get("n_cells"))
        removed = _as_int(r.get("removed"))
        retained = (100 * cells / start) if start and cells is not None else None
        width = 0 if cells is None or cells <= 0 else max(
            1.0, min(100.0, 100 * cells / max_cells))
        filter_text = _html_filter_text(r.get("round"), r.get("filtered_on"),
                                        r.get("filter_source"))
        items.append(
            f"<div class='timeline-row' data-filter='{filter_text}'>"
            f"<div class='timeline-label'><strong>{escape(_fmt(r.get('round')))}</strong>"
            f"<span>{escape(_fmt(cells))} cells</span></div>"
            f"<div class='bar-cell'><div class='bar-track'>"
            f"<span style='width:{width:.2f}%'></span></div>"
            f"<div class='timeline-meta'>"
            f"{escape(_fmt_pct(retained))} retained"
            f" · {escape(_fmt(removed))} removed"
            f" · {escape(_fmt_pct(r.get('removal_pct')))} round removal"
            f"</div></div></div>"
        )
    return "<div class='timeline'>" + "".join(items) + "</div>"


def _render_round_cards(rows: List[Dict[str, Any]]) -> str:
    cards = []
    open_attr = " open" if len(rows) <= 4 else ""
    for r in rows:
        filter_text = _html_filter_text(r.get("round"), r.get("filtered_on"),
                                        r.get("filter_cluster_key"),
                                        r.get("variants"),
                                        r.get("clusters_variant"))
        badges = [_badge(f"{_fmt(r.get('n_cells'))} cells", "good")]
        if r.get("filtered_on"):
            badges.append(_badge(f"filtered on {r.get('filtered_on')}", "neutral"))
        if r.get("filter_source") == "auto":
            badges.append(_badge("auto selection", "warn"))
        if r.get("dirty") is True:
            badges.append(_badge("dirty tree", "warn"))
        fields = [
            ("Updated", _fmt_date(r.get("updated"))),
            ("Genes", _fmt(r.get("n_genes"))),
            ("Model", _fmt(r.get("model_type"))),
            ("Annotation", _fmt(r.get("annotation"))),
            ("Seed", _fmt(r.get("seed"))),
            ("Commit", _fmt(r.get("commit"))),
            ("Clusters", _fmt(r.get("n_clusters"))),
            ("Flagged", _fmt(r.get("n_auto_flagged"))),
            ("Cluster source", _fmt(r.get("clusters_variant"))),
            ("Input cells", _fmt(r.get("input_cells"))),
            ("Output cells", _fmt(r.get("output_cells"))),
            ("Removed", f"{_fmt(r.get('removed'))} ({_fmt_pct(r.get('removal_pct'))})"),
            ("Filter key", _fmt(r.get("filter_cluster_key"))),
            ("scIB best", f"{_fmt(r.get('scib_best'))}{_scib_metric_tag(r.get('scib_metric'))}"
                          f" ({_fmt(r.get('scib_best_variant'))})"),
            ("Selected / parent best", _fmt(r.get("scib_sel_vs_best"))),
        ]
        detail = "".join(
            f"<div><dt>{escape(label)}</dt><dd>{escape(value)}</dd></div>"
            for label, value in fields)
        variants = r.get("variants")
        variants_html = ""
        if variants:
            variants_html = (
                "<div class='variant-list'><dt>Variants</dt>"
                f"<dd>{escape(str(variants))}</dd></div>")
        cards.append(
            f"<details class='round-card' data-filter='{filter_text}'{open_attr}>"
            f"<summary><span><strong>{escape(_fmt(r.get('round')))}</strong>"
            f"<small>{escape(_fmt(r.get('n_clusters')))} clusters"
            f" · {escape(_fmt(r.get('n_variants')))} variants</small></span>"
            f"<span class='badge-row'>{''.join(badges)}</span></summary>"
            f"<dl class='round-fields'>{detail}{variants_html}</dl></details>"
        )
    return "<div class='round-card-grid'>" + "".join(cards) + "</div>"


def _group_ledger(ledger: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for item in ledger:
        grouped.setdefault(str(item.get("round") or "unknown"), []).append(item)
    return grouped


def _render_ledger(ledger: Optional[List[Dict[str, Any]]]) -> str:
    ledger = ledger or []
    if not ledger:
        return "<p class='empty-state'>No keep/remove decisions recorded.</p>"
    groups = _group_ledger(ledger)
    parts = []
    for round_id, actions in groups.items():
        first = actions[0]
        filter_text = _html_filter_text(round_id, first.get("filtered_on"),
                                        first.get("cluster_key"),
                                        *[a.get("action") for a in actions])
        items = "".join(
            f"<li><code>{escape(str(a.get('action')))}</code></li>"
            for a in actions)
        meta = []
        if first.get("filtered_on"):
            meta.append(_badge(first["filtered_on"], "neutral"))
        if first.get("cluster_key"):
            meta.append(_badge(first["cluster_key"], "neutral"))
        parts.append(
            f"<details class='ledger-group' data-filter='{filter_text}' open>"
            f"<summary><strong>{escape(round_id)}</strong>"
            f"<span>{len(actions)} action(s)</span>"
            f"<span class='badge-row'>{''.join(meta)}</span></summary>"
            f"<ul class='action-list'>{items}</ul></details>"
        )
    return "<div class='ledger-list'>" + "".join(parts) + "</div>"


def _render_findings(findings: List[Dict[str, Any]], verified: bool) -> str:
    if not verified:
        return "<p class='empty-state'>Not checked. Run with --verify.</p>"
    if not findings:
        return "<p class='empty-state good'>No integrity issues detected.</p>"
    order = {"error": 0, "warn": 1, "info": 2}
    rows = [
        "<tr><th scope='col'>Severity</th><th scope='col'>Round</th>"
        "<th scope='col'>Check</th><th scope='col'>Detail</th></tr>"
    ]
    for f in sorted(findings, key=lambda x: order.get(x.get("severity"), 9)):
        sev = str(f.get("severity", "info"))
        rows.append(
            f"<tr><td class='{_attr(sev)}'>{escape(sev)}</td>"
            f"<td>{escape(str(f.get('round')))}</td>"
            f"<td>{escape(str(f.get('check')))}</td>"
            f"<td>{escape(str(f.get('message')))}</td></tr>")
    return "<div class='table-wrap'><table class='findings-table'>" + "".join(rows) + "</table></div>"


def render_html(md_table_rows, mermaid_src, findings, verified,
                ledger: Optional[List[Dict[str, Any]]] = None) -> str:
    generated = datetime.now().isoformat(timespec="seconds")
    metrics = _report_metrics(md_table_rows, findings, verified)
    retained = _fmt_pct(metrics["retention_pct"])
    css = """
        :root {
            --bg: #f6f7f2;
            --panel: #ffffff;
            --text: #17211c;
            --muted: #65706a;
            --line: #d9dfd7;
            --head: #153f45;
            --head-2: #4f5f38;
            --accent: #2f7d68;
            --accent-2: #b86f3c;
            --warn: #a86316;
            --risk: #ba3b3b;
            --good-bg: #e6f2e9;
            --warn-bg: #fff3df;
            --risk-bg: #fde8e8;
            --neutral-bg: #edf1f0;
            --shadow: 0 12px 28px rgba(23, 33, 28, 0.09);
        }
        * { box-sizing: border-box; }
        body {
            margin: 0;
            background: var(--bg);
            color: var(--text);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            line-height: 1.45;
        }
        .report-hero {
            background:
                linear-gradient(135deg, rgba(21, 63, 69, 0.98), rgba(79, 95, 56, 0.92)),
                linear-gradient(90deg, var(--head), var(--accent-2));
            color: #fff;
            padding: 34px 24px 30px;
        }
        .hero-inner, main {
            width: min(1180px, calc(100% - 32px));
            margin: 0 auto;
        }
        .eyebrow {
            margin: 0 0 8px;
            color: rgba(255, 255, 255, 0.74);
            font-size: 12px;
            font-weight: 700;
            letter-spacing: 0;
            text-transform: uppercase;
        }
        h1, h2, h3, p { margin-top: 0; }
        h1 {
            margin-bottom: 8px;
            font-size: clamp(32px, 5vw, 54px);
            line-height: 1.02;
            letter-spacing: 0;
        }
        .hero-copy {
            max-width: 760px;
            margin: 0;
            color: rgba(255, 255, 255, 0.82);
            font-size: 16px;
        }
        main { padding: 24px 0 46px; }
        .stat-grid {
            display: grid;
            grid-template-columns: repeat(5, minmax(0, 1fr));
            gap: 12px;
            margin-bottom: 18px;
        }
        .stat-card, .panel, .round-card, .ledger-group {
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 8px;
            box-shadow: var(--shadow);
        }
        .stat-card {
            min-height: 124px;
            padding: 16px;
            border-top: 4px solid var(--accent);
        }
        .stat-card.good { border-top-color: var(--accent); }
        .stat-card.warn, .stat-card.pending { border-top-color: var(--warn); }
        .stat-card.risk { border-top-color: var(--risk); }
        .stat-card span {
            display: block;
            color: var(--muted);
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
        }
        .stat-card strong {
            display: block;
            margin: 10px 0 6px;
            font-size: 28px;
            line-height: 1.1;
        }
        .stat-card p {
            margin: 0;
            color: var(--muted);
            font-size: 13px;
        }
        .panel {
            padding: 18px;
            margin-top: 18px;
        }
        .panel-header {
            display: flex;
            align-items: flex-start;
            justify-content: space-between;
            gap: 14px;
            margin-bottom: 14px;
        }
        .panel-header h2 {
            margin-bottom: 3px;
            font-size: 21px;
        }
        .panel-header p {
            margin: 0;
            color: var(--muted);
            font-size: 13px;
        }
        .filter-input {
            width: min(280px, 100%);
            height: 36px;
            border: 1px solid var(--line);
            border-radius: 6px;
            padding: 0 11px;
            font: inherit;
            background: #fff;
            color: var(--text);
        }
        .table-wrap {
            max-width: 100%;
            overflow: auto;
            border: 1px solid var(--line);
            border-radius: 8px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            min-width: 980px;
            background: #fff;
        }
        th, td {
            padding: 9px 10px;
            border-bottom: 1px solid var(--line);
            font-size: 13px;
            vertical-align: top;
        }
        th {
            position: sticky;
            top: 0;
            z-index: 1;
            background: #eef2ef;
            color: #26332d;
            text-align: right;
            font-weight: 700;
        }
        td.num { text-align: right; white-space: nowrap; }
        td.text, th:first-child, td:first-child { text-align: left; }
        td:first-child {
            font-weight: 700;
            color: #1b423b;
        }
        .warn { color: var(--warn); font-weight: 700; }
        .error { color: var(--risk); font-weight: 700; }
        .info { color: #476474; font-weight: 700; }
        .timeline {
            display: grid;
            gap: 13px;
        }
        .timeline-row {
            display: grid;
            grid-template-columns: minmax(140px, 210px) 1fr;
            gap: 16px;
            align-items: center;
        }
        .timeline-label strong,
        .timeline-label span {
            display: block;
        }
        .timeline-label span,
        .timeline-meta {
            color: var(--muted);
            font-size: 13px;
        }
        .bar-track {
            height: 16px;
            border-radius: 6px;
            background: #e7ebe6;
            overflow: hidden;
            border: 1px solid #d4dcd3;
        }
        .bar-track span {
            display: block;
            height: 100%;
            border-radius: 6px;
            background: linear-gradient(90deg, var(--accent), var(--accent-2));
        }
        .timeline-meta { margin-top: 5px; }
        .round-card-grid,
        .ledger-list {
            display: grid;
            gap: 12px;
        }
        .round-card, .ledger-group {
            box-shadow: none;
        }
        summary {
            cursor: pointer;
            list-style: none;
        }
        summary::-webkit-details-marker { display: none; }
        .round-card summary,
        .ledger-group summary {
            display: flex;
            justify-content: space-between;
            gap: 12px;
            align-items: center;
            padding: 14px 15px;
        }
        .round-card summary strong,
        .ledger-group summary strong {
            display: block;
            font-size: 16px;
        }
        .round-card summary small,
        .ledger-group summary span {
            color: var(--muted);
            font-size: 13px;
        }
        .badge-row {
            display: flex;
            flex-wrap: wrap;
            justify-content: flex-end;
            gap: 6px;
        }
        .badge {
            display: inline-flex;
            align-items: center;
            min-height: 24px;
            padding: 3px 8px;
            border-radius: 999px;
            background: var(--neutral-bg);
            color: #2b3833;
            font-size: 12px;
            font-weight: 700;
            white-space: nowrap;
        }
        .badge.good { background: var(--good-bg); color: #24533d; }
        .badge.warn { background: var(--warn-bg); color: var(--warn); }
        .round-fields {
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 0;
            margin: 0;
            padding: 0 15px 15px;
            border-top: 1px solid var(--line);
        }
        .round-fields div,
        .round-fields .variant-list {
            padding: 12px 12px 10px 0;
            min-width: 0;
        }
        .round-fields dt {
            color: var(--muted);
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
        }
        .round-fields dd {
            margin: 3px 0 0;
            overflow-wrap: anywhere;
        }
        .variant-list {
            grid-column: 1 / -1;
        }
        .lineage-frame {
            overflow: auto;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfcfa;
            padding: 12px;
        }
        pre.mermaid {
            margin: 0;
            min-height: 120px;
            font-size: 13px;
        }
        .action-list {
            margin: 0;
            padding: 0 15px 15px 38px;
            border-top: 1px solid var(--line);
            columns: 2 300px;
        }
        .action-list li {
            margin: 8px 0;
            break-inside: avoid;
        }
        code {
            background: #f1f3ef;
            border: 1px solid #dde4dc;
            border-radius: 5px;
            padding: 2px 5px;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            font-size: 12px;
        }
        .empty-state {
            margin: 0;
            padding: 14px;
            border: 1px solid var(--line);
            border-radius: 8px;
            background: #fbfcfa;
            color: var(--muted);
        }
        .empty-state.good {
            color: #24533d;
            background: var(--good-bg);
            border-color: #c6ddca;
        }
        [hidden] { display: none !important; }
        @media (max-width: 900px) {
            .stat-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
            .round-fields { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        }
        @media (max-width: 680px) {
            .hero-inner, main { width: min(100% - 24px, 1180px); }
            .report-hero { padding: 26px 0 24px; }
            main { padding-top: 14px; }
            .stat-grid { grid-template-columns: 1fr; }
            .panel { padding: 14px; }
            .panel-header,
            .round-card summary,
            .ledger-group summary {
                display: grid;
                justify-items: start;
            }
            .badge-row { justify-content: flex-start; }
            .timeline-row { grid-template-columns: 1fr; gap: 6px; }
            .round-fields { grid-template-columns: 1fr; }
            table { min-width: 860px; }
        }
    """
    js = """
        document.addEventListener('DOMContentLoaded', function () {
            var input = document.getElementById('round-filter');
            if (!input) return;
            var items = Array.prototype.slice.call(
                document.querySelectorAll('[data-filter]'));
            input.addEventListener('input', function () {
                var query = input.value.trim().toLowerCase();
                items.forEach(function (el) {
                    el.hidden = !!query && el.dataset.filter.indexOf(query) === -1;
                });
            });
        });
    """
    head = (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        "<title>Pipeline summary</title>"
        f"<style>{css}</style>"
        "<script src='https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js'>"
        "</script><script>mermaid.initialize({startOnLoad:true,theme:'neutral'});"
        "</script></head>"
    )
    return (
        f"{head}<body><header class='report-hero'><div class='hero-inner'>"
        f"<p class='eyebrow'>Integration rounds report</p>"
        f"<h1>Pipeline summary</h1>"
        f"<p class='hero-copy'>Generated {escape(generated)} · "
        f"{escape(_fmt(len(md_table_rows)))} round(s) · "
        f"{escape(retained)} retained from first to final round</p>"
        f"</div></header><main>"
        f"{_render_overview_cards(md_table_rows, findings, verified)}"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Rounds</h2><p>Per-round integration, filtering, and benchmark summary.</p>"
        f"</div><input class='filter-input' id='round-filter' type='search' "
        f"placeholder='Filter report' aria-label='Filter report'></div>"
        f"{_render_rounds_table(md_table_rows)}</section>"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Retention</h2><p>Cell count trajectory across the lineage.</p>"
        f"</div></div>{_render_retention(md_table_rows)}</section>"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Round details</h2><p>Expanded manifest fields for each round.</p>"
        f"</div></div>{_render_round_cards(md_table_rows)}</section>"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Lineage</h2><p>Parent-child flow with filtering deltas.</p>"
        f"</div></div><div class='lineage-frame'>"
        f"<pre class='mermaid'>{escape(mermaid_src)}</pre></div></section>"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Decisions ledger</h2><p>All recorded keep/remove actions.</p>"
        f"</div></div>{_render_ledger(ledger)}</section>"
        f"<section class='panel'><div class='panel-header'><div>"
        f"<h2>Integrity</h2><p>Verification results and reproducibility checks.</p>"
        f"</div></div>{_render_findings(findings, verified)}</section>"
        f"</main><script>{js}</script></body></html>"
    )


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
    parser.add_argument("--terminal-run-architecture", "--terminal_run_architecture",
                        dest="terminal_run_architecture", default=None,
                        metavar="KEY",
                        help="Architecture to source clusters/flagged from for the "
                             "terminal round(s) — those with no child round to "
                             "declare an acted-on variant. Accepts a variant name "
                             "(e.g. 'harmony') or a cluster key (e.g. "
                             "'leiden_harmony').")
    args = parser.parse_args()
    formats = {f.strip() for f in args.format.split(",") if f.strip()}

    explicit_order = None
    scan_dir = None
    if args.rounds_dir:
        if not os.path.isdir(args.rounds_dir):
            parser.error(f"rounds-dir not found: {args.rounds_dir}")
        manifest_paths = find_manifests(args.rounds_dir)
        manifests = [load_manifest(p) for p in manifest_paths]
        # Also pick up inspection-only sweep rounds — dirs with inspect_*
        # sidecars but no top-level round_manifest.json (e.g. a re-inspection
        # sweep of a prior round). Exclude ordinary rounds and their _sweep
        # siblings so their sidecars aren't counted as standalone rounds.
        exclude = {m["_dir"] for m in manifests}
        exclude |= {m["_dir"] + "_sweep" for m in manifests}
        manifests += discover_inspection_only_rounds(args.rounds_dir, exclude)
        if not manifests:
            parser.error(f"No round_manifest.json found under {args.rounds_dir}")
        scan_dir = args.rounds_dir
        default_out = os.path.join(args.rounds_dir, "summary")
    else:
        manifests = []
        for item in args.rounds:
            if os.path.isdir(item):
                mp = os.path.join(item, "round_manifest.json")
                if os.path.exists(mp):
                    manifests.append(load_manifest(mp))
                    continue
                syn = synthesize_inspection_round(item)
                if syn is None:
                    parser.error(f"No round_manifest.json (or inspect_* sweep "
                                 f"manifests) at: {item}")
                manifests.append(syn)
            else:
                if not os.path.exists(item):
                    parser.error(f"No round_manifest.json at: {item}")
                manifests.append(load_manifest(item))
        explicit_order = [m["_dir"] for m in manifests]  # trust supplied order
        default_out = "round_summary"

    output_dir = args.output_dir or default_out

    parent_map = build_parent_map(manifests, explicit_order)
    # Inspection-only rounds write no h5ad, so nothing's fingerprint points to
    # them and they dangle off their parent. Splice each onto the edge whose
    # filtering it informed before ordering/rendering.
    parent_map = splice_inspection_only_rounds(manifests, parent_map,
                                               explicit_order)
    # Discovery mode: topologically reorder. Explicit list: keep supplied order.
    if explicit_order is None:
        manifests = order_by_lineage(manifests, parent_map)
    print(f"  Found {len(manifests)} round(s): "
          f"{', '.join(m['_round_id'] for m in manifests)}")

    by_dir = {m["_dir"]: m for m in manifests}

    # The variant a round was filtered on for its CHILD is the architecture whose
    # clustering the human reviewed and acted on in the (parent) round. Map each
    # round's _dir → that acted-on variant so extract_row can source the matching
    # inspection manifest for its clusters/flagged counts.
    acted_on: Dict[str, Optional[str]] = {}
    for m in manifests:
        pdir = parent_map.get(m["_dir"])
        v = _g(m, "decisions", "filtered_on", "variant")
        if pdir and v and pdir not in acted_on:
            acted_on[pdir] = v

    # Terminal round(s) have no child to declare an acted-on variant. When the
    # user pins one with --terminal-run-architecture, apply it to the leaves
    # (rounds that are nobody's parent).
    terminal_variant = variant_from_cluster_key(args.terminal_run_architecture)
    if terminal_variant:
        parents = {p for p in parent_map.values() if p}
        for m in manifests:
            if m["_dir"] in parents:
                continue  # not a leaf; its acted-on variant comes from its child
            acted_on[m["_dir"]] = terminal_variant
            available = set(find_inspection_manifests(m["_dir"]))
            if available and terminal_variant not in available:
                print(f"  [WARN] --terminal-run-architecture "
                      f"'{terminal_variant}' not found for {m['_round_id']}; "
                      f"available: {', '.join(sorted(available))}")

    rows = []
    for m in manifests:
        pdir = parent_map.get(m["_dir"])
        rows.append(extract_row(m, by_dir.get(pdir) if pdir else None,
                                acted_on_variant=acted_on.get(m["_dir"])))
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
            h.write(render_html(rows, mermaid_src, findings, args.verify, ledger))
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

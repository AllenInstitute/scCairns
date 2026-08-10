#!/usr/bin/env python3
"""Record provenance for an externally filtered AnnData object.

``cairns record-filter`` does not apply filtering. It audits a completed
interactive/manual filtering step by comparing an input h5ad with an output
h5ad, then writes the normalized decisions block consumed by ``cairns
summarize``.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple

import anndata as ad
import pandas as pd

from .config import (
    fingerprint_file,
    read_parent_provenance,
    update_round_manifest,
    write_yaml_record,
)


def variant_from_embedding_key(key: Optional[str]) -> Optional[str]:
    """Map an embedding key to the variant label used in summary output."""
    if not key:
        return None
    if key == "X_scVI":
        return "default"
    if key == "X_pca_harmony":
        return "harmony"
    if key.startswith("X_scVI_"):
        return key[len("X_scVI_"):]
    return key


def infer_cluster_key(embedding_key: Optional[str]) -> Optional[str]:
    """Infer the likely Leiden cluster key from a latent embedding key."""
    if not embedding_key:
        return None
    if embedding_key == "X_scVI":
        return "leiden"
    if embedding_key == "X_pca_harmony":
        return "leiden_harmony"
    if embedding_key.startswith("X_scVI_"):
        return f"leiden_{embedding_key[len('X_scVI_'):]}"
    return None


def infer_umap_key(embedding_key: Optional[str]) -> Optional[str]:
    """Infer the likely UMAP key from a latent embedding key."""
    if not embedding_key:
        return None
    if embedding_key == "X_scVI":
        return "X_umap"
    if embedding_key == "X_pca_harmony":
        return "X_umap_harmony"
    if embedding_key.startswith("X_scVI_"):
        return f"X_umap_{embedding_key[len('X_scVI_'):]}"
    return None


def _read_h5ad(path: str, label: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} h5ad not found: {path}")
    try:
        return ad.read_h5ad(path, backed="r")
    except Exception as exc:  # pragma: no cover - exact reader errors vary
        raise ValueError(f"{label} h5ad is not readable: {path} ({exc})") from exc


def _close_h5ad(adata) -> None:
    file_obj = getattr(adata, "file", None)
    close = getattr(file_obj, "close", None)
    if callable(close):
        close()


def _relative_to_output(path: str, output_dir: str) -> str:
    abs_path = os.path.abspath(path)
    abs_output = os.path.abspath(output_dir)
    try:
        rel = os.path.relpath(abs_path, abs_output)
    except ValueError:
        return abs_path
    return rel if not rel.startswith("..") else abs_path


def _warn(message: str) -> None:
    print(f"[WARN] {message}", file=sys.stderr)


def _action_summary(action: Dict[str, Any]) -> str:
    if action["type"] == "interactive_filter_cluster":
        return (
            f"interactive_filter: cluster {action.get('cluster', '?')} "
            f"(-{action.get('n_removed', 0)})"
        )
    return f"interactive_filter: {action.get('reason', 'manual filter')} (-{action.get('n_removed', 0)})"


def _cluster_actions(
    input_adata,
    removed_ids: Sequence[str],
    cluster_key: Optional[str],
    filter_used: str,
) -> Tuple[List[Dict[str, Any]], bool]:
    """Return per-cluster action records when a usable cluster key is present."""
    if not cluster_key or cluster_key not in input_adata.obs.columns:
        return [], False
    if not removed_ids:
        return [], True

    removed_obs = input_adata.obs.loc[list(removed_ids)]
    counts = removed_obs[cluster_key].astype(str).value_counts()
    actions = []
    for cluster in sorted(counts.index, key=str):
        actions.append({
            "type": "interactive_filter_cluster",
            "cluster": str(cluster),
            "n_removed": int(counts[cluster]),
            "reason": filter_used,
        })
    return actions, True


def record_filter(
    *,
    input_h5ad: str,
    output_h5ad: str,
    output_dir: str,
    embedding_key: str,
    filter_used: str,
    notes: str = "",
    cluster_key: Optional[str] = None,
    umap_key: Optional[str] = None,
    seed: Optional[int] = None,
    batch_key: Optional[str] = None,
) -> str:
    """Record an already-applied filter by comparing input/output h5ad files."""
    input_adata = _read_h5ad(input_h5ad, "input")
    output_adata = _read_h5ad(output_h5ad, "output")
    try:
        input_names = [str(x) for x in input_adata.obs_names]
        output_names = [str(x) for x in output_adata.obs_names]
        input_set = set(input_names)
        output_set = set(output_names)
        unknown_output = sorted(output_set - input_set)
        if unknown_output:
            preview = ", ".join(unknown_output[:5])
            more = "" if len(unknown_output) <= 5 else f", ... +{len(unknown_output) - 5} more"
            raise ValueError(
                "output h5ad contains cell IDs not present in input h5ad: "
                f"{preview}{more}"
            )

        if embedding_key not in input_adata.obsm:
            _warn(f"embedding key not found in input.obsm: {embedding_key}")

        resolved_cluster_key = cluster_key or infer_cluster_key(embedding_key)
        if resolved_cluster_key and resolved_cluster_key not in input_adata.obs.columns:
            _warn(f"cluster key not found in input.obs: {resolved_cluster_key}")
            usable_cluster_key = None
        else:
            usable_cluster_key = resolved_cluster_key
        if not usable_cluster_key:
            _warn("no usable cluster key found; recording a generic interactive action")

        resolved_umap_key = umap_key or infer_umap_key(embedding_key)
        if resolved_umap_key and resolved_umap_key not in input_adata.obsm:
            _warn(f"UMAP key not found in input.obsm: {resolved_umap_key}")

        removed_ids = [name for name in input_names if name not in output_set]
        kept_ids = output_names
        n_before = int(input_adata.n_obs)
        n_after = int(output_adata.n_obs)
        n_removed = n_before - n_after

        actions, used_clusters = _cluster_actions(
            input_adata, removed_ids, usable_cluster_key, filter_used
        )
        if not used_clusters:
            actions = [{
                "type": "interactive_filter",
                "n_removed": n_removed,
                "reason": filter_used,
            }]

        os.makedirs(output_dir, exist_ok=True)
        cells_to_keep_path = os.path.join(output_dir, "cells_to_keep.csv")
        cells_removed_path = os.path.join(output_dir, "cells_removed.csv")
        decisions_applied_path = os.path.join(output_dir, "decisions_applied.yaml")

        pd.DataFrame({"obs_name": kept_ids}).to_csv(cells_to_keep_path, index=False)
        pd.DataFrame({"obs_name": removed_ids}).to_csv(cells_removed_path, index=False)

        timestamp = datetime.now().isoformat()
        applied = {
            "timestamp": timestamp,
            "n_cells_before": n_before,
            "n_cells_after": n_after,
            "n_removed": n_removed,
            "actions": actions,
            "filter_used": filter_used,
            "notes": notes,
        }
        write_yaml_record(decisions_applied_path, applied)

        manifest = {
            "timestamp": timestamp,
            "input_h5ad": input_h5ad,
            "input_fingerprint": fingerprint_file(input_h5ad),
            "parent_round": read_parent_provenance(input_h5ad),
            "output_h5ad": output_h5ad,
            "output_fingerprint": fingerprint_file(output_h5ad),
            "filter_used": filter_used,
            "notes": notes,
            "cluster_key": resolved_cluster_key,
            "batch_key": batch_key,
            "filtered_on": {
                "variant": variant_from_embedding_key(embedding_key),
                "cluster_key": resolved_cluster_key,
                "latent_key": embedding_key,
                "umap_key": resolved_umap_key,
                "source": "interactive",
            },
            "input_cells": n_before,
            "output_cells": n_after,
            "input_genes": int(input_adata.n_vars),
            "output_genes": int(output_adata.n_vars),
            "removed": n_removed,
            "removal_pct": round(n_removed / n_before * 100, 1) if n_before else 0.0,
            "actions_summary": [_action_summary(a) for a in actions],
            "output_files": {
                "filtered_h5ad": _relative_to_output(output_h5ad, output_dir),
                "cells_to_keep": os.path.basename(cells_to_keep_path),
                "cells_removed": os.path.basename(cells_removed_path),
                "decisions_applied": os.path.basename(decisions_applied_path),
            },
        }

        manifest_path = update_round_manifest(
            output_dir, "decisions", manifest, seed=seed
        )
        print(f"Saved {manifest_path}")
        return manifest_path
    finally:
        _close_h5ad(output_adata)
        _close_h5ad(input_adata)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Record provenance for an externally filtered h5ad pair."
    )
    parser.add_argument("--input-h5ad", required=True, help="Input h5ad before filtering.")
    parser.add_argument("--output-h5ad", required=True, help="Output h5ad after filtering.")
    parser.add_argument("--output-dir", required=True, help="Directory for manifest/sidecars.")
    parser.add_argument("--embedding-key", required=True, help="Embedding used for interactive filtering.")
    parser.add_argument("--filter-used", required=True, help="Human-readable filter description.")
    parser.add_argument("--notes", default="", help="Free-text notes for the filtering step.")
    parser.add_argument("--cluster-key", default=None, help="Cluster key to summarize removals by.")
    parser.add_argument("--umap-key", default=None, help="UMAP key used during filtering.")
    parser.add_argument("--seed", type=int, default=None, help="Seed to stamp into the manifest.")
    parser.add_argument("--batch-key", default=None, help="Batch key recorded for provenance.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        record_filter(
            input_h5ad=args.input_h5ad,
            output_h5ad=args.output_h5ad,
            output_dir=args.output_dir,
            embedding_key=args.embedding_key,
            filter_used=args.filter_used,
            notes=args.notes,
            cluster_key=args.cluster_key,
            umap_key=args.umap_key,
            seed=args.seed,
            batch_key=args.batch_key,
        )
    except (FileNotFoundError, ValueError) as exc:
        parser.exit(1, f"cairns record-filter: error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

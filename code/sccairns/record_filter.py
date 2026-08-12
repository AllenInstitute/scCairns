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
from collections import Counter
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
    if adata is None:
        return
    file_obj = getattr(adata, "file", None)
    close = getattr(file_obj, "close", None)
    if callable(close):
        close()


# Column names a cell-ID export is likely to use. "x" is what R's
# `write.csv(data.frame(x = Cells(obj)))` produces; "" is the unnamed rowname
# column write.csv always prepends.
_ID_COLUMN_NAMES = {
    "x", "obs_name", "obs_names", "cell_id", "cell_ids", "cell", "cells",
    "barcode", "barcodes", "index", "rowname", "row.names",
}


def _looks_like_row_numbers(values: Sequence[str]) -> bool:
    """True for R's implicit 1..N rowname column, which is never a cell ID."""
    return bool(values) and all(v.strip().isdigit() for v in values)


def _pick_id_column(frame: "pd.DataFrame") -> int:
    """Index of the column holding cell IDs.

    Prefers a recognizably-named column, because `write.csv` prepends an unnamed
    column of row numbers — taking column 0 blindly reads those instead of the
    barcodes. Falls back to skipping a leading all-numeric column.
    """
    for position, name in enumerate(frame.columns):
        if str(name).strip().strip('"').lower() in _ID_COLUMN_NAMES:
            return position
    if frame.shape[1] > 1 and _looks_like_row_numbers(
            [str(v) for v in frame.iloc[:, 0]]):
        return 1
    return 0


def _read_cell_ids(path: str, label: str) -> List[str]:
    """Read cell IDs from a CSV/TSV column or a bare ID-per-line list.

    Lets an external tool — Seurat, Loupe, a notebook — hand back only the cells
    it kept, which is all this module needs from the filtered side. Reads this
    module's own ``cells_to_keep.csv``, a bare list, and both common R exports:
    ``write.csv(data.frame(x = Cells(obj)))`` (row numbers then barcodes) and
    ``write.csv(obj@meta.data)`` (barcodes as rownames).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} cell list not found: {path}")
    try:
        probe = pd.read_csv(path, nrows=1, header=None, dtype=str,
                            keep_default_na=False)
    except Exception as exc:
        raise ValueError(f"{label} cell list is not readable: {path} ({exc})") from exc
    if probe.empty:
        raise ValueError(f"{label} cell list is empty: {path}")

    first_row = [str(v).strip().strip('"').lower() for v in probe.iloc[0]]
    has_header = any(v in _ID_COLUMN_NAMES or v == "" for v in first_row)
    frame = pd.read_csv(path, dtype=str, keep_default_na=False,
                        header=0 if has_header else None)
    ids = [str(x).strip() for x in frame.iloc[:, _pick_id_column(frame)]]
    ids = [x for x in ids if x]
    if not ids:
        raise ValueError(f"{label} cell list has no cell IDs: {path}")

    counts = Counter(ids)
    duplicates = sorted(x for x, n in counts.items() if n > 1)
    if duplicates:
        preview = ", ".join(duplicates[:5])
        more = "" if len(duplicates) <= 5 else f", ... +{len(duplicates) - 5} more"
        raise ValueError(
            f"{label} cell list contains duplicate cell IDs: {preview}{more}"
        )
    return ids


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
    input_h5ad: Optional[str] = None,
    output_h5ad: Optional[str] = None,
    input_cell_ids: Optional[str] = None,
    output_cell_ids: Optional[str] = None,
    output_dir: str,
    embedding_key: Optional[str] = None,
    filter_used: str,
    notes: str = "",
    cluster_key: Optional[str] = None,
    umap_key: Optional[str] = None,
    seed: Optional[int] = None,
    batch_key: Optional[str] = None,
) -> str:
    """Record an already-applied filter by comparing the cells before and after.

    Each side may be given as an ``.h5ad`` or as a list of cell IDs — the
    filtered side is only ever read for the IDs that survived, so a tool that
    cannot write AnnData (Seurat, Loupe) can hand back a CSV instead. An h5ad
    input additionally supplies the gene count and the cluster column needed for
    a per-cluster breakdown; with an ID list those degrade to a single action and
    a null gene count rather than failing.

    ``embedding_key`` is optional: a filter applied at ingest, before any
    embedding or clustering exists, records a generic action.
    """
    if bool(input_h5ad) == bool(input_cell_ids):
        raise ValueError("Provide exactly one of input_h5ad or input_cell_ids.")
    if bool(output_h5ad) == bool(output_cell_ids):
        raise ValueError("Provide exactly one of output_h5ad or output_cell_ids.")

    input_adata = _read_h5ad(input_h5ad, "input") if input_h5ad else None
    output_adata = _read_h5ad(output_h5ad, "output") if output_h5ad else None
    try:
        if input_adata is not None:
            input_names = [str(x) for x in input_adata.obs_names]
        else:
            input_names = _read_cell_ids(input_cell_ids, "input")
        if output_adata is not None:
            output_names = [str(x) for x in output_adata.obs_names]
        else:
            output_names = _read_cell_ids(output_cell_ids, "output")

        input_set = set(input_names)
        output_set = set(output_names)
        unknown_output = sorted(output_set - input_set)
        if unknown_output:
            preview = ", ".join(unknown_output[:5])
            more = "" if len(unknown_output) <= 5 else f", ... +{len(unknown_output) - 5} more"
            raise ValueError(
                "filtered side contains cell IDs not present in the input: "
                f"{preview}{more}\n"
                "  Cell IDs must survive the round trip unchanged. Tools that "
                "re-label cells (e.g. Seurat merges adding '_1' suffixes or "
                "sample prefixes) break the match; export the original IDs."
            )

        # An embedding is absent by design at ingest; only complain about a key
        # that was named and then not found.
        if embedding_key and input_adata is not None \
                and embedding_key not in input_adata.obsm:
            _warn(f"embedding key not found in input.obsm: {embedding_key}")

        resolved_cluster_key = cluster_key or infer_cluster_key(embedding_key)
        usable_cluster_key = resolved_cluster_key
        if resolved_cluster_key and input_adata is None:
            _warn("cluster breakdown needs an h5ad input; recording a generic "
                  "interactive action")
            usable_cluster_key = None
        elif resolved_cluster_key and resolved_cluster_key not in input_adata.obs.columns:
            _warn(f"cluster key not found in input.obs: {resolved_cluster_key}")
            usable_cluster_key = None

        resolved_umap_key = umap_key or infer_umap_key(embedding_key)
        if resolved_umap_key and input_adata is not None \
                and resolved_umap_key not in input_adata.obsm:
            _warn(f"UMAP key not found in input.obsm: {resolved_umap_key}")

        removed_ids = [name for name in input_names if name not in output_set]
        kept_ids = output_names
        n_before = len(input_names)
        n_after = len(output_names)
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

        output_files = {
            "cells_to_keep": os.path.basename(cells_to_keep_path),
            "cells_removed": os.path.basename(cells_removed_path),
            "decisions_applied": os.path.basename(decisions_applied_path),
        }
        if output_h5ad:
            output_files["filtered_h5ad"] = _relative_to_output(
                output_h5ad, output_dir)
        else:
            output_files["filtered_cell_ids"] = _relative_to_output(
                output_cell_ids, output_dir)

        manifest = {
            "timestamp": timestamp,
            "input_h5ad": input_h5ad,
            "input_cell_ids": input_cell_ids,
            "input_fingerprint": fingerprint_file(input_h5ad or input_cell_ids),
            # Only an h5ad sits beside a round manifest, so only it can carry a
            # lineage pointer; an ID list makes this the root of the chain.
            "parent_round": read_parent_provenance(input_h5ad),
            "output_h5ad": output_h5ad,
            "output_cell_ids": output_cell_ids,
            "output_fingerprint": fingerprint_file(output_h5ad or output_cell_ids),
            # What was actually compared, so a null gene count or a missing
            # cluster breakdown reads as a consequence rather than a gap.
            "compared": {
                "input": "h5ad" if input_h5ad else "cell_ids",
                "output": "h5ad" if output_h5ad else "cell_ids",
            },
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
            "input_genes": int(input_adata.n_vars) if input_adata is not None else None,
            "output_genes": int(output_adata.n_vars) if output_adata is not None else None,
            "removed": n_removed,
            "removal_pct": round(n_removed / n_before * 100, 1) if n_before else 0.0,
            "actions_summary": [_action_summary(a) for a in actions],
            "output_files": output_files,
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
        description="Record provenance for a filtering step performed outside "
                    "scCairns. Give each side as an .h5ad or as a list of cell "
                    "IDs (one per line, or a column of a CSV) — the latter lets "
                    "Seurat/Loupe hand back only the cells kept.",
    )
    before = parser.add_mutually_exclusive_group(required=True)
    before.add_argument("--input-h5ad", help="Input h5ad before filtering.")
    before.add_argument("--input-cell-ids",
                        help="Cell IDs before filtering (CSV or one per line).")
    after = parser.add_mutually_exclusive_group(required=True)
    after.add_argument("--output-h5ad", help="Output h5ad after filtering.")
    after.add_argument("--output-cell-ids",
                       help="Cell IDs that survived filtering (CSV or one per line).")
    parser.add_argument("--output-dir", required=True, help="Directory for manifest/sidecars.")
    parser.add_argument("--filter-used", required=True, help="Human-readable filter description.")
    parser.add_argument("--embedding-key", default=None,
                        help="Embedding the filtering was done on. Omit for a "
                             "filter applied at ingest, before any embedding "
                             "or clustering exists.")
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
            input_cell_ids=args.input_cell_ids,
            output_cell_ids=args.output_cell_ids,
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

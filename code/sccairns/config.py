"""Shared V1 pipeline configuration and provenance helpers."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Optional


PIPELINE_VERSION = 1


SNS_MARKER_SETS = {
    "Pan-neuronal": ["Snap25", "Tubb3", "Rbfox3", "Elavl4", "Isl1"],
    "Noradrenergic": ["Th", "Dbh", "Ddc", "Slc6a2"],
    "Cholinergic": ["Chat", "Slc18a3", "Slc5a7"],
    "Glutamatergic": ["Slc17a6", "Slc17a7"],
    "GABAergic": ["Slc32a1", "Gad1", "Gad2"],
    "Nitrergic": ["Nos1"],
    "Neuropeptides": ["Npy", "Sst", "Vip", "Pdyn", "Chga"],
    "Transcription": ["Phox2b", "Phox2a", "Sox6", "Shox2"],
    "Satellite glia": ["Sox10", "Fabp7", "S100b"],
    "Immune": ["Ptprc", "Cd68"],
    "Mitochondrial": ["mt-Co1", "mt-Co2", "mt-Cytb"],
}


DEFAULT_CONFIG: Dict[str, Any] = {
    "pipeline_version": PIPELINE_VERSION,
    "reproducibility": {
        # Global RNG seed applied to scvi-tools (which seeds Python random,
        # numpy, and torch via Lightning's seed_everything) and numpy before
        # any stochastic step. Set to null to leave RNG state unseeded
        # (non-reproducible). The resolved seed is recorded in
        # round_manifest.json so every round is traceable to its seed.
        "seed": 0,
    },
    "data": {
        "input_h5ad": None,
        "output_dir": "../results",
        "counts_layer": "counts",
        "batch_key": "data_origin",
        "categorical_covariate_keys": ["tech"],
        "continuous_covariate_keys": [],
        "sample_key": None,
        "technology_key": "tech",
        "species": "mouse",
        "gene_symbol_case": "preserve",
    },
    "qc": {
        "enabled": True,
        "min_genes": 500,
        "min_cells": 3,
        "mt_gene_patterns": ["mt-"],
        "ribo_gene_patterns": ["rps", "rpl"],
        "hb_gene_pattern": "^hb[^(p)]",
        "skip_filter_after_round_1": True,
        "skip_filter": False,
    },
    "integration": {
        "model_type": "scvi",
        "n_hidden": 256,
        "n_layers": 3,
        "n_latent": 32,
        "dispersion": "gene-cell",
        "gene_likelihood": "nb",
        "max_epochs": 200,
        "early_stopping_patience": 20,
        "batch_size": 256,
        "num_workers": 4,
        "hvg": {
            "n_top_genes": 3000,
            "batch_key": "tech",
            "flavor": "seurat_v3",
            "min_batches": None,
        },
        "sweep": [],
        # Optional Harmony integration run alongside scVI. Produces
        # obsm["X_pca_harmony"] (+ X_umap_harmony / leiden_harmony), giving a
        # same-run scVI-vs-Harmony comparison in scIB benchmarking. Runs once
        # per round in both single-model and sweep modes. batch_key defaults to
        # data.batch_key when null; n_pcs is the PCA dimensionality Harmony
        # corrects. hvg_from (sweep mode only) names a sweep entry whose HVG
        # selection Harmony should borrow; null = use the top-level
        # integration.hvg spec (the default, identical to single-model scVI).
        "harmony": {
            "enabled": False,
            "batch_key": None,
            "n_pcs": 30,
            "hvg_from": None,
        },
    },
    "annotation": {
        "enabled": False,
        "method": "scanvi",
        "labels_key": None,
        "unlabeled_category": "Unknown",
        "prediction_key": "scanvi_label",
        "confidence_key": "scanvi_confidence",
        "min_confidence": 0.5,
    },
    "embedding": {
        "n_neighbors": 30,
        "umap_min_dist": 0.4,
        "umap_spread": 3.0,
        "leiden_resolution": 0.3,
    },
    "inspection": {
        "cluster_key": "auto",
        "latent_key": "auto",
        "umap_key": "auto",
        "markers_json": None,
        "marker_sets": copy.deepcopy(SNS_MARKER_SETS),
        "neuronal_markers": [
            "Th",
            "Snap25",
            "Phox2b",
            "Dbh",
            "Chat",
            "Slc18a2",
            "Slc17a6",
        ],
        "neuronal_cutoff": 0.50,
        "marker_threshold": 0.0,
        "pca_color_gene": "Snap25",
        "auto_flag": {
            "mt_threshold": 15.0,
            "min_genes_threshold": 400,
            "min_cells": 20,
            "single_batch_threshold": 0.90,
            "silhouette_threshold": -0.05,
        },
        # Marker-based per-cell contamination flagging (cluster-independent).
        # panels defaults to flag_contamination.DEFAULT_CONTAM_PANELS when unset;
        # override with {label: [gene, ...]} to score custom lineages.
        "contamination": {
            "z_thresh": 2.0,
            "min_genes": 2,
            "layer": None,   # None → adata.X (log-normalized)
            "panels": None,
        },
    },
    "benchmark": {
        "enabled": True,
        "batch_key": "auto",
        "label_key": None,
    },
    "decisions": {
        "ignore_failed_queries": False,
    },
}


class ConfigError(ValueError):
    """Raised when a pipeline config is invalid."""


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    """Return base recursively updated by override without mutating either."""
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if (
            isinstance(value, Mapping)
            and isinstance(merged.get(key), MutableMapping)
        ):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def set_if_provided(config: Dict[str, Any], path: Iterable[str], value: Any) -> None:
    """Set a nested config value if value is not None."""
    if value is None:
        return
    target = config
    parts = list(path)
    for part in parts[:-1]:
        target = target[part]
    target[parts[-1]] = value


def _load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml  # type: ignore[import-untyped]
    except ImportError as exc:
        raise ConfigError(
            "PyYAML is required to read --config YAML files. "
            "Install pyyaml or use the project Docker environment."
        ) from exc

    with open(path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ConfigError(f"Config must contain a mapping at top level: {path}")
    return loaded


def _resolve_path(value: Optional[str], base_dir: Path) -> Optional[str]:
    if value in (None, ""):
        return value
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def _resolve_config_paths(config: Dict[str, Any], config_path: Path) -> None:
    base_dir = config_path.parent
    config["data"]["input_h5ad"] = _resolve_path(
        config["data"].get("input_h5ad"), base_dir
    )
    config["data"]["output_dir"] = _resolve_path(
        config["data"].get("output_dir"), base_dir
    )
    config["inspection"]["markers_json"] = _resolve_path(
        config["inspection"].get("markers_json"), base_dir
    )


def load_pipeline_config(
    config_path: Optional[str] = None,
    *,
    legacy_output_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Load and merge a V1 pipeline config.

    When config_path is omitted, DEFAULT_CONFIG is returned. legacy_output_dir
    lets the current scripts preserve their historical no-config output dirs.
    """
    config = copy.deepcopy(DEFAULT_CONFIG)
    if legacy_output_dir is not None and config_path is None:
        config["data"]["output_dir"] = legacy_output_dir

    if config_path:
        path = Path(config_path).expanduser().resolve()
        loaded = _load_yaml(path)
        version = loaded.get("pipeline_version")
        if version != PIPELINE_VERSION:
            raise ConfigError(
                f"Unsupported pipeline_version {version!r}; expected "
                f"{PIPELINE_VERSION}."
            )
        config = deep_merge(config, loaded)
        _resolve_config_paths(config, path)

    return config


def validate_config(config: Mapping[str, Any], *, mode: str) -> None:
    """Validate common config fields after CLI overrides have been applied."""
    if config.get("pipeline_version") != PIPELINE_VERSION:
        raise ConfigError(
            f"Unsupported pipeline_version {config.get('pipeline_version')!r}; "
            f"expected {PIPELINE_VERSION}."
        )
    if not config.get("data", {}).get("input_h5ad"):
        raise ConfigError("Missing input h5ad. Set data.input_h5ad or --input.")
    if not config.get("data", {}).get("output_dir"):
        raise ConfigError("Missing output directory. Set data.output_dir or --output-dir.")

    model_type = config["integration"].get("model_type")
    if model_type != "scvi":
        raise ConfigError("V1 supports integration.model_type: scvi only.")

    if config["integration"].get("dispersion") not in {
        "gene",
        "gene-batch",
        "gene-label",
        "gene-cell",
    }:
        raise ConfigError("Invalid integration.dispersion.")

    if config["integration"].get("gene_likelihood") not in {"zinb", "nb", "poisson"}:
        raise ConfigError("Invalid integration.gene_likelihood.")

    if config["integration"]["hvg"].get("flavor") not in {
        "seurat_v3",
        "seurat",
        "cell_ranger",
    }:
        raise ConfigError("Invalid integration.hvg.flavor.")

    harmony = config["integration"].get("harmony", {}) or {}
    if harmony.get("enabled"):
        n_pcs = harmony.get("n_pcs", 30)
        if not isinstance(n_pcs, int) or n_pcs < 2:
            raise ConfigError("integration.harmony.n_pcs must be an integer >= 2.")
        # hvg_from is only meaningful in sweep mode (it names a sweep entry);
        # in single-model mode Harmony always uses the top-level integration.hvg.
        hvg_from = harmony.get("hvg_from")
        if hvg_from and not config["integration"].get("sweep"):
            raise ConfigError(
                "integration.harmony.hvg_from names a sweep config entry and is "
                "only valid in sweep mode. Leave it null for single-model runs "
                "(Harmony uses integration.hvg)."
            )

    annotation = config.get("annotation", {})
    if annotation.get("enabled"):
        if annotation.get("method") != "scanvi":
            raise ConfigError("V1 supports annotation.method: scanvi only.")
        if mode == "integration" and config["integration"].get("sweep"):
            raise ConfigError("scANVI annotation is only supported for single-model runs in V1.")
        if not annotation.get("labels_key"):
            raise ConfigError(
                "annotation.labels_key is required when annotation.enabled is true."
            )
        if not annotation.get("unlabeled_category"):
            raise ConfigError(
                "annotation.unlabeled_category is required when annotation.enabled is true."
            )


def write_resolved_config(output_dir: str, config: Mapping[str, Any]) -> str:
    """Write the fully resolved config to run_config_resolved.yml."""
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "run_config_resolved.yml")
    try:
        import yaml  # type: ignore[import-untyped]

        with open(out_path, "w", encoding="utf-8") as handle:
            yaml.safe_dump(dict(config), handle, sort_keys=False)
    except ImportError:
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(config, handle, indent=2)
            handle.write("\n")
    return out_path


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_yaml_record(path: str, data: Mapping[str, Any]) -> str:
    """Serialize a record to a YAML file safely.

    Uses ``yaml.safe_dump`` so values containing quotes, colons, newlines, or
    other YAML-significant characters are escaped correctly. When PyYAML is
    unavailable, falls back to ``json.dump`` — JSON is a subset of YAML, so the
    output remains valid YAML and round-trips through ``yaml.safe_load``. This
    replaces hand-built f-string serialization, which could emit invalid or
    misleading YAML for arbitrary reason/query strings.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    safe = _json_safe(dict(data))
    try:
        import yaml  # type: ignore[import-untyped]

        with open(path, "w", encoding="utf-8") as handle:
            yaml.safe_dump(safe, handle, sort_keys=False, default_flow_style=False)
    except ImportError:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(safe, handle, indent=2)
            handle.write("\n")
    return path


def write_command_args(output_dir: str, stage: str, args: Mapping[str, Any]) -> str:
    """Record command arguments for a stage in command_args.json.

    ``commands[stage]`` holds the most recent invocation of each stage (kept for
    backward compatibility). ``history`` is an append-only list recording every
    invocation in order, so re-running a stage in the same output directory no
    longer silently overwrites the prior record — the full sequence of what was
    run in the round is preserved.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "command_args.json")
    payload: Dict[str, Any] = {}
    if os.path.exists(out_path):
        try:
            with open(out_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except json.JSONDecodeError:
            payload = {}
    entry = {
        "stage": stage,
        "timestamp": datetime.now().isoformat(),
        "argv": sys.argv,
        "args": _json_safe(dict(args)),
    }
    payload.setdefault("commands", {})
    payload["commands"][stage] = entry
    payload.setdefault("history", [])
    payload["history"].append(entry)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    return out_path


def collect_package_versions(packages: Iterable[str]) -> Dict[str, Optional[str]]:
    """Collect installed package versions without importing heavy modules."""
    versions: Dict[str, Optional[str]] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def set_global_seed(seed: Optional[int]) -> Optional[int]:
    """Seed all RNGs used by the pipeline and return the applied seed.

    When seed is None, RNG state is left untouched (non-reproducible) and None
    is returned. Otherwise scvi-tools is seeded via ``scvi.settings.seed``,
    which uses Lightning's ``seed_everything`` to seed Python ``random``,
    NumPy, and Torch. If scvi-tools is unavailable (e.g. config-only tooling),
    NumPy and ``random`` are seeded directly as a best-effort fallback so the
    return value still reflects what was applied.
    """
    if seed is None:
        return None
    seed = int(seed)
    try:
        import scvi  # type: ignore[import-untyped]

        scvi.settings.seed = seed
    except Exception:
        # Fallback: seed the standard RNGs directly. scvi.settings.seed would
        # normally cover these, but we still want deterministic numpy/random
        # behaviour when scvi is not importable.
        import random as _random

        _random.seed(seed)
        try:
            import numpy as _np

            _np.random.seed(seed)
        except Exception:
            pass
    return seed


def collect_code_provenance(repo_dir: Optional[str] = None) -> Dict[str, Any]:
    """Capture the pipeline's own git revision for the running code.

    Returns a dict with the commit SHA, short SHA, branch, a ``dirty`` flag
    (uncommitted tracked changes present), and the commit timestamp. Every
    field is None/False when git metadata is unavailable (e.g. not a checkout
    or git not installed) so the manifest still records the attempt.
    """
    if repo_dir is None:
        repo_dir = os.path.dirname(os.path.abspath(__file__))

    def _git(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(
                ["git", *args],
                cwd=repo_dir,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        if out.returncode != 0:
            return None
        return out.stdout.strip()

    commit = _git("rev-parse", "HEAD")
    provenance: Dict[str, Any] = {
        "commit": commit,
        "commit_short": commit[:12] if commit else None,
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "commit_time": _git("show", "-s", "--format=%cI", "HEAD"),
        "dirty": None,
    }
    if commit is not None:
        status = _git("status", "--porcelain")
        # status is "" when clean; None when the command failed.
        provenance["dirty"] = bool(status) if status is not None else None
    return provenance


def fingerprint_file(
    path: Optional[str], *, max_bytes: Optional[int] = None
) -> Optional[Dict[str, Any]]:
    """Return a content fingerprint for an input file.

    Captures the resolved absolute path, size, modification time, and a SHA-256
    digest of the file contents so a recorded input can be verified rather than
    merely named. Returns None if the path is missing or unreadable. ``max_bytes``,
    when set, hashes only the leading bytes (recorded via ``hashed_bytes``) for
    very large files where a full hash is impractical.
    """
    if not path:
        return None
    try:
        resolved = os.path.abspath(path)
        stat = os.stat(resolved)
    except OSError:
        return None

    digest = hashlib.sha256()
    hashed = 0
    try:
        with open(resolved, "rb") as handle:
            while True:
                if max_bytes is not None and hashed >= max_bytes:
                    break
                chunk_size = 1024 * 1024
                if max_bytes is not None:
                    chunk_size = min(chunk_size, max_bytes - hashed)
                chunk = handle.read(chunk_size)
                if not chunk:
                    break
                digest.update(chunk)
                hashed += len(chunk)
    except OSError:
        return None

    return {
        "path": resolved,
        "size_bytes": stat.st_size,
        "mtime": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "sha256": digest.hexdigest(),
        "hashed_bytes": hashed,
        "partial": max_bytes is not None and hashed < stat.st_size,
    }


def read_parent_provenance(input_path: Optional[str]) -> Optional[Dict[str, Any]]:
    """Find the round_manifest.json beside an input h5ad, if any.

    When a round consumes the output of a previous round (e.g.
    ``round_02/filtered.h5ad``), the producing round's manifest sits in the same
    directory. This returns a compact parent pointer — the manifest path plus its
    recorded round identity — so cross-round lineage is explicit rather than only
    implied by directory naming. Returns None when no sibling manifest exists.
    """
    if not input_path:
        return None
    parent_dir = os.path.dirname(os.path.abspath(input_path))
    manifest_path = os.path.join(parent_dir, "round_manifest.json")
    if not os.path.exists(manifest_path):
        return None
    pointer: Dict[str, Any] = {
        "input_file": os.path.basename(input_path),
        "manifest_path": manifest_path,
        "output_dir": parent_dir,
    }
    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            parent = json.load(handle)
        pointer["created"] = parent.get("created")
        pointer["updated"] = parent.get("updated")
        code_prov = parent.get("code")
        if isinstance(code_prov, Mapping):
            pointer["commit"] = code_prov.get("commit")
        repro = parent.get("reproducibility")
        if isinstance(repro, Mapping):
            pointer["seed"] = repro.get("seed")
    except (json.JSONDecodeError, OSError):
        pass
    return pointer


_SEED_UNSET = object()


def update_round_manifest(
    output_dir: str,
    stage: str,
    payload: Mapping[str, Any],
    *,
    seed: Any = _SEED_UNSET,
) -> str:
    """Merge one stage payload into round_manifest.json.

    Always stamps the manifest with top-level code provenance (git revision of
    the running pipeline). When ``seed`` is provided, records it under a
    top-level ``reproducibility`` block so the round is traceable to the RNG
    seed that produced its stochastic outputs.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "round_manifest.json")
    manifest: Dict[str, Any] = {}
    if os.path.exists(out_path):
        try:
            with open(out_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)
        except json.JSONDecodeError:
            manifest = {}

    manifest.setdefault("pipeline_version", PIPELINE_VERSION)
    manifest.setdefault("created", datetime.now().isoformat())
    manifest["updated"] = datetime.now().isoformat()
    manifest["python"] = {
        "version": platform.python_version(),
        "executable": sys.executable,
    }
    manifest["code"] = collect_code_provenance()
    if seed is not _SEED_UNSET:
        repro = manifest.get("reproducibility")
        if not isinstance(repro, dict):
            repro = {}
        repro["seed"] = seed
        manifest["reproducibility"] = repro
    manifest[stage] = _json_safe(dict(payload))

    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return out_path

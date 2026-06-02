"""Shared V1 pipeline configuration and provenance helpers."""

from __future__ import annotations

import copy
import json
import os
import platform
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


def write_command_args(output_dir: str, stage: str, args: Mapping[str, Any]) -> str:
    """Merge command arguments for a stage into command_args.json."""
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "command_args.json")
    payload: Dict[str, Any] = {}
    if os.path.exists(out_path):
        try:
            with open(out_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except json.JSONDecodeError:
            payload = {}
    payload.setdefault("commands", {})
    payload["commands"][stage] = {
        "timestamp": datetime.now().isoformat(),
        "argv": sys.argv,
        "args": _json_safe(dict(args)),
    }
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


def update_round_manifest(output_dir: str, stage: str, payload: Mapping[str, Any]) -> str:
    """Merge one stage payload into round_manifest.json."""
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
    manifest[stage] = _json_safe(dict(payload))

    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return out_path

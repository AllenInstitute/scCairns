import copy
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from pipeline_config import (  # noqa: E402
    ConfigError,
    DEFAULT_CONFIG,
    deep_merge,
    load_pipeline_config,
    set_if_provided,
    validate_config,
)


def test_deep_merge_preserves_nested_defaults():
    merged = deep_merge(
        {"data": {"batch_key": "a", "counts_layer": "counts"}},
        {"data": {"batch_key": "b"}},
    )

    assert merged["data"]["batch_key"] == "b"
    assert merged["data"]["counts_layer"] == "counts"


def test_set_if_provided_ignores_none():
    config = {"data": {"batch_key": "data_origin"}}

    set_if_provided(config, ["data", "batch_key"], None)
    assert config["data"]["batch_key"] == "data_origin"

    set_if_provided(config, ["data", "batch_key"], "sample")
    assert config["data"]["batch_key"] == "sample"


def test_validate_requires_input():
    config = copy.deepcopy(DEFAULT_CONFIG)

    with pytest.raises(ConfigError, match="Missing input h5ad"):
        validate_config(config, mode="integration")


def test_scanvi_requires_labels_key():
    config = copy.deepcopy(DEFAULT_CONFIG)
    config["data"]["input_h5ad"] = "input.h5ad"
    config["annotation"]["enabled"] = True
    config["annotation"]["labels_key"] = None

    with pytest.raises(ConfigError, match="annotation.labels_key"):
        validate_config(config, mode="integration")


def test_load_yaml_config_resolves_paths(tmp_path):
    pytest.importorskip("yaml")
    cfg_path = tmp_path / "pipeline.yml"
    cfg_path.write_text(
        "\n".join(
            [
                "pipeline_version: 1",
                "data:",
                "  input_h5ad: data/input.h5ad",
                "  output_dir: rounds/round_01",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    config = load_pipeline_config(str(cfg_path))

    assert config["data"]["input_h5ad"] == str(tmp_path / "data/input.h5ad")
    assert config["data"]["output_dir"] == str(tmp_path / "rounds/round_01")

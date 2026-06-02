import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))


anndata = pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from inspect_integration import apply_decisions  # noqa: E402


def _adata():
    adata = anndata.AnnData(np.ones((4, 3)))
    adata.obs_names = [f"cell{i}" for i in range(4)]
    adata.var_names = [f"gene{i}" for i in range(3)]
    adata.obs["leiden"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obs["n_genes_by_counts"] = [600, 700, 300, 800]
    return adata


def test_apply_decisions_filters_valid_query(tmp_path):
    filtered = apply_decisions(
        _adata(),
        {"keep_cells": [{"query": "n_genes_by_counts > 500"}]},
        output_dir=str(tmp_path),
    )

    assert filtered.n_obs == 3
    assert (tmp_path / "filtered.h5ad").exists()
    assert (tmp_path / "round_manifest.json").exists()


def test_apply_decisions_raises_on_invalid_query(tmp_path):
    with pytest.raises(ValueError, match="Keep query failed"):
        apply_decisions(
            _adata(),
            {"keep_cells": [{"query": "missing_column > 0"}]},
            output_dir=str(tmp_path),
        )

"""Make the ``sccairns`` package importable during tests without installing it.

The package lives under ``code/`` (so Code Ocean reproducible runs, which ship only
``code/``, can import it). Putting ``<repo>/code`` on ``sys.path`` lets ``pytest`` run
straight from a checkout — ``import sccairns`` resolves whether or not the package has
been ``pip install``-ed. When installed via ``pip install -e .`` this is a harmless no-op.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "code"))

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def tiny_adata():
    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(np.ones((4, 3)))
    adata.obs_names = [f"cell{i}" for i in range(4)]
    adata.var_names = [f"gene{i}" for i in range(3)]
    adata.obs["leiden"] = pd.Categorical(["0", "0", "1", "1"])
    adata.obs["data_origin"] = pd.Categorical(["a", "a", "b", "b"])
    adata.obs["n_genes_by_counts"] = [600, 700, 300, 800]
    return adata

#!/usr/bin/env python3
"""Backward-compatible shim → :mod:`sccairns.config`.

The shared config/provenance library moved into the ``sccairns`` package. This
re-export shim keeps ``from pipeline_config import ...`` working for any external
notebook or script that added ``code/`` to ``sys.path``. Puts this directory
(``code/``, where the ``sccairns`` package lives) on ``sys.path`` so the package resolves.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sccairns.config import *  # noqa: E402,F401,F403

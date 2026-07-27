#!/usr/bin/env python3
"""Backward-compatible shim → :mod:`sccairns.summarize`.

Kept so the Code Ocean capsule keeps invoking ``./summarize_rounds.py`` from
CWD=``code/`` with no change. Puts the repo root on ``sys.path`` so the package
imports whether or not it has been ``pip install``-ed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sccairns.summarize import main  # noqa: E402

if __name__ == "__main__":
    main()

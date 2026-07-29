#!/usr/bin/env python3
"""Backward-compatible shim → :mod:`sccairns.integrate`.

The implementation moved into the installable ``sccairns`` package. This shim is
kept so the Code Ocean capsule (``code/run`` and the ``integrate_sns_scvi.py``
alias) keeps invoking ``./integrate_scvi.py`` from CWD=``code/`` with no capsule
change. It puts this directory (``code/``, where the ``sccairns`` package lives) on
``sys.path`` so the package imports whether or not it has been ``pip install``-ed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sccairns.integrate import main  # noqa: E402

if __name__ == "__main__":
    main()

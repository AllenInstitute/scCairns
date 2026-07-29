#!/usr/bin/env python3
"""Backward-compatible shim → :mod:`sccairns.contamination`.

Kept so the Code Ocean capsule keeps invoking ``./flag_contamination.py`` from
CWD=``code/`` with no change, and so ``from flag_contamination import ...`` in any
external notebook still resolves. Puts this directory (``code/``, where the ``sccairns``
package lives) on ``sys.path`` so the package imports with or without ``pip install``.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sccairns.contamination import *  # noqa: E402,F401,F403
from sccairns.contamination import (  # noqa: E402,F401
    DEFAULT_CONTAM_PANELS,
    flag_contamination,
    main,
)

if __name__ == "__main__":
    main()

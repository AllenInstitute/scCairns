"""Make the ``sccairns`` package importable during tests without installing it.

The package lives under ``code/`` (so Code Ocean reproducible runs, which ship only
``code/``, can import it). Putting ``<repo>/code`` on ``sys.path`` lets ``pytest`` run
straight from a checkout — ``import sccairns`` resolves whether or not the package has
been ``pip install``-ed. When installed via ``pip install -e .`` this is a harmless no-op.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "code"))

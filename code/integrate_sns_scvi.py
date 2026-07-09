#!/usr/bin/env python3
"""Deprecated alias for ``integrate_scvi.py``.

This script was renamed ``integrate_sns_scvi.py`` -> ``integrate_scvi.py`` (the code
is dataset-agnostic; the SNS framing lives in config/marker defaults, not the engine).
This thin wrapper is kept so existing invocations — including commented Code Ocean
``run`` templates and any external scripts — keep working unchanged. It forwards all
arguments to the renamed module. Prefer calling ``integrate_scvi.py`` directly.
"""

import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    print(
        "[DEPRECATION] integrate_sns_scvi.py has been renamed to integrate_scvi.py; "
        "please update your command. Forwarding to integrate_scvi.py ...",
        file=sys.stderr,
    )
    _target = Path(__file__).with_name("integrate_scvi.py")
    # Run the renamed module as __main__ with the current argv, so argparse sees
    # exactly the flags the caller passed.
    runpy.run_path(str(_target), run_name="__main__")

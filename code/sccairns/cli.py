"""``cairns`` command-line entry point: a thin dispatcher over the stage modules.

Each subcommand forwards to a stage module's ``main()``; the module keeps its own
argparse, so ``cairns integrate --help`` shows exactly that stage's options. The
underlying scripts are unchanged — this only routes ``cairns <cmd> ...`` to them.
"""
from __future__ import annotations

import importlib
import sys
from typing import List, Optional

# subcommand -> module providing main()
COMMANDS = {
    "integrate": "sccairns.integrate",
    "inspect": "sccairns.inspect",
    "summarize": "sccairns.summarize",
    "record-filter": "sccairns.record_filter",
    "flag-contamination": "sccairns.contamination",
}

_USAGE = (
    "usage: cairns <command> [options]\n\n"
    "commands:\n"
    "  integrate           QC, HVG, scVI/scANVI training, UMAP/Leiden, benchmarking\n"
    "  inspect             inspection report, cluster QC, auto-flags, decision filtering\n"
    "  summarize           cross-round summary: table, ledger, lineage, integrity\n"
    "  record-filter       record provenance for an externally filtered h5ad pair\n"
    "  flag-contamination  marker-based per-cell contamination flagging\n\n"
    "Run 'cairns <command> --help' for command-specific options.\n"
    "Run 'cairns --version' to see which scCairns is on PATH.\n"
)


def _version_line() -> str:
    """Identify the running scCairns the same way a round manifest does.

    Answers "which one is installed?" for a pinned deployment, where the version
    — not a commit — is the identifier stamped into every round.
    """
    from .config import collect_code_provenance  # local: keeps --help cheap

    prov = collect_code_provenance()
    version = prov.get("version") or "unknown version"
    if prov.get("source") == "checkout":
        detail = f"checkout {prov.get('commit_short')}"
        if prov.get("repo_root"):
            detail += f" in {prov['repo_root']}"
        if prov.get("dirty"):
            detail += ", dirty"
    else:
        detail = prov.get("source") or "unknown"
    return f"cairns {version} ({detail})\n"


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        sys.stdout.write(_USAGE)
        return 0
    if argv[0] in ("-V", "--version", "version"):
        sys.stdout.write(_version_line())
        return 0
    command = argv[0]
    module_name = COMMANDS.get(command)
    if module_name is None:
        sys.stderr.write(f"cairns: unknown command '{command}'\n\n{_USAGE}")
        return 2
    module = importlib.import_module(module_name)
    # Reshape argv so the stage's argparse sees prog='cairns <command>' + its args.
    sys.argv = [f"cairns {command}"] + argv[1:]
    result = module.main()
    return result if isinstance(result, int) else 0


if __name__ == "__main__":
    raise SystemExit(main())

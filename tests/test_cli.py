"""Dispatcher-level behavior of the ``cairns`` command."""

import importlib

import pytest

from sccairns import cli


def test_help_lists_every_command(capsys):
    assert cli.main([]) == 0
    usage = capsys.readouterr().out
    for command in cli.COMMANDS:
        assert command in usage


@pytest.mark.parametrize("command", sorted(cli.COMMANDS))
def test_every_command_dispatches_to_a_module_main(command):
    """The dispatcher calls ``module.main()``, so every target must export one.

    Listing a command in COMMANDS and in the usage text is not enough: a module
    whose argparse lives only under ``if __name__ == "__main__"`` imports fine
    and prints in --help, then raises AttributeError when actually invoked.
    """
    module = importlib.import_module(cli.COMMANDS[command])
    assert callable(getattr(module, "main", None)), (
        f"{cli.COMMANDS[command]} has no module-level main(); "
        f"'cairns {command}' would fail at dispatch"
    )


def test_unknown_command_is_an_error(capsys):
    assert cli.main(["nope"]) == 2
    assert "unknown command" in capsys.readouterr().err


def test_version_reports_the_running_scCairns(capsys, monkeypatch):
    """`cairns --version` is how you confirm which pinned build is deployed."""
    from sccairns import config as cfg

    monkeypatch.setattr(cfg, "collect_code_provenance", lambda: {
        "version": "0.1.0", "source": "installed", "commit": None,
        "commit_short": None, "repo_root": None, "dirty": None,
    })
    assert cli.main(["--version"]) == 0
    assert capsys.readouterr().out.strip() == "cairns 0.1.0 (installed)"


def test_version_from_a_checkout_names_the_commit(capsys, monkeypatch):
    from sccairns import config as cfg

    monkeypatch.setattr(cfg, "collect_code_provenance", lambda: {
        "version": "0.1.0", "source": "checkout", "commit": "abc123def4567",
        "commit_short": "abc123def456", "repo_root": "scCairns", "dirty": True,
    })
    assert cli.main(["-V"]) == 0
    out = capsys.readouterr().out.strip()
    assert out == "cairns 0.1.0 (checkout abc123def456 in scCairns, dirty)"

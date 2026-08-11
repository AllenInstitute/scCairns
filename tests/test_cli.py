"""Dispatcher-level behavior of the ``cairns`` command."""

from sccairns import cli


def test_help_lists_every_command(capsys):
    assert cli.main([]) == 0
    usage = capsys.readouterr().out
    for command in cli.COMMANDS:
        assert command in usage


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

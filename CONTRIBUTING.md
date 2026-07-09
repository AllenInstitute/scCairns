# Contributing

Notes for developing on `scvi_integration_loops`.

## Setup

Use the local conda/mamba environment (see
[docs/getting-started.md](docs/getting-started.md)):

```bash
mamba create -n scvi-loops -c conda-forge python=3.10 pip scikit-misc -y
mamba activate scvi-loops
python -m pip install -r requirements.txt
```

## Repository layout

| Path | What's there |
|---|---|
| `code/` | The analysis scripts and shared helpers (`pipeline_config.py`, `flag_contamination.py`). |
| `code/run`, `code/run_capsule.py` | Code Ocean capsule entry points — **the GitHub → Code Ocean bridge; keep them working** ([docs/codeocean.md](docs/codeocean.md)). |
| `examples/` | Config templates + the tutorial data generator. |
| `environment/` | Dockerfile and post-install for the capsule image. |
| `docs/` | User documentation. |
| `tests/` | Pytest suite for config, provenance, decisions, harmony, and summary logic. |

## Running the checks

```bash
python -m compileall -q code tests   # everything imports/compiles
python -m pytest -q                   # unit tests pass
```

The tests cover the traceability-critical paths (config merge/validation, provenance
fingerprinting, decision application, sweep-variant selection, and round
summarization). Add or update tests when you touch those.

## Sanity-check against the synthetic data

Before proposing a change to the integrate/inspect/filter flow, run the
[tutorial](docs/tutorial.md) end-to-end on the generated example data — it exercises
the full loop (auto-flagging, contamination flagging, filtering, re-integration,
summary) in a couple of minutes and produces artifacts you can diff.

## Conventions

- **Config is the source of truth.** New behavior should be a config field with a
  validated default in `DEFAULT_CONFIG` (`code/pipeline_config.py`), not a hardcoded
  constant. Bump `PIPELINE_VERSION` and document schema changes if you break
  compatibility.
- **Keep the provenance trail intact.** Anything that changes what a round produces
  should be reflected in `round_manifest.json` so runs stay reproducible.
- **Match the surrounding style** — the scripts use argparse groups, explicit
  `--flag` CLI overrides layered on the YAML config, and readable stdout narration.
- Update `docs/` and `CHANGELOG.md` alongside user-facing changes.

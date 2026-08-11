# Running on Code Ocean

This repo doubles as a [Code Ocean](https://codeocean.com/) capsule. You can develop
and run locally exactly as in the [tutorial](tutorial.md); this page covers the
capsule-specific pieces so the absolute paths and orchestration make sense.

> Local users can ignore everything here. The `/root/capsule/...` and `/data/...`
> paths below only exist inside the capsule container.

## The capsule files

| File | Role | Edit? |
|---|---|---|
| `code/run` | The master bash script Code Ocean executes on "Reproducible Run". Contains the ordered, commented commands for each round. | Yes — this is where you sequence a run. |
| `code/run_capsule.py` | Minimal Python entry point. | Rarely. |
| `environment/Dockerfile` | The container image (CUDA `scvi-tools`, pinned libraries). | Only for dependency changes. |
| `.codeocean/datasets.json` | Attached datasets and their mount points under `/data/`. | Via the Code Ocean UI. |
| `.codeocean/resources.json` | Compute class. | Via the Code Ocean UI. |

These files are the GitHub → Code Ocean bridge: push to GitHub, pull into the capsule,
click Run. Keep their behavior intact.

## How `cairns` gets into the capsule

Two arrangements, depending on whether the capsule *is* this repo.

**Vendored (this repo as a capsule).** The `sccairns` package ships inside `code/` — the
only folder a reproducible run materializes — so `code/run` puts it on `PYTHONPATH` and
runs it as a module. No install, no network:

```bash
HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export PYTHONPATH="${HERE}:${PYTHONPATH:-}"
cairns() { python -m sccairns.cli "$@"; }
```

**Pinned dependency (a project capsule that is not this repo).** The capsule holds only
its own `pipeline.yml`, `run`, and `.codeocean/datasets.json`, and scCairns comes from
the image. Put the install in the capsule's **`postInstall`**:

```bash
pip install --no-cache-dir "git+https://github.com/AllenInstitute/scCairns.git@v0.1.0"
```

Then `run` needs no preamble at all — pip puts a real `cairns` on PATH.

Two things to get right. **The install must happen at image build, never in `run`**: a
reproducible run is offline, so a `pip install` there fails (it works on the Cloud
Workstation, which is why this is easy to miss). And **put it in `postInstall`, not
`environment/Dockerfile`** — Code Ocean generates the Dockerfile from the UI package
list (note its `# hash:sha256:` header) and will overwrite hand edits. Authenticated
access to the private repo is already wired up: the Dockerfile carries `ARG GIT_ASKPASS`
and `COPY git-ask-pass /`.

Pin an exact tag. `@main` means rebuilding the image silently changes the tool, and
since an installed distribution has no git metadata, the pinned version is what each
`round_manifest.json` records as its `code.version`. Confirm what a built image actually
got with `cairns --version`.

## Path conventions inside the capsule

- **Inputs** are mounted read-only under `/data/<mount-name>/` (defined in
  `.codeocean/datasets.json`).
- **Outputs** go under `/root/capsule/results/` — which Code Ocean **overwrites each
  run**. So finished rounds are copied out to durable storage (e.g. an S3 mount)
  between runs, and later referenced by their archived paths.
- `code/run` uses `--config ./pipeline_scvi.yml` (relative to `code/`) and absolute
  `/data` / `/root/capsule/results` paths for inputs and outputs.

## A round in the capsule

The pattern mirrors the local loop, with capsule paths:

```bash
# integrate a mounted dataset
python -u ./integrate_scvi.py --config ./pipeline_scvi.yml \
  --input /data/<mount>/combined.h5ad \
  --output-dir /root/capsule/results/rounds/round_01/

# inspect (report + auto flags, all sweep architectures)
python -u ./inspect_integration.py --config ./pipeline_scvi.yml \
  --input /root/capsule/results/rounds/round_01/integrated.h5ad \
  --output-dir /root/capsule/results/rounds/round_01/ \
  --all-sweep-architectures

# after uploading a decisions.yaml as a new dataset: filter into the next round
python -u ./inspect_integration.py --config ./pipeline_scvi.yml \
  --input /data/<round1-mount>/round_01/integrated.h5ad \
  --cluster-key leiden_<variant> \
  --decisions /data/<decisions-mount>/decisions.yaml \
  --output-dir /root/capsule/results/rounds/round_02
```

Because decisions are reviewed by a human between rounds, the typical capsule rhythm
is: run one round, copy its outputs to S3, review locally, upload a `decisions.yaml`
as a new attached dataset, then run the next round. `code/run` keeps every round's
commands (commented) as a record of the actual analysis history.

## Summarizing archived rounds

Since `/results` is overwritten each run, rounds live at arbitrary archived paths.
`summarize_rounds.py` handles this — pass the rounds explicitly and their **supplied
order is taken as the lineage** when the original run-time paths no longer resolve:

```bash
python -u ./summarize_rounds.py \
  --rounds /data/<r1>/round_01 /data/<r2>/round_02 /data/<r3>/rounds/round_03 \
  --output-dir /root/capsule/results/summary --verify \
  --terminal-run-architecture leiden_harmony
```

`--verify` still hashes each parent's `integrated.h5ad` as it sits on disk now and
checks it against the child's recorded input fingerprint, so lineage continuity is
validated even across archiving. See the
[summarize section of the README](../README.md#summarizing-a-completed-set-of-rounds).

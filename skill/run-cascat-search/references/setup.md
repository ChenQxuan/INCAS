# Environment setup and installation

Use this reference when installation, environment discovery, or installation repair is required. The user's request to install the project authorizes its ordinary local setup steps; respect the host's filesystem and network permissions. Do not turn environment setup into a request to alter unrelated system configuration.

## Choose an interpreter

The checked-out `pyproject.toml` is authoritative. CASCAT 0.4.0 requires Python 3.11 or 3.12. Inspect available interpreters and any project `.venv`, and preserve a compatible environment the user explicitly selected. A Python outside the allowed range is not fixed by removing the project's version constraint.

Prefer a project-local `.venv` for a new installation. Do not install into the system Python or modify a shared conda environment implicitly. On macOS/Linux, use the installed `python3.12` or `python3.11`. For Windows, WSL2 is the project's recommended path; native Windows has not been validated. If no compatible interpreter exists, use the user's available environment/package manager to provision an isolated one when authorized, or explain the exact missing interpreter. Do not require conda if venv is available.

## Install the actual checkout

Commands below assume the working directory is the verified project root and Python 3.12 exists. Adapt the interpreter to the machine; execute the steps, not just print them. Check an existing `.venv/pyvenv.cfg` and its interpreter before deciding to create one. Never replace a non-environment directory or an existing environment merely to reuse the name.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install .
.venv/bin/python -m pip check
.venv/bin/python -m cascat molecules
.venv/bin/python -m cascat run --help
```

Use the absolute environment Python after setup. On native Windows, the executable is `.venv/Scripts/python.exe`. A successful command in a different environment does not validate this installation.

Dependencies are declared in the project: NumPy <2, ASE, ACAT 2.0.3, AutoAdsorbate 0.2.5, and mace-torch 0.3.16, with their transitive dependencies. Do not copy the developer's local conda path or blindly pin the historical torch version from `VALIDATION.md`. Check the installed `cascat` CLI has `--ads`, `--n`, `--parents`, `--screen-keep`, `--relax-keep`, and `--final-keep`. If it exposes the old `ads-search`/`--keep` interface, reinstall from the current checkout into the selected environment.

Dependency installation needs network access unless wheels are already available. MACE uses the `mh-1` model with head `oc20_usemppbe`; its first actual calculation downloads weights if not cached. CLI help and molecule listing do not validate model loading. Describe these as separate milestones.

When an install fails, use the actual error to distinguish interpreter incompatibility, network/permissions, dependency resolution, and missing build prerequisites. Obtain required host tool permissions rather than changing package sources, disabling TLS, using sudo pip, or silently changing dependency constraints. After a successful install, use `pip check` in the isolated environment and rerun the CLI checks. Do not repeatedly retry an unchanged network failure.

## Check or demonstrate use

For “install only”, successful dependency and CLI checks complete the task; report that MACE optimization has not been exercised. For “install and run my structure”, proceed with that structure and its requested settings. Do not insert a full catalogue test suite or a separate expensive demo as an installation requirement.

If the user requests a demo, the repository includes `examples/POSCAR`. A small actual MACE check is:

```bash
.venv/bin/python -m cascat run examples/POSCAR --ads H --n 1 \
  --parents 1 --screen-keep 1 --relax-keep 1 --final-keep 1 \
  --surface 111 --output results/demo-H
```

Choose a fresh output directory if that name already belongs to a different configuration. These small limits validate execution, not exhaustive search quality. `--generate-only` tests candidate generation without downloading/loading MACE; explicitly report that limited scope if used.

## Distribution and skill activation

The repository includes this skill in the visible `skill/run-cascat-search/` directory. Users should open the source checkout and explicitly ask their assistant to read `skill/run-cascat-search/SKILL.md` and follow it. Do not assume this directory is automatically discovered as an installed skill. No home-directory installation is required when the assistant reads the skill directly.

`pip install .` installs the search command, not a global Codex skill. Source distributions retain the skill through `MANIFEST.in`; the Python wheel alone is not a skill installer. Do not copy to a user's global skill directory unless they request global availability. Never invent a public repository URL before the owner publishes one.

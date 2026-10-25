---
name: run-incas-search
description: Install and operate INCAS, the ACAT/MACE incremental coverage-dependent catalytic metal-surface adsorption configuration workflow, from a user's POSCAR or CONTCAR. Use for local setup, choosing search parameters, running or resuming searches, and checking BEST structures and retained results; not for submitting VASP or changing the search algorithm.
---

# Run INCAS search

Carry the user's requested task through actual execution. For “install and run”, install any missing project dependencies in an appropriate isolated environment, run the requested search, and inspect the outputs. For installation only, stop after verifying the installation; for an explanation or status request, inspect and answer without launching a new calculation. Match the user's language.

## Locate the project and environment

Find the checkout containing `pyproject.toml` with project name `incas` and `src/incas/cli.py`. When this skill is inside the checkout, walk up from this skill directory to find that root; it may differ from the current working directory. If the skill was installed separately, inspect the user's named checkout or working directory. Do not assume a personal absolute path or invent a GitHub URL. If no checkout is available, obtain the source location from the user while checking available Python interpreters.

Read the checkout's `README.md` and `pyproject.toml`. For environment creation, installation, or installation errors, read [references/setup.md](references/setup.md) and execute the relevant steps. Use the same verified environment's Python for all subsequent commands, preferably its absolute path with `-m incas`; activation in one tool call may not persist in another.

Verify `python -m incas run --help` and `python -m incas molecules`. Use the installed help as the interface authority if a future version differs from this skill. The current interface is 0.5.0; obsolete `ads-search`, `--keep`, `--count`, and `--adsorbate` commands must not be used.

## Resolve the requested calculation

Use information already supplied before asking questions. Required scientific inputs are a readable clean slab file, the adsorption species, and the target molecular count. Ask for genuinely missing inputs together; do not silently substitute the example Cu slab or an arbitrary molecule for the user's structure. A runnable example is appropriate when the user requests a demo.

Keep these four controls independent:

| User meaning | Flag | Values | Default if omitted |
| --- | --- | --- | --- |
| N=n starts from x final structures at N=n−1; the same x for every coverage | `--parents` | Integer ≥1 or `all` | `1` |
| Initial unique candidates, ranked by single-point energy, entering coarse relaxation | `--screen-keep` | Integer ≥1, percentage, or `all` | `80%` |
| Coarse-relaxed, ranked, graph-unique eligible candidates entering full relaxation | `--relax-keep` | Integer ≥1, percentage, or `all` | `50%` |
| Fully relaxed, ranked, converged unique results to save, including BEST | `--final-keep` | Integer ≥1, percentage, or `all` | `all` |

Percentages use each stage's own eligible unique pool and round up. For example, 11 at 80% yields 9; 9 at 50% yields 5. Counts are upper limits, not guarantees. `all` removes the corresponding selection limit, but not deduplication or validity/convergence checks. Do not apply an additional energy window, hidden limit, or the old 32/128/256 bounds.

`--parents` selects only from the results retained by `--final-keep`; if only 3 were saved, requesting 5 parents uses those 3. State this relationship if the user's settings imply a shortage; do not silently enlarge the final result set. `--parents all` means all retained final results, not all historical candidates. Never collapse these controls into a common value unless the user explicitly supplies that choice.

`--ads CO --n 4` searches CO counts 1 → 2 → 3 → 4. N=1 starts from the clean slab; each later coverage adds one molecule to each selected parent. The three filtering stages occur inside every coverage. A request for N=4 does not bypass the required lower coverages.

Validate the actual POSCAR/CONTCAR through ASE and the project's checks. The current workflow requires a clean metal slab in the xy plane, a z-oriented third cell vector, and at least 6 Å total vacuum. It rejects substrate elements that ACAT treats as adsorbates (H/C/N/O/F/S/X). Do not silently delete atoms, rotate the slab, expand the cell, or remove pre-existing adsorbates to make an unsupported input pass. Explain the concrete issue and obtain the missing scientific choice when a transformation is needed.

Default repeat is `1 1 1`; use `--repeat 2 2 1` only when requested. Default `--surface auto` uses CustomSurface/AutoAdsorbate. Set a numeric FCC label only for a known corresponding FCC surface, not merely because a filename contains 111. Check the actual molecule catalogue; do not map unsupported intermediates to a different molecule or treat the library as an exhaustive CO2RR reaction network.

If the user omits optional settings, state the defaults briefly and proceed within the authorized task. CPU execution defaults to one model worker and two threads; retain that unless resources or user instructions justify changing it. Do not change scientific selection parameters to make a slow or failing search finish.

## Execute and follow the search

For example, using the verified environment's Python:

```bash
python -m incas run POSCAR --ads CO --n 4 --parents 2 \
  --screen-keep 80% --relax-keep 50% --final-keep 10 --output results/CO
```

For fixed stage counts, use `--screen-keep 1024 --relax-keep 256`. Resolve slab and output paths explicitly so changing the working directory cannot change the calculation's inputs. Present the resolved command/settings as an update, then execute when the user requested a run; do not stop at providing a command to copy.

Use a new descriptive output directory for new scientific settings. The root `search.json` contains the input hash, selection settings, versions, and engine version. Reuse an existing output only when it matches; never edit its hash or metadata to bypass a mismatch. If the user asks to resume, reconstruct settings from `search.json` and locate the original slab. The saved `slab.vasp` is already repeated and may have a different byte hash; do not pass it back with the original repeat as a substitute for the source input.

`--generate-only` generates the next unfinished coverage's complete unique candidate pool without MACE. It does not apply energy-ranked percentage/count selection and does not complete all coverages. Use it for a requested generation-only check or when diagnosing site generation; do not substitute it for an authorized full calculation.

Keep the running process attached to the available execution/session tools, follow progress, and report meaningful updates until it completes or fails. Do not launch a second process into the same output directory. If saving a log through a pipeline, preserve the program exit status. A `.running.lock` indicates an active or interrupted run; verify the process has stopped before considering removal. Do not automatically delete a lock on another machine based only on a locally absent PID.

If a run fails, inspect the actual traceback and last stage. Missing dependencies/network access are installation issues; no candidates, invalid geometry, or no converged unique minimum are calculation failures. Do not mark them complete, loosen convergence, switch the model/head, or alter distances/counts without an authorized change. Re-running identical settings resumes completed coverages and recalculates unfinished ones.

## Inspect and hand off

Use the bundled read-only checker with the selected environment:

```bash
python /path/to/run-incas-search/scripts/inspect_results.py /path/to/results/CO
```

Exit 0 means all expected coverages have complete metadata and consistent result files; exit 1 means incomplete, missing, or inconsistent outputs. The checker does not run MACE or repair data. For an ongoing status request, report incomplete stages as progress, not automatically as a failed calculation.

Return the output directory, completed coverages, actual coarse/full/final/parent counts, and links to the user's relevant `BEST.vasp`, `final_results.xyz`, `growth_parents.xyz`, and `ranking.csv`. Final results and parents are distinct: `final_results` contains the retained result set; `growth_parents` contains the subset for N+1. `BEST.vasp` is the lowest retained geometry, not necessarily the only saved result. Full traces may retain unselected structures for audit.

Report verification boundaries candidly: template construction and successful optimization do not prove a global minimum. Energies rank structures within the same species and coverage; do not rank different N by raw total energy or describe MACE results as constant-potential/solvated DFT. Do not submit VASP or publish results merely because this search finished.

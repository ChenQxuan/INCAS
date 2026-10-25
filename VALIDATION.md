# Validation — 2026-10-25

## INCAS 0.5.0 naming

Renamed the project from CASCAT to **INCAS — Incremental Coverage Adsorbate Search**; the CASCAT name now belongs to a separate, extended project. The distribution, import module, CLI, and skill are now `incas`, `incas`, `incas`, and `run-incas-search`, respectively. The search algorithm, options, defaults, and output layout are unchanged from CASCAT 0.4.0, and the result checker accepts both INCAS 0.5.x and CASCAT 0.4.x output trees. Sections below this one describe the former names as they were at the time.

## CASCAT 0.4.0 naming (2026-09-06)

Renamed the project to **CASCAT — Coverage-Aware Search for Catalytic Adsorbate Topologies**. The distribution, import module, CLI, and skill are now `cascat`, `cascat`, `cascat`, and `run-cascat-search`, respectively. The former prototype names remain only in historical validation notes and explicit obsolete-interface warnings.

Built `cascat-0.4.0-py3-none-any.whl`, validated the renamed skill, and verified `python -m cascat molecules` and `python -m cascat run --help`. The complete automated suite was rerun after the rename.

## Repository skill

Added `skill/run-cascat-search/` using the skill-creator structure: execution instructions, installation reference, UI metadata, and a read-only result checker. Skill frontmatter validation passed. Five checker tests passed, covering reordered VASP atoms, a damaged BEST geometry, incomplete status, metadata count mismatch, and absent output paths without writes.

The checker also passed on the actual 0.3.0 two-coverage H results recorded below: N=1 retained 2 final results / 1 parent and N=2 retained 5 / 1. Built a source distribution and verified that it includes SKILL.md, agents/openai.yaml, references/setup.md, and scripts/inspect_results.py through MANIFEST.in. Environment installation instructions reuse the existing project's validated CLI; a fresh online dependency installation and an independent agent execution of the new skill were not performed in this skill-writing task.

## 0.3.0: independent stage and parent retention

The current CLI uses `--parents` (integer or all), independently of `--screen-keep`, `--relax-keep`, and `--final-keep` (integer, percentage, or all). Percentages use each stage's available unique eligible pool and round up. Final results are saved separately from growth parents. All sections below 0.3.0 describe historical interfaces.

All 7 automated tests passed in 46.0 seconds, including independent CLI values, invalid percentages, exact percentage rounding, and parent selection from a separately limited final pool.

Actual MACE integration command:

```bash
cascat run examples/POSCAR --ads H --n 2 --parents 1 \
  --screen-keep 80% --relax-keep 50% --final-keep all \
  --surface 111 --output results/independent-stages
```

| Coverage | Initial unique | Coarse optimized | Coarse graph representatives | Fully optimized | Saved final results | Growth parents |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N=1 | 4 | 4 | 4 | 2 | 2 | 1 |
| N=2 | 11 | 9 | 9 | 5 | 5 | 1 |

Both coverages completed successfully. In N=2, 80% of 11 correctly selected 9 candidates, and 50% of 9 selected 5. All five eligible final results were retained while only one became a growth parent. BEST remains the lowest-energy retained result. This small Cu(111)/H test verifies the independent settings, not global search completeness.

## 0.2.0: unified retention setting

`--keep` now accepts one positive integer (minimum 1) or `all`. The same limit controls single-point-ranked candidates entering short optimization, graph representatives entering full optimization, and final eligible growth parents. Short-optimization graph deduplication retains one lowest-energy representative per graph. There is no additional energy window or candidate-count floor/cap.

The commands in the 0.1.0 sections below are historical and use options removed in 0.2.0. Use the README commands for the current interface.

0.2.0 integration commands (historical flags):

```bash
ads-search run examples/POSCAR --ads H --n 2 --keep 1 \
  --surface 111 --output results/unified-one
ads-search run examples/POSCAR --ads H --n 1 --keep all \
  --surface 111 --output results/unified-all
```

All 6 automated tests passed (69.7 seconds). The added checks cover `--keep 1`, `all`, invalid values, rejection of the former list syntax, shortages, and retention of high-energy candidates with `all`.

Both actual MACE runs completed:

| Setting | Coverage | Initial unique candidates | Short optimizations | Full optimizations | Saved parents |
| --- | ---: | ---: | ---: | ---: | ---: |
| keep=1 | 1 | 4 | 1 | 1 | 1 |
| keep=1 | 2 | 11 | 1 | 1 | 1 |
| keep=all | 1 | 4 | 4 | 4 | 4 |

All final structures in these runs were converged, valid and unique. Repeating the keep=1 command correctly reused both completed coverages. The 0.2.0 wheel also built successfully. These are small Cu(111)/H software integration checks, not proof of global search completeness.

## Environment

Verified on the existing macOS `surface-search` Python 3.12 environment:

| Package | Version |
| --- | --- |
| numpy | 1.26.4 |
| ase | 3.29.0 |
| acat | 2.0.3 |
| autoadsorbate | 0.2.5 |
| mace-torch | 0.3.16 |
| torch | 2.13.0 |

Model: cached MACE-MH-1, `oc20_usemppbe` head, CPU. Float32 energy screening and float64 geometry optimization. No new model download was required.

## 0.1.0 automated checks

`python -m unittest discover -s tests -v`: 4 tests passed (45.3 seconds).

- All 28 catalogue entries generated at least one candidate recognized as an intact molecule on Cu(111).
- CO2, OOH and H2O2 grew from one to two molecules with valid molecular counts.
- Parent selection supports more than two structures, handles shortages, and retains optional legacy energy-gap behavior in the engine.
- Repeating candidate-only generation works; changed slab input is rejected when reusing an output directory.

These checks validate template construction and software behavior, not adsorption stability for every molecule and material.

## 0.1.0 real MACE searches

H on the bundled 36-atom Cu(111) slab:

```bash
ads-search run examples/POSCAR --adsorbate H --count 2 --keep 3 2 \
  --surface 111 --output results/smoke-H \
  --min-relax 2 --max-pre-relax 3 --max-relax 3 --workers 1
```

| Round | Generated | Fully optimized | Eligible unique minima | Saved parents |
| --- | ---: | ---: | ---: | ---: |
| H N=1 | 4 | 3 | 3 | 3 |
| H N=2 | 25 | 3 | 3 | 2 |

Both metadata files report `complete`. Verified the saved BEST structures have the expected atom counts, unchanged slab coordinates and unchanged cells. Repeating the command reused both completed rounds.

CO with the default automatic surface identification:

```bash
ads-search run examples/POSCAR --adsorbate CO --count 1 --keep 2 \
  --output results/smoke-CO \
  --min-relax 2 --max-pre-relax 3 --max-relax 3 --workers 1
```

Generated 3 candidates; all 3 fully optimized structures were converged, valid and unique. Saved 2 parents and BEST.vasp. This checks the CustomSurface/AutoAdsorbate path and a molecular adsorbate through the actual MACE stages.

Small screening caps were used for these integration checks; they do not establish search completeness.

## 0.1.0 packaging

Built `local_adsorption_search-0.1.0-py3-none-any.whl` using `pip wheel --no-deps --no-build-isolation`. Installed the wheel into a separate temporary target and successfully ran its `ads-search molecules` and `ads-search run --help` from outside the source tree.

This validates packaging and the installed entry point using existing dependencies. A clean dependency installation and Windows/Linux execution have not been tested. Complete optimizations for the other 26 templates and larger coverages have not been performed.

## Origin

`engine.py` reuses the scientific functions from the existing `5a_search_H_configurations.py`. The portable wrapper replaces the fixed CP directory input with an explicit VASP file, separates output from input, and records a checked run configuration. Engine changes extend parent selection beyond two and recognize multidentate molecules using their atom sets and ACAT naming conventions. Historical files outside this project were not changed.

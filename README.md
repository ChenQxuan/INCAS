# INCAS

**INCAS — Incremental Coverage Adsorbate Search** searches for adsorbate configurations on catalytic metal surfaces across increasing coverages. Given a clean-metal-surface POSCAR or CONTCAR and an adsorbate, INCAS performs a coverage-by-coverage search locally. The project packages a proven ACAT + MACE workflow as a standalone Python package with an `incas` command; it does not require the original project directory, a user-specific conda path, or a VASP installation.

INCAS was released as CASCAT up to version 0.4.0. The CASCAT name is now used by a separate, extended project; results and commands from CASCAT 0.4.x correspond to INCAS 0.5.0 apart from the name.

## Installation

### Ask an assistant to install and run INCAS

This repository includes the [run-incas-search skill](skill/run-incas-search/SKILL.md). After downloading or cloning the complete repository, open the project directory in your assistant and use a prompt such as:

```text
Read skill/run-incas-search/SKILL.md and follow its workflow to install and run INCAS.
My surface file is /path/to/POSCAR, the adsorbate is CO, and the target is N=4.
At each coverage, start from two configurations retained at the previous coverage.
After single-point screening, keep the top 80%; after coarse relaxation, keep the top 50%;
after full relaxation, save the top 10 results.
```

The skill guides the assistant through environment checks, dependency installation, input validation, execution, progress monitoring, and final-output inspection. It asks only for required information that is missing. You can also request installation only, candidate generation only, a resumed run, or result inspection.

The skill is stored in the visible `skill/` directory. Explicitly ask the assistant to read the `SKILL.md` above instead of relying on automatic discovery. Running `pip install .` alone does not install the skill in a user's global skill directory.

When publishing the project, include the `skill/` directory and `MANIFEST.in`. The skill instructions, installation reference, and inspection script are stored there, and the source distribution includes them as well.

### Manual installation

Use Python 3.11 or 3.12. From the project directory, run the following commands on macOS or Linux; WSL2 is recommended on Windows:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
incas molecules
```

Installing dependencies and running MACE for the first time require network access. INCAS uses the `oc20_usemppbe` head of MACE-MH-1; after it has been cached, calculations can run offline. The defaults are CPU execution, one optimization worker, and two threads per worker. Each worker loads its own model instance; if enough memory is available, use `--workers 4` to run four workers.

## Quick start

```bash
incas run /path/to/POSCAR --ads CO --n 4 --parents 2 \
  --screen-keep 80% --relax-keep 50% --final-keep 10 --output results/CO
```

This command searches for 1, 2, 3, and then 4 CO molecules, starting from the clean slab. `--n` is the target number of molecules and `--ads` is the adsorbate name. **The number of parent configurations and the retention setting at each filtering stage are controlled independently:**

| Option | Meaning | Accepted values | Default |
| --- | --- | --- | --- |
| `--parents` | For N=n, number of lowest-energy configurations taken from the final retained results at N=n−1 | Positive integer or `all` | `1` |
| `--screen-keep` | Number of deduplicated initial candidates sent to coarse relaxation after single-point-energy ranking | Positive integer, percentage, or `all` | `80%` |
| `--relax-keep` | Number of valid graph-unique structures sent to full relaxation after coarse-relaxation ranking | Positive integer, percentage, or `all` | `50%` |
| `--final-keep` | Number of eligible results saved after full-relaxation ranking; the lowest-energy result is BEST | Positive integer, percentage, or `all` | `all` |

In the example, the top 80% of candidates ranked by single-point energy enter coarse relaxation at every N. The top 50% of valid unique coarse-relaxed configurations then enter full relaxation, and at most 10 lowest-energy eligible results are saved. At the next coverage, INCAS adds one molecule to each of at most two lowest-energy parent configurations selected from those saved results. N=1 starts from the clean slab.

A single `--parents` value applies to every coverage and does not affect the three stage-specific retention settings. If only three final results are saved, `--parents 5` can still use only those three. To start from more parent configurations, ensure that `--final-keep` retains enough results as well.

Fixed counts can also differ between stages:

```bash
incas run CONTCAR --ads OH --n 4 --parents 2 \
  --screen-keep 1024 --relax-keep 256 --final-keep 10 --output results/OH
```

Each percentage is calculated independently from the eligible, deduplicated structure count at the current stage and is rounded up. For example, 80% of 11 structures gives 9, while 50% of 7 gives 4. At least one candidate is retained when candidates exist, but no structure is invented when a stage has no eligible candidates. Integer settings start at 1; when fewer structures exist than requested, all available structures are used. Percentages may be fractional and must be greater than 0% and no greater than 100%.

Each option may independently be set to `all`, meaning that no count limit is applied at that point. Deduplication, integrity checks, and final convergence checks still run, with no additional energy window or hidden count cap. `--parents all` uses **all final results retained at the previous coverage**. To make every eligible fully relaxed result a parent, also set `--final-keep all`. Computational cost can grow quickly at higher coverages.

Version 0.3.0 removed the overloaded `--keep` option and replaced it with the four independent options above. Preserve results from older versions and use a new output directory when running with the new rules.

By default, INCAS uses the input cell unchanged and does not create a 2×2 supercell. Specify `--repeat 2 2 1` when a larger cell is required. The surface must lie in the xy plane, with vacuum along z, and adsorption occurs on the top surface. The input should be a clean metal slab without pre-adsorbed species. All substrate atoms are fixed; only adsorbate atoms are relaxed.

By default, `--surface auto` identifies the surface using the AutoAdsorbate / ACAT CustomSurface route from the original workflow. Known FCC surfaces can instead use labels such as `--surface 111`; do not assign an FCC label to a BCC or HCP input. Site-identification algorithms still determine the supported range, so arbitrary inputs are not guaranteed to be recognized correctly.

## Runnable example

`examples/POSCAR` is a 3×3, four-layer Cu(111) example. It can also be regenerated:

```bash
python examples/make_slab.py
incas run examples/POSCAR --ads H --n 2 --parents 1 \
  --screen-keep 80% --relax-keep 50% --final-keep all --surface 111 --output results/demo
```

For a quick trial, set both `--screen-keep` and `--relax-keep` to `1` to reduce the search space.

To check candidate generation without loading MACE:

```bash
incas run examples/POSCAR --ads CO2 --n 4 --parents 2 --surface 111 --output results/CO2 --generate-only
```

This mode generates every unique candidate for the **next incomplete coverage**. It does not truncate the pool with `--screen-keep`, because ranking requires MACE single-point energies. Later coverages require optimized parents from the preceding coverage. Remove `--generate-only` and rerun the same command to perform the full calculation.

## Molecule catalogue

Run `incas molecules` to list the templates and their internal names. The current catalogue contains 28 species covering HER, OER, common CO₂RR C1 pathways, and selected C2 species:

| Scope | Templates |
| --- | --- |
| HER and water-related species | H, H2, H2O, OH |
| OER | O, OH, OOH, O2, H2O, H2O2 |
| CO₂RR C1 | CO2, COOH, HCOO, HCOOH, CO, CHO, COH, CHOH, CH2O, CH2OH, CH3O, CH3OH, C, CH, CH2, CH3, CH4 |
| Selected CO₂RR C2 species | CH2CO, CH3CO, CH3COOH |

These templates are neutral geometric fragments placed on a surface; they do not encode the charge states of solution species. HCOO is a formate geometry template. H⁺, OH⁻, electrons, and the electrolyte are not treated as independent adsorption-search objects. Each search uses a single adsorbate species; mixed adsorption is not supported.

The complete set of CO₂RR species depends on the pathway being studied. This release does not claim to contain every intermediate: complete C–C coupling networks and pathways to ethylene or ethanol, for example, require additional templates and validation. New entries are maintained centrally in `src/incas/molecules.py`. A new template must pass generation and recognition tests; adding only a chemical-formula name is insufficient.

Templates reuse ACAT geometry and binding-order names, such as H2O → OH2 and CH3O → OCH3. Multidentate molecules are counted by their atom sets so that multiple occupied sites are not mistaken for multiple molecules. Recognition uses ACAT connectivity and site analysis. Successful template generation does not guarantee stability after optimization on a particular surface. See the [ACAT building documentation](https://asm-dtu.gitlab.io/acat/build.html) and [coverage-identification documentation](https://asm-dtu.gitlab.io/acat/modules.html).

## Search workflow and output

The main computational stages follow the original scripts: incremental ACAT generation and symmetry deduplication → float32 single-point energies and ranked selection → coarse relaxation → adsorption-graph deduplication and ranked selection → float64 full relaxation → final deduplication and result selection → independent selection of parent configurations for the next coverage. “Short optimization” and “coarse relaxation” refer to the same stage; “full optimization” and “full relaxation” refer to the same stage.

Coarse relaxation uses 10 steps by default. Full relaxation uses up to 120 cumulative steps with a convergence threshold of 0.05 eV/Å. Each stage uses its corresponding retention option. The original H₂/O₂ formation exclusion thresholds are retained for H and O. Run `incas run --help` for all options.

Outputs for each coverage are stored in directories such as `results/job-name/N_1/` and `N_2/`:

- `BEST.vasp`: the lowest-energy eligible configuration for that coverage.
- `final_results.traj` / `.xyz`: results selected by `--final-keep`, sorted by full-relaxation energy; the first frame corresponds to BEST.
- `growth_parents.traj` / `.xyz`: the subset of final results selected by `--parents` for the next coverage.
- `candidate_pool.traj` / `.xyz`: the initial candidate pool.
- `ranking.csv`: complete ranking, convergence, and graph-deduplication information. `retained_final` and `growth_parent` mark the final retained results and growth parents, respectively. Unselected structures remain in the stage trajectories and ranking for auditability.
- `pre_relaxation_ranking.csv` and stage trajectories: intermediate results.
- `metadata.json`: status, actual retained counts, and selection records.

The job root stores `slab.vasp` and `search.json`, including the input-file hash, parameters, and dependency versions. Rerunning the same command reuses coverage rounds whose status is `complete`; incomplete rounds are recalculated. Use a new output directory after changing the molecule, input, retention settings, total number of coverages, or computational parameters. INCAS does not overwrite results created with a different configuration.

Multiple tasks must not use the same output directory concurrently. A normal exit removes `.running.lock`; after forcibly terminating a process, verify that the old process has stopped before removing the lock manually. If the candidate pool cannot grow further, INCAS reports an explicit error. This release does not silently fall back to unselected candidates from earlier stages, thereby preserving the user-specified parent scope.

## Citation

If INCAS contributes to published research, cite the upstream methods actually used. The current workflow uses ACAT to generate and recognize adsorption configurations and the `oc20_usemppbe` head of MACE-MH-1 for single-point energies and geometry optimization. At minimum, cite the following three works:

1. **ACAT** — S. Han, S. Lysgaard, T. Vegge and H. A. Hansen, “Rapid mapping of alloy surface phase diagrams via Bayesian evolutionary multitasking,” *npj Computational Materials* **9**, 139 (2023). [DOI: 10.1038/s41524-023-01087-4](https://doi.org/10.1038/s41524-023-01087-4)
2. **MACE architecture** — I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner and G. Csányi, “MACE: Higher order equivariant message passing neural networks for fast and accurate force fields,” *Advances in Neural Information Processing Systems* **35**, 11423–11436 (2022). [arXiv:2206.07697](https://arxiv.org/abs/2206.07697)
3. **MACE-MH-1 model** — I. Batatia, C. Lin, J. Hart, E. Kasoar, A. M. Elena, S. W. Norwood, T. Wolf and G. Csányi, “Cross Learning between Electronic Structure Theories for Unifying Molecular, Surface, and Inorganic Crystal Foundation Force Fields” (2025). [arXiv:2510.25380](https://arxiv.org/abs/2510.25380)

INCAS also depends on ASE, for which the recommended citation is: A. H. Larsen *et al.*, “The Atomic Simulation Environment—A Python library for working with atoms,” *Journal of Physics: Condensed Matter* **29**, 273002 (2017). [DOI: 10.1088/1361-648X/aa680e](https://doi.org/10.1088/1361-648X/aa680e)

For surfaces not represented by native ACAT FCC labels, INCAS uses [AutoAdsorbate](https://github.com/basf/autoadsorbate) to help identify the surface termination. That project currently provides no formal citation format, so identify AutoAdsorbate by name, version, and repository link in the methods or software section. Publications should also report the INCAS version, MACE-MH-1, the selected head, and the `--parents`, `--screen-keep`, `--relax-keep`, and `--final-keep` settings. Citations do not replace the corresponding software and model licenses; consult the [official MACE-MH-1 model card](https://huggingface.co/mace-foundations/mace-mh-1) for license information.

## Validation and limitations

```bash
python -m unittest discover -s tests -v
```

Tests cover single-molecule candidate generation and intact recognition for every catalogue entry, retention rules with more than two parent configurations, and rejection of an incorrect resume attempt after the input changes. See `VALIDATION.md` for two-coverage tests using MACE and for dependency information.

This release primarily targets the metal surfaces used by the original workflow. Oxides, carbon materials, and other substrates containing elements such as O, C, H, N, S, or F are not supported; INCAS rejects them to avoid confusing substrate atoms with adsorbates. Production searches have not been validated on every operating system, in every fresh environment, or for every molecule.

This is a configurational search over finite candidate and parent sets and does not guarantee the global minimum. MACE energies are used only to rank systems with the same adsorbate and molecule count; they do not represent constant-potential or solvated free energies, nor final DFT results. The original Cu/Pt data and scripts remain in the parent directory. This project can be copied and distributed independently without the historical calculation data.

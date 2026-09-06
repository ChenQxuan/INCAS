from __future__ import annotations

import argparse
import csv
import gc
import itertools
import json
import multiprocessing as mp
import os
import re
import sys
import time
import warnings
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from ase import Atoms
from ase.formula import Formula
from .retention import retention_count, label
from ase.constraints import FixAtoms
from ase.geometry import find_mic
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import FIRE
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
from acat.build.adlayer import SystematicPatternGenerator
from acat.ga.graph_comparators import WLGraphComparator
from acat.settings import CustomSurface, site_heights


MIN_ADSORBATE_DISTANCE = {"H": 1.8, "O": 2.2, "OH": 2.2}
CLEAR_DIMER_DISTANCE = {"H": 0.90, "O": 1.50, "OH": None}
ADSORBATE_STOICHIOMETRY = {
    "H": {"H": 1},
    "O": {"O": 1},
    "OH": {"O": 1, "H": 1},
}
ADSORBATE_ANCHOR = {"H": "H", "O": "O", "OH": "O"}
ACAT_GRAPH_HMAX = 4
DEFAULT_SCREEN_WINDOW_EV = 1.0
DEFAULT_MAX_PRE_RELAX = 256
DEFAULT_PRE_RELAX_STEPS = 10
DEFAULT_PRE_RELAX_KEEP_PER_GRAPH = 2
DEFAULT_MIN_RELAX = 32
DEFAULT_MAX_RELAX = 128
DEFAULT_FMAX = 0.05
DEFAULT_MAX_STEPS = 120
DEFAULT_WORKERS = 4
DEFAULT_THREADS_PER_WORKER = 2

MACE_MODEL = "mh-1"
MACE_HEAD = "oc20_usemppbe"
DEFAULT_SCREEN_DTYPE = "float32"
DEFAULT_MACE_DTYPE = "float64"
DEFAULT_GROWTH_GAP_EV = 0.01
DEFAULT_MAX_GROWTH_PARENTS = 2

BUILTIN_FCC_SURFACES = {
    "111", "100", "110", "211", "221", "311", "322", "331", "332"
}

UNIDENTIFIED_SITE_PATTERN = re.compile(r"Cannot identify site \(([^)]+)\)")


def format_duration(seconds: float) -> str:
    """Format wall time compactly for live terminal progress."""
    seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours:d}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes:d}m{seconds:02d}s"
    return f"{seconds:d}s"


def progress_message(stage: str, completed: int, total: int, started: float) -> str:
    """Return count, percentage, elapsed time and estimated time remaining."""
    elapsed = time.perf_counter() - started
    fraction = completed / total if total else 1.0
    remaining = elapsed * (total - completed) / completed if completed else 0.0
    return (
        f"  [{stage}] {completed}/{total} ({100.0 * fraction:5.1f}%)  "
        f"elapsed={format_duration(elapsed)}  ETA={format_duration(remaining)}"
    )


def unidentified_site_topologies(caught_warnings) -> set[tuple[int, ...]]:
    """Extract atom-index tuples that ACAT detected but could not name."""
    topologies = set()
    for caught in caught_warnings:
        match = UNIDENTIFIED_SITE_PATTERN.fullmatch(str(caught.message))
        if match is None:
            continue
        topologies.add(
            tuple(sorted(int(value.strip()) for value in match.group(1).split(",")))
        )
    return topologies


def wrap_site_position(position, atoms: Atoms) -> np.ndarray:
    """Wrap a site position into the slab's periodic unit cell."""
    fractional = np.asarray(position) @ np.linalg.pinv(atoms.cell.array)
    for axis, periodic in enumerate(atoms.pbc):
        if periodic:
            fractional[axis] %= 1.0
    return fractional @ atoms.cell.array


def append_unknown_site(
    adsorption_sites: SlabAdsorptionSites,
    atoms: Atoms,
    site_type: str,
    position,
    normal,
    indices: tuple[int, ...],
) -> bool:
    """Keep a geometrically valid ACAT site whose morphology is unknown."""
    position = wrap_site_position(position, atoms)
    normal = np.asarray(normal, dtype=float)
    if normal.shape != (3,) or not np.all(np.isfinite(normal)):
        return False
    norm = float(np.linalg.norm(normal))
    if norm == 0.0:
        return False
    normal /= norm
    if normal[2] < 0.0:
        normal *= -1.0

    for existing in adsorption_sites.site_list:
        if existing["site"] != site_type:
            continue
        if periodic_distance(existing["position"], position, atoms) < 1e-5:
            return False

    surface_name = adsorption_sites.site_list[0]["surface"]
    site = adsorption_sites.new_site()
    site.update(
        {
            "site": site_type,
            "surface": surface_name,
            "morphology": "unknown",
            "position": np.round(position, 8),
            "normal": np.round(normal, 8),
            "indices": indices,
            "topology": np.asarray(indices),
        }
    )
    adsorption_sites.site_list.append(site)
    return True


def recover_unclassified_geometric_sites(
    atoms: Atoms,
    adsorption_sites: SlabAdsorptionSites,
    surface,
    metal: str,
    construction_warnings,
) -> dict[str, int]:
    """Restore valid sites that ACAT discards only for missing morphology.

    ACAT 2.0.3 detects these geometric sites but executes ``continue`` when
    its facet-specific morphology table has no matching name.  Morphology is
    not needed for placement, so retain the position/topology/normal and use
    ``unknown`` as the graph label component.
    """
    recovered = {"4fold": 0, "bridge": 0}

    if isinstance(surface, CustomSurface):
        # Disable cell expansion in this probe so warning topologies refer to
        # the original atom indices. The production object above still uses
        # ACAT's normal expanded-cell handling.
        with warnings.catch_warnings(record=True) as probe_warnings:
            warnings.simplefilter("always")
            probe = SlabAdsorptionSites(
                atoms,
                surface=surface,
                allow_6fold=False,
                composition_effect=False,
                both_sides=False,
                label_sites=False,
                surrogate_metal=metal,
                _allow_expand=False,
            )
        unidentified = unidentified_site_topologies(probe_warnings)
        if unidentified:
            _, _, fold4_positions = probe.delaunay_triangulation()
            reference_cell = probe.ref_atoms.cell
            surface_indices = probe.surf_ids
            surface_positions = probe.ref_atoms.positions[surface_indices]
            for reference_position in fold4_positions:
                fractional = reference_position @ np.linalg.pinv(
                    reference_cell.array
                )
                if not (
                    -1e-5 < fractional[0] < 1.0 - 1e-5
                    and -1e-5 < fractional[1] < 1.0 - 1e-5
                ):
                    continue
                _, distances = find_mic(
                    surface_positions - reference_position,
                    reference_cell,
                    pbc=probe.pbc,
                )
                indices = tuple(
                    sorted(
                        surface_indices[index]
                        for index in np.argsort(distances)[:4]
                    )
                )
                if indices not in unidentified:
                    continue
                mapped_fractional = reference_position @ np.linalg.pinv(
                    reference_cell.array
                )
                position = (
                    mapped_fractional
                    + np.average(probe.delta_positions[list(indices)], axis=0)
                ) @ probe.cell.array
                normal = adsorption_sites.get_surface_normal(list(indices[:3]))
                recovered["4fold"] += append_unknown_site(
                    adsorption_sites,
                    atoms,
                    "4fold",
                    position,
                    normal,
                    indices,
                )

    else:
        # Built-in fcc331 also detects a small number of bridge positions that
        # its morphology decision tree cannot name.
        for indices in unidentified_site_topologies(construction_warnings):
            if len(indices) != 2:
                continue
            vector, _ = find_mic(
                atoms.positions[indices[1]] - atoms.positions[indices[0]],
                atoms.cell,
                pbc=atoms.pbc,
            )
            position = atoms.positions[indices[0]] + 0.5 * vector
            shared_normals = [
                np.asarray(site["normal"], dtype=float)
                for site in adsorption_sites.site_list
                if site["site"] in {"3fold", "4fold", "fcc", "hcp"}
                and len(set(site["indices"]).intersection(indices)) == 2
                and np.all(np.isfinite(site["normal"]))
            ]
            if not shared_normals:
                continue
            recovered["bridge"] += append_unknown_site(
                adsorption_sites,
                atoms,
                "bridge",
                position,
                np.average(shared_normals, axis=0),
                indices,
            )

    adsorption_sites.site_list.sort(
        key=lambda site: (site["indices"], site["position"].tolist())
    )
    return recovered


def prepare_site_labels(
    adsorption_sites: SlabAdsorptionSites,
) -> tuple[int, list[str]]:
    """Add deterministic labels for retained sites missing from ACAT tables.

    With ``n_layers=None``, AutoAdsorbate identifies the exposed termination
    and ACAT records its morphology as ``unknown``.  ACAT 2.0.3 has no
    numerical labels for signatures such as ``ontop|unknown`` even though its
    graph-based uniqueness checks require every usable site to have a label.
    Keep ACAT's existing labels unchanged and append labels only for signatures
    absent from its table.  Sites without a finite 3-vector normal cannot be
    used to place an adsorbate and are discarded before generation.
    """
    usable_sites = []
    for site in adsorption_sites.site_list:
        normal = np.asarray(site.get("normal"), dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            continue
        usable_sites.append(site)

    filtered_count = len(adsorption_sites.site_list) - len(usable_sites)
    if not usable_sites:
        raise RuntimeError("AutoAdsorbate found no sites with valid surface normals")
    adsorption_sites.site_list = usable_sites

    def signature(site: dict) -> str:
        morphology = site.get("morphology") or "unknown"
        return f"{site['site']}|{morphology}"

    missing_signatures = sorted(
        {
            signature(site)
            for site in usable_sites
            if signature(site) not in adsorption_sites.label_dict
        }
    )
    next_label = max(adsorption_sites.label_dict.values(), default=0) + 1
    for site_signature in missing_signatures:
        adsorption_sites.label_dict[site_signature] = next_label
        next_label += 1
    for site in usable_sites:
        site["label"] = adsorption_sites.label_dict[signature(site)]

    # Downstream SlabAdsorbateCoverage and SystematicPatternGenerator may now
    # use the completed label dictionary exactly as they do for named sites.
    adsorption_sites.label_sites = True
    return filtered_count, missing_signatures


def build_adsorption_sites(atoms: Atoms, facet: str, metal: str):
    if facet in BUILTIN_FCC_SURFACES:
        surface = f"fcc{facet}"
        mode = surface
        label_sites = True
    else:
        # ACAT's CustomSurface keeps the same site-search algorithm available
        # for Miller indices that do not have a built-in fcc label.
        surface = CustomSurface(atoms.copy(), n_layers=None)
        mode = "CustomSurface"
        # ACAT 2.0.3 cannot label AutoAdsorbate's ``unknown`` morphology during
        # construction.  Build first, then append only the missing labels.
        label_sites = False
    surrogate_metal = (
        metal if metal in {"Ni", "Cu", "Pd", "Ag", "Pt", "Au"} else None
    )
    with warnings.catch_warnings(record=True) as construction_warnings:
        warnings.simplefilter("always")
        sas = SlabAdsorptionSites(
            atoms,
            surface=surface,
            allow_6fold=False,
            composition_effect=False,
            both_sides=False,
            label_sites=label_sites,
            surrogate_metal=surrogate_metal,
        )
    for caught in construction_warnings:
        if UNIDENTIFIED_SITE_PATTERN.fullmatch(str(caught.message)) is None:
            warnings.warn(str(caught.message), caught.category, stacklevel=2)
    recovered_sites = recover_unclassified_geometric_sites(
        atoms,
        sas,
        surface,
        metal,
        construction_warnings,
    )
    if isinstance(surface, CustomSurface) or any(recovered_sites.values()):
        filtered_count, missing_signatures = prepare_site_labels(sas)
    else:
        filtered_count, missing_signatures = 0, []
    recovery_summary = ",".join(
        f"{site_type}:{count}"
        for site_type, count in recovered_sites.items()
        if count
    ) or "none"
    if isinstance(surface, CustomSurface):
        mode = (
            "CustomSurface(AutoAdsorbate; "
            f"fallback_labels={','.join(missing_signatures) or 'none'}; "
            f"recovered_unknown_sites={recovery_summary}; "
            f"filtered_invalid_normals={filtered_count})"
        )
    elif any(recovered_sites.values()):
        mode = f"{surface}(recovered_unknown_sites={recovery_summary})"
    return sas, mode


def periodic_distance(position_a, position_b, atoms: Atoms) -> float:
    vector, _ = find_mic(
        np.asarray(position_a) - np.asarray(position_b), atoms.cell, atoms.pbc
    )
    return float(np.linalg.norm(vector))


def classify_adsorbates(
    atoms: Atoms,
    adsorption_sites: SlabAdsorptionSites,
    adsorbate: str,
) -> dict:
    """Map intact target adsorbates to ACAT sites after generation/relaxation."""
    coverage = SlabAdsorbateCoverage(
        atoms,
        adsorption_sites=adsorption_sites,
        subtract_heights=site_heights,
        label_occupied_sites=True,
        dmax=3.0,
    )
    # ACAT consolidates the two O-H groups in peroxide to O2H2.
    recognized_name = {"OHOH": "O2H2"}.get(adsorbate, adsorbate)
    occupied = [
        site
        for site in coverage.get_sites(occupied=True)
        if site.get("adsorbate")
        and list(Formula(site["adsorbate"])) == list(Formula(recognized_name))
    ]
    return {
        "coverage": coverage,
        "occupied": occupied,
        "site_types": sorted(str(site["site"]) for site in occupied),
        "site_morphologies": sorted(
            str(site.get("morphology", "")) for site in occupied
        ),
    }


def prepare_candidate(
    atoms: Atoms,
    adsorption_sites: SlabAdsorptionSites,
    n_substrate: int,
    adsorbate: str,
) -> Atoms:
    atoms = atoms.copy()
    atoms.calc = None
    atoms.set_constraint(FixAtoms(indices=range(n_substrate)))
    classification = classify_adsorbates(atoms, adsorption_sites, adsorbate)
    atoms.info["site_types"] = classification["site_types"]
    atoms.info["site_morphologies"] = classification["site_morphologies"]
    return atoms


def generate_acat_children(
    parents: list[Atoms],
    adsorption_sites: SlabAdsorptionSites,
    number_adsorbates: int,
    n_substrate: int,
    adsorbate: str,
) -> tuple[list[Atoms], dict]:
    """Add one adsorbate to selected parents with native ACAT dedup."""
    generation_started = time.perf_counter()
    print(
        f"  [ACAT] starting {adsorbate} N={number_adsorbates} from "
        f"{len(parents)} selected parent(s)",
        flush=True,
    )

    with TemporaryDirectory(prefix=f"acat_{adsorbate}_patterns_") as temporary:
        temporary_path = Path(temporary)
        trajectory_path = temporary_path / f"N_{number_adsorbates}.traj"
        log_path = temporary_path / f"N_{number_adsorbates}.log"
        generator = SystematicPatternGenerator(
            [atoms.copy() for atoms in parents],
            adsorbate_species=[adsorbate],
            min_adsorbate_distance=MIN_ADSORBATE_DISTANCE[adsorbate],
            adsorption_sites=adsorption_sites,
            subtract_heights=True,
            smart_skip=False,
            enumerate_orientations=True,
            trajectory=str(trajectory_path),
            logfile=str(log_path),
        )
        generator.run(
            action="add",
            num_act=1,
            unique=True,
            hmax=ACAT_GRAPH_HMAX,
            subsurf_effect=True,
        )
        generator.traj.close()
        generator.logfile.close()
        raw_patterns = read(trajectory_path, index=":") if trajectory_path.stat().st_size else []
        patterns = [
            prepare_candidate(atoms, adsorption_sites, n_substrate, adsorbate)
            for atoms in raw_patterns
            if has_expected_adsorbate_count(
                atoms, n_substrate, adsorbate, number_adsorbates
            )
        ]
    if not patterns:
        raise RuntimeError(
            f"ACAT generated no unique {adsorbate} N={number_adsorbates} patterns"
        )
    metadata = {
        "method": "ACAT_SystematicPatternGenerator_greedy_incremental",
        "action": f"add_one_{adsorbate}",
        "unique": True,
        "hmax": ACAT_GRAPH_HMAX,
        "subsurf_effect": True,
        "subtract_heights": True,
        "smart_skip": False,
        "n_parent_patterns": len(parents),
        "candidate_pool_size": len(patterns),
        "generation_seconds": time.perf_counter() - generation_started,
    }
    print(
        f"  [ACAT] {adsorbate} N={number_adsorbates}: "
        f"{len(parents)} selected parents -> "
        f"{len(patterns)} unique patterns; "
        f"elapsed={format_duration(metadata['generation_seconds'])}",
        flush=True,
    )
    return patterns, metadata


def has_expected_adsorbate_count(
    atoms: Atoms,
    n_substrate: int,
    adsorbate: str,
    number_adsorbates: int,
) -> bool:
    expected = {
        symbol: coefficient * number_adsorbates
        for symbol, coefficient in ADSORBATE_STOICHIOMETRY[adsorbate].items()
    }
    actual = {}
    for atom in atoms[n_substrate:]:
        actual[atom.symbol] = actual.get(atom.symbol, 0) + 1
    return actual == expected


def anchor_minimum(
    atoms: Atoms, n_substrate: int, adsorbate: str
) -> float | None:
    anchor = ADSORBATE_ANCHOR[adsorbate]
    indices = [
        index
        for index in range(n_substrate, len(atoms))
        if atoms[index].symbol == anchor
    ]
    if len(indices) < 2:
        return None
    return min(
        periodic_distance(atoms.positions[i], atoms.positions[j], atoms)
        for i, j in itertools.combinations(indices, 2)
    )


def write_candidate_pool(path: Path, candidates: list[Atoms]) -> None:
    trajectory = Trajectory(path, mode="w")
    try:
        for atoms in candidates:
            trajectory.write(atoms)
    finally:
        trajectory.close()


def write_xyz_pool(path: Path, candidates: list[Atoms]) -> None:
    """Write a human-readable multi-frame extended XYZ structure pool."""
    if not candidates:
        return
    write(path, candidates, format="extxyz")


def select_for_relaxation(
    scored: list[dict],
    energy_key: str,
    energy_window_eV: float,
    minimum: int,
    maximum: int,
) -> tuple[list[dict], dict]:
    """Select a deterministic low-energy window with resource safety bounds."""
    ordered = sorted(scored, key=lambda row: row[energy_key])
    best_energy = ordered[0][energy_key]
    within_window = [
        row
        for row in ordered
        if row[energy_key] <= best_energy + energy_window_eV
    ]
    n_selected = max(len(within_window), min(minimum, len(ordered)))
    window_truncated = bool(maximum and n_selected > maximum)
    if maximum:
        n_selected = min(n_selected, maximum)
    selected = ordered[:n_selected]
    return selected, {
        "energy_key": energy_key,
        "best_energy_eV": best_energy,
        "energy_window_eV": energy_window_eV,
        "n_within_window": len(within_window),
        "minimum": minimum,
        "maximum": maximum,
        "n_selected": len(selected),
        "window_truncated_by_maximum": window_truncated,
    }


def select_lowest_with_cap(
    scored: list[dict], energy_key: str, maximum: int | str | None
) -> tuple[list[dict], dict]:
    """Select the lowest-energy rows with only a hard resource cap."""
    ordered = sorted(scored, key=lambda row: row[energy_key])
    n_selected = retention_count(len(ordered), maximum)
    return ordered[:n_selected], {
        "strategy": "lowest_energy_after_graph_dedup",
        "energy_key": energy_key,
        "requested": label(maximum),
        "n_available": len(ordered),
        "n_selected": n_selected,
        "truncated_by_maximum": n_selected < len(ordered),
    }


WORKER_CALCULATOR = None


def initialize_mace_worker(dtype: str, threads: int) -> None:
    """Load one private MACE model in each spawned worker process."""
    global WORKER_CALCULATOR
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)

    import torch
    from mace.calculators import mace_mp

    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    WORKER_CALCULATOR = mace_mp(
        model=MACE_MODEL,
        device="cpu",
        default_dtype=dtype,
        head=MACE_HEAD,
    )
    # MACECalculator otherwise requests stress for every ASE call.  Geometry
    # optimization needs energy and forces, but no cell derivative because the
    # cell and every substrate atom are fixed.
    WORKER_CALCULATOR.model_type = "MACE_energy_forces_only"


def score_candidates_energy_only(
    candidates: list[Atoms], dtype: str, threads: int
) -> tuple[list[float], float]:
    """Run direct MACE inference with forces and every tensor derivative off."""
    started = time.perf_counter()
    import torch
    from mace.calculators import mace_mp

    torch.set_num_threads(threads)
    calculator = mace_mp(
        model=MACE_MODEL,
        device="cpu",
        default_dtype=dtype,
        head=MACE_HEAD,
    )
    interval = max(1, len(candidates) // 10)
    energies = []
    for completed, atoms in enumerate(candidates, start=1):
        batch = calculator._atoms_to_batch(atoms)
        model_energies = []
        for model in calculator.models:
            model_dtype = next(model.parameters()).dtype
            for key in batch.keys:
                value = batch[key]
                if torch.is_tensor(value) and torch.is_floating_point(value):
                    batch[key] = value.to(dtype=model_dtype)
            with torch.inference_mode():
                output = model(
                    batch.to_dict(),
                    training=False,
                    compute_force=False,
                    compute_virials=False,
                    compute_stress=False,
                    compute_displacement=False,
                    compute_hessian=False,
                    compute_edge_forces=False,
                    compute_atomic_stresses=False,
                )
            model_energies.append(float(output["energy"].reshape(-1)[0]))
        energies.append(
            float(np.mean(model_energies)) * calculator.energy_units_to_eV
        )
        if completed == 1 or completed % interval == 0 or completed == len(candidates):
            print(progress_message("Energy-only", completed, len(candidates), started), flush=True)
    elapsed = time.perf_counter() - started
    del calculator
    gc.collect()
    return energies, elapsed


def relax_candidate_worker(
    task: tuple[int, Atoms, float, int, int],
) -> tuple[int, Atoms, float, bool, int, float]:
    """Relax one fixed-substrate candidate independently in a worker."""
    candidate_id, atoms, fmax, max_steps, n_substrate = task
    atoms.calc = WORKER_CALCULATOR
    optimizer = FIRE(atoms, logfile=None)
    converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
    forces = atoms.get_forces()
    adsorbate_indices = range(n_substrate, len(atoms))
    adsorbate_force_max = max(
        float(np.linalg.norm(forces[index])) for index in adsorbate_indices
    )
    relaxed_energy = float(atoms.get_potential_energy())
    atoms.calc = None
    return (
        candidate_id,
        atoms,
        relaxed_energy,
        converged,
        int(optimizer.nsteps),
        adsorbate_force_max,
    )


def add_geometry_classification(
    row: dict,
    atoms: Atoms,
    adsorption_sites: SlabAdsorptionSites,
    prefix: str,
    n_substrate: int,
    adsorbate: str,
) -> None:
    """Attach ACAT site assignment and simple validity fields to one row."""
    anchor = ADSORBATE_ANCHOR[adsorbate]
    n_adsorbates = (len(atoms) - n_substrate) // sum(ADSORBATE_STOICHIOMETRY[adsorbate].values())
    classification = classify_adsorbates(atoms, adsorption_sites, adsorbate)
    occupied_sites = classification["occupied"]
    molecule_indices = {tuple(site["adsorbate_indices"]) for site in occupied_sites}
    covered_indices = {index for indices in molecule_indices for index in indices}
    valid_geometry = (
        len(molecule_indices) == n_adsorbates
        and covered_indices == set(range(n_substrate, len(atoms)))
        and has_expected_adsorbate_count(atoms, n_substrate, adsorbate, n_adsorbates)
    )
    site_types = "+".join(classification["site_types"])
    site_morphologies = "+".join(classification["site_morphologies"])
    if not valid_geometry:
        site_types = f"unclassified_{len(occupied_sites)}of{n_adsorbates}"
        site_morphologies = ""
    row[f"{prefix}_site_types"] = site_types
    row[f"{prefix}_site_morphologies"] = site_morphologies
    row[f"{prefix}_valid_geometry"] = valid_geometry
    row[f"{prefix}_min_anchor_distance_A"] = anchor_minimum(
        atoms, n_substrate, adsorbate
    )
    row[f"{prefix}_coverage"] = classification["coverage"]


def retain_graph_representatives(
    rows: list[dict],
    energy_key: str,
    valid_key: str,
    coverage_key: str,
    keep_per_graph: int,
    prefix: str,
    converged_key: str | None = None,
) -> list[dict]:
    """Retain low-energy representatives of each ACAT adsorption graph."""
    comparator = WLGraphComparator(hmax=ACAT_GRAPH_HMAX)
    groups: list[dict] = []
    retained = []
    ordering = sorted(
        rows,
        key=lambda item: (
            not item[converged_key] if converged_key else False,
            item[energy_key],
        ),
    )
    for row in ordering:
        row[f"{prefix}_graph_group"] = None
        row[f"{prefix}_graph_representative"] = False
        row[f"{prefix}_duplicate_of_candidate_id"] = None
        if not row[valid_key]:
            # A short relaxation can stop between recognized sites.  Keep such
            # a structure unless a separate physical-invalidity test rejects it.
            row[f"{prefix}_graph_representative"] = True
            retained.append(row)
            continue

        coverage = row[coverage_key]
        labels = coverage.get_occupied_labels(fragmentation=True)
        graph = coverage.get_graph(fragmentation=True, subsurf_effect=True)
        matching_group = None
        for group in groups:
            if labels == group["labels"] and comparator.looks_like(
                graph, group["graph"]
            ):
                matching_group = group
                break
        if matching_group is None:
            matching_group = {
                "index": len(groups) + 1,
                "labels": labels,
                "graph": graph,
                "representative_id": row["candidate_id"],
                "n_retained": 0,
            }
            groups.append(matching_group)

        row[f"{prefix}_graph_group"] = matching_group["index"]
        if matching_group["n_retained"] < keep_per_graph:
            row[f"{prefix}_graph_representative"] = True
            matching_group["n_retained"] += 1
            retained.append(row)
        else:
            row[f"{prefix}_duplicate_of_candidate_id"] = matching_group[
                "representative_id"
            ]
    return retained


def write_pre_relaxation_ranking(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "candidate_id",
        "site_types_initial",
        "site_morphologies_initial",
        "initial_energy_eV",
        "pre_relax_energy_eV",
        "pre_relax_steps",
        "pre_relax_converged",
        "pre_max_adsorbate_force_eV_A",
        "pre_site_types",
        "pre_site_morphologies",
        "pre_valid_geometry",
        "pre_min_anchor_distance_A",
        "pre_clear_dimer",
        "pre_graph_group",
        "pre_graph_representative",
        "pre_duplicate_of_candidate_id",
        "selected_for_full_relaxation",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: item["pre_relax_energy_eV"]):
            writer.writerow({key: row.get(key) for key in fieldnames})


def select_growth_rows(eligible, limit, gap):
    """Select the lowest eligible structures, optionally within an energy gap."""
    best = eligible[0]["relaxed_energy_eV"]
    return [row for row in eligible
            if gap is None or row["relaxed_energy_eV"] - best <= gap][:limit]


def run_mace_search(
    adsorption_sites: SlabAdsorptionSites,
    candidates_in: list[Atoms],
    output_dir: Path,
    n_substrate: int,
    adsorbate: str,
    screen_keep: int | str | None,
    relax_keep: int | str | None,
    final_keep: int | str | None,
    parent_limit: int | None,
    pre_relax_steps: int,
    fmax: float,
    max_steps: int,
    screen_dtype: str,
    dtype: str,
    workers: int,
    threads: int,
) -> tuple[list[dict], Atoms, list[Atoms], dict]:
    search_started = time.perf_counter()
    scored = []
    pool_structures = [atoms.copy() for atoms in candidates_in]
    for candidate_id, atoms in enumerate(pool_structures, start=1):
        atoms.info["candidate_id"] = candidate_id

    print(
        f"  [Screen] one {screen_dtype} MACE model; energy only "
        "(forces/stress/virials/tensors disabled)",
        flush=True,
    )
    score_values, score_seconds = score_candidates_energy_only(
        pool_structures, screen_dtype, threads
    )
    for candidate_id, (atoms, energy) in enumerate(
        zip(pool_structures, score_values), start=1
    ):
        atoms.info["initial_energy_eV"] = energy
        scored.append(
            {
                "candidate_id": candidate_id,
                "site_types_initial": "+".join(atoms.info["site_types"]),
                "site_morphologies_initial": "+".join(
                    atoms.info["site_morphologies"]
                ),
                "initial_energy_eV": energy,
            }
        )
    write_candidate_pool(output_dir / "candidate_pool.traj", pool_structures)
    write_xyz_pool(output_dir / "candidate_pool.xyz", pool_structures)

    context = mp.get_context("spawn")
    print(
        f"  [Optimize] {workers} {dtype} workers x {threads} threads; "
        "energy + forces only",
        flush=True,
    )
    with context.Pool(
        processes=workers,
        initializer=initialize_mace_worker,
        initargs=(dtype, threads),
    ) as process_pool:
        selected_for_pre, initial_selection = select_lowest_with_cap(
            scored,
            "initial_energy_eV",
            screen_keep,
        )
        print(
            f"  [Pre-select] {len(scored)} unique candidates -> "
            f"{len(selected_for_pre)} short relaxations; screen-keep={label(screen_keep)}",
            flush=True,
        )
        pre_selected_by_id = {
            row["candidate_id"]: row for row in selected_for_pre
        }
        pre_tasks = [
            (
                row["candidate_id"],
                pool_structures[row["candidate_id"] - 1].copy(),
                fmax,
                pre_relax_steps,
                n_substrate,
            )
            for row in selected_for_pre
        ]
        pre_results = []
        pre_started = time.perf_counter()
        for completed, result in enumerate(
            process_pool.imap_unordered(
                relax_candidate_worker, pre_tasks, chunksize=1
            ),
            start=1,
        ):
            candidate_id = result[0]
            row = pre_selected_by_id[candidate_id]
            print(
                progress_message(
                    "Pre-relax", completed, len(pre_tasks), pre_started
                )
                + f"  candidate={candidate_id} ({row['site_types_initial']})",
                flush=True,
            )
            pre_results.append(result)
        pre_relax_seconds = time.perf_counter() - pre_started

        pre_rows = []
        pre_structures = []
        for (
            candidate_id,
            atoms,
            pre_energy,
            pre_converged,
            pre_steps_done,
            pre_adsorbate_force_max,
        ) in pre_results:
            row = dict(pre_selected_by_id[candidate_id])
            row.update(
                {
                    "pre_relax_energy_eV": pre_energy,
                    "pre_relax_converged": pre_converged,
                    "pre_relax_steps": pre_steps_done,
                    "pre_max_adsorbate_force_eV_A": pre_adsorbate_force_max,
                    "atoms": atoms,
                }
            )
            add_geometry_classification(
                row, atoms, adsorption_sites, "pre", n_substrate, adsorbate
            )
            minimum_anchor_distance = row["pre_min_anchor_distance_A"]
            dimer_cutoff = CLEAR_DIMER_DISTANCE[adsorbate]
            row["pre_clear_dimer"] = bool(
                dimer_cutoff is not None
                and minimum_anchor_distance is not None
                and minimum_anchor_distance < dimer_cutoff
            )
            pre_rows.append(row)
            pre_structures.append(atoms)

        physically_retained = [
            row for row in pre_rows if not row["pre_clear_dimer"]
        ]
        if not physically_retained:
            raise RuntimeError(
                f"Every short-relaxed {adsorbate} candidate formed a clear dimer"
            )
        graph_representatives = retain_graph_representatives(
            physically_retained,
            energy_key="pre_relax_energy_eV",
            valid_key="pre_valid_geometry",
            coverage_key="pre_coverage",
            keep_per_graph=1,
            prefix="pre",
        )
        write_xyz_pool(
            output_dir / "pre_graph_representatives.xyz",
            [row["atoms"] for row in graph_representatives],
        )
        selected_for_full, full_selection = select_lowest_with_cap(
            graph_representatives,
            "pre_relax_energy_eV",
            relax_keep,
        )
        selected_for_full_ids = {
            row["candidate_id"] for row in selected_for_full
        }
        for row in pre_rows:
            row.setdefault("pre_graph_group", None)
            row.setdefault("pre_graph_representative", False)
            row.setdefault("pre_duplicate_of_candidate_id", None)
            row["selected_for_full_relaxation"] = (
                row["candidate_id"] in selected_for_full_ids
            )

        write_candidate_pool(output_dir / "pre_relaxed_pool.traj", pre_structures)
        write_pre_relaxation_ranking(
            output_dir / "pre_relaxation_ranking.csv", pre_rows
        )
        print(
            f"  [Pre-dedup] {len(pre_rows)} short-relaxed -> "
            f"{len(graph_representatives)} graph representatives -> "
            f"{len(selected_for_full)} full relaxations",
            flush=True,
        )
        if full_selection["truncated_by_maximum"]:
            print(
                f"  [Full-select] {full_selection['n_available']} graph "
                f"representatives available; keeping the lowest-energy "
                f"{len(selected_for_full)}",
                flush=True,
            )

        full_selected_by_id = {
            row["candidate_id"]: row for row in selected_for_full
        }
        full_tasks = [
            (
                row["candidate_id"],
                row["atoms"].copy(),
                fmax,
                max(1, max_steps - row["pre_relax_steps"]),
                n_substrate,
            )
            for row in selected_for_full
        ]
        full_results = []
        full_started = time.perf_counter()
        for completed, result in enumerate(
            process_pool.imap_unordered(
                relax_candidate_worker, full_tasks, chunksize=1
            ),
            start=1,
        ):
            candidate_id = result[0]
            row = full_selected_by_id[candidate_id]
            print(
                progress_message(
                    "Full-relax", completed, len(full_tasks), full_started
                )
                + f"  candidate={candidate_id} ({row['pre_site_types']})",
                flush=True,
            )
            full_results.append(result)
        full_relax_seconds = time.perf_counter() - full_started

    relaxed_rows = []
    relaxed_structures = []
    for (
        candidate_id,
        atoms,
        relaxed_energy,
        converged,
        full_steps_done,
        adsorbate_force_max,
    ) in full_results:
        pre_row = full_selected_by_id[candidate_id]
        row = {
            key: value
            for key, value in pre_row.items()
            if key not in {"atoms", "pre_coverage"}
        }
        row.update(
            {
                "relaxed_energy_eV": relaxed_energy,
                "converged": converged,
                "full_relax_steps": full_steps_done,
                "optimizer_steps": pre_row["pre_relax_steps"] + full_steps_done,
                "max_adsorbate_force_eV_A": adsorbate_force_max,
                "atoms": atoms,
            }
        )
        add_geometry_classification(
            row, atoms, adsorption_sites, "final", n_substrate, adsorbate
        )
        row["site_types_final"] = row.pop("final_site_types")
        row["site_morphologies_final"] = row.pop("final_site_morphologies")
        row["valid_final_geometry"] = row.pop("final_valid_geometry")
        row["min_anchor_distance_A"] = row.pop(
            "final_min_anchor_distance_A"
        )
        row["coverage"] = row.pop("final_coverage")
        atoms.info.update(
            {
                "candidate_id": candidate_id,
                "initial_energy_eV": row["initial_energy_eV"],
                "pre_relax_energy_eV": row["pre_relax_energy_eV"],
                "relaxed_energy_eV": relaxed_energy,
                "converged": converged,
                "optimizer_steps": row["optimizer_steps"],
            }
        )
        relaxed_rows.append(row)
        relaxed_structures.append(atoms)

    relaxed_rows.sort(key=lambda row: row["relaxed_energy_eV"])
    if not relaxed_rows:
        raise RuntimeError("No candidate was relaxed")

    retain_graph_representatives(
        relaxed_rows,
        energy_key="relaxed_energy_eV",
        valid_key="valid_final_geometry",
        coverage_key="coverage",
        keep_per_graph=1,
        prefix="final",
        converged_key="converged",
    )
    for row in relaxed_rows:
        row["unique_minimum"] = bool(
            row.pop("final_graph_representative")
            and row["valid_final_geometry"]
        )
        row["duplicate_of_candidate_id"] = row.pop(
            "final_duplicate_of_candidate_id"
        )
        row.pop("final_graph_group", None)

    final_unique_rows = [row for row in relaxed_rows if row["unique_minimum"]]
    write_xyz_pool(
        output_dir / "final_graph_representatives.xyz",
        [row["atoms"] for row in final_unique_rows],
    )

    eligible = [
        row
        for row in relaxed_rows
        if row["converged"] and row["valid_final_geometry"] and row["unique_minimum"]
    ]
    if not eligible:
        raise RuntimeError("No converged, valid, unique relaxed minimum was found")
    eligible.sort(key=lambda row: row["relaxed_energy_eV"])
    best_row = eligible[0]
    best_energy = best_row["relaxed_energy_eV"]
    retained_rows, final_selection = select_lowest_with_cap(eligible, "relaxed_energy_eV", final_keep)
    retained_structures = [row["atoms"].copy() for row in retained_rows]
    write_candidate_pool(output_dir / "final_results.traj", retained_structures)
    write_xyz_pool(output_dir / "final_results.xyz", retained_structures)
    growth_rows = select_growth_rows(retained_rows, parent_limit, None)
    growth_parents = [row["atoms"].copy() for row in growth_rows]
    write_candidate_pool(output_dir / "growth_parents.traj", growth_parents)
    write_xyz_pool(output_dir / "growth_parents.xyz", growth_parents)
    print(
        f"  [Final-dedup] {len(relaxed_rows)} full-relaxed -> "
        f"{len(final_unique_rows)} unique valid graphs -> "
        f"{len(eligible)} converged eligible structures -> {len(retained_rows)} saved results",
        flush=True,
    )
    print(
        f"  [Growth] saved {len(growth_parents)} parent structure(s) for N+1; "
        f"parents={label(parent_limit)}",
        flush=True,
    )
    eligible_ids = {row["candidate_id"] for row in eligible}
    retained_ids = {row["candidate_id"] for row in retained_rows}
    parent_ids = {row["candidate_id"] for row in growth_rows}
    for row in relaxed_rows:
        row["eligible_for_best"] = row["candidate_id"] in eligible_ids
        row["retained_final"] = row["candidate_id"] in retained_ids
        row["growth_parent"] = row["candidate_id"] in parent_ids
    relaxed_rows.sort(
        key=lambda row: (not row["eligible_for_best"], row["relaxed_energy_eV"])
    )
    for rank, row in enumerate(relaxed_rows, start=1):
        row["rank"] = rank
        row["delta_e_from_best_eV"] = (
            row["relaxed_energy_eV"] - best_energy
            if row["eligible_for_best"]
            else None
        )
        row.pop("coverage", None)

    write_candidate_pool(output_dir / "relaxed_top.traj", relaxed_structures)
    best_atoms = best_row["atoms"]
    output_order = list(range(n_substrate))
    for symbol in ADSORBATE_STOICHIOMETRY[adsorbate]:
        output_order.extend(
            index
            for index in range(n_substrate, len(best_atoms))
            if best_atoms[index].symbol == symbol
        )
    best_atoms_for_vasp = best_atoms[output_order]
    write(
        output_dir / "BEST.vasp",
        best_atoms_for_vasp,
        format="vasp",
        direct=True,
        vasp5=True,
        sort=False,
    )
    selection = {
        "initial_single_point": initial_selection,
        "pre_relaxation": {
            "steps": pre_relax_steps,
            "keep_per_graph": 1,
            "n_short_relaxed": len(pre_rows),
            "n_clear_dimer_rejected": sum(
                row["pre_clear_dimer"] for row in pre_rows
            ),
            "n_graph_representatives": len(graph_representatives),
        },
        "full_relaxation_selection": full_selection,
        "final": {
            "n_relaxed": len(relaxed_rows),
            "n_valid_relaxed": sum(row["valid_final_geometry"] for row in relaxed_rows),
            "n_unique_relaxed_minima": sum(row["unique_minimum"] for row in relaxed_rows),
            "n_eligible_converged_minima": len(eligible),
            "best_candidate_id": best_row["candidate_id"],
            "retention": final_selection,
            "n_retained_results": len(retained_rows),
            "retained_candidate_ids": [row["candidate_id"] for row in retained_rows],
            "parents_requested": label(parent_limit),
            "growth_parent_candidate_ids": [
                row["candidate_id"] for row in growth_rows
            ],
            "n_growth_parents": len(growth_rows),
        },
        "timing_seconds": {
            "single_point": score_seconds,
            "pre_relaxation": pre_relax_seconds,
            "full_relaxation": full_relax_seconds,
            "total_mace_search": time.perf_counter() - search_started,
        },
        "best_candidate_id": best_row["candidate_id"],
    }
    return relaxed_rows, best_atoms, growth_parents, selection


def write_ranking(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "retained_final",
        "growth_parent",
        "rank",
        "candidate_id",
        "site_types_initial",
        "site_morphologies_initial",
        "site_types_final",
        "site_morphologies_final",
        "initial_energy_eV",
        "pre_relax_energy_eV",
        "relaxed_energy_eV",
        "delta_e_from_best_eV",
        "converged",
        "valid_final_geometry",
        "unique_minimum",
        "eligible_for_best",
        "duplicate_of_candidate_id",
        "optimizer_steps",
        "pre_relax_steps",
        "full_relax_steps",
        "max_adsorbate_force_eV_A",
        "min_anchor_distance_A",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def best_geometry_record(
    best_atoms: Atoms, n_substrate: int, adsorbate: str
) -> dict:
    scaled = best_atoms.get_scaled_positions(wrap=True)
    substrate_z = best_atoms.positions[:n_substrate, 2]
    top_z = float(max(substrate_z))
    adsorbate_records = []
    for atom_index in range(n_substrate, len(best_atoms)):
        atom = best_atoms[atom_index]
        adsorbate_records.append(
            {
                "symbol": atom.symbol,
                "fractional_xy": [
                    float(scaled[atom_index, 0]),
                    float(scaled[atom_index, 1]),
                ],
                "height_above_top_A": float(atom.position[2] - top_z),
                "cartesian_reference_A": [float(value) for value in atom.position],
            }
        )
    return {
        "reference_top_z_A": top_z,
        "adsorbate": adsorbate,
        "atoms": adsorbate_records,
    }

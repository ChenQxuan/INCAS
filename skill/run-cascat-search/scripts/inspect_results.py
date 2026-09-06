#!/usr/bin/env python3
"""Read-only validation of CASCAT 0.4.x output trees."""

import argparse
from collections import Counter
import json
from pathlib import Path


def same_geometry(left, right):
    import numpy as np
    from ase.geometry import get_distances
    from scipy.optimize import linear_sum_assignment

    if Counter(left.symbols) != Counter(right.symbols) or not np.allclose(left.cell, right.cell, atol=1e-6):
        return False
    # BEST.vasp groups adsorbate elements; match by element, not atom order.
    for symbol in set(left.get_chemical_symbols()):
        a = left.positions[[atom.symbol == symbol for atom in left]]
        b = right.positions[[atom.symbol == symbol for atom in right]]
        _, distances = get_distances(a, b, cell=left.cell, pbc=left.pbc)
        rows, columns = linear_sum_assignment(distances)
        if not np.all(distances[rows, columns] < 1e-5):
            return False
    return True


def inspect(root):
    import numpy as np
    from ase.io import read

    root = Path(root).expanduser().resolve()
    report = {"output": str(root), "lock_present": (root / ".running.lock").exists(), "rounds": [], "errors": []}
    try:
        config = json.loads((root / "search.json").read_text())
        count = config["count"]
        if type(count) is not int or count < 1:
            raise ValueError("search.json count must be a positive integer")
        if not str(config.get("engine_version", "")).startswith("0.4."):
            raise ValueError("This checker expects CASCAT 0.4.x output with separate final results and parents")
        slab = read(root / "slab.vasp")
    except (OSError, ValueError, KeyError) as error:
        report["errors"].append(str(error))
        report["complete"] = False
        return report
    for n in range(1, count + 1):
        directory = root / f"N_{n}"
        item = {"N": n, "status": "missing", "errors": []}
        report["rounds"].append(item)
        try:
            metadata = json.loads((directory / "metadata.json").read_text())
            item["status"] = metadata.get("status", "unknown")
            if item["status"] != "complete":
                continue
            final = read(directory / "final_results.traj", index=":")
            parents = read(directory / "growth_parents.traj", index=":")
            best = read(directory / "BEST.vasp")
            selection = metadata["selection"]
            item.update(
                coarse=selection["initial_single_point"]["n_selected"],
                full=selection["full_relaxation_selection"]["n_selected"],
                final=len(final), parents=len(parents), best=str(directory / "BEST.vasp"),
            )
            if not final or not parents:
                raise ValueError("Completed coverage has an empty final or parent pool")
            if len(final) != metadata["n_final_results"] or len(parents) != metadata["n_growth_parents"]:
                item["errors"].append("Saved pool lengths disagree with metadata")
            energies = [atoms.info["relaxed_energy_eV"] for atoms in final]
            if not all(np.isfinite(energies)) or energies != sorted(energies):
                item["errors"].append("Final energies are not finite and sorted")
            final_ids = [atoms.info["candidate_id"] for atoms in final]
            parent_ids = [atoms.info["candidate_id"] for atoms in parents]
            if len(set(final_ids)) != len(final_ids):
                item["errors"].append("Final candidate IDs are duplicated")
            if final_ids != selection["final"]["retained_candidate_ids"] or parent_ids != selection["final"]["growth_parent_candidate_ids"]:
                item["errors"].append("Saved candidate IDs disagree with selection metadata")
            if parent_ids != final_ids[:len(parents)] or any(not same_geometry(parent, source) for parent, source in zip(parents, final)):
                item["errors"].append("Parents are not the lowest-energy retained final structures")
            if not same_geometry(final[0], best):
                item["errors"].append("BEST geometry differs from the first final result")
            for atoms in final:
                if not atoms.info.get("converged"):
                    item["errors"].append("A retained result is not marked converged")
                    break
                if len(atoms) <= len(slab) or not same_geometry(atoms[:len(slab)], slab):
                    item["errors"].append("A retained result changes the fixed substrate or lacks adsorbates")
                    break
        except (OSError, ValueError, KeyError, IndexError, TypeError) as error:
            item["errors"].append(str(error))
    report["complete"] = all(item["status"] == "complete" and not item["errors"] for item in report["rounds"])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        report = inspect(args.output)
    except ImportError as error:
        parser.exit(1, f"Use the installed adsorption-search environment: {error}\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

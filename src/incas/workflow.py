"""Portable input handling, round configuration and checked resume."""

import hashlib
import json
import os
from importlib.metadata import version
from pathlib import Path
from .retention import label
from . import __version__


def run(args):
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    from ase.io import read, write
    from acat.settings import adsorbate_elements
    from . import engine
    from .molecules import configure

    slab = read(args.slab, format="vasp").repeat(args.repeat)
    slab.set_constraint()
    # ACAT's existing workflow separates substrate and adsorbates by element.
    overlap = set(slab.get_chemical_symbols()) & set(adsorbate_elements)
    if overlap:
        raise ValueError(f"This version requires a clean metal slab; unsupported substrate elements: {sorted(overlap)}")
    if abs(slab.cell.volume) < 1e-6:
        raise ValueError("The surface must have a nonzero periodic cell")
    if any(abs(slab.cell[i, 2]) > 1e-6 for i in (0, 1)) or any(abs(slab.cell[2, i]) > 1e-6 for i in (0, 1)):
        raise ValueError("Orient the slab in the xy plane with the third cell vector along z")
    if slab.cell[2, 2] - slab.positions[:, 2].ptp() < 6:
        raise ValueError("Provide a slab with at least 6 angstrom vacuum along z")
    name = configure(engine, args.adsorbate, args.min_distance)
    retention = {key: label(getattr(args, key)) for key in ("parents", "screen_keep", "relax_keep", "final_keep")}
    output = args.output.resolve()
    config = {key: value for key, value in vars(args).items()
              if key not in {"command", "slab", "output", "generate_only", "workers", "threads"}}
    config.update(retention)
    config["slab_sha256"] = hashlib.sha256(args.slab.read_bytes()).hexdigest()
    config["engine_version"] = __version__
    config["versions"] = {package: version(package) for package in ("ase", "acat", "autoadsorbate", "mace-torch", "torch", "numpy")}
    output.mkdir(parents=True, exist_ok=True)
    lock = output / ".running.lock"
    try:
        handle = lock.open("x")
    except FileExistsError:
        raise RuntimeError(f"Output is locked: {lock}. If a previous run crashed, confirm it has stopped before removing this lock.") from None
    try:
        with handle:
            handle.write(str(os.getpid()))
        config_path = output / "search.json"
        if config_path.exists():
            if json.loads(config_path.read_text()) != config:
                raise ValueError("Output belongs to a different configuration or dependency version; choose a new --output directory")
        else:
            if any(item.name != lock.name for item in output.iterdir()):
                raise ValueError("Output is not empty and has no search.json; choose a new --output directory")
            config_path.write_text(json.dumps(config, indent=2) + "\n")
            write(output / "slab.vasp", slab, format="vasp", sort=False)
        metal = slab[0].symbol
        sas, mode = engine.build_adsorption_sites(slab, args.surface, metal)
        parents = [slab]
        for number in range(1, args.count + 1):
            directory = output / f"N_{number}"
            metadata_path = directory / "metadata.json"
            if metadata_path.is_file() and json.loads(metadata_path.read_text()).get("status") == "complete":
                if not (directory / "BEST.vasp").is_file():
                    raise RuntimeError(f"Completed round is missing BEST.vasp: {directory}")
                parents = read(directory / "growth_parents.traj", index=":")
                if not parents:
                    raise RuntimeError(f"Completed round has no growth parents: {directory}")
                print(f"[Reuse] N={number}: {len(parents)} parents", flush=True)
                continue
            directory.mkdir(exist_ok=True)
            candidates, generation = engine.generate_acat_children(parents, sas, number, len(slab), name)
            engine.write_candidate_pool(directory / "candidate_pool.traj", candidates)
            engine.write_xyz_pool(directory / "candidate_pool.xyz", candidates)
            metadata = {"status": "candidate_pool_only" if args.generate_only else "running",
                        "number": number, "adsorbate": args.adsorbate, "acat_surface_mode": mode,
                        "retention": retention, "generation": generation}
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            if args.generate_only:
                print(f"[Generated] N={number}: {len(candidates)} candidates in {directory}", flush=True)
                return
            rows, best, parents, selection = engine.run_mace_search(
                sas, candidates, directory, len(slab), name,
                args.screen_keep, args.relax_keep, args.final_keep, args.parents,
                args.pre_relax_steps, args.fmax, args.max_steps,
                "float32", "float64", args.workers, args.threads,
            )
            engine.write_ranking(directory / "ranking.csv", rows)
            metadata.update(status="complete", selection=selection,
                            n_final_results=selection["final"]["n_retained_results"], n_growth_parents=len(parents))
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            print(f"[Complete] N={number}: saved {metadata['n_final_results']} final results, {len(parents)} growth parents; {directory / 'BEST.vasp'}", flush=True)
    finally:
        lock.unlink()

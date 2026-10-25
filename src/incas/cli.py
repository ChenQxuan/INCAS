"""A single command for local deployment; heavy dependencies load only on run."""

import argparse
import math
from pathlib import Path

from .molecules import CATALOGUE
from .retention import parse_retention


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def keep_value(value):
    if value.lower() == "all":
        return None
    return positive(value)


def stage_value(value):
    try:
        return parse_retention(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def parser():
    p = argparse.ArgumentParser(description="INCAS incremental coverage adsorbate search on catalytic metal surfaces")
    commands = p.add_subparsers(dest="command", required=True)
    commands.add_parser("molecules", help="List available molecular templates")
    run = commands.add_parser("run", help="Search from a clean POSCAR or CONTCAR")
    run.add_argument("slab", type=Path)
    run.add_argument("--ads", dest="adsorbate", choices=sorted(CATALOGUE), required=True)
    run.add_argument("--n", dest="count", type=positive, required=True, help="Grow from 1 to this number of molecules")
    run.add_argument("--parents", type=keep_value, default=1, metavar="K|all", help="Previous-coverage final results used to grow each N (default: 1)")
    run.add_argument("--screen-keep", type=stage_value, default="80%", metavar="K|P%|all", help="Single-point-ranked unique candidates entering coarse relaxation (default: 80%%)")
    run.add_argument("--relax-keep", type=stage_value, default="50%", metavar="K|P%|all", help="Coarse-relaxed graph representatives entering full relaxation (default: 50%%)")
    run.add_argument("--final-keep", type=stage_value, default=None, metavar="K|P%|all", help="Converged unique final results to save, including BEST (default: all)")
    run.add_argument("--output", type=Path, default=Path("results"))
    run.add_argument("--surface", choices=["auto", "111", "100", "110", "211", "221", "311", "322", "331", "332"], default="auto", help="auto, or a known FCC facet")
    run.add_argument("--repeat", nargs=3, type=positive, default=[1, 1, 1], metavar=("NX", "NY", "NZ"))
    run.add_argument("--workers", type=positive, default=1)
    run.add_argument("--threads", type=positive, default=2)
    run.add_argument("--pre-relax-steps", type=positive, default=10)
    run.add_argument("--max-steps", type=positive, default=120)
    run.add_argument("--fmax", type=float, default=0.05)
    run.add_argument("--min-distance", type=float)
    run.add_argument("--generate-only", action="store_true", help="Generate the next unfinished round without MACE")
    return p


def validate(args, p):
    if args.pre_relax_steps >= args.max_steps:
        p.error("--pre-relax-steps must be less than --max-steps")
    if not math.isfinite(args.fmax) or args.fmax <= 0:
        p.error("--fmax must be positive and finite")
    if args.min_distance is not None and (not math.isfinite(args.min_distance) or args.min_distance <= 0):
        p.error("--min-distance must be positive and finite")
    if args.repeat[2] != 1:
        p.error("--repeat NZ must be 1 for a surface slab")
    if not args.slab.is_file():
        p.error(f"Slab file not found: {args.slab}")


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.command == "molecules":
        for public, (name, family) in CATALOGUE.items():
            print(f"{public:10s} {family:16s} binding template={name}")
        return 0
    validate(args, p)
    from .workflow import run
    try:
        run(args)
    except (ValueError, RuntimeError, FileNotFoundError) as error:
        p.exit(1, f"Error: {error}\n")
    return 0

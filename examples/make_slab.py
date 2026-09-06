"""Create a small clean Cu(111) POSCAR for the quick start."""
from pathlib import Path
from ase.build import fcc111
from ase.io import write

if __name__ == "__main__":
    output = Path(__file__).with_name("POSCAR")
    write(output, fcc111("Cu", size=(3, 3, 4), vacuum=10.0, periodic=True), format="vasp")
    print(output)

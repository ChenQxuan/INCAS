import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from ase.build import fcc111
from ase.io import write

from cascat import engine
from cascat.cli import main, parser
from cascat.molecules import CATALOGUE, configure
from cascat.retention import retention_count


class SearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.slab = fcc111("Cu", size=(3, 3, 4), vacuum=10, periodic=True)
        cls.sites, _ = engine.build_adsorption_sites(cls.slab, "111", "Cu")

    def test_catalogue_generates_and_recognizes_intact_molecules(self):
        for public in CATALOGUE:
            with self.subTest(molecule=public), contextlib.redirect_stdout(io.StringIO()):
                name = configure(engine, public)
                candidates, _ = engine.generate_acat_children([self.slab], self.sites, 1, len(self.slab), name)
                valid = []
                for candidate in candidates:
                    row = {}
                    engine.add_geometry_classification(row, candidate, self.sites, "test", len(self.slab), name)
                    valid.append(row["test_valid_geometry"])
                self.assertTrue(any(valid), f"No intact recognized candidate for {public}")

    def test_growth_keeps_more_than_two_without_energy_gap(self):
        rows = [{"relaxed_energy_eV": energy} for energy in [-10, -9, -8, -7]]
        self.assertEqual(engine.select_growth_rows(rows, 3, None), rows[:3])
        self.assertEqual(engine.select_growth_rows(rows, 8, None), rows)
        self.assertEqual(engine.select_growth_rows(rows, 3, 0.01), rows[:1])
        self.assertEqual(engine.select_growth_rows(rows, 1, None), rows[:1])
        self.assertEqual(engine.select_growth_rows(rows, None, None), rows)

    def test_stage_retention_for_energy_selection(self):
        # all must retain even candidates far above the former energy window.
        rows = [{"energy": energy} for energy in [20, -10, -5, 100]]
        for limit, expected in ((1, [-10]), (2, [-10, -5]), (10, [-10, -5, 20, 100]), ("50%", [-10, -5]), ("80%", [-10, -5, 20, 100]), (None, [-10, -5, 20, 100])):
            selected, metadata = engine.select_lowest_with_cap(rows, "energy", limit)
            self.assertEqual([row["energy"] for row in selected], expected)
            self.assertEqual(metadata["n_selected"], len(expected))

    def test_independent_retention_cli(self):
        arguments = ["run", "POSCAR", "--ads", "H", "--n", "2"]
        args = parser().parse_args(arguments + ["--parents", "2", "--screen-keep", "80%", "--relax-keep", "256", "--final-keep", "10"])
        self.assertEqual((args.parents, args.screen_keep, args.relax_keep, args.final_keep), (2, "80%", 256, 10))
        self.assertIsNone(parser().parse_args(arguments + ["--parents", "all"]).parents)
        for flag in ("--screen-keep", "--relax-keep", "--final-keep"):
            self.assertIsNone(getattr(parser().parse_args(arguments + [flag, "all"]), flag[2:].replace("-", "_")))
        for values in (["--parents", "50%"], ["--parents", "0"], ["--screen-keep", "0%"], ["--relax-keep", "101%"], ["--final-keep", "NaN%"], ["--final-keep", "1.5"], ["--keep", "2"]):
            with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                parser().parse_args(arguments + values)

    def test_percentage_rounding_and_separate_parent_limit(self):
        self.assertEqual(retention_count(11, "80%"), 9)
        self.assertEqual(retention_count(7, "50%"), 4)
        self.assertEqual(retention_count(100, "7%"), 7)
        self.assertEqual(retention_count(1, "0.1%"), 1)
        self.assertEqual(retention_count(0, "80%"), 0)
        rows = [{"relaxed_energy_eV": energy} for energy in [-10, -9, -8, -7]]
        final, _ = engine.select_lowest_with_cap(rows, "relaxed_energy_eV", "75%")
        self.assertEqual(len(final), 3)
        self.assertEqual(engine.select_growth_rows(final, 1, None), rows[:1])
        self.assertEqual(engine.select_growth_rows(final, None, None), rows[:3])
        self.assertEqual(engine.select_growth_rows(final, 10, None), rows[:3])

    def test_multidentate_two_molecule_growth(self):
        for public in ("CO2", "OOH", "H2O2"):
            with self.subTest(molecule=public), contextlib.redirect_stdout(io.StringIO()):
                name = configure(engine, public)
                first, _ = engine.generate_acat_children([self.slab], self.sites, 1, len(self.slab), name)
                children, _ = engine.generate_acat_children(first[:1], self.sites, 2, len(self.slab), name)
                valid = []
                for candidate in children:
                    row = {}
                    engine.add_geometry_classification(row, candidate, self.sites, "test", len(self.slab), name)
                    valid.append(row["test_valid_geometry"])
                self.assertTrue(any(valid))

    def test_resume_rejects_changed_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            slab = root / "POSCAR"
            write(slab, self.slab, format="vasp")
            arguments = ["run", str(slab), "--ads", "H", "--n", "2", "--surface", "111", "--output", str(root / "run"), "--generate-only"]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(arguments), 0)
                self.assertEqual(main(arguments), 0)
            with slab.open("a") as handle:
                handle.write("\n")
            with self.assertRaises(SystemExit) as error, contextlib.redirect_stderr(io.StringIO()):
                main(arguments)
            self.assertEqual(error.exception.code, 1)
            self.assertFalse((root / "run" / ".running.lock").exists())


if __name__ == "__main__":
    unittest.main()

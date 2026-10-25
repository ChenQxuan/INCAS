"""Check the skill's read-only reporter against complete and damaged outputs."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from ase import Atoms
from ase.io import write


CHECKER = Path(__file__).resolve().parents[1] / "skill/run-incas-search/scripts/inspect_results.py"
spec = importlib.util.spec_from_file_location("skill_result_checker", CHECKER)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class SkillResultTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.round = self.root / "N_1"
        self.round.mkdir()
        (self.root / "search.json").write_text(json.dumps({"count": 1, "engine_version": "0.5.0"}))
        self.atoms = Atoms("CuOH", positions=[[2, 2, 5], [2, 2, 7], [2, 2, 8]], cell=[10, 10, 15], pbc=True)
        self.atoms.info.update(candidate_id=1, relaxed_energy_eV=-10.0, converged=True)
        write(self.root / "slab.vasp", self.atoms[:1], format="vasp")
        write(self.round / "final_results.traj", [self.atoms])
        write(self.round / "growth_parents.traj", [self.atoms])
        # VASP can reorder adsorbate elements relative to the trajectory.
        write(self.round / "BEST.vasp", self.atoms[[0, 2, 1]], format="vasp")
        self.metadata = {"status": "complete", "n_final_results": 1, "n_growth_parents": 1,
                         "selection": {"initial_single_point": {"n_selected": 2},
                                       "full_relaxation_selection": {"n_selected": 1},
                                       "final": {"retained_candidate_ids": [1], "growth_parent_candidate_ids": [1]}}}
        (self.round / "metadata.json").write_text(json.dumps(self.metadata))

    def test_complete_with_reordered_best(self):
        self.assertTrue(checker.inspect(self.root)["complete"])

    def test_changed_best_is_rejected(self):
        changed = self.atoms.copy()
        changed.positions[-1, 0] += 0.5
        write(self.round / "BEST.vasp", changed, format="vasp")
        report = checker.inspect(self.root)
        self.assertFalse(report["complete"])
        self.assertIn("BEST geometry differs from the first final result", report["rounds"][0]["errors"])

    def test_incomplete_is_not_complete(self):
        self.metadata["status"] = "running"
        (self.round / "metadata.json").write_text(json.dumps(self.metadata))
        report = checker.inspect(self.root)
        self.assertFalse(report["complete"])
        self.assertEqual(report["rounds"][0]["status"], "running")

    def test_metadata_count_mismatch(self):
        self.metadata["n_final_results"] = 3
        (self.round / "metadata.json").write_text(json.dumps(self.metadata))
        self.assertFalse(checker.inspect(self.root)["complete"])

    def test_absent_tree_is_not_created(self):
        missing = self.root / "absent"
        self.assertFalse(checker.inspect(missing)["complete"])
        self.assertFalse(missing.exists())

import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from rdkit import Chem

from packages.science import design_docking
from packages.science.dual_e3 import REFERENCE_SOURCE, SOURCE
from packages.science.reference_parents import (
    load_reference_parent,
    reference_parent_catalog,
)
from packages.science.structures import atom_sites


MODULE_ROOT = Path(design_docking.__file__).resolve().parents[2]


class ReferenceParentDockingTests(unittest.TestCase):
    def _registered_ids(self):
        catalog = reference_parent_catalog(REFERENCE_SOURCE)
        return [
            row["id"]
            for row in catalog.get("available", [])
            if row.get("id") != "SMARCA2-FX5"
        ]

    def _registered_parent(self):
        ids = self._registered_ids()
        if not ids:
            self.skipTest(
                "No hash-verified nondefault reference parent is configured"
            )
        parent, metadata = load_reference_parent(ids[0], REFERENCE_SOURCE)
        return ids[0], parent, metadata

    def _protein(self):
        return atom_sites(SOURCE / "6HAZ.cif", ["A"])

    def test_real_registered_parent_passes_fixed_context(self):
        parent_id, parent, _ = self._registered_parent()
        protein = self._protein()

        fixed_parent, fixed_protein, receipt, hashes = (
            design_docking._validated_parent_context(
                parent, protein, parent_id, MODULE_ROOT
            )
        )

        self.assertEqual(fixed_parent.GetNumConformers(), 1)
        self.assertEqual(len(fixed_protein), len(protein))
        self.assertEqual(receipt["original_source_parent_id"], parent_id)
        self.assertEqual(receipt["receptor_frame"], "6HAZ chain A")
        self.assertEqual(
            receipt["technical_role"],
            "exploratory_pose_comparison_not_scientific_approval",
        )
        self.assertTrue(receipt["source_file_hashes"])
        self.assertTrue(hashes)

    def test_registered_parent_rejects_arbitrary_translation(self):
        parent_id, parent, _ = self._registered_parent()
        translated = Chem.Mol(parent)
        conformer = translated.GetConformer()
        for index in range(translated.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            conformer.SetAtomPosition(
                index, (point.x + 0.01, point.y, point.z)
            )

        with self.assertRaisesRegex(ValueError, "coordinate context"):
            design_docking._validated_parent_context(
                translated, self._protein(), parent_id, MODULE_ROOT
            )

    def test_registered_parent_rejects_spoofed_other_parent_id(self):
        ids = self._registered_ids()
        if len(ids) < 2:
            self.skipTest("Two registered nondefault parents are required")
        parent, _ = load_reference_parent(ids[0], REFERENCE_SOURCE)

        with self.assertRaisesRegex(ValueError, "coordinate context"):
            design_docking._validated_parent_context(
                parent, self._protein(), ids[1], MODULE_ROOT
            )

    def test_fake_caller_parent_id_cannot_bypass_registry(self):
        _, parent, _ = self._registered_parent()

        with self.assertRaisesRegex(KeyError, "Unknown reference parent id"):
            design_docking._validated_parent_context(
                parent,
                self._protein(),
                "SMARCA2-CALLER-SPOOF",
                MODULE_ROOT,
            )

    def test_registered_parent_rejects_mismatched_protein_identity(self):
        parent_id, parent, _ = self._registered_parent()
        protein = [dict(row) for row in self._protein()]
        protein[0]["type_symbol"] = (
            "N" if protein[0]["type_symbol"] != "N" else "C"
        )

        with self.assertRaisesRegex(
            ValueError, "DOCKING_PROTEIN_IDENTITY_MISMATCH"
        ):
            design_docking._validated_parent_context(
                parent, protein, parent_id, MODULE_ROOT
            )

    def test_registered_parent_rejects_mismatched_element(self):
        parent_id, parent, _ = self._registered_parent()
        changed = Chem.RWMol(parent)
        atom = changed.GetAtomWithIdx(0)
        atom.SetAtomicNum(7 if atom.GetAtomicNum() != 7 else 6)

        with self.assertRaisesRegex(ValueError, "coordinate context"):
            design_docking._validated_parent_context(
                changed.GetMol(), self._protein(), parent_id, MODULE_ROOT
            )

    def test_coordinate_comparison_rejects_wrong_shape_without_broadcasting(self):
        parent_id, parent, _ = self._registered_parent()
        shortened = Chem.RWMol(parent)
        shortened.RemoveAtom(shortened.GetNumAtoms() - 1)
        shortened_molecule = shortened.GetMol()

        self.assertFalse(
            design_docking._coordinates_match(shortened_molecule, parent)
        )
        with self.assertRaisesRegex(ValueError, "coordinate context"):
            design_docking._validated_parent_context(
                shortened_molecule,
                self._protein(),
                parent_id,
                MODULE_ROOT,
            )

    def test_reference_hash_mutation_is_detected(self):
        parent_id, parent, _ = self._registered_parent()

        with TemporaryDirectory() as temporary_directory:
            copied = Path(temporary_directory) / "reference_parents"
            shutil.copytree(REFERENCE_SOURCE, copied)
            target = copied / "parents" / f"{parent_id}.sdf"
            target.write_bytes(target.read_bytes() + b"\n")

            import packages.science.dual_e3 as dual_e3

            with mock.patch.object(dual_e3, "REFERENCE_SOURCE", copied):
                with self.assertRaisesRegex(
                    ValueError, "Manifest verification failed"
                ):
                    design_docking._validated_parent_context(
                        parent,
                        self._protein(),
                        parent_id,
                        MODULE_ROOT,
                    )

    def test_fx5_fixed_context_behavior_is_retained(self):
        parent = Chem.MolFromMolBlock(
            (SOURCE / "SMARCA2-neutral-design.sdf").read_text()
        )
        self.assertIsNotNone(parent)

        fixed_parent, _, receipt, _ = (
            design_docking._validated_parent_context(
                parent,
                self._protein(),
                "SMARCA2-FX5",
                MODULE_ROOT,
            )
        )

        self.assertEqual(
            design_docking._canonical_heavy(parent),
            design_docking._canonical_heavy(fixed_parent),
        )
        self.assertEqual(
            receipt["source_kind"], "original_neutral_default"
        )

    def test_empty_protected_maps_do_not_invent_a_scientific_core(self):
        parent = Chem.MolFromMolBlock(
            (SOURCE / "SMARCA2-neutral-design.sdf").read_text()
        )
        self.assertIsNotNone(parent)
        ligand = Chem.Mol(parent)
        ligand.RemoveAllConformers()

        self.assertEqual(
            design_docking._validate_inputs(ligand, parent, []), set()
        )

    def test_actual_docking_with_empty_protected_maps_has_no_rmsd_pass(self):
        bundled_vina = (
            MODULE_ROOT / "vendor" / "vina" / "vina_1.2.7_win.exe"
        )
        self.assertTrue(
            bundled_vina.is_file(),
            f"Bundled pinned Vina executable is missing: {bundled_vina}",
        )

        parent = Chem.MolFromMolBlock(
            (SOURCE / "SMARCA2-neutral-design.sdf").read_text()
        )
        self.assertIsNotNone(parent)
        ligand = Chem.Mol(parent)
        ligand.RemoveAllConformers()

        with TemporaryDirectory() as temporary_directory:
            result = design_docking.dock(
                ligand,
                parent,
                self._protein(),
                [],
                Path(temporary_directory),
                exhaustiveness=1,
            )

        self.assertIsInstance(result, dict)
        self.assertTrue(result.get("real_docking"))
        self.assertEqual(result.get("status"), "completed_with_limits")

        results = result.get("results", {})
        preservation = results.get("pose_preservation", {})
        self.assertEqual(
            preservation.get("status"),
            "review_protected_core_not_defined",
        )
        self.assertEqual(preservation.get("docking_pose_preserved"), "review")
        self.assertEqual(preservation.get("protected_maps"), [])

        poses = results.get("poses")
        self.assertIsInstance(poses, list)
        self.assertTrue(poses)
        for pose in poses:
            self.assertIsNone(pose.get("parent_protected_core_rmsd"))

        diagnostics = preservation.get("all_poses_diagnostics")
        self.assertIsInstance(diagnostics, list)
        self.assertTrue(diagnostics)
        for diagnostic in diagnostics:
            self.assertEqual(diagnostic.get("status"), "review")
            self.assertEqual(
                diagnostic.get("reason"), "PROTECTED_CORE_NOT_DEFINED"
            )
            self.assertIsNone(
                diagnostic.get("core_RMSD_A_in_receptor_frame")
            )

    def test_dock_interface_keeps_default_parent_id(self):
        self.assertEqual(
            design_docking.dock.__defaults__[-1], "SMARCA2-FX5"
        )


if __name__ == "__main__":
    unittest.main()
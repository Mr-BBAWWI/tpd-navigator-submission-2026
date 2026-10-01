"""Regression tests for actual immutable zero-based docking pose selection.

These synthetic fixtures exercise service bindings only. They are not expert
review or evidence of scientific performance.
"""
from __future__ import annotations

import copy
import tempfile
import unittest

from rdkit import Chem
from rdkit.Geometry import Point3D

from packages.contracts import ContractError
from packages.platform.scientific_acceptance import ScientificAcceptanceService


class ScientificPoseIndexTests(unittest.TestCase):
    mapped_smiles = "[CH3:1][CH3:2]"

    def setUp(self):
        self.service = ScientificAcceptanceService.__new__(ScientificAcceptanceService)
        self.poses = [self._pose(0.0), self._pose(5.0)]

    def _pose(self, offset):
        mol = Chem.MolFromSmiles(self.mapped_smiles)
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for index in range(mol.GetNumAtoms()):
            conformer.SetAtomPosition(index, Point3D(offset + index, 0.0, 0.0))
        mol.AddConformer(conformer)
        return mol

    def _row(self, index, status):
        return {
            "alignment_applied": False,
            "caller_supplied_score": -20.0 if index == 0 else -1.0,
            "common_atom_maps": [1, 2],
            "docking_pose_preserved": status == "pass",
            "pose_index_zero_based": index,
            "pose_origin": "caller supplied and not inferred or verified",
            "score_used_for_pass": False,
            "status": status,
        }

    def _source(self, rows):
        return {
            "real_docking": True,
            "results": {
                "pose_count": 2,
                "pose_preservation": {
                    "status": "pass",
                    "docking_pose_preserved": True,
                    "all_poses_diagnostics": rows,
                },
                "poses": [
                    {
                        "mode": 1,
                        "mapping_ambiguous": False,
                        "mapping_count": 1,
                        "selected_mapping": "retained_mapping",
                    },
                    {
                        "mode": 2,
                        "mapping_ambiguous": False,
                        "mapping_count": 1,
                        "selected_mapping": "retained_mapping",
                    },
                ],
            },
        }

    def _real_multimodel_pdbqt(self):
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mapped = "[CH3:1][CH2:5000][CH2:5001][OH:5002]"
        blocks = []
        for mode, offset in enumerate((0.0, 5.0), start=1):
            molecule = Chem.MolFromSmiles(mapped)
            conformer = Chem.Conformer(molecule.GetNumAtoms())
            for atom_index in range(molecule.GetNumAtoms()):
                conformer.SetAtomPosition(
                    atom_index, Point3D(offset + atom_index, 1.0, 2.0)
                )
            molecule.AddConformer(conformer)
            molecule = Chem.AddHs(molecule, addCoords=True)
            setups = MoleculePreparation().prepare(molecule)
            self.assertEqual(len(setups), 1)
            pdbqt, success, error = PDBQTWriterLegacy.write_string(setups[0])
            self.assertTrue(success, error)
            blocks.append(f"MODEL {mode}\n{pdbqt.rstrip()}\nENDMDL\n")
        return mapped, "".join(blocks).encode("utf-8")

    def test_real_meeko_multimodel_parser_restores_each_model_once(self):
        mapped, pdbqt = self._real_multimodel_pdbqt()
        restored = self.service._pdbqt_poses(pdbqt)

        self.assertEqual(len(restored), 2)
        self.assertTrue(all(mol.GetNumConformers() == 1 for mol in restored))
        for molecule in restored:
            maps = {atom.GetAtomMapNum() for atom in molecule.GetAtoms()}
            self.assertTrue({1, 5000, 5001, 5002}.issubset(maps))
        first = restored[0].GetConformer().GetAtomPosition(0)
        second = restored[1].GetConformer().GetAtomPosition(0)
        self.assertGreater(abs(second.x - first.x), 4.0)
        self.assertIn("5000", mapped)

    def test_multimodel_parser_rejects_invalid_records_and_sequences(self):
        _, pdbqt = self._real_multimodel_pdbqt()
        text = pdbqt.decode("utf-8")
        invalid_documents = {
            "duplicate": text.replace("MODEL 2", "MODEL 1", 1),
            "gap": text.replace("MODEL 2", "MODEL 3", 1),
            "mixed_identifier": text.replace("MODEL 2", "MODEL two", 1),
            "extra_model_fields": text.replace("MODEL 2", "MODEL 2 score", 1),
            "unmatched_model": text.replace("ENDMDL", "", 1),
            "extra_record": "REMARK outside model\n" + text,
        }
        for name, document in invalid_documents.items():
            with self.subTest(name=name):
                with self.assertRaisesRegex(
                    ContractError, "SCIENTIFIC_INTERACTION_PDBQT_RESTORATION"
                ):
                    self.service._pdbqt_poses(document.encode("utf-8"))
        with self.assertRaisesRegex(
            ContractError, "SCIENTIFIC_INTERACTION_PDBQT_RESTORATION"
        ):
            self.service._pdbqt_poses(b"MODEL 1\n\xff\nENDMDL\n")

    def test_actual_shape_without_row_mode_selects_zero_based_zero(self):
        rows = [self._row(0, "pass"), self._row(1, "fail")]
        index, chosen, declared, pose, receipt = self.service._actual_passing_pose(
            self._source(rows), self.poses, {"mapped_smiles": self.mapped_smiles}
        )
        self.assertEqual(index, 0)
        self.assertNotIn("mode", chosen)
        self.assertEqual(declared["mode"], 1)
        self.assertIs(pose, self.poses[0])
        self.assertIsNone(receipt)

    def test_first_failed_pose_does_not_trigger_score_or_default_selection(self):
        rows = [self._row(0, "fail"), self._row(1, "pass")]
        index, chosen, declared, pose, receipt = self.service._actual_passing_pose(
            self._source(rows), self.poses, {"mapped_smiles": self.mapped_smiles}
        )
        self.assertEqual(index, 1)
        self.assertEqual(chosen["caller_supplied_score"], -1.0)
        self.assertEqual(declared["mode"], 2)

    def test_negative_bool_and_out_of_bounds_indices_are_rejected(self):
        for invalid in (-1, True, 2):
            with self.subTest(index=invalid):
                row = self._row(0, "pass")
                row["pose_index_zero_based"] = invalid
                with self.assertRaises(ContractError):
                    self.service._actual_passing_pose(
                        self._source([row, self._row(1, "fail")]),
                        self.poses,
                        {"mapped_smiles": self.mapped_smiles},
                    )

    def test_mapping_mode_and_graph_mismatches_are_not_fallbacks(self):
        rows = [self._row(1, "pass")]
        source = self._source(rows)
        source["results"]["poses"][1]["mode"] = 1
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_INTERACTION_POSE_MAPPING"):
            self.service._actual_passing_pose(
                source, self.poses, {"mapped_smiles": self.mapped_smiles}
            )

        wrong_graph = Chem.MolFromSmiles("[CH3:1][OH:2]")
        conformer = Chem.Conformer(wrong_graph.GetNumAtoms())
        for index in range(wrong_graph.GetNumAtoms()):
            conformer.SetAtomPosition(index, Point3D(index, 0.0, 0.0))
        wrong_graph.AddConformer(conformer)
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_INTERACTION_POSE_GRAPH"):
            self.service._actual_passing_pose(
                self._source([self._row(0, "pass")]),
                [wrong_graph, self.poses[1]],
                {"mapped_smiles": self.mapped_smiles},
            )

    def test_sdf_records_select_first_pass_at_zero_based_one_without_mutation(self):
        rows = [self._row(0, "fail"), self._row(1, "pass")]
        source = self._source(rows)
        original_source = copy.deepcopy(source)

        with tempfile.TemporaryDirectory() as directory:
            sdf_path = f"{directory}/poses.sdf"
            writer = Chem.SDWriter(sdf_path)
            for pose in self.poses:
                writer.write(pose)
            writer.close()

            poses = [pose for pose in Chem.SDMolSupplier(sdf_path, removeHs=False)
                     if pose is not None]
            self.assertEqual(len(poses), 2)
            index, chosen, declared, pose, receipt = self.service._actual_passing_pose(
                source, poses, {"mapped_smiles": self.mapped_smiles}
            )

            self.assertEqual(index, 1)
            self.assertEqual(chosen["caller_supplied_score"], -1.0)
            self.assertEqual(declared["mode"], 2)
            self.assertIsNone(receipt)
            coordinates = pose.GetConformer().GetPositions()
            self.assertAlmostEqual(coordinates[0][0], 5.0)
            self.assertAlmostEqual(coordinates[1][0], 6.0)

        self.assertEqual(source, original_source)

    def test_v2000_truncated_maps_restore_from_real_selected_pdbqt_pose(self):
        _, pdbqt = self._real_multimodel_pdbqt()
        restored = self.service._pdbqt_poses(pdbqt)
        mapped = Chem.MolToSmiles(restored[1], canonical=True)

        raw_for_v2000 = Chem.Mol(restored[1])
        for atom in raw_for_v2000.GetAtoms():
            value = atom.GetAtomMapNum()
            if value >= 1000:
                atom.SetAtomMapNum(int(str(value)[:3]))
        v2000 = Chem.MolToMolBlock(raw_for_v2000, forceV3000=False)
        raw = Chem.MolFromMolBlock(v2000, removeHs=False)
        self.assertIsNotNone(raw)
        original_coordinates = raw.GetConformer().GetPositions().copy()

        rows = [self._row(0, "fail"), self._row(1, "pass")]
        rows[0]["common_atom_maps"] = [1, 5000, 5001, 5002]
        rows[1]["common_atom_maps"] = [1, 5000, 5001, 5002]
        source = self._source(rows)

        index, chosen, declared, normalized, receipt = (
            self.service._actual_passing_pose(
                source,
                [Chem.Mol(restored[0]), raw],
                {"mapped_smiles": mapped},
                pdbqt_loader=lambda: pdbqt,
            )
        )

        self.assertEqual(index, 1)
        self.assertEqual(declared["mode"], 2)
        self.assertEqual(chosen["pose_index_zero_based"], 1)
        normalized_maps = {
            atom.GetAtomMapNum() for atom in normalized.GetAtoms()
        }
        self.assertTrue({1, 5000, 5001, 5002}.issubset(normalized_maps))
        raw_maps = {atom.GetAtomMapNum() for atom in raw.GetAtoms()}
        self.assertTrue({1, 500}.issubset(raw_maps))
        normalized_coordinates = normalized.GetConformer().GetPositions()
        for atom_index in range(raw.GetNumAtoms()):
            for axis in range(3):
                self.assertEqual(
                    normalized_coordinates[atom_index][axis],
                    original_coordinates[atom_index][axis],
                )
        self.assertEqual(receipt["status"], "explicit_atom_map_restoration")
        self.assertFalse(receipt["scientific_approval"])
        self.assertFalse(
            receipt["coordinate_preservation"]["alignment_applied"]
        )

    def test_truncated_maps_require_registered_pdbqt_source(self):
        mapped = "[CH3:1][CH2:5000][OH:5001]"
        raw = Chem.MolFromSmiles(mapped)
        conformer = Chem.Conformer(raw.GetNumAtoms())
        for atom_index in range(raw.GetNumAtoms()):
            conformer.SetAtomPosition(atom_index, Point3D(atom_index, 0.0, 0.0))
        raw.AddConformer(conformer)
        for atom in raw.GetAtoms():
            if atom.GetAtomMapNum() >= 1000:
                atom.SetAtomMapNum(500)
        row = self._row(0, "pass")
        row["common_atom_maps"] = [1, 5000, 5001]
        source = self._source([row])
        source["results"]["pose_count"] = 1
        source["results"]["poses"] = [source["results"]["poses"][0]]
        with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_INTERACTION_PDBQT_SOURCE_REQUIRED"):
            self.service._actual_passing_pose(
                source, [raw], {"mapped_smiles": mapped}
            )

    def test_unexpected_colliding_raw_maps_are_not_repaired(self):
        mapped = "[CH3:1][CH2:5000][CH2:5001][OH:5002]"
        restored = Chem.MolFromSmiles(mapped)
        conformer = Chem.Conformer(restored.GetNumAtoms())
        for atom_index in range(restored.GetNumAtoms()):
            conformer.SetAtomPosition(atom_index, Point3D(atom_index, 0.0, 0.0))
        restored.AddConformer(conformer)
        raw = Chem.Mol(restored)
        raw.GetAtomWithIdx(1).SetAtomMapNum(499)
        raw.GetAtomWithIdx(2).SetAtomMapNum(499)
        raw.GetAtomWithIdx(3).SetAtomMapNum(499)
        with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_INTERACTION_MAP_RESTORATION_AMBIGUOUS"):
            self.service._repair_actual_pose_maps(
                raw, restored, Chem.MolFromSmiles(mapped)
            )

    def test_coordinate_change_beyond_v2000_rounding_tolerance_is_rejected(self):
        mapped = "[CH3:1][CH2:5000][OH:5001]"
        restored = Chem.MolFromSmiles(mapped)
        conformer = Chem.Conformer(restored.GetNumAtoms())
        for atom_index in range(restored.GetNumAtoms()):
            conformer.SetAtomPosition(atom_index, Point3D(atom_index, 0.0, 0.0))
        restored.AddConformer(conformer)
        raw = Chem.Mol(restored)
        raw.GetConformer().SetAtomPosition(1, Point3D(1.0002, 0.0, 0.0))
        for atom in raw.GetAtoms():
            if atom.GetAtomMapNum() >= 1000:
                atom.SetAtomMapNum(500)
        with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_INTERACTION_MAP_RESTORATION_AMBIGUOUS"):
            self.service._repair_actual_pose_maps(
                raw, restored, Chem.MolFromSmiles(mapped)
            )


if __name__ == "__main__":
    unittest.main()

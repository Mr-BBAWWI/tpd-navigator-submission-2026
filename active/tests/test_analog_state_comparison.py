from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.science.analog_state_comparison import (
    _verify_heavy_graph_and_stereo,
    build_state_set,
)
from scripts.export_analog_states import _poses, _repair_pose_maps


SMILES = {
    "W-c2afc5e73c1a": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][OH:5001])[CH2:8][CH2:11]3)[c:6]([NH2:17])"
        "[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
    "W-4c0a639c0a41": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][CH2:5001][CH2:5002][OH:5003])[CH2:8][CH2:11]3)"
        "[c:6]([NH2:17])[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
    "W-80f8f4a11b5d": (
        "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12]"
        "[N:19]([CH2:5000][NH2:5001])[CH2:8][CH2:11]3)[c:6]([NH2:17])"
        "[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
    ),
}


def conformer(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    status = AllChem.EmbedMolecule(mol, randomSeed=71)
    if status != 0:
        raise AssertionError("test conformer embedding failed")
    return mol


def parent_from_pose(pose):
    rw = Chem.RWMol(Chem.Mol(pose))
    remove = sorted(
        [atom.GetIdx() for atom in rw.GetAtoms() if atom.GetAtomMapNum() >= 5000],
        reverse=True,
    )
    for index in remove:
        rw.RemoveAtom(index)
    parent = rw.GetMol()
    Chem.SanitizeMol(parent)
    return parent


def protein():
    return [{
        "label_asym_id": "A",
        "auth_seq_id": "1408",
        "label_comp_id": "VAL",
        "label_atom_id": "CG1",
        "type_symbol": "C",
        "xyz": [100.0, 100.0, 100.0],
    }]


def recursive_keys(value):
    if isinstance(value, dict):
        result = set(value)
        for item in value.values():
            result |= recursive_keys(item)
        return result
    if isinstance(value, list):
        result = set()
        for item in value:
            result |= recursive_keys(item)
        return result
    return set()


class AnalogStateComparisonTests(unittest.TestCase):
    def build(self, analog_id):
        pose = conformer(SMILES[analog_id])
        parent = parent_from_pose(pose)
        return pose, parent, build_state_set(analog_id, pose, parent, protein())

    def test_exact_state_counts_and_charges(self):
        expected = {
            "W-c2afc5e73c1a": [0, 1],
            "W-4c0a639c0a41": [0, 1],
            "W-80f8f4a11b5d": [0, 1, 1, 2],
        }
        for analog_id, charges in expected.items():
            with self.subTest(analog_id=analog_id):
                _, _, report = self.build(analog_id)
                self.assertEqual(len(report["states"]), len(charges))
                self.assertEqual(
                    [row["formal_charge"] for row in report["states"]], charges
                )
                self.assertIsNone(report["selected_state"])
                self.assertEqual(report["microstate_gate"], "pending")
                self.assertFalse(report["tautomer_enumeration_performed"])
                json.dumps(report, allow_nan=False)

    def test_only_requested_nitrogens_change(self):
        _, _, report = self.build("W-80f8f4a11b5d")
        source = report["source_mapped_nitrogens"]
        for state in report["states"]:
            changed = set(state["protonated_atom_maps"])
            current = state["mapped_nitrogens"]
            for atom_map in (7, 16, 17, 18, 19, 5001):
                key = str(atom_map)
                if atom_map not in changed:
                    self.assertEqual(current[key], source[key])
            for atom_map in changed:
                key = str(atom_map)
                self.assertEqual(
                    current[key]["formal_charge"],
                    source[key]["formal_charge"] + 1,
                )
                self.assertEqual(
                    current[key]["hydrogen_count"],
                    source[key]["hydrogen_count"] + 1,
                )
        self.assertEqual(source["7"], report["states"][-1]["mapped_nitrogens"]["7"])
        self.assertEqual(source["16"], report["states"][-1]["mapped_nitrogens"]["16"])
        self.assertEqual(source["18"], report["states"][-1]["mapped_nitrogens"]["18"])

    def test_heavy_coordinates_are_preserved_without_alignment(self):
        _, _, report = self.build("W-c2afc5e73c1a")
        self.assertEqual(report["coordinate_frame"],
                         "shared input frame; no alignment or realignment applied")
        for state in report["states"]:
            self.assertIn(state["status"], {"computed", "optimization_review"})
            self.assertTrue(state["heavy_coordinates_exactly_preserved"])
            self.assertEqual(state["heavy_coordinate_max_displacement_A"], 0.0)

    def test_unknown_analog_and_invalid_maps_are_rejected(self):
        pose = conformer(SMILES["W-c2afc5e73c1a"])
        parent = parent_from_pose(pose)
        with self.assertRaisesRegex(ValueError, "unknown analog"):
            build_state_set("not-an-analog", pose, parent, protein())

        invalid = Chem.Mol(pose)
        atoms = [atom for atom in invalid.GetAtoms() if atom.GetAtomicNum() > 1]
        atoms[1].SetAtomMapNum(atoms[0].GetAtomMapNum())
        with self.assertRaisesRegex(ValueError, "unique positive"):
            build_state_set("W-c2afc5e73c1a", invalid, parent, protein())

    def test_source_pose_graph_mismatch_is_rejected(self):
        pose = conformer(SMILES["W-c2afc5e73c1a"])
        parent = parent_from_pose(pose)
        rw = Chem.RWMol(pose)
        mapping = {atom.GetAtomMapNum(): atom.GetIdx() for atom in rw.GetAtoms()}
        bond = rw.GetBondBetweenAtoms(mapping[5000], mapping[5001])
        rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        mismatch = rw.GetMol()
        with self.assertRaisesRegex(ValueError, "graph does not match"):
            build_state_set("W-c2afc5e73c1a", mismatch, parent, protein())

    def test_parent_core_topology_mismatch_is_rejected(self):
        pose = conformer(SMILES["W-4c0a639c0a41"])
        parent = parent_from_pose(pose)
        rw = Chem.RWMol(parent)
        mapping = {atom.GetAtomMapNum(): atom.GetIdx() for atom in rw.GetAtoms()}
        bond = rw.GetBondBetweenAtoms(mapping[8], mapping[19])
        rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        broken = rw.GetMol()
        with self.assertRaisesRegex(ValueError, "parent-core bond topology"):
            build_state_set("W-4c0a639c0a41", pose, broken, protein())

    def test_no_fake_population_or_pka_fields(self):
        _, _, report = self.build("W-c2afc5e73c1a")
        keys = recursive_keys(report)
        self.assertNotIn("population", keys)
        self.assertNotIn("pKa", keys)
        self.assertFalse(report["population_prediction_performed"])
        self.assertFalse(report["pKa_prediction_performed"])
        self.assertFalse(report["scientific_approved"])

    def test_v3000_state_exports_preserve_maps_graph_and_coordinates(self):
        expected_counts = {
            "W-c2afc5e73c1a": 2,
            "W-4c0a639c0a41": 2,
            "W-80f8f4a11b5d": 4,
        }
        with tempfile.TemporaryDirectory() as temporary:
            for analog_id, expected_count in expected_counts.items():
                with self.subTest(analog_id=analog_id):
                    pose = conformer(SMILES[analog_id])
                    parent = parent_from_pose(pose)
                    output = Path(temporary) / analog_id
                    report = build_state_set(
                        analog_id, pose, parent, protein(), output_dir=output
                    )
                    self.assertEqual(len(report["states"]), expected_count)
                    self.assertEqual(
                        sum(len(state["artifacts"]) for state in report["states"]),
                        expected_count * 2,
                    )
                    for state in report["states"]:
                        self.assertEqual(
                            set(state["artifacts"]), {"before_sdf", "after_sdf"}
                        )
                        for relative in state["artifacts"].values():
                            raw = (output / relative).read_bytes()
                            self.assertIn(b"V3000", raw)
                            loaded = list(Chem.ForwardSDMolSupplier(
                                io.BytesIO(raw), removeHs=False,
                                sanitize=True, strictParsing=True,
                            ))
                            self.assertEqual(len(loaded), 1)
                            self.assertIsNotNone(loaded[0])
                            maps = {
                                atom.GetAtomMapNum() for atom in loaded[0].GetAtoms()
                                if atom.GetAtomicNum() > 1
                            }
                            self.assertIn(5000, maps)
                            self.assertIn(5001, maps)
                            xyz = np.asarray(
                                loaded[0].GetConformer().GetPositions(), dtype=float
                            )
                            self.assertTrue(np.isfinite(xyz).all())

    def test_real_carbon_stereo_inversion_is_rejected(self):
        source = Chem.MolFromSmiles("[F:1][C@@:2]([Cl:3])([Br:4])[CH3:5]")
        self.assertIsNotNone(source)
        inverted = Chem.Mol(source)
        center = next(
            atom for atom in inverted.GetAtoms() if atom.GetAtomMapNum() == 2
        )
        center.SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            if center.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
            else Chem.ChiralType.CHI_TETRAHEDRAL_CW
        )
        with self.assertRaisesRegex(ValueError, "explicit source stereo"):
            _verify_heavy_graph_and_stereo(source, inverted, "regression")

    def test_same_index_service_helper_restores_v2000_maps(self):
        expected = conformer(SMILES["W-c2afc5e73c1a"])
        block = Chem.MolToMolBlock(expected, forceV3000=False) + "\n$$$$\n"
        raw_poses = _poses(block.encode("utf-8") * 5)
        self.assertIn(500, {
            atom.GetAtomMapNum() for atom in raw_poses[0].GetAtoms()
            if atom.GetAtomicNum() > 1
        })
        restored = [Chem.Mol(expected) for _ in range(5)]
        service = object.__new__(ScientificAcceptanceService)
        repaired, receipts = _repair_pose_maps(
            service, raw_poses, restored,
            Chem.MolToSmiles(expected, canonical=True, isomericSmiles=True),
        )
        self.assertEqual(len(repaired), 5)
        self.assertEqual(len(receipts), 5)
        for index, (pose, receipt) in enumerate(zip(repaired, receipts)):
            maps = [
                atom.GetAtomMapNum() for atom in Chem.RemoveHs(pose).GetAtoms()
                if atom.GetAtomicNum() > 1
            ]
            self.assertEqual(len(maps), len(set(maps)))
            self.assertIn(5000, maps)
            self.assertIn(5001, maps)
            self.assertEqual(receipt["pose_index_zero_based"], index)
            self.assertEqual(
                receipt["correspondence"],
                "same_zero_based_sdf_record_and_pdbqt_model",
            )
            self.assertFalse(receipt["nearest_pose_guess_used"])
            self.assertFalse(receipt["alignment_applied"])
            self.assertTrue(any(
                row["restoration"] == "v2000_first_three_digits"
                for row in receipt["full_atom_map_mapping"]
            ))


if __name__ == "__main__":
    unittest.main()

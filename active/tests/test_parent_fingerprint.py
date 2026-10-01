import io
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from packages.science.parent_fingerprint import (
    _contacts,
    _protein_index,
    _svg,
    build_parent_fingerprint,
    transfer_heavy_coordinates,
)


def _mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for index, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetAtomMapNum(index)
    embedded = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(embedded, randomSeed=7) != 0:
        raise AssertionError("fixture embedding failed")
    return Chem.RemoveHs(embedded)


def _protein():
    return [
        {"label_asym_id": "A", "auth_seq_id": "10", "label_comp_id": "ALA", "label_atom_id": "CB", "type_symbol": "C", "xyz": [0.0, 0.0, 0.0]},
    ]


def _read_sdf(path):
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(Path(path).read_bytes()), removeHs=False)
    mol = next(iter(supplier))
    if mol is None:
        raise AssertionError(f"failed to parse {path}")
    return mol


def _heavy_by_map(mol):
    conf = mol.GetConformer()
    return {
        atom.GetAtomMapNum(): np.asarray(conf.GetAtomPosition(atom.GetIdx()))
        for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    }


class ParentFingerprintTests(unittest.TestCase):
    def test_transfer_preserves_shared_geometry_for_charged_and_neutral_nitrogen(self):
        charged = _mol("C[NH2+]C")
        neutral = Chem.MolFromSmiles("CNC")
        for index, atom in enumerate(neutral.GetAtoms(), 1):
            atom.SetAtomMapNum(index)
        transferred = transfer_heavy_coordinates(neutral, charged)
        expected = _heavy_by_map(charged)
        actual = _heavy_by_map(transferred)
        self.assertEqual(set(expected), set(actual))
        for number in expected:
            self.assertTrue(np.array_equal(expected[number], actual[number]))
        self.assertEqual(Chem.GetFormalCharge(charged), 1)
        self.assertEqual(Chem.GetFormalCharge(transferred), 0)

    def test_invalid_duplicate_map_and_absent_3d_are_rejected(self):
        source = _mol("CCO")
        source.GetAtomWithIdx(1).SetAtomMapNum(source.GetAtomWithIdx(0).GetAtomMapNum())
        state = Chem.MolFromSmiles("CCO")
        for index, atom in enumerate(state.GetAtoms(), 1):
            atom.SetAtomMapNum(index)
        with self.assertRaisesRegex(ValueError, "unique positive"):
            transfer_heavy_coordinates(state, source)
        no_3d = Chem.MolFromSmiles("CCO")
        for index, atom in enumerate(no_3d.GetAtoms(), 1):
            atom.SetAtomMapNum(index)
        with self.assertRaisesRegex(ValueError, "actual 3D"):
            build_parent_fingerprint(no_3d, _mol("CCO"), _protein(), max_states=1)

    def test_svg_is_xml_and_does_not_change_original_xyz(self):
        mol = _mol("c1ccccc1O")
        before = np.asarray(mol.GetConformer().GetPositions()).copy()
        svg = _svg(mol)
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))
        self.assertTrue(np.array_equal(before, np.asarray(mol.GetConformer().GetPositions())))
        self.assertTrue(mol.GetConformer().Is3D())

    def test_unicode_export_raw_and_optimized_sdf_preserve_heavy_coordinates(self):
        charged, neutral = _mol("C[NH2+]C"), _mol("CN(C)C")
        expected = {"source_charged_parent": _heavy_by_map(charged), "neutral_design_parent": _heavy_by_map(neutral)}
        with tempfile.TemporaryDirectory() as base:
            destination = Path(base) / "한글 패킷"
            result = build_parent_fingerprint(charged, neutral, _protein(), output_dir=destination, max_states=1)
            for parent in result["parents"]:
                state = parent["states"][0]
                raw = _read_sdf(destination / state["files"]["raw_sdf"])
                optimized = _read_sdf(destination / state["files"]["hydrogen_optimized_sdf"])
                for exported in (raw, optimized):
                    observed = _heavy_by_map(exported)
                    self.assertEqual(set(expected[parent["parent_role"]]), set(observed))
                    for number, point in expected[parent["parent_role"]].items():
                        self.assertTrue(np.allclose(point, observed[number], atol=5.1e-5))
            ET.parse(destination / "charged-parent.svg")
            ET.parse(destination / "neutral-parent.svg")

    def test_aromatic_residue_expands_only_known_ring_atoms_within_radius(self):
        mol = _mol("c1ccccc1")
        ligand_map = mol.GetAtomWithIdx(0).GetAtomMapNum()
        point = np.asarray(mol.GetConformer().GetAtomPosition(0))
        records = [
            {"label_asym_id": "A", "auth_seq_id": "5", "label_comp_id": "PHE", "label_atom_id": "CG", "xyz": point.tolist()},
            {"label_asym_id": "A", "auth_seq_id": "5", "label_comp_id": "PHE", "label_atom_id": "CZ", "xyz": (point + [5.0, 0.0, 0.0]).tolist()},
            {"label_asym_id": "A", "auth_seq_id": "5", "label_comp_id": "PHE", "label_atom_id": "CA", "xyz": point.tolist()},
            {"label_asym_id": "A", "auth_seq_id": "5", "label_comp_id": "PHE", "label_atom_id": "CB", "xyz": point.tolist()},
        ]
        profile = {"interactions": [{"kind": "aromatic_contact_geometry", "participants": [f"L:{ligand_map}", "A:5:PHE"], "centroid_distance_A": 3.2}]}
        rows = _contacts(profile, mol, _protein_index(records), False)
        self.assertEqual([row["protein_atom"] for row in rows], ["CG"])
        self.assertTrue(rows[0]["chemical_role_ambiguity"])
        self.assertFalse(rows[0]["source_interaction"].get("confirmed_geometry", False))

    def test_no_protein_hydrogen_prevents_protein_donor_confirmation(self):
        mol = _mol("CO")
        ligand_map = mol.GetAtomWithIdx(1).GetAtomMapNum()
        point = np.asarray(mol.GetConformer().GetAtomPosition(1))
        protein = _protein_index([{"label_asym_id": "A", "auth_seq_id": "7", "label_comp_id": "SER", "label_atom_id": "OG", "xyz": point.tolist()}])
        item = {"kind": "directional_hbond", "participants": [f"L:{ligand_map}", "A:7:SER:OG"], "donor": "A:7:SER:OG", "distance_HA_A": 1.8, "angle_DHA_deg": 170.0}
        row = _contacts({"interactions": [item]}, mol, protein, False)[0]
        self.assertFalse(row["directional_hydrogen"]["confirmed_geometry"])
        self.assertEqual(row["directional_hydrogen"]["ambiguity"], "protein donor hydrogen absent")
        self.assertTrue(row["chemical_role_ambiguity"])

    def test_ligand_donor_geometry_can_be_reported_without_protein_hydrogen(self):
        mol = _mol("CO")
        ligand_map = mol.GetAtomWithIdx(1).GetAtomMapNum()
        point = np.asarray(mol.GetConformer().GetAtomPosition(1))
        protein = _protein_index([{"label_asym_id": "A", "auth_seq_id": "8", "label_comp_id": "ASP", "label_atom_id": "OD1", "xyz": point.tolist()}])
        item = {"kind": "directional_hbond", "participants": [f"L:{ligand_map}", "A:8:ASP:OD1"], "donor": f"L:{ligand_map}", "distance_HA_A": 1.9, "angle_DHA_deg": 165.0}
        row = _contacts({"interactions": [item]}, mol, protein, False)[0]
        self.assertTrue(row["directional_hydrogen"]["confirmed_geometry"])
        self.assertIsNone(row["directional_hydrogen"]["ambiguity"])

    def test_packet_distinguishes_parent_charge_and_disclaims_state_certainty(self):
        charged, neutral = _mol("C[NH2+]C"), _mol("CN(C)C")
        result = build_parent_fingerprint(charged, neutral, _protein(), max_states=2)
        self.assertEqual([parent["source_formal_charge"] for parent in result["parents"]], [1, 0])
        self.assertEqual(result["pH_context"], 7.4)
        self.assertFalse(result["pH_is_selection_evidence"])
        self.assertFalse(result["pKa_prediction_performed"])
        self.assertFalse(result["population_prediction_performed"])
        self.assertFalse(result["scientific_approved"])
        self.assertEqual(result["enumeration"], "truncated_nonexhaustive")
        for parent in result["parents"]:
            self.assertEqual(parent["proposed_protected_interaction_set"]["status"], "calculated_proposal_not_expert_approved")
            for state in parent["states"]:
                self.assertEqual(state["pKa"], "unknown")
                self.assertEqual(state["population"], "unknown")
                self.assertFalse(state["enumeration_exhaustive"])

    def test_saved_raw_is_optimizer_input_structure(self):
        source = _mol("CCO")
        captured = []

        def fake_optimize(mol, seed=23):
            raw = Chem.Mol(mol)
            raw.SetProp("raw_marker", "optimizer-input")
            optimized = Chem.AddHs(Chem.Mol(mol), addCoords=True)
            return {"input": raw, "hydrogens_initial": Chem.Mol(optimized), "hydrogens_optimized": optimized}, {"heavy_coordinates_exactly_preserved": True}

        def fake_profile(mol, protein_atoms, protein_hydrogens=None):
            return {"interactions": []}

        with tempfile.TemporaryDirectory() as base, patch("packages.science.parent_fingerprint.optimize_ligand_hydrogens", side_effect=fake_optimize), patch("packages.science.parent_fingerprint.interaction_profile", side_effect=fake_profile):
            result = build_parent_fingerprint(source, Chem.Mol(source), _protein(), output_dir=Path(base) / "한글", max_states=1)
            for parent in result["parents"]:
                raw_path = Path(base) / "한글" / parent["states"][0]["files"]["raw_sdf"]
                captured.append(_read_sdf(raw_path).GetProp("raw_marker"))
        self.assertEqual(captured, ["optimizer-input", "optimizer-input"])


if __name__ == "__main__":
    unittest.main()

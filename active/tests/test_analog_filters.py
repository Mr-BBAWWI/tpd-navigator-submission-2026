from __future__ import annotations

import unittest

from rdkit import Chem

from packages.science.analog_filters import cheap_filter


class CheapFilterBoundaryTests(unittest.TestCase):
    def test_invalid_valence_returns_invalid_without_exception(self) -> None:
        parent = Chem.MolFromSmiles("[CH4:1]")
        self.assertIsNotNone(parent)
        result = cheap_filter(parent, {"mapped_smiles": "[CH5:1]"})
        self.assertFalse(result["valid"])
        self.assertTrue(any(reason.startswith("PRODUCT_SANITIZATION_OR_VALENCE_FAILED")
                            for reason in result["hard_reasons"]))

    def test_missing_product_atom_map_is_rejected(self) -> None:
        parent = Chem.MolFromSmiles("[CH4:1]")
        self.assertIsNotNone(parent)
        result = cheap_filter(parent, {
            "mapped_smiles": "C",
            "removed_atom_maps": [1],
            "added_atom_maps": [],
        })
        self.assertFalse(result["valid"])
        self.assertIn("PRODUCT_ATOM_MAP_MISSING_OR_ZERO", result["hard_reasons"])

    def test_protected_boundary_bond_order_and_hydrogens_are_rejected(self) -> None:
        parent = Chem.MolFromSmiles("[CH3:1][CH3:2]")
        result = cheap_filter(parent, {
            "mapped_smiles": "[CH2:1]=[CH2:2]",
            "removed_atom_maps": [], "added_atom_maps": [],
            "protected_atom_maps": [1],
        })
        self.assertFalse(result["valid"])
        self.assertTrue(any(reason.startswith("PROTECTED_") for reason in result["hard_reasons"]))

    def test_protected_boundary_bond_removal_is_rejected(self) -> None:
        parent = Chem.MolFromSmiles("[CH3:1][CH3:2]")
        result = cheap_filter(parent, {
            "mapped_smiles": "[CH4:2]",
            "removed_atom_maps": [1], "added_atom_maps": [],
            "protected_atom_maps": [1],
        })
        self.assertFalse(result["valid"])
        self.assertIn("PROTECTED_MAP_REMOVED:1", result["hard_reasons"])

    def test_protected_boundary_stereo_change_is_rejected(self) -> None:
        parent = Chem.MolFromSmiles("[F:1]/[CH:2]=[CH:3]/[F:4]")
        result = cheap_filter(parent, {
            "mapped_smiles": "[F:1]/[CH:2]=[CH:3]\\[F:4]",
            "removed_atom_maps": [], "added_atom_maps": [],
            "protected_atom_maps": [1],
        })
        self.assertFalse(result["valid"])
        self.assertTrue(any("STEREOBOND" in reason or "PROTECTED_BOND" in reason
                            for reason in result["hard_reasons"]))

    def test_unprotected_edit_is_allowed(self) -> None:
        parent = Chem.MolFromSmiles("[CH3:1][CH3:2]")
        result = cheap_filter(parent, {
            "mapped_smiles": "[CH2:1]=[CH2:2]",
            "removed_atom_maps": [], "added_atom_maps": [],
            "protected_atom_maps": [],
        })
        self.assertTrue(result["valid"], result["hard_reasons"])

    def test_stable_protected_roundtrip_is_allowed(self) -> None:
        parent = Chem.MolFromSmiles("[CH3:1][CH3:2]")
        result = cheap_filter(parent, {
            "mapped_smiles": "[CH3:1][CH3:2]",
            "removed_atom_maps": [], "added_atom_maps": [],
            "protected_atom_maps": [1],
        })
        self.assertTrue(result["valid"], result["hard_reasons"])


if __name__ == "__main__":
    unittest.main()

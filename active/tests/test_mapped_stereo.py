from __future__ import annotations

import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from packages.science.analog_filters import cheap_filter
from packages.science.mapped_stereo import (
    MappedStereoError,
    mapped_tetrahedral_parity,
    same_mapped_tetrahedral_stereo,
)


class MappedStereoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        path = (
            Path(__file__).resolve().parents[1]
            / "cases/reference_parents/parents/SMARCA2-9D12-A1A1P.sdf"
        )
        if not path.is_file():
            raise AssertionError(f"required fixture is missing: {path}")
        mol = Chem.MolFromMolBlock(path.read_text(), removeHs=True)
        if mol is None:
            raise AssertionError(f"failed to load required fixture: {path}")
        if not any(atom.GetAtomMapNum() for atom in mol.GetAtoms()):
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
        maps = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
        if any(atom_map <= 0 for atom_map in maps) or len(set(maps)) != len(maps):
            raise AssertionError("9D12 fixture must have complete unique atom maps")
        centers = [
            atom.GetAtomMapNum()
            for atom in mol.GetAtoms()
            if atom.GetChiralTag() in {
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            }
        ]
        if not centers:
            raise AssertionError("9D12 fixture has no defined tetrahedral centers")
        if 26 not in centers:
            raise AssertionError("9D12 fixture does not define mapped stereocenter 26")
        cls.parent_9d12 = mol

    def _record(
        self, mol: Chem.Mol, protected: list[int] | None = None
    ) -> dict[str, object]:
        return {
            "mapped_smiles": Chem.MolToSmiles(
                mol, canonical=True, isomericSmiles=True
            ),
            "removed_atom_maps": [],
            "added_atom_maps": [],
            "protected_atom_maps": protected or [],
        }

    @staticmethod
    def _atom(mol: Chem.Mol, atom_map: int) -> Chem.Atom:
        return next(
            atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == atom_map
        )

    def test_synthetic_regression_survives_renumber_and_smiles_roundtrip(self) -> None:
        parent = Chem.MolFromSmiles("[F:1][C@H:26]([Cl:2])[CH2:18][OH:3]")
        self.assertIsNotNone(parent)
        expected = mapped_tetrahedral_parity(self._atom(parent, 26))
        reversed_mol = Chem.RenumberAtoms(
            parent, list(reversed(range(parent.GetNumAtoms())))
        )
        roundtrip = Chem.MolFromSmiles(Chem.MolToSmiles(
            reversed_mol, canonical=True, isomericSmiles=True
        ))
        self.assertIsNotNone(roundtrip)
        self.assertEqual(expected, mapped_tetrahedral_parity(self._atom(roundtrip, 26)))

    def test_synthetic_regression_is_valid_for_cheap_filter(self) -> None:
        parent = Chem.MolFromSmiles("[F:1][C@H:26]([Cl:2])[CH2:18][OH:3]")
        self.assertIsNotNone(parent)
        result = cheap_filter(parent, self._record(parent, [18, 26]))
        self.assertTrue(result["valid"], result["hard_reasons"])

    def test_shipped_9d12_canonical_and_reverse_renumbering_preserve_identity(self) -> None:
        parent = Chem.Mol(self.parent_9d12)
        canonical = Chem.MolFromSmiles(Chem.MolToSmiles(
            parent, canonical=True, isomericSmiles=True
        ))
        self.assertIsNotNone(canonical)
        reversed_mol = Chem.RenumberAtoms(
            parent, list(reversed(range(parent.GetNumAtoms())))
        )
        expected = mapped_tetrahedral_parity(self._atom(parent, 26))
        for product in (canonical, reversed_mol):
            with self.subTest(order=tuple(atom.GetAtomMapNum() for atom in product.GetAtoms())):
                self.assertEqual(
                    expected, mapped_tetrahedral_parity(self._atom(product, 26))
                )
                result = cheap_filter(parent, self._record(product, [26]))
                self.assertTrue(result["valid"], result["hard_reasons"])

    def test_shipped_9d12_inversion_at_map26_is_rejected(self) -> None:
        parent = Chem.Mol(self.parent_9d12)
        product = Chem.Mol(parent)
        self._atom(product, 26).InvertChirality()
        result = cheap_filter(parent, self._record(product))
        self.assertFalse(result["valid"])
        self.assertIn("PARENT_CARBON_STEREO_CHANGED:26", result["hard_reasons"])

    def test_shipped_9d12_stereo_loss_at_map26_is_rejected(self) -> None:
        parent = Chem.Mol(self.parent_9d12)
        product = Chem.Mol(parent)
        self._atom(product, 26).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        result = cheap_filter(parent, self._record(product))
        self.assertFalse(result["valid"])
        self.assertIn("PARENT_CARBON_STEREO_CHANGED:26", result["hard_reasons"])

    def test_implicit_and_addhs_removehs_parity_match(self) -> None:
        implicit = Chem.MolFromSmiles("[C@H:26]([F:1])([Cl:2])[Br:3]")
        self.assertIsNotNone(implicit)
        explicit = Chem.AddHs(implicit)
        removed = Chem.RemoveHs(explicit)
        expected = mapped_tetrahedral_parity(self._atom(implicit, 26))
        self.assertEqual(expected, mapped_tetrahedral_parity(self._atom(explicit, 26)))
        self.assertEqual(expected, mapped_tetrahedral_parity(self._atom(removed, 26)))

    def test_non_tetrahedral_stereo_fails_closed(self) -> None:
        mol = Chem.MolFromSmiles("[C:26]([F:1])([Cl:2])([Br:3])[I:4]")
        self.assertIsNotNone(mol)
        center = self._atom(mol, 26)
        center.SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
        with self.assertRaises(MappedStereoError):
            mapped_tetrahedral_parity(center)
        self.assertFalse(same_mapped_tetrahedral_stereo(center, center))

    def test_cip_change_preserves_mapped_parity(self) -> None:
        parent = Chem.MolFromSmiles(
            "[C@H:26]([CH2:2][F:5])([CH2:3][Cl:6])[Br:4]"
        )
        self.assertIsNotNone(parent)
        product = Chem.Mol(parent)
        self._atom(product, 5).SetAtomicNum(53)
        self._atom(product, 6).SetAtomicNum(9)
        Chem.SanitizeMol(product)
        before, after = self._atom(parent, 26), self._atom(product, 26)

        # Atom maps are identifiers, not chemical properties, and legacy CIP
        # ranking can include them. Calculate chemically correct labels on
        # unmapped copies with the accurate labeler and compare by atom index.
        # CIP alone must not be compared for stereo preservation: changing
        # ligand priorities can flip R/S while mapped physical parity is kept.
        cip_parent, cip_product = Chem.Mol(parent), Chem.Mol(product)
        for mol in (cip_parent, cip_product):
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
                if atom.HasProp("_CIPCode"):
                    atom.ClearProp("_CIPCode")
            rdCIPLabeler.AssignCIPLabels(mol)
        before_cip = cip_parent.GetAtomWithIdx(before.GetIdx()).GetProp("_CIPCode")
        after_cip = cip_product.GetAtomWithIdx(after.GetIdx()).GetProp("_CIPCode")
        self.assertNotEqual(before_cip, after_cip)
        self.assertEqual(
            mapped_tetrahedral_parity(before), mapped_tetrahedral_parity(after)
        )
        result = cheap_filter(parent, self._record(product, [26]))
        self.assertTrue(result["valid"], result["hard_reasons"])

    def test_unassigned_carbon_roundtrip_remains_valid(self) -> None:
        parent = Chem.MolFromSmiles("[CH3:1][CH2:2][OH:3]")
        self.assertIsNotNone(parent)
        product = Chem.RenumberAtoms(parent, [2, 0, 1])
        result = cheap_filter(parent, self._record(product, [1, 2, 3]))
        self.assertTrue(result["valid"], result["hard_reasons"])


if __name__ == "__main__":
    unittest.main()

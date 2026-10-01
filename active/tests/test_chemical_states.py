import json
import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from packages.science.chemical_states import (
    before_after_report,
    enumerate_microstates,
    interaction_profile,
    optimize_ligand_hydrogens,
)


class ChemicalStateTests(unittest.TestCase):
    def mapped(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        return mol

    def with_conformer(self, mol, points):
        self.assertEqual(mol.GetNumAtoms(), len(points))
        copy = Chem.Mol(mol)
        copy.RemoveAllConformers()
        conf = Chem.Conformer(copy.GetNumAtoms())
        conf.Set3D(True)
        for index, point in enumerate(points):
            conf.SetAtomPosition(index, point)
        copy.AddConformer(conf)
        return copy

    def protein_atom(self, xyz, residue="ASP", atom="OD1", element="O", **extra):
        result = {
            "xyz": xyz,
            "type_symbol": element,
            "label_comp_id": residue,
            "auth_seq_id": "10",
            "label_atom_id": atom,
            "label_asym_id": "A",
        }
        result.update(extra)
        return result

    def test_microstates_retain_source_and_reject_impossible_quaternary_protonation(self):
        source = self.mapped("C[N+](C)(C)C")
        source_smiles = Chem.MolToSmiles(source, isomericSmiles=True)
        states = enumerate_microstates(source, max_states=12)
        self.assertGreaterEqual(len(states), 1)
        self.assertEqual(states[0].GetProp("chemical_state_origin"), "source_ionization_state_retained")
        self.assertEqual(Chem.MolToSmiles(source, isomericSmiles=True), source_smiles)
        self.assertTrue(all(Chem.GetFormalCharge(state) == 1 for state in states))
        expected_maps = {1, 2, 3, 4, 5}
        for state in states:
            self.assertEqual(
                {a.GetAtomMapNum() for a in state.GetAtoms() if a.GetAtomicNum() > 1},
                expected_maps,
            )
            self.assertEqual(state.GetProp("chemical_state_population"), "unknown")
            self.assertEqual(state.GetProp("chemical_state_pKa"), "unknown")

    def test_microstates_include_explicit_acid_base_conjugate_without_population_claim(self):
        acid = self.mapped("CC(=O)O")
        states = enumerate_microstates(acid, max_states=16, pH=6.8)
        charges = {Chem.GetFormalCharge(state) for state in states}
        self.assertIn(0, charges)
        self.assertIn(-1, charges)
        self.assertTrue(all(state.GetProp("chemical_state_pH_context") == "6.8" for state in states))

    def test_amine_conjugates_change_formula_by_exactly_one_hydrogen(self):
        methylamine = self.mapped("CN")
        states = enumerate_microstates(methylamine, max_states=16)
        protonated = [s for s in states if Chem.GetFormalCharge(s) == 1]
        self.assertTrue(protonated)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(protonated[0]), "CH6N+")

        ammonium = self.mapped("C[NH3+]")
        states = enumerate_microstates(ammonium, max_states=16)
        neutral = [s for s in states if Chem.GetFormalCharge(s) == 0]
        self.assertTrue(neutral)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(neutral[0]), "CH5N")

    def test_explicit_hydrogen_before_donor_index_is_removed_safely(self):
        ammonium = Chem.AddHs(self.mapped("C[NH3+]"))
        order = list(range(ammonium.GetNumAtoms()))
        hydrogen = next(a.GetIdx() for a in ammonium.GetAtoms() if a.GetAtomicNum() == 1)
        order.remove(hydrogen)
        renumbered = Chem.RenumberAtoms(ammonium, [hydrogen] + order)
        states = enumerate_microstates(renumbered, max_states=16)
        neutral = [s for s in states if Chem.GetFormalCharge(s) == 0]
        self.assertTrue(neutral)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(neutral[0]), "CH5N")

    def test_neutral_amide_is_not_enumerated_as_an_anion(self):
        amide = self.mapped("CC(=O)N")
        states = enumerate_microstates(amide, max_states=16)
        self.assertNotIn(-1, {Chem.GetFormalCharge(state) for state in states})

    def test_bad_coordinate_input_is_rejected(self):
        mol = self.with_conformer(self.mapped("CO"), [(0, 0, 0), (float("nan"), 0, 0)])
        with self.assertRaises(ValueError):
            optimize_ligand_hydrogens(mol)
        with self.assertRaises(ValueError):
            interaction_profile(
                self.with_conformer(self.mapped("CO"), [(0, 0, 0), (1.4, 0, 0)]),
                [self.protein_atom([float("inf"), 0, 0])],
            )

    def test_hydrogen_optimization_preserves_heavy_coordinates_exactly(self):
        mol = self.mapped("CCO")
        status = AllChem.EmbedMolecule(mol, randomSeed=7)
        self.assertEqual(status, 0)
        before = np.asarray(mol.GetConformer().GetPositions(), dtype=float).copy()
        structures, receipt = optimize_ligand_hydrogens(mol)
        optimized = structures["hydrogens_optimized"]
        heavy = [a.GetIdx() for a in optimized.GetAtoms() if a.GetAtomicNum() > 1]
        after = np.asarray(optimized.GetConformer().GetPositions(), dtype=float)[heavy]
        self.assertTrue(np.array_equal(before, after))
        self.assertEqual(receipt["final_heavy_max_displacement_A"], 0.0)
        self.assertTrue(receipt["heavy_coordinates_exactly_preserved"])
        json.dumps(receipt, ensure_ascii=False, allow_nan=False)

    def test_missing_hydrogen_is_only_possible_contact(self):
        ligand = self.with_conformer(self.mapped("CO"), [(0, 0, 0), (1.4, 0, 0)])
        protein = [self.protein_atom([3.9, 0, 0])]
        profile = interaction_profile(ligand, protein)
        kinds = {item["kind"] for item in profile["interactions"]}
        self.assertIn("possible_hbond_contact", kinds)
        self.assertNotIn("directional_hbond", kinds)

    def explicit_methylamine(self):
        ligand = Chem.AddHs(self.mapped("CN"))
        conf = Chem.Conformer(ligand.GetNumAtoms())
        nitrogen = next(a.GetIdx() for a in ligand.GetAtoms() if a.GetAtomicNum() == 7)
        conf.SetAtomPosition(nitrogen, (0.0, 0.0, 0.0))
        for atom in ligand.GetAtoms():
            if atom.GetIdx() == nitrogen:
                continue
            if atom.GetAtomicNum() == 1 and atom.GetNeighbors()[0].GetIdx() == nitrogen:
                conf.SetAtomPosition(atom.GetIdx(), (1.0, 0.0, 0.0))
            else:
                conf.SetAtomPosition(atom.GetIdx(), (-1.5, atom.GetIdx() * 0.1, 0.0))
        ligand.RemoveAllConformers()
        ligand.AddConformer(conf)
        return ligand

    def test_explicit_hydrogen_bad_angle_rejects_directional_hbond(self):
        ligand = self.explicit_methylamine()
        protein = [self.protein_atom([1.0, 2.5, 0.0])]
        profile = interaction_profile(ligand, protein)
        kinds = [item["kind"] for item in profile["interactions"]]
        self.assertIn("directional_hbond_rejected", kinds)
        self.assertNotIn("directional_hbond", kinds)

    def test_explicit_hydrogen_good_angle_allows_directional_geometry(self):
        ligand = self.explicit_methylamine()
        protein = [self.protein_atom([2.7, 0.0, 0.0])]
        profile = interaction_profile(ligand, protein)
        directional = [i for i in profile["interactions"] if i["kind"] == "directional_hbond"]
        self.assertTrue(directional)
        self.assertGreaterEqual(directional[0]["angle_DHA_deg"], 120.0)

    def test_unspecified_histidine_role_is_marked_uncertain(self):
        ligand = self.with_conformer(self.mapped("C=O"), [(0, 0, 0), (1.2, 0, 0)])
        protein = [self.protein_atom([3.5, 0, 0], residue="HIS", atom="NE2", element="N")]
        profile = interaction_profile(ligand, protein)
        uncertain = [u for u in profile["uncertainties"] if u["kind"] == "uncertain_protein_chemical_role"]
        self.assertTrue(uncertain)
        possible = [i for i in profile["interactions"] if i["kind"] == "possible_hbond_contact"]
        self.assertTrue(possible)
        self.assertTrue(all(i["requires_review"] for i in possible))

    def test_before_after_uses_real_identifiers_without_efficacy_score(self):
        ligand = self.with_conformer(self.mapped("CO"), [(0, 0, 0), (1.4, 0, 0)])
        near = interaction_profile(ligand, [self.protein_atom([3.9, 0, 0])])
        far = interaction_profile(ligand, [self.protein_atom([20.0, 0, 0])])
        report = before_after_report(near, far)
        self.assertTrue(report["lost_identifiers"])
        self.assertEqual(report["gained_identifiers"], [])
        self.assertIsNone(report["quantitative_efficacy_score"])


if __name__ == "__main__":
    unittest.main()

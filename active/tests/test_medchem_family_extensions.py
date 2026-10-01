import unittest

from rdkit import Chem

from packages.science.analog_generation import generate


EXTENSION_RULES = {
    "BI_CARBONYL_O_TO_S",
    "BI_CARBONYL_TO_SULFONYL",
    "EV_RELOCATE_OH",
    "EV_RELOCATE_NH2",
    "HB_OH_TO_NH2",
    "HB_NH2_TO_OH",
    "HB_ETHER_TO_NH",
}


def mapped_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise AssertionError("test SMILES did not parse")
    for atom_map, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetAtomMapNum(atom_map)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def policy_sites(mol, protected_maps=(), allowed_rules=EXTENSION_RULES):
    protected_maps = set(protected_maps)
    sites = []
    for atom in mol.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        protected = atom_map in protected_maps
        sites.append({
            "atom_map": atom_map,
            "state": "PROTECTED" if protected else "MODIFIABLE",
            "evidence": {
                "SAR": {
                    "allowed_transformations": [] if protected else sorted(allowed_rules),
                },
                "scientific_policy": {
                    "rationale": "Unit-test policy fixture; not scientific approval.",
                },
            },
        })
    return sites


def constitutional_smiles(mapped_smiles):
    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        raise AssertionError("generated mapped SMILES did not parse")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


class MedicinalChemistryFamilyExtensionTests(unittest.TestCase):
    def test_phenol_and_ketone_generate_genuine_new_family_graphs(self):
        parent = mapped_molecule("CCCc1ccc(O)cc1C(=O)C")
        protected_map = parent.GetAtomWithIdx(0).GetAtomMapNum()
        result = generate(
            parent,
            policy_sites(parent, protected_maps={protected_map}),
            exploratory=False,
        )

        by_rule = {}
        for record in result["analogs"]:
            if record["rule_id"] in {
                "BI_CARBONYL_O_TO_S",
                "EV_RELOCATE_OH",
            }:
                by_rule.setdefault(record["rule_id"], record)

        self.assertEqual(
            set(by_rule),
            {"BI_CARBONYL_O_TO_S", "EV_RELOCATE_OH"},
        )
        graphs = {
            rule_id: constitutional_smiles(record["mapped_smiles"])
            for rule_id, record in by_rule.items()
        }
        self.assertEqual(len(set(graphs.values())), 2)

        for record in by_rule.values():
            self.assertTrue(record["protected_graph_preserved"])
            self.assertIn(protected_map, record["protected_atom_maps"])
            self.assertTrue(record["review_required"])
            self.assertTrue(record["qualified_for_counts"])
            self.assertEqual(
                record["original_site_map"][str(protected_map)],
                "PROTECTED",
            )

        qualified = result["summary"]["qualified"]
        self.assertGreaterEqual(qualified["constitutional_graph_count"], 2)
        self.assertIn("bioisosteric_replacement", qualified["transformation_families"])
        self.assertIn("exit_vector_relocation", qualified["transformation_families"])

    def test_aromatic_oh_relocation_preserves_roles_bonds_stereo_and_counts(self):
        families = [
            "Oc1ccccc1C(=O)C",
            "CCc1cc(O)ccc1C(=O)OC",
            "C[C@H](F)c1ccc(O)cc1C(=O)N",
        ]
        for smiles in families:
            with self.subTest(smiles=smiles):
                parent = mapped_molecule(smiles)
                parent_by_map = {
                    atom.GetAtomMapNum(): atom for atom in parent.GetAtoms()
                }
                hydroxyl = next(
                    atom for atom in parent.GetAtoms()
                    if atom.GetAtomicNum() == 8
                    and atom.GetDegree() == 1
                    and atom.GetTotalNumHs(includeNeighbors=True) == 1
                )
                source = hydroxyl.GetNeighbors()[0]
                result = generate(
                    parent,
                    policy_sites(parent, allowed_rules={"EV_RELOCATE_OH"}),
                    exploratory=False,
                )
                relocated = [
                    record for record in result["analogs"]
                    if record["rule_id"] == "EV_RELOCATE_OH"
                ]
                self.assertTrue(relocated)
                self.assertEqual(
                    len({record["constitutional_smiles"] for record in relocated}),
                    len(relocated),
                )

                for record in relocated:
                    product = Chem.MolFromSmiles(record["mapped_smiles"])
                    self.assertIsNotNone(product)
                    product_by_map = {
                        atom.GetAtomMapNum(): atom for atom in product.GetAtoms()
                    }
                    touched = set(record["modified_atom_maps"])
                    self.assertEqual(len(touched), 3)
                    self.assertIn(hydroxyl.GetAtomMapNum(), touched)
                    self.assertIn(source.GetAtomMapNum(), touched)
                    dest_maps = touched - {
                        hydroxyl.GetAtomMapNum(), source.GetAtomMapNum()
                    }
                    self.assertEqual(len(dest_maps), 1)
                    dest_map = next(iter(dest_maps))
                    original_dest = parent_by_map[dest_map]
                    self.assertTrue(original_dest.GetIsAromatic())
                    self.assertEqual(original_dest.GetAtomicNum(), 6)
                    self.assertEqual(original_dest.GetTotalNumHs(includeNeighbors=True), 1)

                    product_hetero = product_by_map[hydroxyl.GetAtomMapNum()]
                    product_source = product_by_map[source.GetAtomMapNum()]
                    product_dest = product_by_map[dest_map]
                    self.assertIsNone(product.GetBondBetweenAtoms(
                        product_hetero.GetIdx(), product_source.GetIdx()
                    ))
                    moved_bond = product.GetBondBetweenAtoms(
                        product_hetero.GetIdx(), product_dest.GetIdx()
                    )
                    self.assertIsNotNone(moved_bond)
                    self.assertEqual(moved_bond.GetBondType(), Chem.BondType.SINGLE)
                    self.assertEqual(product_source.GetTotalNumHs(includeNeighbors=True), 1)
                    self.assertEqual(product_dest.GetTotalNumHs(includeNeighbors=True), 0)
                    self.assertTrue(product_source.GetIsAromatic())
                    self.assertTrue(product_dest.GetIsAromatic())

                    for bond in parent.GetBonds():
                        begin_map = bond.GetBeginAtom().GetAtomMapNum()
                        end_map = bond.GetEndAtom().GetAtomMapNum()
                        if hydroxyl.GetAtomMapNum() in {begin_map, end_map}:
                            continue
                        product_bond = product.GetBondBetweenAtoms(
                            product_by_map[begin_map].GetIdx(),
                            product_by_map[end_map].GetIdx(),
                        )
                        self.assertIsNotNone(product_bond)
                        self.assertEqual(product_bond.GetBondType(), bond.GetBondType())
                        self.assertEqual(product_bond.GetStereo(), bond.GetStereo())

    def test_relocation_exact_rule_permission_does_not_authorize_unknown_sites(self):
        parent = mapped_molecule("Oc1ccccc1C(=O)C")
        sites = policy_sites(parent, allowed_rules={"EV_RELOCATE_OH"})
        destination = next(
            atom for atom in parent.GetAtoms()
            if atom.GetIsAromatic()
            and atom.GetAtomicNum() == 6
            and atom.GetTotalNumHs(includeNeighbors=True) == 1
        )
        sites[destination.GetIdx()]["state"] = "UNKNOWN"
        result = generate(parent, sites, exploratory=False)
        self.assertFalse(any(
            row["rule_id"] == "EV_RELOCATE_OH"
            and destination.GetAtomMapNum() in row["modified_atom_maps"]
            for row in result["analogs"]
        ))
        self.assertTrue(any(
            row["rule_id"] == "EV_RELOCATE_OH"
            and row["reason"] == "UNKNOWN_REQUIRES_EXPLORATORY"
            and destination.GetAtomMapNum() in row["unknown_atom_maps"]
            for row in result["rejections"]
        ))

    def test_acyclic_ether_rewiring_is_real_and_requires_exact_rule_permission(self):
        parent = mapped_molecule("CCOCC")
        protected_map = parent.GetAtomWithIdx(0).GetAtomMapNum()

        approved = generate(
            parent,
            policy_sites(parent, protected_maps={protected_map}),
            exploratory=False,
        )
        rewired = [
            record for record in approved["analogs"]
            if record["rule_id"] == "HB_ETHER_TO_NH"
        ]
        self.assertTrue(rewired)
        record = rewired[0]
        self.assertNotEqual(
            constitutional_smiles(record["mapped_smiles"]),
            constitutional_smiles(Chem.MolToSmiles(parent, isomericSmiles=True)),
        )
        self.assertTrue(record["protected_graph_preserved"])
        self.assertTrue(record["qualified_for_counts"])

        class_only_sites = policy_sites(
            parent,
            protected_maps={protected_map},
            allowed_rules={"hbond_rewiring"},
        )
        not_exactly_approved = generate(parent, class_only_sites, exploratory=False)
        self.assertFalse(any(
            row["rule_id"] == "HB_ETHER_TO_NH"
            for row in not_exactly_approved["analogs"]
        ))
        self.assertTrue(any(
            row["rule_id"] == "HB_ETHER_TO_NH"
            and row["reason"] == "UNKNOWN_REQUIRES_EXPLORATORY"
            for row in not_exactly_approved["rejections"]
        ))

    def test_outputs_have_unique_maps_and_no_stereo_only_count_inflation(self):
        parent = mapped_molecule("C[C@H](O)c1ccc(O)cc1C(=O)C")
        result = generate(parent, policy_sites(parent), exploratory=True)

        constitutional = []
        for record in result["analogs"]:
            product = Chem.MolFromSmiles(record["mapped_smiles"])
            self.assertIsNotNone(product)
            maps = [atom.GetAtomMapNum() for atom in product.GetAtoms()]
            self.assertTrue(all(atom_map > 0 for atom_map in maps))
            self.assertEqual(len(maps), len(set(maps)))
            constitutional.append(constitutional_smiles(record["mapped_smiles"]))

        self.assertEqual(len(constitutional), len(set(constitutional)))
        self.assertEqual(
            result["summary"]["all_passed"]["constitutional_graph_count"],
            len(result["analogs"]),
        )

    def test_unknown_exploration_does_not_become_qualified(self):
        parent = mapped_molecule("CCOCC")
        sites = [
            {"atom_map": atom.GetAtomMapNum(), "state": "UNKNOWN", "evidence": None}
            for atom in parent.GetAtoms()
        ]
        result = generate(parent, sites, exploratory=True)
        rewired = [
            record for record in result["analogs"]
            if record["rule_id"] == "HB_ETHER_TO_NH"
        ]
        self.assertTrue(rewired)
        self.assertTrue(all(record["scientific_pending"] for record in rewired))
        self.assertTrue(all(not record["qualified_for_counts"] for record in rewired))
        self.assertNotIn(
            "hbond_rewiring",
            result["summary"]["qualified"]["transformation_families"],
        )


if __name__ == "__main__":
    unittest.main()

import unittest
from unittest.mock import patch

from rdkit import Chem

from packages.science import medchem_sar
from packages.science.medchem_sar import link_selected_parent


BASE = "ClC1=CC=CC(N23)=C1C(N=C3C4(CCCCC4)C5=C2C=CC(C6CCNCC6)=C5)=O"


def mapped(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for index, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetAtomMapNum(index)
    return mol


def source(smiles=BASE, measurements=None, pairs=None):
    return {
        "status": "configured_hash_verified",
        "data": {
            "parent_catalog": {"records": [{
                "compound_id": "SMI-6080",
                "record_class": "measured_parent_ligand",
                "identity": {"source_smiles": smiles},
                "measurements": measurements if measurements is not None else [{
                    "target": "SMARCA2", "endpoint": "IC50", "relation": "=",
                    "value": 43.0, "unit": "nM", "conditions": {"buffer": "unknown"},
                }],
                "source": {"doi": "10.1021/acs.jmedchem.4c01903", "record_number_1_based_including_header": 5},
            }]},
            "sar_pair_evidence": pairs or [],
        },
    }


class MedchemSarTests(unittest.TestCase):
    def test_exact_graph_preserves_reported_measurement(self):
        result = link_selected_parent(mapped(BASE), "SMARCA2-9D12-A1A1P", source())
        self.assertEqual(result["status"], "exact_measured_parent_join")
        self.assertEqual(result["matched_record"]["compound_id"], "SMI-6080")
        self.assertEqual(result["matched_record"]["identity_status"], "connectivity_matched_stereochemistry_not_established")
        self.assertEqual(result["measurements"][0]["value"], 43.0)
        self.assertEqual(result["measurements"][0]["unit"], "nM")
        self.assertEqual(result["new_analog_measurement_status"], "not_proved_for_any_new_analog")

    def test_charge_mismatch_does_not_join(self):
        charged = mapped(BASE)
        nitrogen = next(atom for atom in charged.GetAtoms() if atom.GetSymbol() == "N" and not atom.GetIsAromatic())
        nitrogen.SetFormalCharge(1)
        nitrogen.SetNumExplicitHs(1)
        Chem.SanitizeMol(charged)
        result = link_selected_parent(charged, "charged", source())
        self.assertEqual(result["status"], "not_joined")
        self.assertEqual(result["measurements"], [])

    def test_different_graph_does_not_join(self):
        result = link_selected_parent(mapped(BASE.replace("Cl", "Br", 1)), "other", source())
        self.assertEqual(result["status"], "not_joined")
        self.assertIsNone(result["matched_record"])

    def test_no_measurement_is_invented(self):
        result = link_selected_parent(mapped(BASE), "parent", source(measurements=[]))
        self.assertEqual(result["status"], "exact_measured_parent_join")
        self.assertEqual(result["measurements"], [])
        self.assertEqual(result["new_analog_measurement_status"], "not_proved_for_any_new_analog")

    def test_neutral_amine_sar_sites_preserve_distinct_source_and_parent_maps(self):
        smiles = "CCN"
        source_probe = Chem.MolFromSmiles(smiles)
        ranks = list(Chem.CanonicalRankAtoms(
            source_probe, breakTies=True, includeChirality=True
        ))
        order = sorted(
            range(source_probe.GetNumAtoms()),
            key=lambda index: (ranks[index], index),
        )
        old_to_map = {old: new + 1 for new, old in enumerate(order)}
        source_nitrogen_map = old_to_map[next(
            atom.GetIdx() for atom in source_probe.GetAtoms()
            if atom.GetSymbol() == "N"
        )]

        parent = Chem.AddHs(Chem.MolFromSmiles(smiles))
        heavy = [atom.GetIdx() for atom in parent.GetAtoms() if atom.GetAtomicNum() > 1]
        hydrogens = [atom.GetIdx() for atom in parent.GetAtoms() if atom.GetAtomicNum() == 1]
        parent = Chem.RenumberAtoms(
            parent,
            [heavy[0], hydrogens[0], heavy[1], hydrogens[1], heavy[2]]
            + hydrogens[2:],
        )
        for atom in parent.GetAtoms():
            if atom.GetAtomicNum() > 1:
                atom.SetAtomMapNum(907 if atom.GetSymbol() == "N" else 800 + atom.GetIdx())

        pairs = [
            {
                "pair_id": "selected-left",
                "left_compound_id": "SMI-6080",
                "right_compound_id": "other-right",
                "mapping_possibilities": [{
                    "left_affected_parent_atom_maps": [source_nitrogen_map],
                }],
            },
            {
                "pair_id": "selected-right",
                "left_compound_id": "other-left",
                "right_compound_id": "SMI-6080",
                "mapping_possibilities": [{
                    "right_affected_atom_maps": [source_nitrogen_map],
                }],
            },
        ]

        result = link_selected_parent(parent, "neutral-amine", source(
            smiles=smiles, pairs=pairs
        ))
        self.assertEqual(result["status"], "exact_measured_parent_join")
        self.assertNotEqual(source_nitrogen_map, 907)
        by_pair = {entry["pair_id"]: entry for entry in result["source_sar"]}
        self.assertEqual(set(by_pair), {"selected-left", "selected-right"})
        for entry in by_pair.values():
            self.assertEqual(entry["source_affected_atom_maps"], [source_nitrogen_map])
            self.assertEqual(entry["selected_parent_affected_atom_maps"], [907])
            self.assertTrue(entry["source_site_consensus"])
            self.assertEqual(entry["status"], "informational_exact_observed_pair_only")
            self.assertFalse(entry["protected_atom_activation"])
        self.assertFalse(result["changes_design_policy"])
        self.assertEqual(
            result["matched_record"]["identity_status"],
            "connectivity_matched_stereochemistry_not_established",
        )

    def test_isomorphism_enumeration_cap_fails_closed(self):
        with patch.object(medchem_sar, "MAX_FULL_ISOMORPHISMS", 1):
            result = link_selected_parent(mapped(BASE), "parent", source())
        self.assertEqual(result["status"], "not_joined")
        self.assertEqual(result["reason"], "FULL_GRAPH_ISOMORPHISM_LIMIT_REACHED")
        self.assertIsNone(result["matched_record"])
        self.assertEqual(result["measurements"], [])
        self.assertEqual(result["source_sar"], [])

    def test_symmetric_mapping_uncertain_site_is_skipped(self):
        pair = {
            "pair_id": "x", "left_compound_id": "SMI-6080", "right_compound_id": "other",
            "mapping_possibilities": [
                {"left_affected_parent_atom_maps": [19]},
                {"left_affected_parent_atom_maps": [23]},
            ],
        }
        result = link_selected_parent(mapped(BASE), "parent", source(pairs=[pair]))
        self.assertEqual(result["source_sar"], [])
        self.assertFalse(result["changes_design_policy"])


if __name__ == "__main__":
    unittest.main()

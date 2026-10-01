from __future__ import annotations

import copy
import math
import unittest
import xml.etree.ElementTree as ET

from rdkit import Chem

from packages.science import synthesis_proposals as synthesis


MAPPED_SMILES = "[cH:1]1[c:2](-[c:4]2[cH:5][c:15]([N:18]3[CH2:9][CH2:12][N:19]([CH2:5000][O:5001][CH2:2001][CH2:2002][CH2:2003][CH2:2004][CH2:2005][CH2:2006][O:4037][c:4009]4[cH:4010][c:4011](-[c:4027]5[c:4013]([CH3:4030])[n:4028][cH:4014][s:4029]5)[cH:4026][cH:4012][c:4025]4[CH2:4024][NH:4023][C:4008]([C@H:4006]4[N:4017]([C:4002]([C@@H:4015]([NH:4001][C:4031]([C:4032]5([F:4034])[CH2:4035][CH2:4036]5)=[O:4033])[C:4016]([CH3:4003])([CH3:4004])[CH3:4005])=[O:4018])[CH2:4007][C@H:4020]([OH:4021])[CH2:4019]4)=[O:4022])[CH2:8][CH2:11]3)[c:6]([NH2:17])[n:16][n:7]2)[c:3]([OH:20])[cH:10][cH:13][cH:14]1"
CANONICAL_SMILES = "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)C2(F)CC2)C(C)(C)C)c(OCCCCCCOCN2CCN(c3cc(-c4ccccc4O)nnc3N)CC2)c1"


def actual_candidate() -> dict:
    return {
        "candidate_id": "D-bf12802fe7d0",
        "canonical_smiles": CANONICAL_SMILES,
        "mapped_smiles": MAPPED_SMILES,
        "atom_roles": {
            "warhead_maps": list(range(1, 21)) + [5000, 5001],
            "linker_maps": list(range(2001, 2007)),
            "recruiter_maps": list(range(4001, 4038)),
        },
        "attachment_metadata": {
            "warhead_linker_bond": {
                "attachment_atom_map": 5001,
                "partner_atom_map": 2001,
                "bond_type": "SINGLE",
            },
            "recruiter_linker_bond": {
                "attachment_atom_map": 4037,
                "partner_atom_map": 2006,
                "bond_type": "SINGLE",
            },
        },
        "route_evidence": {
            "reusable_precedent_records": [
                {
                    "record_id": "original-route-C01",
                    "source_record": {
                        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6600871/#SD1",
                        "locator": "Synthetic methods for the preparation of PROTAC 1 (2)",
                    },
                    "scope": "source compound only; not this candidate graph",
                }
            ]
        },
    }


def ester_candidate() -> dict:
    mapped = "[CH3:1][O:2][C:3](=[O:4])[CH2:5][CH2:6][NH2:7]"
    mol = Chem.MolFromSmiles(mapped)
    assert mol is not None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return {
        "candidate_id": "ester-fixture",
        "mapped_smiles": mapped,
        "canonical_smiles": canonical,
        "atom_roles": {
            "warhead_maps": [1, 2],
            "linker_maps": [3, 4, 5],
            "recruiter_maps": [6, 7],
        },
        "attachment_metadata": {
            "warhead_linker_bond": {
                "attachment_atom_map": 2,
                "partner_atom_map": 3,
                "bond_type": "SINGLE",
            },
            "recruiter_linker_bond": {
                "attachment_atom_map": 6,
                "partner_atom_map": 5,
                "bond_type": "SINGLE",
            },
        },
    }


class TestSynthesisProposals(unittest.TestCase):
    def test_actual_candidate_cuts_produce_exact_role_partitions(self):
        candidate = actual_candidate()
        result = synthesis.propose_synthesis(candidate)
        fragments = result["actual_fragment_products"]
        self.assertEqual(len(fragments), 3)
        self.assertEqual(
            {frozenset(row["atom_maps"]) for row in fragments},
            {
                frozenset(candidate["atom_roles"]["warhead_maps"]),
                frozenset(candidate["atom_roles"]["linker_maps"]),
                frozenset(candidate["atom_roles"]["recruiter_maps"]),
            },
        )
        self.assertEqual(
            {row["name"] for row in result["exact_attachment_cuts"]},
            {"warhead_linker_bond", "recruiter_linker_bond"},
        )

    def test_fragmentation_does_not_mutate_actual_chiral_source(self):
        candidate = actual_candidate()
        original = copy.deepcopy(candidate)
        result = synthesis.propose_synthesis(candidate)
        self.assertEqual(candidate, original)
        stereo = [
            item
            for fragment in result["actual_fragment_products"]
            for item in fragment["stereochemistry"]
        ]
        self.assertTrue(stereo)
        self.assertTrue(all(item["CIP"] in {"R", "S"} for item in stereo))
        self.assertTrue(
            any(
                "@" in fragment["mapped_smiles_with_attachment_dummies"]
                for fragment in result["actual_fragment_products"]
            )
        )

    def test_actual_competing_sites_and_ester_precursor_policy_are_flagged(self):
        actual = synthesis.propose_synthesis(actual_candidate())
        actual_flags = {
            row["flag"]
            for row in actual["functional_group_and_competing_site_checks"]
        }
        self.assertIn("phenol", actual_flags)
        self.assertIn("free_amine", actual_flags)
        self.assertIn(
            "multiple_competing_phenol_amine_or_alcohol_sites", actual_flags
        )

        ester = synthesis.propose_synthesis(ester_candidate())
        ester_flags = {
            row["flag"]
            for row in ester["functional_group_and_competing_site_checks"]
        }
        self.assertIn("ester", ester_flags)
        self.assertIn("ester_handle_policy", ester_flags)

    def test_absent_exact_route_remains_unverified_not_unsynthesizable(self):
        result = synthesis.propose_synthesis(
            actual_candidate(), exact_route_records=None
        )
        self.assertEqual(
            result["independent_synthesis_disposition"],
            "unverified_proposal_requires_review",
        )
        self.assertIs(result["synthesis_unverified"], True)
        self.assertIs(result["unsynthesizable"], False)
        self.assertIs(result["scientific_approved"], False)
        self.assertIs(result["expert_approved"], False)

    def test_forged_compatibility_boolean_cannot_approve_source(self):
        candidate = actual_candidate()
        candidate["route_evidence"] = {
            "reusable_precedent_records": [
                {
                    "candidate_specific_chemical_compatibility": True,
                    "source_record": {
                        "doi": "10.1000/forged",
                        "locator": "claimed exact paragraph",
                        "candidate_specific_chemical_compatibility": True,
                    },
                }
            ]
        }
        result = synthesis.propose_synthesis(candidate)
        sources = result["source_template_suggestions"]
        self.assertEqual(
            sources["status"],
            "unknown_no_graph_confirmed_chemically_compatible_template",
        )
        self.assertEqual(sources["suggestions"], [])

    def test_invalid_graph_roles_cut_maps_bond_types_and_canonical_are_rejected(self):
        mutations = []

        duplicate = actual_candidate()
        duplicate["mapped_smiles"] = duplicate["mapped_smiles"].replace(
            "[cH:1]", "[cH:2]", 1
        )
        mutations.append(duplicate)

        incomplete_roles = actual_candidate()
        incomplete_roles["atom_roles"]["warhead_maps"].remove(1)
        mutations.append(incomplete_roles)

        absent_cut_map = actual_candidate()
        absent_cut_map["attachment_metadata"]["warhead_linker_bond"][
            "attachment_atom_map"
        ] = 999999
        mutations.append(absent_cut_map)

        non_single = actual_candidate()
        non_single["attachment_metadata"]["warhead_linker_bond"][
            "bond_type"
        ] = "DOUBLE"
        mutations.append(non_single)

        wrong_canonical = actual_candidate()
        wrong_canonical["canonical_smiles"] = "CC"
        mutations.append(wrong_canonical)

        for position, candidate in enumerate(mutations):
            with self.subTest(position=position):
                with self.assertRaises(ValueError):
                    synthesis.propose_synthesis(candidate)

    def test_generated_scheme_is_valid_svg_with_three_components_and_map_labels(self):
        svg = synthesis.scheme_svg(actual_candidate())
        root = ET.fromstring(svg)
        self.assertEqual(root.tag.rsplit("}", 1)[-1], "svg")
        metadata = next(
            element
            for element in root.iter()
            if element.tag.rsplit("}", 1)[-1] == "metadata"
            and element.attrib.get("id") == "attachment-map-labels"
        )
        self.assertEqual(metadata.attrib["component-count"], "3")
        self.assertTrue(metadata.text)
        self.assertIn("5001", metadata.text)
        self.assertIn("2001", metadata.text)
        self.assertIn("4037", metadata.text)
        self.assertIn("2006", metadata.text)

    def test_sa_score_is_finite_when_available_and_never_approval(self):
        result = synthesis.propose_synthesis(actual_candidate())
        score = result["SA_Score"]
        if score["status"] == "available":
            self.assertIsInstance(score["value"], float)
            self.assertTrue(math.isfinite(score["value"]))
            self.assertIn("screening", score["interpretation"])
        else:
            self.assertEqual(score["status"], "unavailable")
            self.assertIsNone(score["value"])
        self.assertIs(result["scientific_approved"], False)
        self.assertIs(result["expert_approved"], False)


if __name__ == "__main__":
    unittest.main()

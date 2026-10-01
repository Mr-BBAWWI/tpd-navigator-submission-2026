"""Tests for the isolated native parent expansion v7."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

from packages.science.native_parent_expansion import (
    FRAGMENTS,
    NEW_FRAGMENTS,
    _actual_counts,
    _execute_case,
    audit_rule_applicability,
    enumerate_expansion_candidates,
)


def _parent() -> Chem.Mol:
    mol = Chem.MolFromSmiles("[CH3:1][NH:3][CH2:2][C@@H:4]([F:5])[Cl:6]")
    if mol is None:
        raise AssertionError("Test parent failed to parse")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _stereo(parent: Chem.Mol, atom_map: int) -> str:
    atom = next(atom for atom in parent.GetAtoms() if atom.GetAtomMapNum() == atom_map)
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else ""


class ExpansionGraphTests(unittest.TestCase):
    def test_twenty_unique_constitutional_graphs(self):
        records = enumerate_expansion_candidates(_parent())
        self.assertEqual(20, len(FRAGMENTS))
        self.assertEqual(8, len(NEW_FRAGMENTS))
        self.assertEqual(20, len(records))
        graphs = set()
        for record in records:
            mol = Chem.MolFromSmiles(record["mapped_smiles"])
            self.assertIsNotNone(mol)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            graphs.add(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False))
        self.assertEqual(20, len(graphs))

    def test_parent_maps_bonds_and_stereo_are_preserved(self):
        parent = _parent()
        records = enumerate_expansion_candidates(parent)
        original_atoms = {atom.GetAtomMapNum(): atom.GetAtomicNum() for atom in parent.GetAtoms()}
        original_bonds = {
            tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))): bond.GetBondType()
            for bond in parent.GetBonds()
        }
        expected_stereo = _stereo(parent, 4)
        for record in records:
            candidate = Chem.MolFromSmiles(record["mapped_smiles"])
            self.assertIsNotNone(candidate)
            Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
            by_map = {atom.GetAtomMapNum(): atom for atom in candidate.GetAtoms()}
            for atom_map, atomic_number in original_atoms.items():
                self.assertIn(atom_map, by_map)
                self.assertEqual(atomic_number, by_map[atom_map].GetAtomicNum())
            for pair, bond_type in original_bonds.items():
                bond = candidate.GetBondBetweenAtoms(by_map[pair[0]].GetIdx(), by_map[pair[1]].GetIdx())
                self.assertIsNotNone(bond)
                self.assertEqual(bond_type, bond.GetBondType())
            self.assertEqual(expected_stereo, _stereo(candidate, 4))

    def test_all_new_endpoints_are_valid_and_family_is_not_inflated(self):
        records = enumerate_expansion_candidates(_parent())
        families = {record["broad_family"] for record in records}
        self.assertEqual({"linker_handle_introduction"}, families)
        for record in records:
            mol = Chem.MolFromSmiles(record["mapped_smiles"])
            self.assertIsNotNone(mol)
            endpoint_map = max(record["added_atom_maps"])
            endpoint = next(atom for atom in mol.GetAtoms() if atom.GetAtomMapNum() == endpoint_map)
            self.assertIn(endpoint.GetAtomicNum(), (7, 8))
            self.assertEqual(0, endpoint.GetFormalCharge())
            self.assertGreaterEqual(endpoint.GetTotalNumHs(), 1)
            self.assertFalse(record["fragment_motif_counts_as_new_family"])
            self.assertFalse(record["strict_qualified_for_counts"])


class RuleAuditTests(unittest.TestCase):
    def test_catalog_availability_does_not_inflate_actual_counts(self):
        parent = Chem.MolFromSmiles("[CH3:1][OH:2]")
        sites = [
            {"atom_map": 1, "state": "UNKNOWN", "evidence": None},
            {"atom_map": 2, "state": "PROTECTED", "evidence": None},
        ]
        calls = []

        def generated(_parent, supplied_sites, exploratory=False):
            calls.append((exploratory, json.dumps(supplied_sites, sort_keys=True)))
            analogs = []
            if exploratory:
                analogs.append({
                    "candidate_id": "fixture-pass",
                    "mapped_smiles": "[CH3:1][O:2][CH3:5000]",
                    "rule_id": "LH_HYDROXYETHYL",
                    "transformation_class": "linker_handle_introduction",
                    "attachment_site_atom_maps": [1],
                    "modified_atom_maps": [2],
                    "removed_atom_maps": [3],
                    "touched_atom_maps": [4],
                    "added_atom_maps": [5000],
                    "protected_atom_maps": [],
                    "qualified_for_counts": False,
                })
                analogs.append({
                    "candidate_id": "fixture-fail",
                    "mapped_smiles": "[CH3:1][O:2][CH2:5001][OH:5002]",
                    "rule_id": "LH_HYDROXYETHYL",
                    "transformation_class": "linker_handle_introduction",
                    "attachment_site_atom_maps": [1],
                    "modified_atom_maps": [2],
                    "removed_atom_maps": [3],
                    "touched_atom_maps": [4],
                    "protected_atom_maps": [],
                    "qualified_for_counts": False,
                })
            return {
                "analogs": analogs,
                "rejections": [
                    {"rule_id": "RE_INSERT_CH2", "reason": "NO_APPLICABLE_SITE", "candidate": None},
                    {"rule_id": "RC_REMOVE_CH2", "reason": "PROTECTED_TOUCHED", "touched_atom_maps": [2], "protected_atom_maps": [2]},
                ],
            }

        def passed(_parent, record):
            if record.get("candidate_id") == "fixture-fail":
                return {
                    "valid": False,
                    "hard_reasons": ["EXACT_FIXTURE_REASON"],
                    "soft_flags": ["fixture diagnostic"],
                }
            return {"valid": True, "hard_reasons": [], "soft_flags": []}

        audit = audit_rule_applicability(
            parent, sites, generate_fn=generated, cheap_filter_fn=passed
        )
        self.assertEqual([False, True], [row[0] for row in calls])
        self.assertEqual(0, audit["strict"]["actual_broad_family_count"])
        self.assertEqual(1, audit["exploratory"]["actual_broad_family_count"])
        self.assertGreater(audit["catalog_rule_count"], 6)
        self.assertTrue(audit["catalog_availability_does_not_count_as_actual"])
        self.assertEqual("ring_modification", audit["ring_expansion_and_contraction_merged_as"])
        ring_rows = {
            row["rule_id"]: row for row in audit["strict"]["rules"]
            if row["rule_id"] in {"RE_INSERT_CH2", "RC_REMOVE_CH2"}
        }
        self.assertEqual("ring_modification", ring_rows["RE_INSERT_CH2"]["broad_family"])
        self.assertIn("chemical_pattern_missing", ring_rows["RE_INSERT_CH2"]["rejection_reason_counts"])
        self.assertIn("protected", ring_rows["RC_REMOVE_CH2"]["rejection_reason_counts"])
        linker_row = next(
            row for row in audit["exploratory"]["rules"]
            if row["rule_id"] == "LH_HYDROXYETHYL"
        )
        self.assertEqual([1, 2, 3, 4], linker_row["exact_touched_maps"])
        self.assertEqual([{
            "candidate_id": "fixture-fail",
            "reasons": ["EXACT_FIXTURE_REASON"],
            "modified_atom_maps": [2],
            "diagnostics": {
                "valid": False,
                "hard_reasons": ["EXACT_FIXTURE_REASON"],
                "soft_flags": ["fixture diagnostic"],
            },
        }], linker_row["cheap_filter_failures"])


class FailurePreservationTests(unittest.TestCase):
    def test_safe_mocked_worker_failure_is_retained_without_real_claim(self):
        parent = _parent()
        record = enumerate_expansion_candidates(parent)[0]

        def failed_docker(*args, **kwargs):
            return {
                "status": "failed_docking",
                "real_docking": False,
                "results": {"error": "fixture failure"},
            }

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "case"
            outcome = _execute_case(
                parent=parent,
                protein=[{"fixture": True}],
                protected_maps=[1, 2, 4, 5, 6],
                receptor_pdbqt=Path(directory) / "unused.pdbqt",
                record=record,
                case_dir=root,
                audit_only=False,
                docker=failed_docker,
            )
            self.assertEqual("docking_failed_not_replaced", outcome["status"])
            self.assertEqual("raw-dock-receipt.json", outcome["raw_dock_receipt"])
            self.assertFalse(outcome["docking"]["real_docking"])
            self.assertTrue((root / "raw-candidate.sdf").is_file())
            self.assertTrue((root / "raw-dock-receipt.json").is_file())
            self.assertTrue((root / "outcome.json").is_file())
            self.assertFalse(outcome["claims"]["binding"])

    def test_raised_worker_counts_as_one_failed_attempt(self):
        parent = _parent()
        record = enumerate_expansion_candidates(parent)[0]

        def raised_docker(*args, **kwargs):
            raise RuntimeError("fixture worker exception")

        with tempfile.TemporaryDirectory() as directory:
            outcome = _execute_case(
                parent, [], [], Path(directory) / "unused", record,
                Path(directory) / "case", audit_only=False, docker=raised_docker,
            )
            counts = _actual_counts([outcome], [record])
            self.assertEqual(1, counts["dock_attempt"])
            self.assertEqual(0, counts["dock_completed"])
            self.assertEqual(1, counts["dock_failed"])

    def test_completed_flag_with_false_real_docking_counts_zero_completed(self):
        parent = _parent()
        record = enumerate_expansion_candidates(parent)[0]

        def mocked_completion(*args, **kwargs):
            return {"status": "completed_with_limits", "real_docking": False}

        with tempfile.TemporaryDirectory() as directory:
            outcome = _execute_case(
                parent, [], [], Path(directory) / "unused", record,
                Path(directory) / "case", audit_only=False, docker=mocked_completion,
            )
            counts = _actual_counts([outcome], [record])
            self.assertEqual("docking_failed_not_replaced", outcome["status"])
            self.assertEqual(1, counts["dock_attempt"])
            self.assertEqual(0, counts["dock_completed"])
            self.assertEqual(1, counts["dock_failed"])

    def test_audit_only_never_calls_docker(self):
        parent = _parent()
        record = enumerate_expansion_candidates(parent)[0]

        def forbidden(*args, **kwargs):
            raise AssertionError("Docker must not run in audit-only mode")

        with tempfile.TemporaryDirectory() as directory:
            outcome = _execute_case(
                parent, [], [], Path(directory) / "unused", record,
                Path(directory) / "case", audit_only=True, docker=forbidden,
            )
            self.assertEqual("audit_only_no_docking", outcome["status"])
            self.assertIsNone(outcome["docking"])


if __name__ == "__main__":
    unittest.main()

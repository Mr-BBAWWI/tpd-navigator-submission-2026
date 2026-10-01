"""Tests for the isolated exploratory native N3 attachment experiment."""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

from packages.science.native_attachment_experiment import (
    FRAGMENTS,
    _pose_metrics,
    _preservation_pass_count,
    enumerate_n3_candidates,
    find_source_linkage,
    verify_native_probe,
)


def _parent() -> Chem.Mol:
    mol = Chem.MolFromSmiles("[CH3:1][NH:3][CH2:2][C@@H:4]([F:5])[Cl:6]")
    if mol is None:
        raise AssertionError("test parent failed to parse")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class NativeAttachmentGraphTests(unittest.TestCase):
    def test_twelve_distinct_graphs_alter_n3_only(self):
        parent = _parent()
        records = enumerate_n3_candidates(parent)
        self.assertEqual(12, len(records))
        self.assertEqual(12, len({record["candidate_id"] for record in records}))
        parent_maps = {
            atom.GetAtomMapNum(): atom.GetAtomicNum() for atom in parent.GetAtoms()
        }
        parent_bonds = {
            tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))):
            str(bond.GetBondType())
            for bond in parent.GetBonds()
        }
        for record in records:
            candidate = Chem.MolFromSmiles(record["mapped_smiles"])
            self.assertIsNotNone(candidate)
            candidate_atoms = {
                atom.GetAtomMapNum(): atom for atom in candidate.GetAtoms()
            }
            for atom_map, atomic_number in parent_maps.items():
                self.assertIn(atom_map, candidate_atoms)
                self.assertEqual(atomic_number, candidate_atoms[atom_map].GetAtomicNum())
            for pair, bond_type in parent_bonds.items():
                bond = candidate.GetBondBetweenAtoms(
                    candidate_atoms[pair[0]].GetIdx(), candidate_atoms[pair[1]].GetIdx()
                )
                self.assertIsNotNone(bond)
                self.assertEqual(bond_type, str(bond.GetBondType()))
            external = []
            original = set(parent_maps)
            for bond in candidate.GetBonds():
                left = bond.GetBeginAtom().GetAtomMapNum()
                right = bond.GetEndAtom().GetAtomMapNum()
                if (left in original) != (right in original):
                    external.append({left, right})
            self.assertEqual(1, len(external))
            self.assertIn(3, external[0])
            self.assertEqual(0, candidate_atoms[3].GetFormalCharge())
            self.assertEqual(0, candidate_atoms[3].GetTotalNumHs())
            self.assertEqual("UNKNOWN", record["site_status"][0]["state"])
            self.assertTrue(record["site_status"][0]["scientific_pending"])
            self.assertFalse(record["strict_qualified_for_counts"])

    def test_terminal_handlers_and_added_maps(self):
        records = enumerate_n3_candidates(_parent())
        self.assertEqual(12, len(FRAGMENTS))
        for record in records:
            self.assertTrue(record["added_atom_maps"])
            self.assertEqual(
                len(record["added_atom_maps"]), len(set(record["added_atom_maps"]))
            )
            self.assertTrue(all(atom_map >= 5000 for atom_map in record["added_atom_maps"]))
            mol = Chem.MolFromSmiles(record["mapped_smiles"])
            endpoint = next(
                atom for atom in mol.GetAtoms()
                if atom.GetAtomMapNum() == max(record["added_atom_maps"])
            )
            self.assertIn(endpoint.GetAtomicNum(), (7, 8))
            self.assertGreaterEqual(endpoint.GetTotalNumHs(), 1)
            self.assertTrue(record["protected_graph_preserved"])

    def test_no_approval_or_bioisostere_flags(self):
        encoded = json.dumps(enumerate_n3_candidates(_parent()), sort_keys=True).lower()
        self.assertNotIn("bioisostere", encoded)
        self.assertNotIn('"approved": true', encoded)
        self.assertNotIn('"approval_claim": true', encoded)
        self.assertNotIn('"scientific_pending": false', encoded)


class SourceTests(unittest.TestCase):
    def _fixtures(self):
        locators = [
            {"compound_id_column": "Ligands", "doi": "10.1021/acs.jmedchem.4c01903", "filename": "supporting/extra/01_jm4c01903_si_006.csv", "record_number_1_based_including_header": 5},
            {"compound_id_column": "Ligands", "doi": "10.1021/acs.jmedchem.4c01903", "filename": "supporting/extra/01_jm4c01903_si_006.csv"},
        ]
        selected_row = {
            "pair_id": "SMI-6080__SMD-6087",
            "requested_change_label": "terminal_N_derivatization_to_PROTAC",
            "selected_parent_pair_side": "left",
            "other_compound_id": "SMD-6087",
            "source_row_locators": locators,
            "source_affected_atom_maps": [19],
            "selected_parent_affected_atom_maps": [3],
            "source_site_consensus": True,
        }
        selected = {
            "status": "exact_measured_parent_join",
            "matched_record": {"compound_id": "SMI-6080"},
            "source_sar": [selected_row],
        }
        raw_row = {
            "pair_id": "SMI-6080__SMD-6087",
            "left_compound_id": "SMI-6080",
            "right_compound_id": "SMD-6087",
            "requested_change_label": "terminal_N_derivatization_to_PROTAC",
            "source_row_locators": locators,
            "mapping_possibilities": [
                {"left_affected_parent_atom_maps": [19]},
                {"left_affected_parent_atom_maps": [19]},
            ],
        }
        return selected, {"sar_pair_evidence": [raw_row]}

    def test_realistic_selected_and_raw_map_spaces_join(self):
        selected, source = self._fixtures()
        result = find_source_linkage(selected, source)
        self.assertEqual([3], result["selected_parent_atom_maps"])
        self.assertEqual([19], result["source_parent_atom_maps"])
        self.assertEqual(
            "result.selected_parent_measured_evidence.source_sar[0]",
            result["source_full_pointer"],
        )

    def test_source_is_read_only(self):
        selected, source = self._fixtures()
        before = json.dumps([selected, source], sort_keys=True)
        result = find_source_linkage(selected, source)
        result["source_record"]["pair_id"] = "changed-copy"
        self.assertEqual(before, json.dumps([selected, source], sort_keys=True))

    def test_missing_duplicate_and_unrelated_pairs_fail(self):
        selected, source = self._fixtures()
        with self.assertRaises(ValueError):
            find_source_linkage(selected, {"sar_pair_evidence": []})
        duplicate = json.loads(json.dumps(source))
        duplicate["sar_pair_evidence"].append(duplicate["sar_pair_evidence"][0])
        with self.assertRaises(ValueError):
            find_source_linkage(selected, duplicate)
        unrelated = json.loads(json.dumps(source))
        unrelated["sar_pair_evidence"][0]["pair_id"] = "SMI-6080__OTHER"
        with self.assertRaises(ValueError):
            find_source_linkage(selected, unrelated)

    def test_unrelated_or_missing_locators_do_not_fallback(self):
        selected, source = self._fixtures()
        selected["source_sar"][0]["source_row_locators"] = []
        source["unrelated_locator"] = {
            "doi": "10.1021/acs.jmedchem.4c01903",
            "filename": "unrelated.csv",
            "compound_id_column": "Ligands",
        }
        with self.assertRaises(ValueError):
            find_source_linkage(selected, source)

    def test_raw_source_map_mismatch_fails(self):
        selected, source = self._fixtures()
        source["sar_pair_evidence"][0]["mapping_possibilities"][1]["left_affected_parent_atom_maps"] = [3]
        with self.assertRaises(ValueError):
            find_source_linkage(selected, source)


class PoseDiagnosticTests(unittest.TestCase):
    def test_exact_booleans_and_status_only_count(self):
        run = {"results": {"pose_preservation": {"passing_pose_count": 99, "all_poses_diagnostics": [
            {"docking_pose_preserved": True, "status": "pass"},
            {"docking_pose_preserved": True, "status": "fail"},
            {"docking_pose_preserved": 1, "status": "pass"},
            {"docking_pose_preserved": False, "status": "pass", "pass": True},
        ]}}}
        self.assertEqual(1, _preservation_pass_count(run))

    def test_best_geometry_retention_comes_from_same_pose(self):
        run = {"results": {
            "poses": [
                {"score": {"affinity": -7.0}},
                {"score": {"affinity": -10.0}},
            ],
            "pose_preservation": {"all_poses_diagnostics": [
                {
                    "pose_index_zero_based": 0,
                    "caller_supplied_score": -7.0,
                    "core_RMSD_A_in_receptor_frame": 0.5,
                    "protected_atom_near_residue_preservation": {"retention_fraction": 0.6},
                    "docking_pose_preserved": True,
                    "status": "pass",
                },
                {
                    "pose_index_zero_based": 1,
                    "caller_supplied_score": -10.0,
                    "core_RMSD_A_in_receptor_frame": 4.0,
                    "protected_atom_near_residue_preservation": {"retention_fraction": 0.95},
                    "docking_pose_preserved": False,
                    "status": "fail",
                },
            ]},
        }}
        metrics = _pose_metrics(run)
        self.assertEqual(0, metrics["best_geometry_pose_index_zero_based"])
        self.assertEqual(0.5, metrics["best_geometry_core_RMSD_A"])
        self.assertEqual(0.6, metrics["best_geometry_contact_retention"])
        self.assertEqual(-10.0, metrics["global_best_raw_score"])


class NativeManifestTests(unittest.TestCase):
    def _probe(self, root: Path) -> None:
        protocol = root / "protocol.json"
        result = root / "native-probe-result.json"
        protocol.write_text("{}\n", encoding="utf-8")
        result.write_text("{}\n", encoding="utf-8")
        rows = []
        for path in (result, protocol):
            rows.append({
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            })
        (root / "manifest.json").write_text(
            json.dumps({"files": rows}, sort_keys=True) + "\n", encoding="utf-8"
        )

    def test_verified_manifest_accepts_all_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._probe(root)
            receipt = verify_native_probe(root)
            self.assertIn("manifest_sha256", receipt)

    def test_tampered_manifest_artifact_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._probe(root)
            (root / "native-probe-result.json").write_text("{\"tampered\":true}\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_native_probe(root)

    def test_unlisted_file_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._probe(root)
            (root / "extra.txt").write_text("not declared", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_native_probe(root)

    def test_duplicate_manifest_path_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._probe(root)
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
            manifest["files"].append(dict(manifest["files"][0]))
            (root / "manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_native_probe(root)

    def test_protocol_and_result_must_both_be_declared(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._probe(root)
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
            manifest["files"] = [
                row for row in manifest["files"] if row["path"] != "protocol.json"
            ]
            (root / "manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_native_probe(root)


if __name__ == "__main__":
    unittest.main()

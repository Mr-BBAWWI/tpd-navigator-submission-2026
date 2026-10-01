from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

from packages.science.parent_sar_probe import (
    EvidenceError,
    PARENT_ID,
    build_exact_halogen_probe,
    summarise_exports,
    validate_predeclared_seeds,
)


def _bytes(value) -> bytes:
    return (json.dumps(value, sort_keys=True) + "\n").encode()


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _make_export(root: Path, index: int, result: dict | None = None) -> Path:
    parent = PARENT_ID if index == 0 else f"PARENT-{index}"
    directory = root / f"export-{index}"
    blobs = directory / "blobs"
    blobs.mkdir(parents=True)
    if result is None:
        result = {
            "parameters": {"parent_id": parent},
            "analogs": [
                {
                    "candidate_id": f"C-{index}-{number}",
                    "transformation_class": "ring_expansion" if number == 0 else "heteroatom_swap",
                    "mapped_smiles": "[CH3:1][OH:2]",
                    "cheap_filter": {"valid": True},
                    "qualified_for_counts": number == 0,
                }
                for number in range(6)
            ],
            "rule_catalog": ["ring_contraction", "bioisosteric_replacement"],
            "summary": {"cheap_passed_count": 999},
            "rejections": [{"reason": "ATTACHMENT_SITE_EXCLUDED"}],
        }
    artifact_id = f"a-{index}"
    artifact_data = _bytes({"provenance": index})
    artifact_sha = _sha(artifact_data)
    (blobs / artifact_id).write_bytes(artifact_data)
    ref = {"artifact_id": artifact_id, "version": 1, "sha256": artifact_sha}
    receipt = {
        "id": f"job-{index}",
        "state": "completed",
        "stage": "finished",
        "parameters": {"parent_id": parent},
        "runtime": {"artifacts": {"provenance": ref}},
        "result": {
            **result,
            "archived_refs": {"provenance": ref},
            "source_metadata": {"sha256": "f" * 64},
        },
        "outputs": [],
    }
    receipt_data = _bytes(receipt)
    receipt_sha = _sha(receipt_data)
    (directory / "job-receipt.json").write_bytes(receipt_data)
    manifest = {
        "job_id": f"job-{index}",
        "parent_id": parent,
        "raw_receipt_sha256": receipt_sha,
        "artifacts": [{"artifact_id": artifact_id, "version": 1, "sha256": artifact_sha}],
    }
    (directory / "manifest.json").write_bytes(_bytes(manifest))
    return directory


def _mapped_parent() -> Chem.Mol:
    # Deliberately tiny test-only fixture: map8 terminal Cl on aromatic map9.
    mol = Chem.MolFromSmiles("[Cl:8][c:9]1[cH:1][cH:2][cH:3][cH:4][cH:5]1")
    assert mol is not None
    return mol


def _source(parent_smiles="Clc1ccccc1", product_smiles="Brc1ccccc1"):
    data = {
        "freeze_provenance": {"original_ref": "original.json", "original_content_sha256": "b" * 64},
        "parent_catalog": {"records": [
            {"compound_id": "SMI-6080", "identity": {"source_smiles": parent_smiles}},
            {"compound_id": "SMI-6085", "identity": {"source_smiles": product_smiles}},
        ]},
        "sar_pair_evidence": [{
            "pair_id": "SMI-6080__SMI-6085",
            "requested_change_label": "chlorine_to_bromine",
            "source_row_locators": [{"filename": "source.csv", "record_number_1_based_including_header": 5}],
        }],
    }
    return {"status": "configured_hash_verified", "relative_path": "medchem.json", "sha256": "a" * 64, "data": data}


def _linkage(protected=True):
    return {
        "status": "exact_measured_parent_join",
        "selected_parent_id": PARENT_ID,
        "matched_record": {"compound_id": "SMI-6080"},
        "source_sar": [{
            "pair_id": "SMI-6080__SMI-6085",
            "requested_change_label": "chlorine_to_bromine",
            "other_compound_id": "SMI-6085",
            "selected_parent_affected_atom_maps": [9],
            "source_row_locators": [],
        }],
        "protected_atom_maps": [8] if protected else [],
    }


class ExportAuditTests(unittest.TestCase):
    def test_all_nine_and_aggregate_not_strict(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for index in range(9):
                _make_export(root, index)
            audit = summarise_exports(root)
            self.assertEqual(audit["parent_count"], 9)
            self.assertFalse(audit["pooled_union_satisfies_per_parent_gate"])
            self.assertIn("ring_modification", audit["pooled_cheap_passed_family_union"])
            self.assertEqual(audit["parents"][0]["selected_candidate_ids"], [])
            self.assertFalse(audit["parents"][0]["per_parent_at_least_6_gate"])
            self.assertTrue(audit["parents"][0]["stored_summary_count_mismatches"])

    def test_corrupted_receipt_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directories = [_make_export(root, index) for index in range(9)]
            (directories[0] / "job-receipt.json").write_text("tampered")
            with self.assertRaises(EvidenceError):
                summarise_exports(root)

    def test_corrupted_blob_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directories = [_make_export(root, index) for index in range(9)]
            blob = next((directories[0] / "blobs").iterdir())
            blob.write_bytes(blob.read_bytes() + b"x")
            with self.assertRaises(EvidenceError):
                summarise_exports(root)

    def test_missing_nested_archived_ref_declaration_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = _make_export(root, 0)
            receipt_path = directory / "job-receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["result"]["archived_refs"]["missing"] = {"artifact_id": "a-dead", "version": 1, "sha256": "0" * 64}
            receipt_path.write_bytes(_bytes(receipt))
            manifest_path = directory / "manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["raw_receipt_sha256"] = _sha(receipt_path.read_bytes())
            manifest_path.write_bytes(_bytes(manifest))
            with self.assertRaises(EvidenceError):
                summarise_exports(root, expected_count=1)

    def test_three_families_across_one_hundred_analogs_is_not_six_family_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {"parameters": {"parent_id": PARENT_ID}, "analogs": [
                {"candidate_id": f"C-{n}", "transformation_class": f"family-{n % 3}",
                 "cheap_filter": {"valid": True}, "qualified_for_counts": True}
                for n in range(100)
            ]}
            _make_export(root, 0, result)
            self.assertFalse(summarise_exports(root, expected_count=1)["parents"][0]["per_parent_at_least_6_gate"])

    def test_truthy_string_booleans_do_not_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {"parameters": {"parent_id": PARENT_ID}, "analogs": [{
                "candidate_id": "truthy", "transformation_class": "swap",
                "cheap_filter": {"valid": "true"}, "qualified_for_counts": "true",
                "selected": "true", "pipeline_status": "qualified",
            }]}
            _make_export(root, 0, result)
            row = summarise_exports(root, expected_count=1)["parents"][0]
            self.assertEqual(row["cheap_passed_count"], 0)
            self.assertEqual(row["selected_candidate_ids"], [])

    def test_actual_embedded_result_pose_shape_is_counted_without_invented_redocking_block(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {"parameters": {"parent_id": PARENT_ID}, "analogs": [{
                "candidate_id": "actual-shape",
                "pipeline_status": "qualified",
                "selected": True,
                "qualified_for_counts": False,
                "docking": {
                    "status": "completed_with_limits",
                    "pose_preserved": True,
                    "error": None,
                },
                "parent_redocking_supported": True,
            }]}
            _make_export(root, 0, result)
            row = summarise_exports(root, expected_count=1)["parents"][0]
            self.assertEqual(row["pose_qualified_candidate_ids"], ["actual-shape"])
            self.assertEqual(row["selected_candidate_ids"], ["actual-shape"])
            self.assertEqual(row["strict_qualified_candidate_ids"], [])

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {"parameters": {"parent_id": PARENT_ID}, "analogs": [{
                "candidate_id": "truthy", "transformation_class": "swap",
                "cheap_filter": {"valid": "true"}, "qualified_for_counts": "true",
                "selected": "true", "pipeline_status": "qualified",
            }]}
            _make_export(root, 0, result)
            row = summarise_exports(root, expected_count=1)["parents"][0]
            self.assertEqual(row["cheap_passed_count"], 0)
            self.assertEqual(row["selected_candidate_ids"], [])

    def test_fake_positive_summary_boolean_does_not_qualify(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = {
                "parameters": {"parent_id": PARENT_ID},
                "summary": {"qualified_for_counts": True, "meets_all_acceptance_requirements": True},
                "analogs": [{
                    "candidate_id": "fake", "transformation_class": "heteroatom_swap",
                    "mapped_smiles": "[CH3:1][OH:2]", "cheap_filter": {"valid": False},
                    "qualified_for_counts": True,
                }],
            }
            _make_export(root, 0, result)
            audit = summarise_exports(root, expected_count=1)
            row = audit["parents"][0]
            self.assertEqual(row["strict_qualified_candidate_ids"], [])
            self.assertFalse(row["per_parent_at_least_6_gate"])


class ExactPairTests(unittest.TestCase):
    def test_exact_pair_and_protected_conflict_visible_no_policy_activation(self):
        proposal = build_exact_halogen_probe(_mapped_parent(), _linkage(True), _source())
        self.assertTrue(proposal["graph_verified_against_source_SMI_6085"])
        self.assertTrue(proposal["protected_map_conflict"])
        self.assertEqual(proposal["protected_conflict_atom_maps"], [8])
        self.assertNotIn(8, proposal["candidate_diagnostic_mask"])
        self.assertFalse(proposal["qualified_for_counts"])
        self.assertFalse(proposal["changes_design_policy"])
        self.assertFalse(proposal["scientific_approval"])

    def test_missing_source_record_refused(self):
        source = _source()
        source["data"]["parent_catalog"]["records"].pop()
        with self.assertRaises(EvidenceError):
            build_exact_halogen_probe(_mapped_parent(), _linkage(), source)

    def test_wrong_element_or_mapping_refused(self):
        wrong = Chem.MolFromSmiles("[Br:8][c:9]1[cH:1][cH:2][cH:3][cH:4][cH:5]1")
        with self.assertRaises(EvidenceError):
            build_exact_halogen_probe(wrong, _linkage(), _source())
        wrong_map = Chem.MolFromSmiles("[Cl:7][c:9]1[cH:1][cH:2][cH:3][cH:4][cH:5]1")
        with self.assertRaises(EvidenceError):
            build_exact_halogen_probe(wrong_map, _linkage(), _source())

    def test_unexpected_product_graph_refused(self):
        with self.assertRaises(EvidenceError):
            build_exact_halogen_probe(_mapped_parent(), _linkage(), _source(product_smiles="Brc1ccncc1"))

    def test_repeated_seeds_must_all_be_retained(self):
        rows = [
            {"ligand": ligand, "seed": seed, "result": {"status": "failed"}}
            for ligand in ("parent", "observed_br") for seed in (23, 41, 61)
        ]
        validate_predeclared_seeds([23, 41, 61], rows)
        with self.assertRaises(EvidenceError):
            validate_predeclared_seeds([23, 41, 61], rows[:-1])
        with self.assertRaises(ValueError):
            validate_predeclared_seeds([23, 23], rows)
        with self.assertRaises(ValueError):
            validate_predeclared_seeds([True], rows)
        with self.assertRaises(ValueError):
            validate_predeclared_seeds([-1], rows)
        with self.assertRaises(EvidenceError):
            validate_predeclared_seeds([23, 41, 61], rows + [rows[0]])


if __name__ == "__main__":
    unittest.main()

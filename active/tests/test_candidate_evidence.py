"""Stored real source consumption and explicitly synthetic literature exchange tests."""
import copy
from pathlib import Path
import unittest

from packages.science.boltz_collection import empty_collection
from packages.science.candidate_evidence import build_evidence, markdown, validate_evidence
from packages.science.evidence_common import seal
from packages.science.handoff import parse

ACTIVE = Path(__file__).resolve().parents[1]


class CandidateEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.material = parse((ACTIVE / "outputs/b_handoff_20260923/material.json").read_bytes())
        self.collection = empty_collection(self.material)

    def exchange(self):
        source = self.material["h1"]["observations"][0]
        return {"format": "tpd-b-literature-exchange/0.1.0-draft", "producer": "B synthetic adapter test; not teammate output",
                "material_digest": self.material["material_digest"], "data_mode": "synthetic_test",
                "records": [{"record_id": "sample-1", "query": {"doi": source["citation"]["doi"], "paper_name": "PROTAC 1", "compound_number": "2"},
                             "observation": copy.deepcopy(source["original"]),
                             "source": {"doi": source["citation"]["doi"], "sha256": source["citation"]["source_sha256"], "locator": source["citation"]["locator"], "review_status": "pending"}}]}

    def test_actual_known_candidates_keep_missing_observation_and_no_predictions(self):
        before = copy.deepcopy(self.material)
        result = build_evidence(self.material)
        self.assertEqual([c["compound"]["compound_id"] for c in result["candidates"]], ["C01", "C02"])
        self.assertIsNone(result["candidates"][1]["literature_observations"][0]["original"]["value"])
        self.assertIsNone(result["candidates"][0]["literature_observations"][0]["original"]["exposure_time_h"])
        self.assertTrue(all(c["prediction"]["status"] == "not_run" and not c["prediction"]["samples"] for c in result["candidates"]))
        self.assertEqual(result["literature_reconciliation"]["status"], "not_provided")
        self.assertEqual(self.material, before)

    def test_cpu_properties_and_start_geometry_not_relabelled_as_prediction(self):
        result = build_evidence(self.material)
        self.assertEqual(result["starting_ligand"]["compound_id"], "START")
        self.assertEqual(result["shared_hypothesis"]["descriptive_geometry"]["attachment_ccd_atom"], "N18")
        for c in result["candidates"]:
            self.assertEqual(c["computed_properties"]["kind"], "computed")
            self.assertEqual(c["computed_properties"]["values"], c["compound"]["identity"]["properties"])

    def test_projection_cannot_silently_change_value_or_status(self):
        for mutate in (lambda b: b["candidates"][1]["literature_observations"][0]["original"].update(value=300),
                       lambda b: b["candidates"][0]["prediction"].update(status="available_unreviewed")):
            result = build_evidence(self.material)
            mutate(result)
            result = seal({k: v for k, v in result.items() if k != "digest"})
            with self.assertRaisesRegex(ValueError, "PROJECTION"):
                validate_evidence(result, self.material, self.collection)

    def test_approval_and_public_ready_not_promoted(self):
        result = build_evidence(self.material)
        for key in ("dispatch_authorized", "approval_record_created", "public_release_ready"):
            copy_result = copy.deepcopy(result)
            copy_result["authority"][key] = True
            with self.assertRaisesRegex(ValueError, "SCHEMA"):
                validate_evidence(copy_result, self.material, self.collection)

    def test_synthetic_exchange_requires_explicit_test_option(self):
        with self.assertRaisesRegex(ValueError, "SYNTHETIC_LITERATURE"):
            build_evidence(self.material, exchange=self.exchange())
        result = build_evidence(self.material, exchange=self.exchange(), allow_synthetic=True)
        self.assertEqual(result["data_mode"], "synthetic_test")
        self.assertEqual(result["literature_reconciliation"]["records"][0]["comparison"]["status"], "same_record_pending_review")

    def test_incoming_conditions_and_values_never_replace_reference(self):
        exchange = self.exchange()
        exchange["records"][0]["observation"]["value"] = 100
        result = build_evidence(self.material, exchange=exchange, allow_synthetic=True)
        row = result["literature_reconciliation"]["records"][0]
        self.assertEqual(row["comparison"]["incoming"]["value"], 100)
        self.assertEqual(row["comparison"]["status"], "difference_requires_review")
        self.assertEqual(result["candidates"][0]["literature_observations"][0]["original"]["value"], 300)

    def test_unresolved_identity_and_source_version_are_visible(self):
        exchange = self.exchange()
        exchange["records"][0]["query"]["paper_name"] = "ACBI1"
        exchange["records"][0]["source"]["sha256"] = "0" * 64
        result = build_evidence(self.material, exchange=exchange, allow_synthetic=True)
        row = result["literature_reconciliation"]["records"][0]
        self.assertEqual(row["comparison"]["status"], "identity_unresolved")
        self.assertIn("INCOMING_SOURCE_VERSION_NOT_VERIFIED", row["source_flags"])

    def test_stale_and_duplicate_literature_records_rejected(self):
        exchange = self.exchange()
        exchange["material_digest"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "STALE_MATERIAL"):
            build_evidence(self.material, exchange=exchange, allow_synthetic=True)
        exchange = self.exchange()
        exchange["records"].append(copy.deepcopy(exchange["records"][0]))
        with self.assertRaisesRegex(ValueError, "DUPLICATE_LITERATURE"):
            build_evidence(self.material, exchange=exchange, allow_synthetic=True)

    def test_report_has_real_missing_states_and_no_invented_gpu_metrics(self):
        report = markdown(build_evidence(self.material)).decode()
        for phrase in ("미실행", "미추출", "START 실험 구조", "동료 A의 실제 출력은 아직", "= 300 nM"):
            self.assertIn(phrase, report)
        self.assertNotIn("confidence |", report)

    def test_changed_collection_version_rejected(self):
        collection = copy.deepcopy(self.collection)
        collection["material_digest"] = "0" * 64
        collection = seal({k: v for k, v in collection.items() if k != "digest"})
        with self.assertRaisesRegex(ValueError, "STALE_MATERIAL"):
            build_evidence(self.material, collection)

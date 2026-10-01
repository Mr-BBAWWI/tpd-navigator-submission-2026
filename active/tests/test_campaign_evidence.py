import hashlib
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

from packages.agents.campaign_evidence import (
    EvidenceValidationError,
    build_campaign_evidence,
)


CRITERION_IDS = [
    "parent_funnel",
    "expert_parent_selection",
    "modifiable_sites",
    "distinct_constitutional_graphs",
    "actual_broad_families",
    "qualified_panel",
    "both_e3_assembly",
    "core_interaction_preservation",
    "microstates_h_direction",
    "known_crbn_calibration",
    "novel_ternary_repeats",
    "novel_ternary_geometry",
    "exact_synthesis_review",
    "formal_expert_decision",
]


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def artifact_id(number):
    return "a-" + format(number, "032x")


class CampaignEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.assessment = self.root / "assessment.json"
        self.write_assessment()

    def tearDown(self):
        self.temp.cleanup()

    def write_assessment(self, criteria=None, **extra):
        if criteria is None:
            criteria = [
                {"id": criterion_id, "title": criterion_id, "status": "failed"}
                for criterion_id in CRITERION_IDS
            ]
        value = {
            "id": "assessment-1",
            "job_id": "assessment-job",
            "revision": 1,
            "scientific_accepted": True,
            "criteria": criteria,
        }
        value.update(extra)
        self.assessment.write_bytes(encoded(value))

    def analog(
        self,
        analog_id="analog-1",
        family="heteroatom_swap",
        status="completed_with_limits",
        selected=True,
        qualified=False,
        exploratory=True,
    ):
        return {
            "id": analog_id,
            "canonical_smiles": "C(C)O",
            "transformation_class": family,
            "attachment_site_atom_maps": [1],
            "cheap_filter": {"valid": True},
            "qualified_for_counts": qualified,
            "selected": selected,
            "assembly_eligible": True,
            "parent_redocking_supported": True,
            "docking": {
                "status": status,
                "pose_preserved": status.startswith("completed"),
                "passing_pose_count": 1 if status.startswith("completed") else 0,
            },
        }

    def candidate(self, candidate_id="candidate-1", analog_id="analog-1", **updates):
        value = {
            "candidate_id": candidate_id,
            "canonical_smiles": "NCCO",
            "warhead_analog_id": analog_id,
            "e3_type": "CRBN",
            "assembly_mode": "pose_supported_hypothesis",
            "status": "hypothesis_pending_review",
            "final_candidate": False,
            "synthetic_route_status": "not_assessed",
        }
        value.update(updates)
        return value

    def full_ref(self, aid, value):
        raw = encoded(value)
        return {
            "artifact_id": aid,
            "version": 1,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "media_type": "application/json",
            "schema_id": "research-artifact/1",
            "provenance": {"producer": "test-fixture"},
        }

    def make_export(
        self,
        name,
        job_id,
        parent_id,
        analogs=None,
        candidates=None,
        exploratory=True,
        artifacts=None,
        result_extra=None,
    ):
        root = self.root / name
        blobs = root / "blobs"
        blobs.mkdir(parents=True)
        result = {
            "analogs": analogs if analogs is not None else [self.analog()],
            "candidates": candidates if candidates is not None else [self.candidate()],
            "sites": [{"state": "MODIFIABLE", "atom_map": 1}],
        }
        if result_extra:
            result.update(result_extra)
        receipt = {
            "job_id": job_id,
            "parent_id": parent_id,
            "parameters": {"exploratory": exploratory},
            "result": result,
        }
        records = []
        for number, value in enumerate(artifacts or [], 1):
            aid = artifact_id(number)
            raw = encoded(value)
            (blobs / aid).write_bytes(raw)
            records.append(
                {
                    "artifact_id": aid,
                    "version": 1,
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "media_type": "application/json",
                    "schema_id": "research-artifact/1",
                    "provenance": {"producer": "test-fixture"},
                }
            )
        receipt_raw = encoded(receipt)
        (root / "job-receipt.json").write_bytes(receipt_raw)
        manifest = {
            "job_id": job_id,
            "parent_id": parent_id,
            "raw_receipt_sha256": hashlib.sha256(receipt_raw).hexdigest(),
            "artifacts": records,
        }
        (root / "manifest.json").write_bytes(encoded(manifest))
        return root

    def build(self, exports):
        return build_campaign_evidence(self.assessment, exports)

    def test_roundtrip_two_parents_scopes_same_candidate_and_digest_is_deterministic(self):
        first = self.make_export(
            "first", "job-one", "SMARCA2-9QAD-A1I45", candidates=[self.candidate("same")]
        )
        second = self.make_export(
            "second", "job-two", "SMARCA2-9D12-A1A1P", candidates=[self.candidate("same")]
        )
        one = self.build([second, first])
        two = self.build([first, second])
        self.assertEqual(one["digest"], two["digest"])
        self.assertEqual(
            ["SMARCA2-9D12-A1A1P/job-two/same", "SMARCA2-9QAD-A1I45/job-one/same"],
            one["retrieval_indexes"]["candidate_keys"],
        )
        self.assertEqual("CCO", one["analogs"][0]["canonical_smiles"])
        self.assertFalse(one["approval"])
        self.assertFalse(one["scientific_accepted"])

    def test_family_status_and_assembly_counts_are_recomputed(self):
        analogs = [
            self.analog("a1", "linker_handle_introduction", "completed_with_limits"),
            self.analog("a2", "heteroatom_swap", "not_eligible"),
            self.analog("a3", "ring_expansion", "failed_docking"),
            self.analog("a4", "conformational_restriction", "mystery"),
        ]
        candidates = [
            self.candidate("c1", "a1"),
            self.candidate("c2", "a1", assembly_mode="preview"),
            self.candidate("c3", "absent"),
            self.candidate("c4", "a1", status="failed"),
            self.candidate("c5", "a1", canonical_smiles="invalid"),
        ]
        export = self.make_export(
            "counts", "job-counts", "parent-counts", analogs, candidates
        )
        bundle = self.build([export])
        counts = bundle["campaign_counts_recomputed_from_records"]
        self.assertEqual(4, counts["actual_generated_family_count"])
        self.assertEqual(
            [
                "conformational_restriction",
                "heteroatom_swap",
                "linker_handle_introduction",
                "ring_modification",
            ],
            counts["actual_generated_families"],
        )
        self.assertEqual(0, counts["supported_site_family_count"])
        self.assertEqual(2, counts["docking_attempts"])
        self.assertEqual(1, counts["docking_completed"])
        self.assertEqual(1, counts["failed_docking"])
        self.assertEqual(1, counts["assemblies"])
        parent_counts = bundle["parents"][0]["counts_recomputed_from_records"]
        self.assertEqual(1, parent_counts["docking_skipped"])
        self.assertEqual({"mystery": 1}, parent_counts["docking_unknown_statuses"])

    def test_unsupported_or_exploratory_family_never_strict_passes(self):
        exploratory = self.make_export(
            "exploratory",
            "job-exploratory",
            "parent-exploratory",
            [self.analog(qualified=True)],
            exploratory=True,
        )
        bundle = self.build([exploratory])
        self.assertFalse(bundle["analogs"][0]["strict_eligible"])
        self.assertFalse(bundle["candidates"][0]["strict_eligible"])

        unsupported = self.make_export(
            "unsupported",
            "job-unsupported",
            "parent-unsupported",
            [self.analog(qualified=True)],
            exploratory=False,
            result_extra={"sites": []},
        )
        bundle = self.build([unsupported])
        self.assertFalse(bundle["analogs"][0]["strict_eligible"])

    def test_missing_explicit_exploratory_parameter_is_rejected(self):
        export = self.make_export("missing-mode", "job-mode", "parent-mode")
        receipt_path = export / "job-receipt.json"
        receipt = json.loads(receipt_path.read_text())
        del receipt["parameters"]["exploratory"]
        raw = encoded(receipt)
        receipt_path.write_bytes(raw)
        manifest = json.loads((export / "manifest.json").read_text())
        manifest["raw_receipt_sha256"] = hashlib.sha256(raw).hexdigest()
        (export / "manifest.json").write_bytes(encoded(manifest))
        with self.assertRaises(EvidenceValidationError):
            self.build([export])

    def test_nested_diagnostics_do_not_supply_rows_but_closure_is_verified(self):
        nested = {
            "diagnostic": {
                "analogs": [self.analog("blob-injected")],
                "candidates": [self.candidate("blob-injected", "blob-injected")],
            }
        }
        export = self.make_export(
            "diagnostics",
            "job-diagnostics",
            "parent-diagnostics",
            artifacts=[nested],
            result_extra={
                "diagnostic": {
                    "analogs": [self.analog("result-injected")],
                    "candidates": [self.candidate("result-injected", "result-injected")],
                }
            },
        )
        bundle = self.build([export])
        self.assertEqual(1, len(bundle["analogs"]))
        self.assertEqual(1, len(bundle["candidates"]))
        self.assertEqual(1, bundle["evidence_sources"][1]["artifact_count"])

    def test_missing_nested_blob_reference_is_rejected(self):
        missing = artifact_id(99)
        export = self.make_export(
            "missing-ref",
            "job-missing-ref",
            "parent-missing-ref",
            artifacts=[{"nested": self.full_ref(missing, {"missing": True})}],
        )
        with self.assertRaisesRegex(EvidenceValidationError, "absent from manifest"):
            self.build([export])

    def test_reference_with_only_artifact_id_is_invalid(self):
        target = artifact_id(2)
        export = self.make_export(
            "closure",
            "job-closure",
            "parent-closure",
            artifacts=[{"nested": {"artifact_id": target}}, {"value": 1}],
        )
        with self.assertRaisesRegex(EvidenceValidationError, "version"):
            self.build([export])

    def test_complete_artifact_reference_closes_against_manifest(self):
        target_value = {"value": 1}
        target = artifact_id(2)
        export = self.make_export(
            "complete-closure",
            "job-complete-closure",
            "parent-complete-closure",
            artifacts=[{"nested": self.full_ref(target, target_value)}, target_value],
        )
        self.build([export])

    def test_realistic_sites_object_and_schema_id_refs_strictly_qualify(self):
        target_value = {"value": 1}
        target = artifact_id(2)
        analog = self.analog(qualified=True)
        analog["attachment_site_atom_maps"] = [19]
        export = self.make_export(
            "realistic-sites",
            "job-realistic-sites",
            "parent-realistic-sites",
            analogs=[analog],
            exploratory=False,
            artifacts=[{"nested": self.full_ref(target, target_value)}, target_value],
            result_extra={
                "sites": {
                    "atoms": [{"atom_map": 19, "state": "MODIFIABLE"}],
                    "method": "expert_annotation",
                    "diagnostic": {
                        "atoms": [{"atom_map": 1, "state": "MODIFIABLE"}]
                    },
                }
            },
        )
        bundle = self.build([export])
        self.assertTrue(bundle["analogs"][0]["supported_site"])
        self.assertTrue(bundle["analogs"][0]["strict_eligible"])
        self.assertTrue(bundle["candidates"][0]["strict_eligible"])
        self.assertEqual(
            {"19": 1},
            bundle["parents"][0]["counts_recomputed_from_records"][
                "distinct_graphs_by_supported_site"
            ],
        )

    def test_blob_tamper_is_rejected(self):
        export = self.make_export(
            "tamper", "job-tamper", "parent-tamper", artifacts=[{"value": 1}]
        )
        (export / "blobs" / artifact_id(1)).write_bytes(b"{}")
        with self.assertRaisesRegex(EvidenceValidationError, "hash mismatch"):
            self.build([export])

    def test_artifact_traversal_ads_and_duplicate_version_are_rejected(self):
        for index, bad_id in enumerate(("../blob", "a:" + "0" * 32), 1):
            export = self.make_export(
                f"unsafe-{index}", f"job-unsafe-{index}", f"parent-unsafe-{index}"
            )
            manifest = json.loads((export / "manifest.json").read_text())
            manifest["artifacts"] = [{
                "artifact_id": bad_id,
                "version": 1,
                "sha256": "0" * 64,
                "media_type": "application/json",
            }]
            (export / "manifest.json").write_bytes(encoded(manifest))
            with self.assertRaises(EvidenceValidationError):
                self.build([export])

        export = self.make_export(
            "duplicate-artifact",
            "job-duplicate-artifact",
            "parent-duplicate-artifact",
            artifacts=[{"value": 1}],
        )
        manifest = json.loads((export / "manifest.json").read_text())
        manifest["artifacts"].append(dict(manifest["artifacts"][0]))
        (export / "manifest.json").write_bytes(encoded(manifest))
        with self.assertRaisesRegex(EvidenceValidationError, "artifact_id/version"):
            self.build([export])

    def test_duplicate_job_and_parent_are_rejected(self):
        one = self.make_export("one", "job-duplicate", "parent-one")
        two = self.make_export("two", "job-duplicate", "parent-two")
        with self.assertRaisesRegex(EvidenceValidationError, "duplicate job_id"):
            self.build([one, two])

        three = self.make_export("three", "job-three", "parent-duplicate")
        four = self.make_export("four", "job-four", "parent-duplicate")
        with self.assertRaisesRegex(EvidenceValidationError, "duplicate parent_id"):
            self.build([three, four])

    def test_exact_fourteen_unique_known_criteria_are_required(self):
        export = self.make_export("criteria", "job-criteria", "parent-criteria")
        self.write_assessment(
            [{"id": value, "status": "failed"} for value in CRITERION_IDS[:-1]]
        )
        with self.assertRaisesRegex(EvidenceValidationError, "exactly 14"):
            self.build([export])

        duplicate = list(CRITERION_IDS)
        duplicate[-1] = duplicate[0]
        self.write_assessment([{"id": value, "status": "failed"} for value in duplicate])
        with self.assertRaisesRegex(EvidenceValidationError, "unique known ids"):
            self.build([export])

        self.write_assessment(
            [{"id": value, "status": "invented"} for value in CRITERION_IDS]
        )
        with self.assertRaisesRegex(EvidenceValidationError, "criterion status"):
            self.build([export])

        self.write_assessment(
            [{"id": value, "status": "fail"} for value in CRITERION_IDS]
        )
        with self.assertRaisesRegex(EvidenceValidationError, "criterion status"):
            self.build([export])

    def test_failed_status_and_realistic_parent_id_are_preserved(self):
        export = self.make_export("realistic", "job-realistic", "SMARCA2-9QAD-A1I45")
        bundle = self.build([export])
        self.assertEqual("SMARCA2-9QAD-A1I45", bundle["parents"][0]["parent_id"])
        self.assertEqual(
            ["failed"] * 14,
            [row["status"] for row in bundle["assessment_snapshot"]["criteria14"]],
        )

    def test_duplicate_final_result_rows_and_nonlists_are_rejected(self):
        duplicate = self.analog("duplicate")
        export = self.make_export(
            "duplicate-rows",
            "job-duplicate-rows",
            "parent-duplicate-rows",
            analogs=[duplicate, dict(duplicate)],
        )
        with self.assertRaisesRegex(EvidenceValidationError, "duplicate analog record"):
            self.build([export])

        export = self.make_export(
            "invalid-list",
            "job-invalid-list",
            "parent-invalid-list",
            result_extra={"analogs": {"id": "not-a-list"}},
        )
        with self.assertRaisesRegex(EvidenceValidationError, "analogs must be a list"):
            self.build([export])

    def test_nonfinite_json_is_rejected(self):
        self.assessment.write_text(
            '{"criteria":[],"value":NaN}', encoding="utf-8"
        )
        export = self.make_export("nan", "job-nan", "parent-nan")
        with self.assertRaisesRegex(EvidenceValidationError, "nonfinite"):
            self.build([export])

    def test_required_datatypes_are_enforced(self):
        export = self.make_export("types", "job-types", "parent-types")
        receipt_path = export / "job-receipt.json"
        receipt = json.loads(receipt_path.read_text())
        receipt["parameters"]["exploratory"] = 1
        raw = encoded(receipt)
        receipt_path.write_bytes(raw)
        manifest = json.loads((export / "manifest.json").read_text())
        manifest["raw_receipt_sha256"] = hashlib.sha256(raw).hexdigest()
        (export / "manifest.json").write_bytes(encoded(manifest))
        with self.assertRaisesRegex(EvidenceValidationError, "explicitly boolean"):
            self.build([export])

    def test_symlink_is_rejected_when_supported(self):
        export = self.make_export(
            "symlink", "job-symlink", "parent-symlink", artifacts=[{"value": 1}]
        )
        target = export / "blobs" / artifact_id(1)
        replacement = export / "outside.json"
        replacement.write_bytes(target.read_bytes())
        target.unlink()
        try:
            target.symlink_to(replacement)
        except (OSError, NotImplementedError):
            self.skipTest("symlinks unavailable")
        with self.assertRaisesRegex(EvidenceValidationError, "symlink"):
            self.build([export])


if __name__ == "__main__":
    unittest.main()

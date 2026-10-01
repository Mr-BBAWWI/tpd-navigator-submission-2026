"""Service and CPU-ledger contracts for protein-hydrogen evidence."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from packages.contracts import ContractError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.science import scientific_assessment


class ProteinHydrogenPlatformContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(Path(self.temp.name))
        self.service = ScientificAcceptanceService(self.store, "protein-h-test")

    @staticmethod
    def _binding(job_id="job-a"):
        return {
            "project": "protein-h-test",
            "job_id": job_id,
            "input_sha256": "1" * 64,
            "result_sha256": "2" * 64,
            "result_binding_kind": "exact_archived_design_json_bytes",
            "parent_id": "SMARCA2-FX5",
        }

    def _register(self, job_id="job-a", value=None):
        paths = []
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            return self.service._register_evidence(
                db, paths, job_id, "protein_h", "SMARCA2-FX5",
                self._binding(job_id), value or {
                    "status": "computed_diagnostic_pending_human_review",
                    "protein_hydrogens": [],
                },
            )

    def test_review_cannot_fall_back_to_evidence_from_another_job(self):
        registered = self._register("job-a")
        data = {
            "evidence_ref": registered["ref"],
            "source_ref": registered["ref"],
            "accepted": True,
            "rationale": "Reviewed only for job A.",
        }
        result = {"files": {"json": {"sha256": "2" * 64}}}
        with self.store.db() as db:
            with self.assertRaisesRegex(
                    ContractError, "SCIENTIFIC_PROTEIN_H_REVIEW_EVIDENCE"):
                self.service._validate_policy_data(
                    db, "job-b", "protein_hydrogen_review", data, result
                )

    def test_exact_registry_reimport_is_idempotent(self):
        first = self._register()
        second = self._register()
        self.assertTrue(first["created"])
        self.assertFalse(second["created"])
        self.assertEqual(first["id"], second["id"])
        self.assertEqual(first["ref"], second["ref"])

    def test_registered_evidence_blob_tampering_is_rejected(self):
        registered = self._register()
        blob = self.store.root / "blobs" / registered["ref"]["artifact_id"]
        blob.write_bytes(b'{"tampered":true}')
        with self.assertRaisesRegex(ContractError, "STORED_ARTIFACT_CORRUPTED"):
            self.service.port.json(registered["ref"])

    def test_arbitrary_effective_status_cannot_clear_review(self):
        receipt = {
            "requires_review": True,
            "effective_human_review_status": "accepted",
            "original_state_flags": {
                "added_hydrogen_orientation_requires_review": True,
            },
        }
        self.assertTrue(scientific_assessment._receipt_requires_review(
            {"protein_hydrogen_receipt": receipt}
        ))

    def test_authenticated_active_evidence_binding_clears_only_protein_review(self):
        evidence_ref = {
            "artifact_id": "a-evidence", "version": 1,
            "sha256": "3" * 64, "media_type": "application/json",
            "schema_id": "urn:test", "provenance": "computed",
        }
        review = {
            "accepted": True,
            "evidence_ref": evidence_ref,
            "decision_id": "decision-1",
            "actor": {"id": "reviewer-1"},
        }
        receipt = {
            "requires_review": True,
            "effective_human_review_status": "accepted",
            "evidence_id": "evidence-1",
            "evidence_ref": evidence_ref,
            "human_review": review,
            "platform_authenticated_human_review_binding": {
                "platform_authenticated": True,
                "project_id": "protein-h-test",
                "job_id": "job-a",
                "evidence_id": "evidence-1",
                "evidence_ref": evidence_ref,
                "decision_id": "decision-1",
                "actor_id": "reviewer-1",
            },
        }
        self.assertFalse(scientific_assessment._receipt_requires_review(
            {"protein_hydrogen_receipt": receipt}
        ))
        self.assertTrue(scientific_assessment._receipt_requires_review({
            "protein_hydrogen_receipt": receipt,
            "pose_hydrogen_receipt": {
                "requires_review": True,
                "effective_human_review_status": "accepted",
            },
        }))

    def test_cpu_interaction_coverage_requires_authenticated_h_review(self):
        base = {
            "_validated_by_platform": True,
            "job_id": "job-a",
            "current_job_id": "job-a",
            "policy_digest": "policy",
            "analog_id": "analog-1",
            "report_id": "interaction-analog-1",
            "requirements": [{
                "id": "required-1", "required": True,
                "status": "preserved", "computed_pass": True,
                "pending_missing_protein_hydrogen": False,
            }],
            "failures": [],
        }
        arbitrary = dict(base, protein_hydrogen_receipt={
            "requires_review": True,
            "effective_human_review_status": "accepted",
        })
        status, evidence, _, _ = scientific_assessment._covered_interactions(
            [arbitrary], ["analog-1"], {}, "policy", {"job_id": "job-a"}
        )
        self.assertEqual(status, "pending")
        self.assertTrue(evidence["pending"])


if __name__ == "__main__":
    unittest.main()

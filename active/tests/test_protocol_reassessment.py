"""Synthetic integration tests for immutable protocol reassessment.

Every protocol record in this module is a synthetic contract fixture and makes no
claim about scientific performance.
"""
from __future__ import annotations

import copy
import json
import math
import subprocess
import sys
import unittest
import uuid
from unittest import mock

from packages.contracts import ContractError
from packages.platform import protocol_evidence, protocol_reassessment
from packages.platform.protocol_evidence import TABLE_SQL
from tests.test_scientific_acceptance_service import ScientificAcceptanceServiceTests


class ProtocolReassessmentTests(ScientificAcceptanceServiceTests):
    def test_scientific_acceptance_imports_in_fresh_process(self):
        completed = subprocess.run(
            [
                sys.executable,
                "-c",
                "from packages.platform.scientific_acceptance import "
                "ScientificAcceptanceService; print(ScientificAcceptanceService.__name__)",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout.strip(), "ScientificAcceptanceService")

    def _known_summary(self, selected_pass=2):
        seeds = [23, 41, 61, 79, 97]
        selected = [{
            "seed": seed,
            "selected_model_index": 0,
            "selected_e3_CA_RMSD_after_target_alignment_A": float(index + 1),
            "selected_pass": index < selected_pass,
        } for index, seed in enumerate(seeds)]
        return {
            "benchmark_scope":
                "BRD4-CRBN 6BOY related diagnostic; not the SMARCA2 novel target",
            "ligand_scope": "Synthetic RN6-scoped fixture only.",
            "raw_pass_count": 16,
            "raw_model_denominator": 25,
            "selected_pass_count": selected_pass,
            "selected_seed_denominator": 5,
            "claimed_protocol_pass": selected_pass >= 3,
            "failure_count": 0,
            "selected": selected,
            "preserved_original_baseline": f"{selected_pass}_of_5",
        }

    def _register_known(self, protocol_id="calibration:synthetic", *,
                        selected_pass=2, job_id=None, project=None,
                        engine_hashes=None):
        job_id = job_id or self.job_id
        project = project or self.project
        with self.store.db() as db:
            db.execute(TABLE_SQL)
            _, _result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=False
            )
        binding = dict(binding, project=project, job_id=job_id)
        archive_ref = self.store.scope(project).put_raw(
            b"synthetic protocol archive", "application/zip", "computed"
        )
        identifier = "pdiag-" + uuid.uuid4().hex
        created_at = "2026-10-03T00:00:00+00:00"
        record = {
            "id": identifier,
            "project_id": project,
            "job_id": job_id,
            "protocol_id": protocol_id,
            "kind": protocol_reassessment.KNOWN_KIND,
            "created_at": created_at,
            "source_binding": binding,
            "assessment_link": None,
            "summary": self._known_summary(selected_pass),
            "measurement_verification": {
                "verification_level": "synthetic strict fixture",
                "geometry_recomputed": True,
                "stored_metrics_only": False,
                "claimed_pass_reported_separately": True,
            },
            "pack_manifest_sha256": "1" * 64,
            "compact_sha256": "2" * 64,
            "archive_ref": archive_ref,
            "verification_engine": {
                "version": 1,
                "source_hashes": engine_hashes or
                    protocol_reassessment.current_verification_source_hashes(),
            },
            "scientific_approved": False,
            "formal_acceptance": False,
            "gates_affected": False,
            "status": protocol_evidence.STATUS,
            "operator_kind": protocol_evidence.OPERATOR_KIND,
        }
        raw = protocol_evidence._canonical(record)
        record_ref = self.store.scope(project).put_raw(
            raw, "application/json", "computed"
        )
        with self.store.db() as db:
            db.execute(
                """INSERT INTO protocol_diagnostics(
                    id,project,job_id,protocol_id,kind,pack_manifest_sha256,
                    binding,record_ref,archive_ref,created_at)
                    VALUES(?,?,?,?,?,?,?,?,?,?)""",
                (identifier, project, job_id, protocol_id,
                 protocol_reassessment.KNOWN_KIND, "1" * 64,
                 json.dumps(binding, sort_keys=True),
                 json.dumps(record_ref, sort_keys=True),
                 json.dumps(archive_ref, sort_keys=True), created_at),
            )
        return identifier, record_ref, archive_ref

    def test_no_selection_preserves_create_behavior(self):
        assessment, created = self.service.create(self.job_id)
        repeated, created_again = self.service.create(self.job_id)
        self.assertTrue(created)
        self.assertFalse(created_again)
        self.assertEqual(assessment["id"], repeated["id"])
        self.assertNotIn("protocol_evaluation", assessment)

    def test_exact_two_of_five_fails_and_old_rows_are_immutable(self):
        diagnostic_id, record_ref, archive_ref = self._register_known(selected_pass=2)
        with self.store.db() as db:
            before = dict(db.execute(
                "SELECT * FROM protocol_diagnostics WHERE id=?", (diagnostic_id,)
            ).fetchone())
        assessment, created = self.service.create(
            self.job_id, protocol_ids=["calibration:synthetic"]
        )
        self.assertTrue(created)
        criterion = {row["id"]: row for row in assessment["criteria"]}[
            "known_crbn_calibration"]
        current = criterion["observed"]["current_protocol"]
        self.assertEqual(criterion["status"], "failed")
        self.assertEqual(
            current["machine_subchecks"]["selected_current_quality"]["observed_pass"], 2
        )
        self.assertEqual(
            current["machine_subchecks"]["raw_model_accounting"]["observed_pass"], 16
        )
        self.assertIn("CURRENT selected top-confidence result is 2/5", criterion["reason"])
        self.assertEqual(assessment["protocol_evaluation"]["records"][0]["record_ref"],
                         record_ref)
        self.assertEqual(assessment["protocol_evaluation"]["records"][0]["archive_ref"],
                         archive_ref)
        with self.store.db() as db:
            after = dict(db.execute(
                "SELECT * FROM protocol_diagnostics WHERE id=?", (diagnostic_id,)
            ).fetchone())
        self.assertEqual(before, after)

    def test_exact_three_of_five_passes_numerically_but_scope_remains_pending(self):
        self._register_known(
            "calibration:three-of-five", selected_pass=3
        )
        assessment, created = self.service.create(
            self.job_id, protocol_ids=["calibration:three-of-five"]
        )
        self.assertTrue(created)
        criterion = {row["id"]: row for row in assessment["criteria"]}[
            "known_crbn_calibration"
        ]
        current = criterion["observed"]["current_protocol"]
        self.assertEqual(current["numerical_status"], "pass")
        self.assertEqual(current["scope_status"],
                         "scope_review_not_implemented")
        self.assertEqual(current["effective_status"], "pending")
        self.assertEqual(criterion["status"], "pending")
        quality = current["machine_subchecks"]["selected_current_quality"]
        self.assertEqual(quality["observed_pass"], 3)
        self.assertEqual(quality["observed_denominator"], 5)
        self.assertEqual(quality["required_pass"], 3)
        self.assertEqual(quality["trusted_expert_cutoff"],
                         "at_least_3_of_5")
        original_claim = current["machine_subchecks"][
            "original_record_protocol_claim"
        ]
        self.assertTrue(original_claim["claimed_protocol_pass"])
        self.assertEqual(original_claim["claim_rule_required_pass"], 3)
        self.assertIn("최소 3/5", current["reason"])

    def test_same_choice_is_idempotent_and_different_choice_is_new_revision(self):
        self._register_known("calibration:first")
        first, made = self.service.create(
            self.job_id, protocol_ids=["calibration:first"]
        )
        repeated, repeated_made = self.service.create(
            self.job_id, protocol_ids=["calibration:first"]
        )
        self.assertTrue(made)
        self.assertFalse(repeated_made)
        self.assertEqual(first["id"], repeated["id"])
        with self.store.db() as db:
            db.execute("DELETE FROM protocol_diagnostics WHERE protocol_id=?",
                       ("calibration:first",))
        self._register_known("calibration:second", selected_pass=5)
        second, second_made = self.service.create(
            self.job_id, protocol_ids=["calibration:second"]
        )
        self.assertTrue(second_made)
        self.assertNotEqual(first["id"], second["id"])
        self.assertGreater(second["revision"], first["revision"])
        criterion = {row["id"]: row for row in second["criteria"]}[
            "known_crbn_calibration"]
        self.assertEqual(criterion["status"], "pending")
        self.assertEqual(criterion["observed"]["prior_baseline"]["status"], "pending")
        self.assertEqual(
            criterion["observed"]["current_protocol"]["scope_status"],
            "scope_review_not_implemented",
        )

    def test_multiple_same_kind_is_rejected(self):
        self._register_known("calibration:a")
        self._register_known("calibration:b")
        with self.assertRaisesRegex(
            ContractError, "PROTOCOL_REASSESSMENT_AMBIGUOUS_PROTOCOL_KIND"
        ):
            self.service.create(
                self.job_id, protocol_ids=["calibration:a", "calibration:b"]
            )

    def test_missing_wrong_job_cross_project_and_stale_engine_fail(self):
        with self.store.db() as db:
            db.execute(TABLE_SQL)
        with self.assertRaisesRegex(ContractError, "PROTOCOL_REASSESSMENT_PROTOCOL_NOT_FOUND"):
            self.service.create(self.job_id, protocol_ids=["missing"])
        self._register_known("calibration:wrong-job", job_id="wrong-job")
        with self.assertRaisesRegex(ContractError, "PROTOCOL_REASSESSMENT_PROTOCOL_NOT_FOUND"):
            self.service.create(self.job_id, protocol_ids=["calibration:wrong-job"])
        self._register_known("calibration:foreign", project="foreign-project")
        with self.assertRaisesRegex(ContractError, "PROTOCOL_REASSESSMENT_PROTOCOL_NOT_FOUND"):
            self.service.create(self.job_id, protocol_ids=["calibration:foreign"])
        self._register_known("calibration:stale", engine_hashes={
            **protocol_reassessment.current_verification_source_hashes(),
            "packages/science/protocol_evidence.py": "0" * 64,
        })
        with self.assertRaisesRegex(
            ContractError, "PROTOCOL_REASSESSMENT_STALE_VERIFICATION_ENGINE"
        ):
            self.service.create(self.job_id, protocol_ids=["calibration:stale"])

    def test_nonfinite_bool_duplicate_and_mismatched_counts_fail(self):
        cases = []
        value = self._known_summary()
        value["selected"][0]["selected_e3_CA_RMSD_after_target_alignment_A"] = math.inf
        cases.append(value)
        value = self._known_summary()
        value["selected_pass_count"] = True
        cases.append(value)
        value = self._known_summary()
        value["selected"][1]["seed"] = value["selected"][0]["seed"]
        cases.append(value)
        value = self._known_summary()
        value["selected_pass_count"] = 3
        cases.append(value)
        for summary in cases:
            record = {
                "summary": summary,
                "measurement_verification": {
                    "geometry_recomputed": True,
                    "stored_metrics_only": False,
                    "claimed_pass_reported_separately": True,
                },
            }
            with self.assertRaises(ContractError):
                protocol_reassessment._known_evaluation(record)

    def test_pending_action_categories_distinguish_development_and_evidence(self):
        assessment = {
            "criteria": [
                {"id": "parent_funnel", "status": "failed",
                 "needs_expert": True, "reason": "parent evidence missing"},
                {"id": "broad_chemotype_families", "status": "pending",
                 "needs_expert": True, "reason": "family evidence missing"},
                {"id": "qualified_panel_size", "status": "failed",
                 "needs_expert": True, "reason": "panel evidence missing"},
                {"id": "core_interaction_requirements", "status": "pending",
                 "needs_expert": True, "reason": "interaction evidence missing"},
                {"id": "microstate_population", "status": "pending",
                 "needs_expert": True, "reason": "population evidence missing"},
            ]
        }
        protocol_reassessment._recompute_summary(assessment)
        categories = {
            row["criterion_id"]: row["category"]
            for row in assessment["pending_actions"]
        }
        self.assertEqual(categories["parent_funnel"],
                         "additional_development_needed")
        self.assertEqual(categories["broad_chemotype_families"],
                         "additional_development_needed")
        self.assertEqual(categories["qualified_panel_size"],
                         "additional_development_needed")
        self.assertEqual(categories["core_interaction_requirements"],
                         "evidence_missing")
        self.assertEqual(categories["microstate_population"],
                         "evidence_missing")

    def test_non_object_observed_values_do_not_break_pending_action_classification(self):
        for observed in (None, [], "not-structured"):
            with self.subTest(observed=observed):
                assessment = {
                    "criteria": [{
                        "id": "formal_expert_decision",
                        "status": "pending",
                        "needs_expert": True,
                        "observed": observed,
                        "reason": "Formal expert decision remains pending.",
                    }]
                }
                protocol_reassessment._recompute_summary(assessment)
                self.assertEqual(len(assessment["pending_actions"]), 1)
                self.assertEqual(
                    assessment["pending_actions"][0]["category"],
                    "human_review_required",
                )
                self.assertFalse(assessment["scientific_accepted"])

    def test_explicit_human_reject_is_preserved_but_nested_numeric_failure_is_not(self):
        self._register_known("calibration:better", selected_pass=5)
        original = protocol_reassessment._human_disposition
        with mock.patch.object(
            protocol_reassessment, "_human_disposition",
            side_effect=original,
        ):
            assessment, _ = self.service.create(
                self.job_id, protocol_ids=["calibration:better"]
            )
        criterion = {row["id"]: row for row in assessment["criteria"]}[
            "known_crbn_calibration"]
        self.assertEqual(criterion["status"], "pending")
        self.assertEqual(criterion["observed"]["prior_baseline"]["status"],
                         "pending")

        closure = {
            "records": [{
                "kind": protocol_reassessment.KNOWN_KIND,
                "evaluation": copy.deepcopy(
                    criterion["observed"]["current_protocol"]),
            }]
        }
        numeric_failure = {
            "criteria": [{
                "id": "known_crbn_calibration",
                "status": "failed",
                "reason": "Synthetic numerical baseline failed.",
                "observed": {
                    "current_protocol": {"numerical_status": "failed"},
                },
            }]
        }
        protocol_reassessment.apply(numeric_failure, closure)
        updated = numeric_failure["criteria"][0]
        self.assertEqual(updated["observed"]["prior_baseline"]["status"],
                         "failed")
        self.assertEqual(updated["status"], "pending")

        explicit_reject = {
            "criteria": [{
                "id": "known_crbn_calibration",
                "status": "failed",
                "reason": "Expert geometry rejected.",
                "observed": {"disposition": "reject"},
                "evidence": {},
            }]
        }
        protocol_reassessment.apply(explicit_reject, closure)
        self.assertEqual(explicit_reject["criteria"][0]["status"], "failed")
        self.assertEqual(explicit_reject["criteria"][0]["reason"],
                         "Expert geometry rejected.")

        unexpected = {
            "status": "failed",
            "reason": "Unexpected structured value.",
            "observed": {"disposition": []},
        }
        self.assertEqual(protocol_reassessment._human_disposition(unexpected),
                         (None, None))

    def test_view_freshness_does_not_rerun_pack_geometry_and_tamper_fails_closed(self):
        _, record_ref, archive_ref = self._register_known("calibration:view")
        assessment, _ = self.service.create(
            self.job_id, protocol_ids=["calibration:view"]
        )
        with mock.patch(
            "packages.platform.protocol_evidence.verify_pack",
            side_effect=AssertionError("geometry verifier must not run on view"),
        ):
            viewed = self.service.view(assessment["id"])
        self.assertTrue(viewed["protocol_evidence_current"])
        blob = self.store.root / "blobs" / archive_ref["artifact_id"]
        original = blob.read_bytes()
        blob.write_bytes(original + b"tamper")
        viewed = self.service.view(assessment["id"])
        self.assertFalse(viewed["protocol_evidence_current"])
        self.assertFalse(viewed["scientific_accepted"])
        criterion = {row["id"]: row for row in viewed["criteria"]}[
            "known_crbn_calibration"]
        self.assertIn(criterion["status"], {"failed", "pending"})

    def test_missing_registry_table_has_contract_error(self):
        with self.store.db() as db:
            db.execute("DROP TABLE IF EXISTS protocol_diagnostics")
        with self.assertRaisesRegex(
            ContractError, "PROTOCOL_REASSESSMENT_REGISTRY_MISSING"
        ):
            self.service.create(self.job_id, protocol_ids=["calibration:none"])


if __name__ == "__main__":
    unittest.main()

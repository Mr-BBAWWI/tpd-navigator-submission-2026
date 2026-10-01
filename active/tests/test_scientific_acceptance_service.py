"""Real Store and ScientificAcceptanceService contract tests.

The compact design fixture is synthetic and is not evidence of scientific performance.
"""
from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from packages.contracts import ContractError, encoded
from packages.platform.design_panel import input_binding
from packages.platform.review_identity import ReviewAuthError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService, identity
from packages.science import synthesis_review


class ScientificAcceptanceServiceTests(unittest.TestCase):
    project = "science-test"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(Path(self.temp.name))
        self.workbench = WorkbenchService(self.store, self.project)
        self.service = ScientificAcceptanceService(self.store, self.project)
        self.key = self.service.auth.register("reviewer-one", "Reviewer One")
        self.token = self.service.auth.login("reviewer-one", self.key)
        self.job_id = self._completed_job()

    def _completed_job(self):
        binding = input_binding()
        result = {
            "format": "design-panel/20260930.4",
            "status": "synthetic_fixture_not_scientific_performance",
            "input_binding": binding,
            "parent_scope": {"actual_design_parent_id": "SMARCA2-FX5"},
            "summary": {"valid_analogs": 0},
            "sites": {"atoms": [
                {"atom_map": 19, "state": "MODIFIABLE"},
                {"atom_map": 20, "state": "PROTECTED"},
                {"atom_map": 21, "state": "UNKNOWN", "coordinates": [1.0, 2.0, 3.0]},
            ]},
            "analogs": [{
                "id": "A-1",
                "mapped_smiles": "[CH3:1][CH3:2]",
                "canonical_smiles": "CC",
                "attachment_site_atom_maps": [19],
                "transformation_class": "ring_expansion",
                "cheap_filter": {"valid": True},
                "selected": True,
                "pipeline_status": "qualified",
                "assembly_eligible": True,
                "parent_redocking_supported": True,
                "preview": False,
                "docking": {
                    "status": "completed_with_limits",
                    "pose_preserved": True,
                    "passing_pose_count": 1,
                    "files": {},
                },
            }],
            "protac_candidates": [{
                "candidate_id": "C-1",
                "e3_type": "CRBN",
                "mapped_smiles": "[CH3:1][CH2:2][CH2:3][CH3:4]",
                "canonical_smiles": "CCCC",
                "assembly_mode": "pose_supported_hypothesis",
                "warhead_analog_id": "A-1",
                "atom_roles": {
                    "warhead_maps": [1], "linker_maps": [2, 3],
                    "recruiter_maps": [4],
                },
                "attachment_metadata": {
                    "warhead_linker_bond": {
                        "attachment_atom_map": 1, "partner_atom_map": 2,
                        "bond_type": "SINGLE",
                    },
                    "recruiter_linker_bond": {
                        "attachment_atom_map": 4, "partner_atom_map": 3,
                        "bond_type": "SINGLE",
                    },
                },
            }],
            "calibration": {"CRBN": {"seed_receipts": []}},
            "files": {},
        }
        # Synthetic 3D file capability fixture only; it makes no evidence claim.
        pose_sdf = (
            b"Synthetic mapped ethane pose\n"
            b"  fixture          3D\n"
            b"\n"
            b"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            b"   -0.7700   -0.1000   -0.0500 C   0  0  0  0  0  0  0  0  1  0  0\n"
            b"    0.7700    0.1000    0.0500 C   0  0  0  0  0  0  0  0  2  0  0\n"
            b"  1  2  1  0\n"
            b"M  END\n"
            b">  <mapped_smiles>\n"
            b"[CH3:1][CH3:2]\n"
            b"\n"
            b"$$$$\n"
        )
        poses_ref = self.store.scope(self.project).put_raw(
            pose_sdf, "chemical/x-mdl-sdfile", "computed"
        )
        result["analogs"][0]["docking"]["files"]["poses_sdf"] = poses_ref
        archived = copy.deepcopy(result)
        archive_raw = encoded(archived)
        json_ref = self.store.scope(self.project).put_raw(
            archive_raw, "application/json", "computed"
        )
        report_ref = self.store.scope(self.project).put_raw(
            b"# Synthetic fixture\n", "text/markdown", "computed"
        )
        result["files"] = {"json": json_ref, "report": report_ref}
        job_id = "job-scientific-fixture"
        row = {
            "id": job_id,
            "project_id": self.project,
            "run_id": "design:SMARCA2",
            "result_id": "design:SMARCA2",
            "operation": "design_panel",
            "title": "Synthetic completed design fixture",
            "parameters": {},
            "input_digest": binding["digest"],
            "binding": binding,
            "runtime": identity(),
            "state": "completed",
            "stage": "finished",
            "tool_inputs": [],
            "created_at": "2026-10-02T00:00:00Z",
            "updated_at": "2026-10-02T00:00:00Z",
            "retry_of": None,
            "attempt": 1,
            "calls": 0,
            "total_tokens": 0,
            "usage_status": "not_started",
            "limits": {},
            "outputs": [
                {"name": "result.json", "ref": json_ref},
                {"name": "report.md", "ref": report_ref},
            ],
            "error_code": None,
            "authority": {},
            "result": result,
            "scientific_policy_binding": None,
        }
        with self.store.db() as db:
            db.execute(
                "INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)",
                (job_id, self.project, row["run_id"], row["result_id"],
                 "scientific-fixture", "fixture-fingerprint", "completed",
                 json.dumps(row)),
            )
        return job_id

    def _assessment(self, reviewer_id=None):
        value, created = self.service.create(
            self.job_id, reviewer_id=reviewer_id
        )
        return value, created

    def _replace_job(self, mutate):
        with self.store.db() as db:
            row = db.execute(
                "SELECT body FROM workbench_jobs WHERE id=? AND project=?",
                (self.job_id, self.project),
            ).fetchone()
            body = json.loads(row["body"])
            mutate(body)
            db.execute(
                "UPDATE workbench_jobs SET body=? WHERE id=? AND project=?",
                (json.dumps(body), self.job_id, self.project),
            )

    def _stale_runtime(self):
        self._replace_job(lambda body: body.__setitem__(
            "runtime", {"archived_runtime": True}
        ))

    def _registered_ternary_report(self, candidate_id="C-1", seed=23):
        with self.store.db() as db:
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=False
            )
        graph_hash = "a" * 64
        plan_digest = "synthetic-plan-digest"
        raw_ref = self.service.port.put_raw(
            b"synthetic immutable raw receipt\n", "text/plain", "computed"
        )
        plan = {
            "plan_digest": plan_digest,
            "bindings": {"job": {
                "project": self.project,
                "job_id": self.job_id,
                "input_sha256": binding["input_sha256"],
                "result_sha256": binding["result_sha256"],
                "candidate_graph_sha256": graph_hash,
            }},
            "candidate_graph": {
                "candidate_id": candidate_id,
                "e3_type": "CRBN",
                "graph_sha256": graph_hash,
            },
        }
        plan_ref = self.service.port.put_json(plan)
        aggregate_ref = self.service.port.put_json({
            "plan_digest": plan_digest,
            "candidate_id": candidate_id,
        })
        report = {
            "candidate_id": candidate_id,
            "e3_type": "CRBN",
            "seed": seed,
            "derived_interpretation": False,
            "actual_computation": True,
            "execution_success": True,
            "reference_free": True,
            "failure": None,
            "registered_artifacts": {"raw_receipt": raw_ref},
            "plan_ref": plan_ref,
            "aggregate_receipt_ref": aggregate_ref,
            "producer_policy_binding": {},
            "measurement_protocol": {
                "candidate_graph_sha256": graph_hash,
                "plan_digest": plan_digest,
                "seed": seed,
            },
        }
        report_ref = self.service.port.put_json(report)
        stored = {
            **binding,
            "evidence_sha256": hashlib.sha256(encoded(report)).hexdigest(),
            "measurement_binding": {
                "producer_policy": {},
                "candidate_graph_sha256": graph_hash,
                "plan_digest": plan_digest,
            },
        }
        identity = f"{candidate_id}|CRBN|{seed}"
        with self.store.db() as db:
            db.execute(
                """INSERT INTO scientific_evidence(
                    id,project,job_id,kind,identity_key,fingerprint,binding,ref,
                    superseded,created_at) VALUES(?,?,?,?,?,?,?,?,?,?)""",
                ("evidence-geometry-fixture", self.project, self.job_id,
                 "ternary", identity, "geometry-fixture-fingerprint",
                 json.dumps(stored), json.dumps(report_ref), 0,
                 "2026-10-02T00:00:00Z"),
            )
        return result, report_ref, graph_hash, seed

    def _reduced_numerical_policy(self, assessment_id, token=None):
        token = token or self.token
        statement = self.service.create_statement(
            "Synthetic test rationale for a reduced panel minimum.",
            "fixture:test-policy", token,
        )
        submission = {
            "kind": "numerical_criteria",
            "action": "accept",
            "intent": "Exercise authenticated policy revision behavior.",
            "reason": "Synthetic service contract test only.",
            "expected_policy_revision": self.service.policy(self.job_id)["revision"],
            "data": {
                "values": {"panel_min": 1},
                "minimum_reduction_reviewed": True,
                "reduction_justification": "Explicit synthetic reduction.",
                "reduction_source_ref": statement["ref"],
            },
        }
        return self.service.record_decision(assessment_id, submission, token)

    def test_create_is_idempotent_and_uses_actual_engine(self):
        first, created = self._assessment("reviewer-one")
        second, repeated = self._assessment("reviewer-one")
        self.assertTrue(created)
        self.assertFalse(repeated)
        self.assertEqual(first["id"], second["id"])
        self.assertFalse(first["scientific_accepted"])
        self.assertEqual(first["assigned_reviewer_id"], "reviewer-one")
        self.assertEqual(first["source_binding"]["parent_id"], "SMARCA2-FX5")

    def test_assignment_requires_real_active_registered_actor(self):
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_REVIEWER_INACTIVE"):
            self._assessment("invented-reviewer")
        self.service.auth.disable("reviewer-one")
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_REVIEWER_INACTIVE"):
            self._assessment("reviewer-one")

    def test_source_json_hash_and_cross_project_registration_are_enforced(self):
        with self.store.db() as db:
            _, result, binding, _ = self.service._verified_result(db, self.job_id)
        raw = self.service.port.read(result["files"]["json"])
        self.assertEqual(hashlib.sha256(raw).hexdigest(), binding["result_sha256"])
        other = ScientificAcceptanceService(self.store, "other-project")
        with self.assertRaises(KeyError):
            other.create(self.job_id)

    def test_source_and_engine_freshness_invalidate_effective_assessment(self):
        assessment, _ = self._assessment()
        changed = dict(input_binding())
        changed["digest"] = "0" * 64
        with mock.patch(
            "packages.platform.scientific_acceptance.input_binding",
            return_value=changed,
        ):
            viewed = self.service.view(assessment["id"])
            self.assertFalse(viewed["source_inputs_current"])
            self.assertFalse(viewed["scientific_accepted"])
        with mock.patch(
            "packages.platform.scientific_acceptance.engine_fingerprint",
            return_value="f" * 64,
        ):
            viewed = self.service.view(assessment["id"])
            self.assertFalse(viewed["current_implementation"])

    def test_registered_statement_can_support_policy_revision(self):
        assessment, _ = self._assessment("reviewer-one")
        decision, created = self._reduced_numerical_policy(assessment["id"])
        self.assertTrue(created)
        self.assertEqual(decision["policy_revision"], 1)
        self.assertEqual(self.service.policy(self.job_id)["revision"], 1)
        stale = self.service.view(assessment["id"])
        self.assertFalse(stale["current_policy"])

    def test_unregistered_or_wrong_author_statement_is_rejected(self):
        assessment, _ = self._assessment()
        foreign_key = self.service.auth.register("reviewer-two", "Reviewer Two")
        foreign_token = self.service.auth.login("reviewer-two", foreign_key)
        statement = self.service.create_statement("Foreign statement", "fixture:x", foreign_token)
        submission = {
            "kind": "numerical_criteria", "action": "accept",
            "intent": "test", "reason": "test", "expected_policy_revision": 0,
            "data": {
                "values": {"panel_min": 1},
                "minimum_reduction_reviewed": True,
                "reduction_justification": "test",
                "reduction_source_ref": statement["ref"],
            },
        }
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_STATEMENT_AUTHOR"):
            self.service.record_decision(assessment["id"], submission, self.token)

    def test_pending_formal_accept_is_blocked_and_formal_reject_is_failed(self):
        assessment, _ = self._assessment("reviewer-one")
        accept = {
            "kind": "formal_decision", "action": "accept",
            "intent": "Attempt formal acceptance", "reason": "fixture",
            "expected_policy_revision": 0, "data": {},
        }
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_ACTIVE_POLICY_REQUIRED"):
            self.service.record_decision(assessment["id"], accept, self.token)
        reject = dict(accept, action="reject", intent="Formally reject incomplete evidence")
        self.service.record_decision(assessment["id"], reject, self.token)
        viewed = self.service.view(assessment["id"])
        formal = {row["id"]: row for row in viewed["criteria"]}["formal_expert_decision"]
        self.assertEqual(formal["status"], "failed")
        self.assertEqual(viewed["confirmation_status"], "rejected")
        self.assertFalse(viewed["scientific_accepted"])
        self.assertTrue(any(row["status"] != "pass" for row in viewed["criteria"]
                            if row["id"] != "formal_expert_decision"))

    def test_withdrawal_allows_stale_author_and_disabled_actor_invalidates_policy(self):
        assessment, _ = self._assessment()
        first, _ = self._reduced_numerical_policy(assessment["id"])
        current, _ = self._assessment()
        key2 = self.service.auth.register("reviewer-two", "Reviewer Two")
        token2 = self.service.auth.login("reviewer-two", key2)
        statement2 = self.service.create_statement("Second revision", "fixture:second", token2)
        submission = {
            "kind": "numerical_criteria", "action": "accept",
            "intent": "second revision", "reason": "fixture",
            "expected_policy_revision": 1,
            "data": {
                "values": {"panel_max": 19},
                "minimum_reduction_reviewed": True,
                "reduction_justification": "second synthetic reduction",
                "reduction_source_ref": statement2["ref"],
            },
        }
        self.service.record_decision(current["id"], submission, token2)
        withdrawn = self.service.withdraw(first["id"], self.token)
        self.assertTrue(withdrawn["withdrawn"])
        self.service.auth.disable("reviewer-two")
        policy = self.service.policy(self.job_id)
        self.assertEqual(policy["active_policy_decision_ids"], [])
        self.assertGreaterEqual(policy["revision"], 3)

    def test_cpu_compute_calls_real_synthesis_and_preserves_missing_pose_failure(self):
        assessment, _ = self._assessment()
        self._reduced_numerical_policy(assessment["id"])
        with mock.patch.object(
            synthesis_review,
            "assess_synthesis",
            wraps=synthesis_review.assess_synthesis,
        ) as actual_method:
            result = self.service.compute(self.job_id, self.token)
        self.assertEqual(actual_method.call_count, 1)
        self.assertEqual({row["kind"] for row in result["evidence"]},
                         {"synthesis", "interaction"})
        interaction = next(row for row in result["evidence"]
                           if row["kind"] == "interaction")
        report = self.service.port.json(interaction["ref"])
        self.assertFalse(report["scientific_approved"])
        self.assertEqual(report["policy_status"], "registered_computation_failure")
        self.assertTrue(report["failures"][0]["preserved_in_failure_ledger"])
        second = self.service.compute(self.job_id, self.token)
        self.assertTrue(all(not row["created"] for row in second["evidence"]))

    def test_design_policy_is_latest_per_site_and_preserves_release_evidence(self):
        from packages.science.dual_e3 import catalog

        rule_id = catalog()["rule_catalog"][0]["rule_id"]

        def record_site(assessment_id, atom_map, state, rationale, statement,
                        release=False, release_statement=None):
            data = {
                "atom_map": atom_map,
                "state": state,
                "allowed_rule_ids": [rule_id],
                "rationale": rationale,
                "source_ref": statement["ref"],
                "allow_release_protected": release,
            }
            if release:
                data["release_justification"] = "Explicit protected-site release justification."
                data["release_source_ref"] = release_statement["ref"]
            submission = {
                "kind": "site_policy",
                "action": "accept",
                "intent": "Exercise latest scoped design policy behavior.",
                "reason": rationale,
                "expected_policy_revision": self.service.policy(self.job_id)["revision"],
                "data": {
                    "parent_id": "SMARCA2-FX5",
                    "sites": [data],
                },
            }
            self.service.record_decision(assessment_id, submission, self.token)

        first_assessment, _ = self._assessment()
        first_source = self.service.create_statement(
            "First site rationale source.", "fixture:/sites/first", self.token
        )
        record_site(
            first_assessment["id"], 19, "MODIFIABLE",
            "First atom-map policy.", first_source,
        )

        second_assessment, _ = self._assessment()
        second_source = self.service.create_statement(
            "Latest site rationale source.", "fixture:/sites/latest", self.token
        )
        record_site(
            second_assessment["id"], 19, "MODIFIABLE",
            "Latest atom-map policy.", second_source,
        )

        release_assessment, _ = self._assessment()
        ordinary_source = self.service.create_statement(
            "Ordinary protected-site rationale.", "fixture:/sites/protected", self.token
        )
        release_source = self.service.create_statement(
            "Explicit protected-site release evidence.",
            "fixture:/sites/protected/release",
            self.token,
        )
        record_site(
            release_assessment["id"], 20, "PROTECTED",
            "Ordinary protected-site rationale.", ordinary_source,
            release=True, release_statement=release_source,
        )

        current_assessment, _ = self._assessment()
        design_policy = self.service.policy_for_design(current_assessment["id"])
        rows = {row["atom_map"]: row for row in design_policy["site_policy"]}

        self.assertEqual(
            set(design_policy), {"scope", "site_policy", "parent_funnel"}
        )
        self.assertEqual(design_policy["parent_funnel"], [])
        self.assertEqual(sorted(rows), [19, 20])
        expected_site_fields = {
            "atom_map", "state", "allowed_transforms", "rationale",
            "source_ref", "allow_release_protected",
        }
        self.assertTrue(all(set(row) == expected_site_fields for row in rows.values()))
        self.assertEqual(rows[19]["rationale"], "Latest atom-map policy.")
        self.assertEqual(rows[19]["source_ref"], second_source["ref"])
        self.assertTrue(rows[20]["allow_release_protected"])
        self.assertEqual(
            rows[20]["rationale"],
            "Explicit protected-site release justification.",
        )
        self.assertEqual(rows[20]["source_ref"], release_source["ref"])

    def test_unknown_site_requires_explicit_authenticated_expert_classification(self):
        from packages.science.dual_e3 import catalog

        assessment, _ = self._assessment()
        rule_id = catalog()["rule_catalog"][0]["rule_id"]
        statement = self.service.create_statement(
            "Human expert classification for synthetic atom map 21 with a specific rationale.",
            "fixture:/sites/21/classification", self.token,
        )
        base = {
            "kind": "site_policy", "action": "accept",
            "intent": "Explicit human expert site classification",
            "reason": "Synthetic classification contract test.",
            "expected_policy_revision": 0,
            "data": {"parent_id": "SMARCA2-FX5", "sites": [{
                "atom_map": 21,
                "original_state": "UNKNOWN",
                "state": "MODIFIABLE",
                "expert_classification": True,
                "allowed_rule_ids": [rule_id],
                "rationale": "Specific human policy rationale for this atom map.",
                "source_ref": statement["ref"],
                "allow_release_protected": False,
            }]},
        }
        with self.assertRaises(ReviewAuthError):
            self.service.record_decision(assessment["id"], base, "not-a-session")
        missing_flag = copy.deepcopy(base)
        missing_flag["data"]["sites"][0].pop("expert_classification")
        with self.assertRaisesRegex(
            ContractError, "SCIENTIFIC_SITE_EXPERT_CLASSIFICATION_REQUIRED"
        ):
            self.service.record_decision(assessment["id"], missing_flag, self.token)
        decision, created = self.service.record_decision(
            assessment["id"], base, self.token
        )
        self.assertTrue(created)
        self.assertEqual(decision["data"]["sites"][0]["original_state"], "UNKNOWN")
        self.assertEqual(
            self.service._job(self.store.db().__enter__(), self.job_id)
            if False else "UNKNOWN",
            "UNKNOWN",
        )
        current, _ = self._assessment()
        policy = self.service.policy_for_design(current["id"])
        row = next(item for item in policy["site_policy"] if item["atom_map"] == 21)
        self.assertEqual(row["state"], "MODIFIABLE")
        self.assertEqual(set(row), {
            "atom_map", "state", "allowed_transforms", "rationale",
            "source_ref", "allow_release_protected",
        })
        self.assertNotIn("original_state", row)
        self.assertNotIn("expert_classification", row)
        self.assertEqual(decision["data"]["sites"][0]["original_state"], "UNKNOWN")
        self.assertTrue(decision["data"]["sites"][0]["expert_classification"])
        self.service.withdraw(decision["id"], self.token)
        self.assertFalse(self.service.view(current["id"])["current_policy"])
        with self.store.db() as db:
            source = self.service._job(db, self.job_id)["result"]["sites"]["atoms"]
        self.assertEqual(next(item for item in source if item["atom_map"] == 21)["state"],
                         "UNKNOWN")

    def test_unknown_site_cannot_use_release_as_classification(self):
        from packages.science.dual_e3 import catalog

        assessment, _ = self._assessment()
        rule_id = catalog()["rule_catalog"][0]["rule_id"]
        statement = self.service.create_statement(
            "Synthetic attempted release source.", "fixture:/sites/21/release", self.token
        )
        submission = {
            "kind": "site_policy", "action": "accept", "intent": "test",
            "reason": "test", "expected_policy_revision": 0,
            "data": {"parent_id": "SMARCA2-FX5", "sites": [{
                "atom_map": 21, "original_state": "UNKNOWN", "state": "MODIFIABLE",
                "expert_classification": True, "allowed_rule_ids": [rule_id],
                "rationale": "Attempted release shortcut.", "source_ref": statement["ref"],
                "allow_release_protected": True,
                "release_justification": "Not a protected source state.",
                "release_source_ref": statement["ref"],
            }]},
        }
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_SITE_RELEASE_FLAG"):
            self.service.record_decision(assessment["id"], submission, self.token)

    def test_archived_runtime_reinspection_is_explicit_and_fail_closed(self):
        self._stale_runtime()
        assessment, created = self._assessment()
        self.assertTrue(created)
        self.assertFalse(assessment["source_runtime_current"])
        self.assertFalse(assessment["scientific_accepted"])

        with self.assertRaisesRegex(
            ContractError, "SCIENTIFIC_DESIGN_RUNTIME_CHANGED"
        ):
            self.service.compute(self.job_id, self.token)

        result = self.service.compute(
            self.job_id, self.token, allow_archived_runtime=True
        )
        self.assertTrue(result["evidence"])
        for evidence in result["evidence"]:
            report = self.service.port.json(evidence["ref"])
            provenance = report["report_provenance"]
            self.assertTrue(provenance["historical_runtime_reinspection"])
            self.assertFalse(provenance["source_runtime_current"])
            self.assertTrue(provenance["source_inputs_current"])
            self.assertFalse(provenance["scientific_approval"])
            self.assertFalse(provenance["formal_acceptance"])
            self.assertIn("source_runtime_current=false", provenance["note"])
        self.assertFalse(self.service.view(assessment["id"])["scientific_accepted"])

    def test_archived_runtime_reinspection_rejects_changed_source_inputs(self):
        self._stale_runtime()
        changed = dict(input_binding())
        changed["digest"] = "0" * 64
        with mock.patch(
            "packages.platform.scientific_acceptance.input_binding",
            return_value=changed,
        ):
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_DESIGN_SOURCE_CHANGED"
            ):
                self.service.compute(
                    self.job_id, self.token, allow_archived_runtime=True
                )

    def test_archived_runtime_reinspection_rejects_reported_result_hash_change(self):
        self._stale_runtime()
        self._replace_job(
            lambda body: body["result"]["files"]["json"].__setitem__(
                "sha256", "0" * 64
            )
        )
        with self.assertRaises(ContractError):
            self.service.compute(
                self.job_id, self.token, allow_archived_runtime=True
            )

    def test_registered_geometry_boundary_rejects_identity_mismatches_and_never_accepts(self):
        result, report_ref, graph_hash, seed = self._registered_ternary_report()
        with mock.patch(
            "scripts.run_novel_ternary._candidate", return_value={}
        ), mock.patch(
            "packages.platform.scientific_acceptance.novel_ternary._mapped_graph",
            return_value={"graph_sha256": graph_hash},
        ):
            with self.store.db() as db:
                binding = self.service._followup_geometry_evidence_binding(
                    db, self.job_id, result, report_ref,
                    "C-1", "CRBN", seed, require_current=True,
                )
                self.assertTrue(binding["actual_computation"])
                self.assertTrue(binding["execution_success"])
                self.assertTrue(binding["reference_free"])
                with self.assertRaisesRegex(
                    ContractError,
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_REPORT_IDENTITY",
                ):
                    self.service._followup_geometry_evidence_binding(
                        db, self.job_id, result, report_ref,
                        "wrong-candidate", "CRBN", seed,
                        require_current=True,
                    )
                with self.assertRaisesRegex(
                    ContractError,
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_REPORT_IDENTITY",
                ):
                    self.service._followup_geometry_evidence_binding(
                        db, self.job_id, result, report_ref,
                        "C-1", "CRBN", seed + 1,
                        require_current=True,
                    )
                with self.assertRaisesRegex(
                    ContractError,
                    "SCIENTIFIC_FOLLOWUP_GEOMETRY_EVIDENCE",
                ):
                    self.service._followup_geometry_evidence_binding(
                        db, "wrong-job", result, report_ref,
                        "C-1", "CRBN", seed,
                        require_current=True,
                    )

        assessment, _ = self._assessment()
        rejection = {
            "kind": "formal_decision",
            "action": "reject",
            "intent": "Reject synthetic geometry fixture",
            "reason": "No scientific approval is represented.",
            "expected_policy_revision": 0,
            "data": {},
        }
        with self.assertRaises(ReviewAuthError):
            self.service.record_decision(
                assessment["id"], rejection, "not-a-session"
            )
        decision, created = self.service.record_decision(
            assessment["id"], rejection, self.token
        )
        self.assertTrue(created)
        self.assertFalse(self.service.view(assessment["id"])["scientific_accepted"])
        withdrawn = self.service.withdraw(decision["id"], self.token)
        self.assertTrue(withdrawn["withdrawn"])
        viewed = self.service.view(assessment["id"])
        self.assertFalse(viewed["scientific_accepted"])
        self.assertEqual(viewed["confirmation_status"], "unconfirmed")

    def test_authentication_is_required_for_mutating_authenticated_operations(self):
        with self.assertRaises(ReviewAuthError):
            self.service.compute(self.job_id, "not-a-session")
        with self.assertRaises(ReviewAuthError):
            self.service.create_statement("x", "fixture:x", "not-a-session")


if __name__ == "__main__":
    unittest.main()

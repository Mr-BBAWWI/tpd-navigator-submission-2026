"""Synthetic contract tests for replaying registered actual-computation records.

These fixtures exercise service behavior only and are not evidence of scientific
performance.
"""
from __future__ import annotations

import copy
import json
import unittest

from packages.contracts import ContractError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.science import novel_ternary

try:
    import tests.test_scientific_acceptance_service as fixture_module
except ImportError:
    import test_scientific_acceptance_service as fixture_module


class ScientificTernaryReplayTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture_module.ScientificAcceptanceServiceTests(
            methodName="runTest"
        )
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.service = self.fixture.service
        self.store = self.fixture.store
        self.job_id = self.fixture.job_id
        self.token = self.fixture.token
        self.report_provenance = {
            "kind": "authenticated_immutable_archived_result_reinspection",
            "operator_id": "synthetic-test-operator",
            "source_inputs_current": True,
            "source_runtime_current": False,
            "historical_runtime_reinspection": True,
            "note": (
                "reinspection of immutable archived result; "
                "source_runtime_current=false"
            ),
            "scientific_approval": False,
            "formal_acceptance": False,
        }
        self.original_report_provenance = copy.deepcopy(self.report_provenance)

    def _register_actual_sources(self):
        paths = []
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            self.service._authenticate(self.token, db)
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=True
            )
            from scripts.run_novel_ternary import _candidate as qualified_candidate
            graph = novel_ternary._mapped_graph(qualified_candidate(result, "C-1"))
            producer_policy = {"revision": 0, "digest": "producer-policy-digest"}
            plan = {
                "plan_digest": "fixture-plan-digest",
                "bindings": {
                    "job": {
                        "project": self.fixture.project,
                        "job_id": self.job_id,
                        "input_sha256": binding["input_sha256"],
                        "result_sha256": binding["result_sha256"],
                        "candidate_graph_sha256": graph["graph_sha256"],
                    },
                    "policy": producer_policy,
                },
                "candidate_graph": {
                    "candidate_id": "C-1",
                    "e3_type": "CRBN",
                    "graph_sha256": graph["graph_sha256"],
                },
            }
            plan_ref = self.service._write(db, paths, plan)
            aggregate_ref = self.service._write(
                db, paths, {"plan_digest": plan["plan_digest"]}
            )
            raw_ref = self.service._write(
                db, paths, b"immutable registered GPU fixture bytes",
                "application/octet-stream", "computed"
            )
            base_binding = {
                **binding,
                "evaluation_policy_digest": producer_policy["digest"],
                "evaluation_policy": producer_policy,
                "measurement_binding": {
                    "producer_policy": producer_policy,
                    "plan_digest": plan["plan_digest"],
                    "candidate_graph_sha256": graph["graph_sha256"],
                    "msa_mode": "single_sequence",
                },
            }

            created = []
            for seed, success, score in ((11, True, 0.4), (12, False, 0.1)):
                inspection = {
                    "status": "computed_hypothesis" if success else "failed",
                    "actual_computation": success,
                    "descriptive_metrics": {"geometry": {"score": score}},
                }
                value = {
                    "format": novel_ternary.RECEIPT_FORMAT,
                    "candidate_id": "C-1",
                    "e3_type": "CRBN",
                    "seed": seed,
                    "reference_free": True,
                    "actual_computation": success,
                    "execution_success": success,
                    "failure": None if success else {
                        "status": "failed",
                        "process_state": "failed",
                        "exit_code": 7,
                        "reason": "fixture seed failed",
                    },
                    "inspection": inspection,
                    "inspection_summary": {
                        "compatibility": "unknown",
                        "metrics": {"geometry.score": score},
                        "criterion_metric": None,
                        "criterion_observed": None,
                    },
                    "producer_policy_binding": producer_policy,
                    "measurement_protocol": {
                        "plan_digest": plan["plan_digest"],
                        "candidate_graph_sha256": graph["graph_sha256"],
                        "msa_mode": "single_sequence",
                        "boltz_input_sha256": "fixture-input",
                        "seed": seed,
                    },
                    "evaluation_policy": producer_policy,
                    "policy_incompatibilities": [
                        "current geometry metric is absent from measured output"
                    ],
                    "scientific_status": (
                        "computed_hypothesis" if success else "not_computed"
                    ),
                    "registered_artifacts": {"raw": raw_ref},
                    "plan_ref": plan_ref,
                    "aggregate_receipt_ref": aggregate_ref,
                }
                identity = f"C-1|CRBN|{seed}"
                created.append(self.service._register_evidence(
                    db, paths, self.job_id, "ternary", identity,
                    base_binding, value
                ))
        return created, raw_ref

    def test_replays_current_policy_without_changing_measurements(self):
        original, raw_ref = self._register_actual_sources()
        raw_before = self.service.port.read(raw_ref)
        policy = {
            "revision": 1,
            "decisions": {
                "ternary_geometry_criterion": {
                    "metric": "geometry.score",
                    "threshold": 0.5,
                    "comparison": "lte",
                }
            },
        }
        paths = []
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            self.service._authenticate(self.token, db)
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=True
            )
            replayed = self.service._replay_actual_ternary(
                db, paths, self.job_id, result, binding, policy,
                self.report_provenance,
            )

        self.assertEqual(len(replayed), 2)
        self.assertEqual(
            self.report_provenance, self.original_report_provenance
        )
        reports = {
            self.service.port.json(row["ref"])["seed"]:
            self.service.port.json(row["ref"])
            for row in replayed
        }
        self.assertEqual(
            reports[11]["inspection_summary"]["compatibility"], "compatible"
        )
        self.assertEqual(
            reports[12]["inspection_summary"]["compatibility"], "unknown"
        )
        self.assertFalse(reports[12]["execution_success"])
        self.assertEqual(reports[12]["failure"]["exit_code"], 7)
        self.assertEqual(
            reports[11]["original_inspection_summary"]["compatibility"], "unknown"
        )
        self.assertEqual(reports[11]["evaluation_policy"]["revision"], 1)
        for report in reports.values():
            self.assertEqual(
                report["report_provenance"], self.original_report_provenance
            )
        self.assertFalse(reports[11]["scientific_approved"])
        self.assertFalse(reports[11]["efficacy_approved"])
        self.assertEqual(self.service.port.read(raw_ref), raw_before)

        restarted = ScientificAcceptanceService(self.store, self.fixture.project)
        for report in reports.values():
            with self.store.db() as db:
                restarted._verify_refs(db, report)

        with self.store.db() as db:
            source_id = original[0]["id"]
            row = db.execute(
                "SELECT binding FROM scientific_evidence WHERE id=?", (source_id,)
            ).fetchone()
            tampered = json.loads(row["binding"])
            tampered["result_sha256"] = "0" * 64
            db.execute(
                "UPDATE scientific_evidence SET binding=? WHERE id=?",
                (json.dumps(tampered), source_id),
            )
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=True
            )
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_TERNARY_REPLAY_SOURCE_BINDING"
            ):
                self.service._replay_actual_ternary(
                    db, [], self.job_id, result, binding, copy.deepcopy(policy),
                    self.report_provenance,
                )

    def test_withdrawn_current_criterion_replays_as_unknown_and_preserves_raw(self):
        _, raw_ref = self._register_actual_sources()
        raw_before = self.service.port.read(raw_ref)
        policy = {"revision": 2, "decisions": {}}

        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            self.service._authenticate(self.token, db)
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=True
            )
            replayed = self.service._replay_actual_ternary(
                db, [], self.job_id, result, binding, policy,
                self.report_provenance,
            )

        reports = {
            self.service.port.json(row["ref"])["seed"]:
            self.service.port.json(row["ref"])
            for row in replayed
        }
        successful = reports[11]
        self.assertEqual(
            successful["report_provenance"], self.original_report_provenance
        )
        self.assertEqual(
            self.report_provenance, self.original_report_provenance
        )
        self.assertTrue(successful["execution_success"])
        self.assertTrue(successful["derived_interpretation"])
        self.assertEqual(
            successful["inspection_summary"]["compatibility"], "unknown"
        )
        self.assertIsNone(successful["inspection_summary"]["criterion_metric"])
        self.assertIsNone(successful["inspection_summary"]["criterion_observed"])
        self.assertIsNone(successful["inspection_summary"]["criterion_comparison"])
        self.assertIn(
            "current geometry criterion is absent or invalid",
            successful["policy_incompatibilities"],
        )
        self.assertEqual(successful["evaluation_policy"]["revision"], 2)
        self.assertEqual(successful["registered_artifacts"]["raw"], raw_ref)
        self.assertEqual(self.service.port.read(raw_ref), raw_before)

    def test_missing_current_metric_cannot_default_to_compatible(self):
        _, raw_ref = self._register_actual_sources()
        raw_before = self.service.port.read(raw_ref)
        policy = {
            "revision": 3,
            "decisions": {
                "ternary_geometry_criterion": {
                    "metric": "geometry.metric_not_measured",
                    "threshold": 1.0,
                    "comparison": "lte",
                }
            },
        }

        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            self.service._authenticate(self.token, db)
            _, result, binding, _ = self.service._verified_result(
                db, self.job_id, require_current=True
            )
            replayed = self.service._replay_actual_ternary(
                db, [], self.job_id, result, binding, policy,
                self.report_provenance,
            )

        reports = {
            self.service.port.json(row["ref"])["seed"]:
            self.service.port.json(row["ref"])
            for row in replayed
        }
        successful = reports[11]
        self.assertEqual(
            successful["report_provenance"], self.original_report_provenance
        )
        self.assertEqual(
            self.report_provenance, self.original_report_provenance
        )
        summary = successful["inspection_summary"]
        self.assertTrue(successful["execution_success"])
        self.assertEqual(summary["metrics"]["geometry.score"], 0.4)
        self.assertEqual(summary["criterion_metric"], "geometry.metric_not_measured")
        self.assertIsNone(summary["criterion_observed"])
        self.assertEqual(summary["compatibility"], "unknown")
        self.assertNotEqual(summary["compatibility"], "compatible")
        self.assertIn(
            "current geometry metric is absent from measured output",
            successful["policy_incompatibilities"],
        )
        self.assertEqual(successful["registered_artifacts"]["raw"], raw_ref)
        self.assertEqual(self.service.port.read(raw_ref), raw_before)


if __name__ == "__main__":
    unittest.main()

"""Regression tests for transaction-local scientific source verification."""
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
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService, identity


class ScientificTransactionReadTests(unittest.TestCase):
    project = "scientific-transaction-test"
    job_id = "job-uncommitted-scientific-fixture"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(Path(self.temp.name))

        # Initialize every required service table before the import transaction.
        self.workbench = WorkbenchService(self.store, self.project)
        self.service = ScientificAcceptanceService(self.store, self.project)
        self.created_paths = []

    def _artifact_ref(self, db, identifier, raw, media_type, *, stored_raw=None):
        reference = {
            "artifact_id": identifier,
            "version": 1,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "media_type": media_type,
            "schema_id": "urn:tpd-navigator:raw:1",
            "provenance": "computed",
        }
        path = self.store.root / "blobs" / identifier
        path.write_bytes(raw if stored_raw is None else stored_raw)
        self.created_paths.append(path)
        db.execute(
            "INSERT INTO artifacts VALUES(?,?,?)",
            (identifier, self.project, json.dumps(reference)),
        )
        return reference

    def _stage_uncommitted_fixture(self, db, *, corrupt_report=False):
        binding = input_binding()
        archived = {
            "format": "design-panel/20260930.4",
            "status": "synthetic_transaction_fixture",
            "input_binding": binding,
            "parent_scope": {"actual_design_parent_id": "SMARCA2-FX5"},
            "summary": {"valid_analogs": 0},
            "sites": {"atoms": []},
            "analogs": [],
            "protac_candidates": [],
            "calibration": {"CRBN": {"seed_receipts": []}},
            "files": {},
        }
        json_raw = encoded(archived)
        json_ref = self._artifact_ref(
            db,
            "a-uncommitted-result-json",
            json_raw,
            "application/json",
        )
        report_raw = b"# Uncommitted synthetic report\n"
        report_ref = self._artifact_ref(
            db,
            "a-uncommitted-result-report",
            report_raw,
            "text/markdown",
            stored_raw=(b"corrupt report bytes\n" if corrupt_report else None),
        )

        result = copy.deepcopy(archived)
        result["files"] = {"json": json_ref, "report": report_ref}
        job = {
            "id": self.job_id,
            "project_id": self.project,
            "run_id": "design:transaction-fixture",
            "result_id": "design:transaction-fixture",
            "operation": "design_panel",
            "title": "Uncommitted transaction fixture",
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
        db.execute(
            "INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)",
            (
                self.job_id,
                self.project,
                job["run_id"],
                job["result_id"],
                "transaction-fixture",
                "transaction-fixture-fingerprint",
                "completed",
                json.dumps(job),
            ),
        )
        return result, json_ref, report_ref

    def _remove_created_blobs(self):
        for path in self.created_paths:
            path.unlink(missing_ok=True)
        self.created_paths.clear()

    def _assert_nothing_persisted(self):
        with self.store.db() as db:
            jobs = db.execute(
                "SELECT COUNT(*) FROM workbench_jobs WHERE project=?",
                (self.project,),
            ).fetchone()[0]
            artifacts = db.execute(
                "SELECT COUNT(*) FROM artifacts WHERE project=?",
                (self.project,),
            ).fetchone()[0]
        self.assertEqual(jobs, 0)
        self.assertEqual(artifacts, 0)

    def test_verified_result_reads_real_uncommitted_refs_on_existing_connection(self):
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                result, json_ref, report_ref = self._stage_uncommitted_fixture(db)

                # A distinct real SQLite connection cannot see these rows yet.
                with self.store.db() as observer:
                    self.assertIsNone(observer.execute(
                        "SELECT 1 FROM workbench_jobs WHERE id=?",
                        (self.job_id,),
                    ).fetchone())
                    self.assertEqual(observer.execute(
                        "SELECT COUNT(*) FROM artifacts WHERE project=?",
                        (self.project,),
                    ).fetchone()[0], 0)

                with mock.patch.object(
                    self.service.port,
                    "read",
                    side_effect=AssertionError(
                        "transactional verification must not open the port reader"
                    ),
                ):
                    job, verified, source_binding, source = (
                        self.service._verified_result(
                            db, self.job_id, require_current=False
                        )
                    )

                self.assertEqual(job["id"], self.job_id)
                self.assertEqual(verified, result)
                self.assertEqual(source_binding, {
                    "project": self.project,
                    "job_id": self.job_id,
                    "input_sha256": job["input_digest"],
                    "result_sha256": json_ref["sha256"],
                    "result_binding_kind": "exact_archived_design_json_bytes",
                    "parent_id": "SMARCA2-FX5",
                })
                self.assertEqual(verified["files"]["report"], report_ref)
                self.assertTrue(source["source_inputs_current"])
                self.assertTrue(source["source_runtime_current"])
                db.rollback()
        finally:
            self._remove_created_blobs()

        self._assert_nothing_persisted()

    def test_corrupt_uncommitted_report_rolls_back_job_and_artifacts(self):
        try:
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_ARTIFACT_HASH"
            ):
                with self.store.db() as db:
                    db.execute("BEGIN IMMEDIATE")
                    self._stage_uncommitted_fixture(db, corrupt_report=True)
                    with mock.patch.object(
                        self.service.port,
                        "read",
                        side_effect=AssertionError(
                            "transactional verification must not open the port reader"
                        ),
                    ):
                        self.service._verified_result(
                            db, self.job_id, require_current=False
                        )
        finally:
            self._remove_created_blobs()

        self._assert_nothing_persisted()


if __name__ == "__main__":
    unittest.main()

"""Transactional scientific artifact-reference integration tests."""
from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from packages.contracts import ContractError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store


class ScientificTransactionReferenceTests(unittest.TestCase):
    project = "scientific-transaction-test"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(Path(self.temp.name))
        self.service = ScientificAcceptanceService(self.store, self.project)

    def test_nested_raw_references_registered_in_same_transaction_succeed(self):
        paths = []
        binding = {
            "project": self.project,
            "job_id": "job-h-import-fixture",
            "input_sha256": "1" * 64,
            "result_sha256": "2" * 64,
            "result_binding_kind": "exact_archived_design_json_bytes",
            "parent_id": "SMARCA2-FX5",
        }
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            source_ref = self.service._write(
                db, paths, b"ATOM fixture\n", "chemical/x-pdb", "source"
            )
            diagnostic_ref = self.service._write(
                db, paths, b"diagnostic fixture\n", "text/plain", "computed"
            )
            value = {
                "status": "computed_diagnostic_pending_human_review",
                "protein_hydrogens": [],
                "registered_artifacts": {
                    "source/source.pdb": source_ref,
                    "pdb2pqr/stdout.log": diagnostic_ref,
                },
                "submitted_evidence_ref": source_ref,
                "state_flags": {
                    "computed_diagnostic": True,
                    "formal_scientific_approval": False,
                },
            }
            registered = self.service._register_evidence(
                db, paths, "job-h-import-fixture", "protein_h",
                "SMARCA2-FX5", binding, value,
            )

        self.assertTrue(registered["created"])
        archived = self.service.port.json(registered["ref"])
        self.assertEqual(archived["registered_artifacts"]["source/source.pdb"], source_ref)
        self.assertEqual(self.service.port.read(source_ref), b"ATOM fixture\n")
        with self.store.db() as db:
            row = db.execute(
                "SELECT kind,ref FROM scientific_evidence WHERE id=? AND project=?",
                (registered["id"], self.project),
            ).fetchone()
        self.assertIsNotNone(row)
        self.assertEqual(row["kind"], "protein_h")
        self.assertEqual(json.loads(row["ref"]), registered["ref"])

    def test_tampered_reference_metadata_and_blob_are_rejected(self):
        ref = self.service.port.put_raw(b"immutable fixture", "text/plain", "source")
        changed = copy.deepcopy(ref)
        changed["sha256"] = "f" * 64
        with self.store.db() as db:
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_ARTIFACT_METADATA"
            ):
                self.service._verify_refs(db, {"nested": [changed]})

        blob = self.store.root / "blobs" / ref["artifact_id"]
        blob.write_bytes(b"tampered fixture")
        with self.store.db() as db:
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_ARTIFACT_HASH"
            ):
                self.service._verify_refs(db, ref)

    def test_cross_project_and_registered_path_escape_are_rejected(self):
        ref = self.service.port.put_raw(b"project-bound", "text/plain", "source")
        other = ScientificAcceptanceService(self.store, "other-project")
        with self.store.db() as db:
            with self.assertRaisesRegex(
                ContractError, "SCIENTIFIC_ARTIFACT_PROJECT"
            ):
                other._verify_refs(db, ref)

        escaping = {
            "artifact_id": "../outside",
            "version": 1,
            "sha256": "0" * 64,
            "media_type": "application/octet-stream",
            "schema_id": "urn:tpd-navigator:raw:1",
            "provenance": "source",
        }
        with self.store.db() as db:
            db.execute(
                "INSERT INTO artifacts(id,project,metadata) VALUES(?,?,?)",
                (escaping["artifact_id"], self.project, json.dumps(escaping)),
            )
        with self.store.db() as db:
            with self.assertRaisesRegex(
                ContractError, "(?:SCHEMA:ArtifactRef|SCIENTIFIC_ARTIFACT_PATH)"
            ):
                self.service._verify_refs(db, escaping)

    def test_failed_transaction_removes_only_its_new_blob(self):
        key = self.service.auth.register("reviewer", "Reviewer")
        token = self.service.auth.login("reviewer", key)
        existing = self.service.port.put_raw(b"existing", "text/plain", "source")
        existing_path = self.store.root / "blobs" / existing["artifact_id"]
        before = {path.name for path in (self.store.root / "blobs").iterdir()}
        original = self.service.saved._write

        def write_then_fail(db, paths, raw, media="application/json",
                            provenance="computed"):
            original(db, paths, raw, media, provenance)
            raise RuntimeError("forced transaction failure")

        with mock.patch.object(self.service.saved, "_write", side_effect=write_then_fail):
            with self.assertRaisesRegex(RuntimeError, "forced transaction failure"):
                self.service.create_statement("rollback fixture", "fixture:rollback", token)

        after = {path.name for path in (self.store.root / "blobs").iterdir()}
        self.assertEqual(after, before)
        self.assertTrue(existing_path.is_file())
        self.assertEqual(existing_path.read_bytes(), b"existing")
        with self.store.db() as db:
            count = db.execute(
                "SELECT COUNT(*) FROM scientific_statements WHERE project=?",
                (self.project,),
            ).fetchone()[0]
        self.assertEqual(count, 0)


if __name__ == "__main__":
    unittest.main()

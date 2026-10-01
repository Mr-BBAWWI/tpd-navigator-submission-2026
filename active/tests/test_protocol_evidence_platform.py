"""Platform registration tests use only synthetic packs and a mocked verifier."""
from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from packages.contracts import ContractError
from packages.platform.protocol_evidence import ProtocolEvidenceService
from packages.platform.store import Store


class ProtocolEvidencePlatformTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.store = Store(self.temp.name)
        self.service = ProtocolEvidenceService(self.store, "alpha")
        with self.store.db() as db:
            db.execute("CREATE TABLE IF NOT EXISTS workbench_jobs(id TEXT, project TEXT)")
            db.execute("INSERT INTO workbench_jobs(id,project) VALUES(?,?)", ("job-1", "alpha"))
        self.pack = Path(self.temp.name) / "pack"
        self.pack.mkdir()
        self.manifest_raw = b'{"files":[]}'
        (self.pack / "output-manifest.json").write_bytes(self.manifest_raw)
        self.manifest_sha = hashlib.sha256(self.manifest_raw).hexdigest()
        self.binding = {
            "project": "alpha", "job_id": "job-1", "input_sha256": "1" * 64,
            "result_sha256": "2" * 64,
            "result_binding_kind": "exact_archived_design_json_bytes", "parent_id": "parent-1"
        }
        self.result = {"protac_candidates": []}
        self.bundle = {
            "pack_manifest_sha256": self.manifest_sha,
            "compact_sha256": "3" * 64,
            "files": [{"path": "output-manifest.json", "sha256": self.manifest_sha,
                       "bytes": len(self.manifest_raw)}],
            "protocols": [
                {"protocol_id": "protocol-a", "kind": "core", "summary": {"row_count": 40},
                 "measurement_verification": {"verified": True},
                 "scientific_approved": False, "gates_affected": False},
                {"protocol_id": "protocol-b", "kind": "ternary", "summary": {"seed_count": 10},
                 "measurement_verification": {"verified": True},
                 "scientific_approved": False, "gates_affected": False},
            ],
        }

    def tearDown(self):
        self.temp.cleanup()

    def source(self, db, job_id, assessment_id):
        if job_id != "job-1":
            raise KeyError(job_id)
        link = None if assessment_id is None else {
            "assessment_id": assessment_id,
            "revision": 1,
            "original_status_summary": "computed",
        }
        return {}, self.result, dict(self.binding), {}, link

    def snapshots(self):
        with self.store.db() as db:
            result = {}
            for table in ("scientific_assessments", "scientific_decisions",
                          "scientific_policy_events", "scientific_evidence"):
                exists = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                                    (table,)).fetchone()
                result[table] = [] if exists is None else [tuple(row) for row in db.execute(
                    "SELECT * FROM " + table + " ORDER BY rowid").fetchall()]
            return result

    def register(self, bundle=None, pack=None, manifest_sha=None, assessment_id=None):
        bundle = self.bundle if bundle is None else bundle
        pack = self.pack if pack is None else pack
        manifest_sha = self.manifest_sha if manifest_sha is None else manifest_sha
        with patch.object(self.service, "_verified_source", side_effect=self.source), \
             patch("packages.platform.protocol_evidence.verify_pack", return_value=bundle):
            return self.service.register_pack(
                "job-1", pack, expected_manifest_sha256=manifest_sha,
                expected_assessment_id=assessment_id)

    def independent_pack(self, protocols):
        pack = Path(self.temp.name) / ("pack-" + str(len(list(Path(self.temp.name).glob("pack-*")))))
        pack.mkdir()
        raw = json.dumps({"files": [], "nonce": pack.name}, sort_keys=True).encode()
        (pack / "output-manifest.json").write_bytes(raw)
        digest = hashlib.sha256(raw).hexdigest()
        bundle = {
            "pack_manifest_sha256": digest,
            "compact_sha256": "4" * 64,
            "files": [{"path": "output-manifest.json", "sha256": digest, "bytes": len(raw)}],
            "protocols": copy.deepcopy(protocols),
        }
        return pack, digest, bundle

    def test_registers_both_and_is_idempotent_without_approval_mutation(self):
        before = self.snapshots()
        first = self.register()
        second = self.register()
        self.assertTrue(first["created"])
        self.assertFalse(second["created"])
        self.assertEqual(2, len(first["records"]))
        self.assertEqual(before, self.snapshots())
        for record in first["records"]:
            self.assertFalse(record["scientific_approved"])
            self.assertFalse(record["formal_acceptance"])
            self.assertFalse(record["gates_affected"])
            self.assertEqual("computed", record["archive_ref"]["provenance"])
            self.assertEqual(1, record["verification_engine"]["version"])
        self.assertEqual(2, len(self.service.list("job-1")))

    def test_existing_job_with_no_records_and_absent_job(self):
        self.assertEqual([], self.service.list("job-1"))
        with self.assertRaises(KeyError):
            self.service.list("missing-job")

    def test_missing_workbench_table_is_not_found(self):
        with self.store.db() as db:
            db.execute("DROP TABLE workbench_jobs")
        with self.assertRaises(KeyError):
            self.service.list("job-1")

    def test_cross_project_not_found(self):
        record = self.register()["records"][0]
        other = ProtocolEvidenceService(self.store, "beta")
        with self.assertRaises(KeyError):
            other.list("job-1")
        with self.assertRaises(KeyError):
            other.view(record["id"])

    def test_independent_pack_reuses_identical_protocol_and_adds_new_protocol(self):
        first = self.register()
        protocol_a = copy.deepcopy(self.bundle["protocols"][0])
        protocol_c = copy.deepcopy(self.bundle["protocols"][1])
        protocol_c["protocol_id"] = "protocol-c"
        protocol_c["summary"] = {"seed_count": 20}
        pack, digest, bundle = self.independent_pack([protocol_a, protocol_c])

        second = self.register(bundle, pack, digest)

        self.assertTrue(second["created"])
        self.assertEqual(["protocol-a", "protocol-c"],
                         [record["protocol_id"] for record in second["records"]])
        self.assertEqual(first["records"][0]["id"], second["records"][0]["id"])
        self.assertEqual(first["records"][0]["archive_ref"],
                         second["records"][0]["archive_ref"])
        self.assertNotEqual(second["records"][0]["archive_ref"],
                            second["records"][1]["archive_ref"])
        self.assertEqual(3, len(self.service.list("job-1")))

    def test_changed_protocol_summary_conflicts_without_overwrite(self):
        first = self.register()
        changed = copy.deepcopy(self.bundle)
        changed["protocols"][0]["summary"] = {"row_count": 41}
        with self.assertRaises(ContractError):
            self.register(changed)
        stored = self.service.view(first["records"][0]["id"])
        self.assertEqual({"row_count": 40}, stored["summary"])
        self.assertEqual(2, len(self.service.list("job-1")))

    def test_changed_assessment_link_conflicts(self):
        self.register()
        with self.assertRaises(ContractError):
            self.register(assessment_id="assessment-2")
        self.assertEqual(2, len(self.service.list("job-1")))

    def test_failure_rolls_back_registry_and_new_blobs(self):
        original = self.service._write_artifact
        calls = 0
        def failing(db, paths, raw, media, provenance):
            nonlocal calls
            calls += 1
            ref = original(db, paths, raw, media, provenance)
            if calls == 2:
                raise RuntimeError("synthetic failure")
            return ref
        with patch.object(self.service, "_verified_source", side_effect=self.source), \
             patch("packages.platform.protocol_evidence.verify_pack", return_value=self.bundle), \
             patch.object(self.service, "_write_artifact", side_effect=failing):
            with self.assertRaises(RuntimeError):
                self.service.register_pack("job-1", self.pack,
                    expected_manifest_sha256=self.manifest_sha)
        with self.store.db() as db:
            self.assertEqual(0, db.execute("SELECT COUNT(*) FROM protocol_diagnostics").fetchone()[0])
            self.assertEqual(0, db.execute("SELECT COUNT(*) FROM artifacts").fetchone()[0])
        self.assertEqual([], list((Path(self.temp.name) / "blobs").iterdir()))

    def test_corrupt_record_blob_fails_closed(self):
        record = self.register()["records"][0]
        with self.store.db() as db:
            row = db.execute("SELECT record_ref FROM protocol_diagnostics WHERE id=?",
                             (record["id"],)).fetchone()
        ref = json.loads(row["record_ref"])
        (Path(self.temp.name) / "blobs" / ref["artifact_id"]).write_bytes(b"{}")
        with self.assertRaises(ContractError):
            self.service.view(record["id"])
        with self.assertRaises(ContractError):
            self.service.list("job-1")


if __name__ == "__main__":
    unittest.main()

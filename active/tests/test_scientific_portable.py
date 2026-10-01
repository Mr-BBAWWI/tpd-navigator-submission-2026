from __future__ import annotations

import hashlib
import importlib.util
import json
import sqlite3
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


SOURCE = Path(__file__).resolve().parents[1] / "reviewer_delivery" / "scientific_portable.py"
BUILDER_SOURCE = Path(__file__).resolve().parents[1] / "reviewer_delivery" / "build_scientific_package.py"


def load_portable(root):
    program = root / "program"
    program.mkdir(parents=True)
    design = types.ModuleType("design_portable")
    design._immutable_files = lambda: {}
    design.dump = lambda path, value: path.write_text(json.dumps(value), encoding="utf-8")
    design.check = lambda **kwargs: None
    design.prepare = lambda **kwargs: None
    design.start = lambda **kwargs: None
    design.stop = lambda: None
    design.ping = lambda port: {}
    old = sys.modules.get("design_portable")
    sys.modules["design_portable"] = design
    try:
        spec = importlib.util.spec_from_file_location("scientific_portable_test", SOURCE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if old is None:
            sys.modules.pop("design_portable", None)
        else:
            sys.modules["design_portable"] = old
    module.ROOT = root
    module.APP = program
    module.STATE = root / "design-data"
    module.SEED = root / "evidence" / "scientific-seed"
    module.design = design
    return module


class ScientificPortableTests(unittest.TestCase):
    def test_review_url_uses_scientific_route(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            manifest = {"release": "fixture"}
            (module.ROOT / "PACKAGE-MANIFEST.json").write_text(
                json.dumps(manifest), encoding="utf-8")
            with mock.patch.object(
                    module, "_verified_v5_release", return_value="v7") as verified, \
                    mock.patch.object(
                        module, "_seed_data",
                        return_value={"job": {"id": "job/a b"}}) as seeded:
                self.assertEqual(module._review_url(8790), "http://127.0.0.1:8790/scientific-review/job%2Fa%20b")
            verified.assert_called_once_with(manifest)
            seeded.assert_called_once_with("v7")

    def test_recursive_json_blob_closure_includes_nested_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            module.SEED.mkdir(parents=True)
            blobs = module.SEED / "blobs"
            blobs.mkdir()
            leaf_raw = b"actual-3d-pose"
            leaf = {"artifact_id": "a-leaf", "version": 1, "sha256": hashlib.sha256(leaf_raw).hexdigest(),
                    "media_type": "chemical/x-sdf", "schema_id": "raw", "provenance": "computed"}
            report_raw = json.dumps({"diagnostic": {"pose_ref": leaf}}).encode()
            report = {"artifact_id": "a-report", "version": 1, "sha256": hashlib.sha256(report_raw).hexdigest(),
                      "media_type": "application/json", "schema_id": "raw", "provenance": "computed"}
            (blobs / "a-report").write_bytes(report_raw)
            (blobs / "a-leaf").write_bytes(leaf_raw)
            data = {"job": {"report": report}, "calls": [], "assessments": [], "evidence": [],
                    "artifacts": [report, leaf]}
            refs, raw = module._load_blob_closure(data, lambda kind, value: None)
            self.assertEqual(set(refs), {"a-report", "a-leaf"})
            self.assertEqual(raw["a-leaf"], leaf_raw)

    def test_corrupt_nested_blob_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            blobs = module.SEED / "blobs"
            blobs.mkdir(parents=True)
            leaf = {"artifact_id": "a-leaf", "version": 1, "sha256": hashlib.sha256(b"correct").hexdigest(),
                    "media_type": "application/octet-stream", "schema_id": "raw", "provenance": "computed"}
            report_raw = json.dumps({"nested": leaf}).encode()
            report = {"artifact_id": "a-report", "version": 1, "sha256": hashlib.sha256(report_raw).hexdigest(),
                      "media_type": "application/json", "schema_id": "raw", "provenance": "computed"}
            (blobs / "a-report").write_bytes(report_raw)
            (blobs / "a-leaf").write_bytes(b"tampered")
            data = {"job": {"report": report}, "calls": [], "assessments": [], "evidence": [],
                    "artifacts": [report, leaf]}
            with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_BLOB_HASH_INVALID"):
                module._load_blob_closure(data, lambda kind, value: None)

    def test_repeat_verification_preserves_human_and_superseded_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            module.STATE.mkdir()
            (module.STATE / "blobs").mkdir()
            db_path = module.STATE / "index.sqlite3"
            connection = sqlite3.connect(db_path)
            connection.row_factory = sqlite3.Row
            connection.executescript("""
                CREATE TABLE artifacts(id TEXT PRIMARY KEY,project TEXT,metadata TEXT);
                CREATE TABLE workbench_jobs(id TEXT PRIMARY KEY,project TEXT,run_id TEXT,result_id TEXT,request_key TEXT,fingerprint TEXT,state TEXT,body TEXT);
                CREATE TABLE workbench_calls(job_id TEXT,ordinal INTEGER,body TEXT);
                CREATE TABLE scientific_decisions(id TEXT PRIMARY KEY,project TEXT,job_id TEXT,body TEXT);
                CREATE TABLE scientific_evidence(id TEXT PRIMARY KEY,project TEXT,job_id TEXT,superseded INTEGER);
            """)
            job = {"id": "job-1", "run_id": "run", "result_id": "design:SMARCA2", "state": "completed"}
            binding = {"request_key": "request", "fingerprint": "fingerprint"}
            ref = {"artifact_id": "a-one", "version": 1, "sha256": hashlib.sha256(b"blob").hexdigest(),
                   "media_type": "application/octet-stream", "schema_id": "raw", "provenance": "computed"}
            connection.execute("INSERT INTO artifacts VALUES(?,?,?)", ("a-one", "project", json.dumps(ref)))
            connection.execute("INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)",
                               ("job-1", "project", "run", "design:SMARCA2", "request", "fingerprint", "completed", json.dumps(job)))
            connection.execute("INSERT INTO scientific_decisions VALUES(?,?,?,?)", ("decision-local", "project", "job-1", "{}"))
            connection.execute("INSERT INTO scientific_evidence VALUES(?,?,?,?)", ("evidence-local", "project", "job-1", 1))
            connection.commit()
            (module.STATE / "blobs" / "a-one").write_bytes(b"blob")
            data = {"job": job, "calls": [], "source_binding": binding}
            module._verify_imported_base(connection, data, "project", {"a-one": ref}, {"a-one": b"blob"})
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM scientific_decisions").fetchone()[0], 1)
            self.assertEqual(connection.execute("SELECT superseded FROM scientific_evidence").fetchone()[0], 1)
            connection.close()

    def test_lifecycle_adapter_accepts_both_base_keywords_without_recursion(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            original_check = module.design.check
            original_prepare = module.design.prepare
            with mock.patch.object(module, "check_manifest", return_value=True) as checked:
                with module._design_lifecycle_adapter({"seed": True}):
                    self.assertTrue(module.design.check(require_seed=True))
                    self.assertEqual(module.design.prepare(required=True), {"seed": True})
                    self.assertEqual(module.design.prepare(require_seed=True), {"seed": True})
            self.assertIs(module.design.check, original_check)
            self.assertIs(module.design.prepare, original_prepare)
            checked.assert_called()


if __name__ == "__main__":
    unittest.main()

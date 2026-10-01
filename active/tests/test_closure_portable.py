"""Synthetic v7 portable-registry tests; no fake readiness receipt is produced."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
DELIVERY = ROOT / "reviewer_delivery"
if str(DELIVERY) not in sys.path:
    sys.path.insert(0, str(DELIVERY))

import scientific_portable as portable


def ref(identifier, raw, media="application/json"):
    return {
        "artifact_id": identifier,
        "version": 1,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "media_type": media,
        "schema_id": "urn:tpd-navigator:raw:1",
        "provenance": "computed",
    }


class ClosurePortableTests(unittest.TestCase):
    def test_strict_json_rejects_duplicate_and_nonfinite(self):
        for raw in (b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":Infinity}'):
            with self.subTest(raw=raw), self.assertRaises(RuntimeError):
                portable._strict_loads(raw, "SYNTHETIC_INVALID")

    def test_only_v7_accepts_protocol_seed_field(self):
        data = {
            "format": portable.SEED_FORMAT,
            "job": {}, "calls": [], "assessments": [], "evidence": [],
            "artifacts": [], "source_binding": {},
            "row_hashes": {}, "human_decisions": [],
            "scientific_approval": False, "protocol_diagnostics": [],
        }
        for name in ("job", "calls", "assessments", "evidence", "protocol_diagnostics"):
            data["row_hashes"][name] = portable.digest(data[name])
        data["seed_sha256"] = portable.digest(data)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed = root / "seed.json"
            seed.write_text(json.dumps(data), encoding="utf-8")
            with mock.patch.object(portable, "SEED", root):
                with self.assertRaisesRegex(RuntimeError, "V7_REQUIRED"):
                    portable._seed_data("v6")
                self.assertEqual(portable._seed_data("v7")["protocol_diagnostics"], [])

    def test_review_url_verifies_actual_v7_metadata_and_seed_without_import_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed_root = root / "evidence" / "scientific-seed"
            seed_root.mkdir(parents=True)
            state = root / "design-data"
            release = {
                "release_version": "v7",
                "artifact": "final-closure-package",
                "archive_name": "TPD_closure_v7.zip",
                "package_schema": portable.MANIFEST_FORMAT,
                "seed_schema": portable.SEED_FORMAT,
                "stage_name": "TPD",
                "source_followup": portable.V5_SOURCE_FOLLOWUP,
                "source_followup_rules_installed": True,
                "scientific_accepted": False,
                "scientific_approval": False,
                "automatic_approval": False,
                "all_pass_claimed": False,
                "scientific_status": "pending",
                "portable_runtime_included": True,
                "source_included": True,
                "verified_source_bound_seed_required_for_archive": True,
                "runtime_policy": {
                    "archived_job": "read-only historical evidence under current implementation policy",
                    "original_runtime_distinct": True,
                },
            }
            release_raw = json.dumps(release).encode("utf-8")
            (root / portable.V7_RELEASE_METADATA_NAME).write_bytes(release_raw)
            manifest = {
                "format": portable.MANIFEST_FORMAT,
                "files": {
                    portable.V7_RELEASE_METADATA_NAME: hashlib.sha256(release_raw).hexdigest(),
                },
            }
            (root / "PACKAGE-MANIFEST.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            seed = {
                "format": portable.SEED_FORMAT,
                "job": {"id": "job/v7"},
                "calls": [],
                "assessments": [],
                "evidence": [],
                "protocol_diagnostics": [],
                "artifacts": [],
                "source_binding": {},
                "row_hashes": {},
                "human_decisions": [],
                "scientific_approval": False,
            }
            for name in (
                "job", "calls", "assessments", "evidence", "protocol_diagnostics"
            ):
                seed["row_hashes"][name] = portable.digest(seed[name])
            seed["seed_sha256"] = portable.digest(seed)
            (seed_root / "seed.json").write_text(json.dumps(seed), encoding="utf-8")
            with mock.patch.object(portable, "ROOT", root), \
                    mock.patch.object(portable, "SEED", seed_root), \
                    mock.patch.object(portable, "STATE", state):
                self.assertEqual(
                    portable._review_url(8123),
                    "http://127.0.0.1:8123/scientific-review/job%2Fv7",
                )
            self.assertFalse(state.exists())

    def test_metadata_ambiguity_and_v7_version_tamper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = {"format": portable.MANIFEST_FORMAT, "files": {}}
            for name in (portable.V6_RELEASE_METADATA_NAME, portable.V7_RELEASE_METADATA_NAME):
                (root / name).write_text("{}", encoding="utf-8")
                manifest["files"][name] = hashlib.sha256(b"{}").hexdigest()
            with mock.patch.object(portable, "ROOT", root):
                with self.assertRaisesRegex(RuntimeError, "AMBIGUOUS"):
                    portable._verified_v5_release(manifest)

            (root / portable.V6_RELEASE_METADATA_NAME).unlink()
            manifest["files"].pop(portable.V6_RELEASE_METADATA_NAME)
            value = {
                "release_version": "v6",
                "artifact": "final-closure-package",
                "archive_name": "TPD_closure_v7.zip",
                "package_schema": portable.MANIFEST_FORMAT,
                "seed_schema": portable.SEED_FORMAT,
                "stage_name": "TPD",
                "source_followup": portable.V5_SOURCE_FOLLOWUP,
                "source_followup_rules_installed": True,
                "scientific_accepted": False,
                "scientific_approval": False,
                "automatic_approval": False,
                "all_pass_claimed": False,
                "scientific_status": "pending",
                "portable_runtime_included": True,
                "source_included": True,
                "verified_source_bound_seed_required_for_archive": True,
                "runtime_policy": {
                    "archived_job": "read-only historical evidence under current implementation policy",
                    "original_runtime_distinct": True,
                },
            }
            raw = json.dumps(value).encode()
            (root / portable.V7_RELEASE_METADATA_NAME).write_bytes(raw)
            manifest["files"][portable.V7_RELEASE_METADATA_NAME] = hashlib.sha256(raw).hexdigest()
            with mock.patch.object(portable, "ROOT", root):
                with self.assertRaisesRegex(RuntimeError, "V7_RELEASE_METADATA_INVALID"):
                    portable._verified_v5_release(manifest)

    def test_protocol_plan_rejects_duplicate_id_and_bad_ref(self):
        binding = {"project": "p", "job_id": "j", "input_sha256": "1" * 64,
                   "result_sha256": "2" * 64}
        archive_raw = b"archive"
        archive = ref("a-archive", archive_raw, "application/zip")
        record = {
            "id": "pdiag-stable", "project_id": "p", "job_id": "j",
            "protocol_id": "proto", "kind": "synthetic", "created_at": "t",
            "source_binding": binding,
            "assessment_link": {
                "assessment_id": "assessment-1",
                "revision": 1,
                "original_status_summary": None,
            },
            "archive_ref": archive,
            "pack_manifest_sha256": "3" * 64,
            "verification_engine": {"version": 1, "source_hashes": {}},
        }
        record_raw = portable.canonical(record)
        record_ref = ref("a-record", record_raw)
        row = {
            "id": "pdiag-stable", "project": "p", "job_id": "j",
            "protocol_id": "proto", "kind": "synthetic",
            "pack_manifest_sha256": "3" * 64,
            "binding": json.dumps(binding), "record_ref": json.dumps(record_ref),
            "archive_ref": json.dumps(archive), "created_at": "t",
        }
        item = {"row": row, "record": record,
                "record_ref": record_ref, "archive_ref": archive}
        data = {"job": {"id": "j", "project_id": "p"},
                "protocol_diagnostics": [item, copy.deepcopy(item)]}
        service = mock.Mock()
        service.verification_engine = record["verification_engine"]
        service._validate_record.return_value = None
        refs = {"a-record": record_ref, "a-archive": archive}
        blobs = {"a-record": record_raw, "a-archive": archive_raw}
        with self.assertRaisesRegex(RuntimeError, "ID_INVALID"):
            portable._validate_protocol_registry(data, "p", refs, blobs, service)
        data["protocol_diagnostics"] = [item]
        bad = copy.deepcopy(refs)
        bad["a-record"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(RuntimeError, "REF|RECORD"):
            portable._validate_protocol_registry(data, "p", bad, blobs, service)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
DELIVERY = ROOT / "reviewer_delivery"
if str(DELIVERY) not in sys.path:
    sys.path.insert(0, str(DELIVERY))

import build_sprint_package as sprint
from tests.test_scientific_portable import load_portable


def canonical(value):
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")


def seed_value(job_id="job-one", assessments=None):
    value = {
        "format": "scientific-seed/3",
        "job": {"id": job_id, "project_id": "project-one"},
        "assessments": assessments if assessments is not None else [
            {"id": "assessment-2", "revision": 2, "job_id": job_id, "project": "project-one"},
            {"id": "assessment-7", "revision": 7, "job_id": job_id, "project": "project-one"},
        ],
    }
    value["seed_sha256"] = hashlib.sha256(canonical(value)).hexdigest()
    return value


def write_builder_fixture(stage, seed=None):
    seed = seed or seed_value()
    seed_path = stage / "evidence/scientific-seed/seed.json"
    seed_path.parent.mkdir(parents=True)
    seed_path.write_text(json.dumps(seed), encoding="utf-8")
    pack = stage / sprint.PROTOCOL_PACK_RELATIVE_PATH
    pack.mkdir(parents=True)
    output = pack / "output-manifest.json"
    output.write_bytes(b'{"protocols":[]}')
    return seed, pack, output


class ProtocolDescriptorBuilderTests(unittest.TestCase):
    def test_descriptor_pins_highest_seed_revision_and_pack_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "TPD"
            stage.mkdir()
            seed, pack, output = write_builder_fixture(stage)
            before = {path.relative_to(pack): path.read_bytes() for path in pack.rglob("*") if path.is_file()}
            value = sprint._write_protocol_import(stage, "job-one")
            stored = json.loads((stage / sprint.PROTOCOL_IMPORT_DESCRIPTOR).read_text(encoding="utf-8"))
            self.assertEqual(value, stored)
            self.assertEqual(stored["format"], "sprint-protocol-import/1")
            self.assertEqual(stored["expected_assessment_id"], "assessment-7")
            self.assertEqual(stored["seed_sha256"], seed["seed_sha256"])
            self.assertEqual(stored["expected_manifest_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
            self.assertFalse(stored["scientific_approval"])
            self.assertEqual(before, {path.relative_to(pack): path.read_bytes() for path in pack.rglob("*") if path.is_file()})
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_DESCRIPTOR_EXISTS"):
                sprint._write_protocol_import(stage, "job-one")

    def test_absent_pack_does_not_create_descriptor(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "TPD"
            stage.mkdir()
            self.assertIsNone(sprint._write_protocol_import(stage, "job-one"))
            self.assertFalse((stage / sprint.PROTOCOL_IMPORT_DESCRIPTOR).exists())

    def test_wrong_job_and_duplicate_revision_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "wrong" / "TPD"
            stage.mkdir(parents=True)
            write_builder_fixture(stage)
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_JOB_INVALID"):
                sprint._write_protocol_import(stage, "another-job")
        duplicate = [
            {"id": "a", "revision": 3, "job_id": "job-one", "project": "project-one"},
            {"id": "b", "revision": 3, "job_id": "job-one", "project": "project-one"},
        ]
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "duplicate" / "TPD"
            stage.mkdir(parents=True)
            write_builder_fixture(stage, seed_value(assessments=duplicate))
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_ASSESSMENT_INVALID"):
                sprint._write_protocol_import(stage, "job-one")


class ProtocolPortableImportTests(unittest.TestCase):
    def _setup(self, root):
        module = load_portable(root)
        data = seed_value()
        pack = root / module.PROTOCOL_PACK_RELATIVE_PATH
        pack.mkdir(parents=True)
        output = pack / "output-manifest.json"
        output.write_bytes(b'{"format":"protocol-pack/1","protocols":[]}')
        payload = sprint._release_payload()
        release = root / module.V6_RELEASE_METADATA_NAME
        release.write_text(json.dumps(payload), encoding="utf-8")
        descriptor = {
            "format": module.PROTOCOL_IMPORT_FORMAT,
            "job_id": "job-one",
            "expected_assessment_id": "assessment-7",
            "pack_relative_path": module.PROTOCOL_PACK_RELATIVE_PATH,
            "expected_manifest_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            "seed_sha256": data["seed_sha256"],
            "scientific_approval": False,
        }
        descriptor_path = root / module.PROTOCOL_IMPORT_DESCRIPTOR
        descriptor_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor_path.write_text(json.dumps(descriptor), encoding="utf-8")
        files = {}
        for path in (release, descriptor_path, output):
            files[path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest = {
            "format": module.MANIFEST_FORMAT,
            "scientific_seed_job_id": "job-one",
            "files": files,
        }
        (root / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
        return module, data, descriptor, descriptor_path, output

    def _manifest(self, root):
        return json.loads((root / "PACKAGE-MANIFEST.json").read_text(encoding="utf-8"))

    def _save_manifest(self, root, manifest):
        (root / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")

    def _rewrite_descriptor(self, root, module, descriptor, update_hash=True):
        path = root / module.PROTOCOL_IMPORT_DESCRIPTOR
        path.write_text(json.dumps(descriptor), encoding="utf-8")
        if update_hash:
            manifest = self._manifest(root)
            manifest["files"][module.PROTOCOL_IMPORT_DESCRIPTOR] = hashlib.sha256(path.read_bytes()).hexdigest()
            self._save_manifest(root, manifest)

    def _service_module(self, calls, result=None):
        fake = types.ModuleType("packages.platform.protocol_evidence")
        returned = result if result is not None else {"created": True, "records": [{"id": "record"}]}

        class ProtocolEvidenceService:
            def __init__(self, store, project):
                calls.append(("init", store, project))

            def register_pack(self, job_id, pack_root, **kwargs):
                calls.append(("register", job_id, pack_root, kwargs))
                return returned

        fake.ProtocolEvidenceService = ProtocolEvidenceService
        return fake

    def test_valid_v6_descriptor_calls_only_new_service_with_exact_binding(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, descriptor, _, _ = self._setup(root)
            calls = []
            store = object()
            with mock.patch.dict(sys.modules, {
                "packages.platform.protocol_evidence": self._service_module(calls)
            }):
                result = module._import_protocol_pack(data, store, "project-one")
            self.assertTrue(result["created"])
            self.assertEqual(calls[0], ("init", store, "project-one"))
            self.assertEqual(calls[1][1], "job-one")
            self.assertEqual(calls[1][2], root / module.PROTOCOL_PACK_RELATIVE_PATH)
            self.assertEqual(calls[1][3], {
                "expected_manifest_sha256": descriptor["expected_manifest_sha256"],
                "expected_assessment_id": "assessment-7",
            })

    def test_descriptor_hash_tamper_fails_before_service(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, descriptor, _, _ = self._setup(root)
            descriptor["job_id"] = "tampered"
            self._rewrite_descriptor(root, module, descriptor, update_hash=False)
            calls = []
            with mock.patch.dict(sys.modules, {
                "packages.platform.protocol_evidence": self._service_module(calls)
            }), self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_DESCRIPTOR_HASH_INVALID"):
                module._import_protocol_pack(data, object(), "project-one")
            self.assertEqual(calls, [])

    def test_semantic_descriptor_tampering_fails_before_service(self):
        mutations = {
            "wrong-job": (lambda value: value.__setitem__("job_id", "job-other"), "JOB_MISMATCH"),
            "wrong-assessment": (lambda value: value.__setitem__("expected_assessment_id", "assessment-2"), "ASSESSMENT_MISMATCH"),
            "wrong-seed": (lambda value: value.__setitem__("seed_sha256", "0" * 64), "SEED_HASH_MISMATCH"),
            "approval": (lambda value: value.__setitem__("scientific_approval", True), "DESCRIPTOR_INVALID"),
            "traversal": (lambda value: value.__setitem__("pack_relative_path", "../outside"), "DESCRIPTOR_INVALID"),
        }
        for label, (mutate, error) in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                module, data, descriptor, _, _ = self._setup(root)
                mutate(descriptor)
                self._rewrite_descriptor(root, module, descriptor)
                calls = []
                with mock.patch.dict(sys.modules, {
                    "packages.platform.protocol_evidence": self._service_module(calls)
                }), self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_" + error):
                    module._import_protocol_pack(data, object(), "project-one")
                self.assertEqual(calls, [])

    def test_pack_hash_and_unmanifested_file_fail_before_service(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, _, _, output = self._setup(root)
            output.write_bytes(b"tampered")
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_OUTPUT_MANIFEST_INVALID"):
                module._import_protocol_pack(data, object(), "project-one")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, _, _, output = self._setup(root)
            (output.parent / "undeclared.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_PACK_UNMANIFESTED"):
                module._import_protocol_pack(data, object(), "project-one")

    def test_manifested_missing_pack_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, _, _, _ = self._setup(root)
            manifest = self._manifest(root)
            manifest["files"][module.PROTOCOL_PACK_RELATIVE_PATH + "/missing.json"] = "1" * 64
            self._save_manifest(root, manifest)
            with self.assertRaisesRegex(RuntimeError, "PROTOCOL_IMPORT_PACK_MANIFEST_MISMATCH"):
                module._import_protocol_pack(data, object(), "project-one")

    def test_no_descriptor_preserves_v5_without_service_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module = load_portable(root)
            payload = {
                "release_version": "v5",
                "source_followup": module.V5_SOURCE_FOLLOWUP,
                "source_followup_rules_installed": True,
                "scientific_accepted": False,
                "scientific_approval": False,
                "automatic_approval": False,
                "all_pass_claimed": False,
            }
            release = root / module.V5_RELEASE_METADATA_NAME
            release.write_text(json.dumps(payload), encoding="utf-8")
            manifest = {"files": {
                module.V5_RELEASE_METADATA_NAME: hashlib.sha256(release.read_bytes()).hexdigest()
            }}
            (root / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
            self.assertIsNone(module._import_protocol_pack(seed_value(), object(), "project-one"))

    def test_repeat_service_result_does_not_mutate_seed_or_human_sentinel(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module, data, _, _, _ = self._setup(root)
            seed_before = copy.deepcopy(data)
            sentinel = root / "design-data/human-sentinel.json"
            sentinel.parent.mkdir()
            sentinel.write_text('{"approved":false}', encoding="utf-8")
            sentinel_before = sentinel.read_bytes()
            calls = []
            fake = self._service_module(calls, {"created": False, "records": [{"id": "same"}]})
            with mock.patch.dict(sys.modules, {"packages.platform.protocol_evidence": fake}):
                first = module._import_protocol_pack(data, object(), "project-one")
                second = module._import_protocol_pack(data, object(), "project-one")
            self.assertFalse(first["created"])
            self.assertEqual(first, second)
            self.assertEqual(data, seed_before)
            self.assertEqual(sentinel.read_bytes(), sentinel_before)
            self.assertEqual(len([call for call in calls if call[0] == "register"]), 2)


if __name__ == "__main__":
    unittest.main()

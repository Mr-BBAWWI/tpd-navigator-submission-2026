from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from tests import test_v5_archived_runtime as archived_helpers

archived_pair = archived_helpers.archived_pair
load_portable = archived_helpers.load_portable

import build_sprint_package as sprint


class V6ReleaseCompatibilityTests(unittest.TestCase):
    def _write_release(self, root, module, payload=None, raw=None, declared=True, manifest_format=True):
        if raw is None:
            raw = json.dumps(payload if payload is not None else sprint._release_payload()).encode("utf-8")
        path = root / module.V6_RELEASE_METADATA_NAME
        path.write_bytes(raw)
        files = {}
        if declared:
            files[module.V6_RELEASE_METADATA_NAME] = hashlib.sha256(raw).hexdigest()
        manifest = {"files": files}
        if manifest_format:
            manifest["format"] = module.MANIFEST_FORMAT
        (root / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
        return manifest, path

    def _write_v5(self, root, module):
        payload = {
            "release_version": "v5",
            "source_followup": module.V5_SOURCE_FOLLOWUP,
            "source_followup_rules_installed": True,
            "scientific_accepted": False,
            "scientific_approval": False,
            "automatic_approval": False,
            "all_pass_claimed": False,
        }
        raw = json.dumps(payload).encode("utf-8")
        path = root / module.V5_RELEASE_METADATA_NAME
        path.write_bytes(raw)
        return {"files": {module.V5_RELEASE_METADATA_NAME: hashlib.sha256(raw).hexdigest()}}

    def test_real_builder_v6_payload_is_verified_without_v5_file(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            manifest, _ = self._write_release(root, module)
            self.assertFalse((root / module.V5_RELEASE_METADATA_NAME).exists())
            self.assertTrue(module._verified_v5_release(manifest))

    def test_existing_minimal_v5_contract_still_passes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            self.assertTrue(module._verified_v5_release(self._write_v5(root, module)))

    def test_missing_or_unmanifested_release_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            self.assertFalse(module._verified_v5_release({"format": module.MANIFEST_FORMAT, "files": {}}))
            manifest, path = self._write_release(root, module, declared=False)
            self.assertFalse(module._verified_v5_release(manifest))
            path.unlink()
            manifest["files"][module.V6_RELEASE_METADATA_NAME] = "1" * 64
            self.assertFalse(module._verified_v5_release(manifest))

    def test_v6_hash_errors_are_precise(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            manifest, _ = self._write_release(root, module)
            for invalid_hash in ("not-a-hash", "g" * 64, "0" * 64):
                changed = copy.deepcopy(manifest)
                changed["files"][module.V6_RELEASE_METADATA_NAME] = invalid_hash
                with self.subTest(invalid_hash=invalid_hash), self.assertRaisesRegex(
                    RuntimeError, "^V6_RELEASE_METADATA_HASH_INVALID$"
                ):
                    module._verified_v5_release(changed)

    def test_v6_contract_rejects_wrong_fields_and_nonliteral_boolean(self):
        mutations = {
            "version": lambda value: value.__setitem__("release_version", "v5"),
            "artifact": lambda value: value.__setitem__("artifact", "final-expert-package"),
            "runtime_policy": lambda value: value.__setitem__("runtime_policy", {
                "archived_job": "execute archived runtime",
                "original_runtime_distinct": True,
            }),
            "runtime_policy_bool_one": lambda value: value["runtime_policy"].__setitem__(
                "original_runtime_distinct", 1
            ),
            "claim_true": lambda value: value.__setitem__("all_pass_claimed", True),
            "claim_zero": lambda value: value.__setitem__("all_pass_claimed", 0),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                module = load_portable(root)
                payload = copy.deepcopy(sprint._release_payload())
                mutate(payload)
                manifest, _ = self._write_release(root, module, payload=payload)
                with self.assertRaisesRegex(RuntimeError, "^V6_RELEASE_METADATA_INVALID$"):
                    module._verified_v5_release(manifest)

    def test_v6_malformed_json_wrong_type_and_manifest_format_are_rejected(self):
        cases = ((b"{", True), (b"[]", True), (json.dumps(sprint._release_payload()).encode("utf-8"), False))
        for raw, with_format in cases:
            with self.subTest(raw=raw[:8], with_format=with_format), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                module = load_portable(root)
                manifest, _ = self._write_release(root, module, raw=raw, manifest_format=with_format)
                with self.assertRaisesRegex(RuntimeError, "^V6_RELEASE_METADATA_INVALID$"):
                    module._verified_v5_release(manifest)

    def test_both_known_metadata_versions_are_always_ambiguous(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            manifest, _ = self._write_release(root, module)
            v5 = self._write_v5(root, module)
            manifest["files"].update(v5["files"])
            with self.assertRaisesRegex(RuntimeError, "^RELEASE_METADATA_AMBIGUOUS$"):
                module._verified_v5_release(manifest)

    def test_declared_v5_and_on_disk_v6_are_ambiguous(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            self._write_release(root, module)
            manifest = {"files": {module.V5_RELEASE_METADATA_NAME: "0" * 64}}
            with self.assertRaisesRegex(RuntimeError, "^RELEASE_METADATA_AMBIGUOUS$"):
                module._verified_v5_release(manifest)


class V6ArchivedRuntimeIntegrationTests(unittest.TestCase):
    def _setup(self, root):
        module = load_portable(root)
        payload = sprint._release_payload()
        raw = json.dumps(payload).encode("utf-8")
        release = root / module.V6_RELEASE_METADATA_NAME
        release.write_bytes(raw)
        manifest = {
            "format": module.MANIFEST_FORMAT,
            "files": {module.V6_RELEASE_METADATA_NAME: hashlib.sha256(raw).hexdigest()},
        }
        (root / "PACKAGE-MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
        archived, current = archived_pair()
        helper = archived_helpers.ArchivedRuntimePortableTests()
        data = helper._seed(module, archived, current)
        return module, helper, data, archived, current, release

    def test_real_v6_metadata_accepts_exact_archived_runtime_marker(self):
        with tempfile.TemporaryDirectory() as temp:
            module, helper, data, _, current, _ = self._setup(Path(temp))
            with helper._install_identity_modules(current):
                module._verify_current(data)

    def test_removing_v6_metadata_denies_stale_archived_runtime(self):
        with tempfile.TemporaryDirectory() as temp:
            module, helper, data, _, current, release = self._setup(Path(temp))
            release.unlink()
            with helper._install_identity_modules(current), self.assertRaisesRegex(
                RuntimeError, "SCIENTIFIC_SEED_WORKBENCH_ENGINE_CHANGED"
            ):
                module._verify_current(data)

    def test_v6_does_not_waive_workbench_dependency_or_unallowed_file_changes(self):
        changes = {
            "workbench": lambda current: current["files"].__setitem__(
                "packages/platform/workbench.py", "e" * 64
            ),
            "dependency": lambda current: current["dependencies"].__setitem__("numpy", "2.0"),
            "unallowed_file": lambda current: current["files"].__setitem__("packages/science/new.py", "f" * 64),
        }
        for label, change in changes.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                module, helper, data, _, current, _ = self._setup(Path(temp))
                changed = copy.deepcopy(current)
                change(changed)
                with helper._install_identity_modules(changed), self.assertRaisesRegex(
                    RuntimeError, "SCIENTIFIC_SEED_ARCHIVED_RUNTIME_INVALID"
                ):
                    module._verify_current(data)


if __name__ == "__main__":
    unittest.main()

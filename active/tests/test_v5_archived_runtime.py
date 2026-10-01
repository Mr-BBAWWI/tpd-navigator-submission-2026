from __future__ import annotations

import copy
import hashlib
import importlib.util
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

import build_final_expert_package as final_expert
import build_scientific_package as scientific

PORTABLE_SOURCE = DELIVERY / "scientific_portable.py"
ALLOWED = final_expert.ARCHIVED_RUNTIME_ALLOWED_FILES


def runtime(files=None, dependencies=None, data_sources=None, version="runtime-v1", prompt="prompt-v1"):
    return {
        "version": version,
        "prompt": prompt,
        "files": files or {
            "packages/platform/workbench.py": "1" * 64,
            "packages/platform/design_panel.py": "2" * 64,
            "packages/science/design_docking.py": "3" * 64,
        },
        "dependencies": dependencies or {"numpy": "1.0"},
        "data_sources": data_sources or {"source": "4" * 64},
    }


def archived_pair():
    current = runtime()
    archived = copy.deepcopy(current)
    archived["files"]["packages/platform/design_panel.py"] = "a" * 64
    archived["files"]["packages/science/design_docking.py"] = "b" * 64
    return archived, current


def load_portable(root: Path):
    program = root / "program"
    program.mkdir(parents=True)
    design = types.ModuleType("design_portable")
    design._immutable_files = lambda: {
        path.relative_to(root).as_posix(): path
        for path in root.rglob("*")
        if path.is_file() and path.name != "PACKAGE-MANIFEST.json"
    }
    design.dump = lambda path, value: path.write_text(json.dumps(value), encoding="utf-8")
    design.check = lambda **kwargs: None
    design.prepare = lambda **kwargs: None
    design.start = lambda **kwargs: None
    design.stop = lambda: None
    design.ping = lambda port: {}
    old = sys.modules.get("design_portable")
    sys.modules["design_portable"] = design
    try:
        spec = importlib.util.spec_from_file_location("scientific_portable_v5_test", PORTABLE_SOURCE)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
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


class ArchivedRuntimeBuilderTests(unittest.TestCase):
    def test_marker_records_exact_actual_diff_and_hash(self):
        archived, current = archived_pair()
        marker = final_expert._archived_runtime_review(archived, current)
        self.assertEqual(set(marker["changed_files"]), ALLOWED)
        self.assertEqual(marker["archived_runtime_sha256"], scientific.sha256_bytes(scientific.canonical(archived)))
        self.assertEqual(marker["current_runtime_sha256"], scientific.sha256_bytes(scientific.canonical(current)))
        self.assertFalse(marker["scientific_approval"])

    def test_empty_extra_and_dependency_diffs_are_rejected(self):
        current = runtime()
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_DIFF_INVALID"):
            final_expert._archived_runtime_review(copy.deepcopy(current), current)
        archived, current = archived_pair()
        archived["files"]["packages/science/extra.py"] = "c" * 64
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_DIFF_INVALID"):
            final_expert._archived_runtime_review(archived, current)
        archived, current = archived_pair()
        archived["dependencies"]["numpy"] = "0.9"
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_IDENTITY_MISMATCH"):
            final_expert._archived_runtime_review(archived, current)
        archived, current = archived_pair()
        archived["previously_unreviewed_key"] = {"value": True}
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_IDENTITY_MISMATCH"):
            final_expert._archived_runtime_review(archived, current)

    def test_wrong_marker_hash_is_rejected(self):
        archived, current = archived_pair()
        marker = final_expert._archived_runtime_review(archived, current)
        marker["changed_files"]["packages/platform/design_panel.py"]["archived_sha256"] = "f" * 64
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_MARKER_INVALID"):
            final_expert._validate_archived_runtime_review(marker, archived, current)

    def test_freshness_accepts_only_verified_marker(self):
        archived, current = archived_pair()
        job = {"runtime": archived}
        view = {"freshness": {"current": False}}
        marker = final_expert._archived_runtime_review(archived, current)
        final_expert._require_export_runtime(job, view, current, marker)
        self.assertEqual(job["runtime"], archived)
        with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_SOURCE_NOT_CURRENT"):
            final_expert._require_export_runtime(job, view, current, None)
        job = {"runtime": current}
        view["freshness"]["current"] = True
        with self.assertRaisesRegex(RuntimeError, "ARCHIVED_RUNTIME_MARKER_NOT_ALLOWED_FOR_CURRENT_JOB"):
            final_expert._require_export_runtime(job, view, current, marker)


class ArchivedRuntimePortableTests(unittest.TestCase):
    def _seed(self, module, archived, current):
        marker = final_expert._archived_runtime_review(archived, current)
        job = {
            "id": "job-1", "project_id": "project", "operation": "design_panel",
            "state": "completed", "runtime": archived, "input_digest": "input-now",
        }
        return {
            "job": job,
            "source_binding": {
                "project_id": "project", "job_id": "job-1", "operation": "design_panel",
                "request_key": "request", "fingerprint": "fingerprint",
                "workbench_runtime": archived, "scientific_engine_fingerprint": {"engine": "now"},
            },
            "archived_runtime_review": marker,
        }

    def _install_identity_modules(self, current):
        apps = types.ModuleType("apps")
        api = types.ModuleType("apps.api")
        main = types.ModuleType("apps.api.main")
        main.PROJECT = "project"
        packages = types.ModuleType("packages")
        platform = types.ModuleType("packages.platform")
        workbench = types.ModuleType("packages.platform.workbench")
        workbench.identity = lambda: current
        acceptance = types.ModuleType("packages.platform.scientific_acceptance")
        acceptance.engine_fingerprint = lambda: {"engine": "now"}
        design_panel = types.ModuleType("packages.platform.design_panel")
        design_panel.input_binding = lambda: {"digest": "input-now"}
        return mock.patch.dict(sys.modules, {
            "apps": apps, "apps.api": api, "apps.api.main": main,
            "packages": packages, "packages.platform": platform,
            "packages.platform.workbench": workbench,
            "packages.platform.scientific_acceptance": acceptance,
            "packages.platform.design_panel": design_panel,
        })

    def test_marked_v5_archived_runtime_is_accepted(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            archived, current = archived_pair()
            data = self._seed(module, archived, current)
            with mock.patch.object(module, "_verified_v5_release", return_value=True), self._install_identity_modules(current):
                module._verify_current(data)

    def test_unmarked_stale_runtime_retains_original_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            archived, current = archived_pair()
            data = self._seed(module, archived, current)
            data.pop("archived_runtime_review")
            with self._install_identity_modules(current):
                with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_WORKBENCH_ENGINE_CHANGED"):
                    module._verify_current(data)

    def test_marker_tamper_and_new_source_runtime_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            archived, current = archived_pair()
            data = self._seed(module, archived, current)
            data["archived_runtime_review"]["current_runtime_sha256"] = "0" * 64
            with mock.patch.object(module, "_verified_v5_release", return_value=True), self._install_identity_modules(current):
                with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_ARCHIVED_RUNTIME_INVALID"):
                    module._verify_current(data)

    def test_archived_runtime_does_not_exempt_input_drift(self):
        with tempfile.TemporaryDirectory() as temp:
            module = load_portable(Path(temp))
            archived, current = archived_pair()
            data = self._seed(module, archived, current)
            with mock.patch.object(module, "_verified_v5_release", return_value=True), self._install_identity_modules(current):
                from packages.platform import design_panel
                design_panel.input_binding = lambda: {"digest": "different-input"}
                with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_SOURCE_INPUTS_CHANGED"):
                    module._verify_current(data)
            data = self._seed(module, archived, current)
            changed_current = copy.deepcopy(current)
            changed_current["files"]["packages/platform/workbench.py"] = "e" * 64
            with mock.patch.object(module, "_verified_v5_release", return_value=True), self._install_identity_modules(changed_current):
                with self.assertRaisesRegex(RuntimeError, "SCIENTIFIC_SEED_ARCHIVED_RUNTIME_INVALID"):
                    module._verify_current(data)

    def test_v5_release_requires_false_claims_followup_and_manifested_marker(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            module = load_portable(root)
            release = {
                "release_version": "v5", "source_followup": "cases/expert_opinions/20261001",
                "source_followup_rules_installed": True, "scientific_accepted": False,
                "scientific_approval": False, "automatic_approval": False, "all_pass_claimed": False,
            }
            path = root / module.V5_RELEASE_METADATA_NAME
            path.write_text(json.dumps(release), encoding="utf-8")
            manifest = {"files": {module.V5_RELEASE_METADATA_NAME: hashlib.sha256(path.read_bytes()).hexdigest()}}
            self.assertTrue(module._verified_v5_release(manifest))
            release["all_pass_claimed"] = True
            path.write_text(json.dumps(release), encoding="utf-8")
            manifest["files"][module.V5_RELEASE_METADATA_NAME] = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "V5_RELEASE_METADATA_INVALID"):
                module._verified_v5_release(manifest)


if __name__ == "__main__":
    unittest.main()

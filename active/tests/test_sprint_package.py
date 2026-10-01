"""Packaging-boundary tests; fixtures do not claim real runtime science."""
from __future__ import annotations

import errno
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
DELIVERY = ROOT / "reviewer_delivery"
if str(DELIVERY) not in sys.path:
    sys.path.insert(0, str(DELIVERY))

import build_sprint_package as sprint


def _spec(directory: Path, source: Path, destination="evidence/sprint/run") -> Path:
    path = directory / "extras.json"
    path.write_text(json.dumps({
        "format": "sprint-extras/1",
        "trees": [{"source": str(source.absolute()), "destination": destination}],
    }), encoding="utf-8")
    return path


def _record(path: Path) -> dict:
    raw = path.read_bytes()
    return {
        "path": path,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
    }


class SprintPackageTests(unittest.TestCase):
    def test_wrapper_globals_restored_on_exception(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            before = sprint._wrapper_globals()
            with mock.patch.object(
                sprint.v5, "stage_package", side_effect=RuntimeError("fixture")
            ):
                with self.assertRaisesRegex(RuntimeError, "fixture"):
                    sprint.stage_package(root, False)
            self.assertEqual(sprint._wrapper_globals(), before)
            self.assertEqual(sprint.v5.RELEASE_VERSION, "v5")
            self.assertEqual(sprint.v5.ARCHIVE_NAME, "TPD_expert_v5.zip")

    def test_release_identity_and_readme_are_truthful(self):
        payload = sprint._release_payload()
        text = sprint._readme().decode()
        self.assertEqual(payload["release_version"], "v6")
        self.assertFalse(payload["scientific_approval"])
        self.assertFalse(payload["all_pass_claimed"])
        self.assertIn("127.0.0.1", text)
        self.assertIn("DACON_API_KEY", text)
        self.assertIn("10_연구비교실_열기.cmd", text)
        self.assertIn("3개 query 버튼", text)

    def test_destination_raw_normalization_and_overlap_blocked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            source = workspace / "chosen"
            source.mkdir(parents=True)
            (source / "result.json").write_text("{}", encoding="utf-8")
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                for destination in (
                    "evidence/sprint/../escape",
                    "evidence//sprint/run",
                    "evidence/sprint/./run",
                ):
                    with self.subTest(destination=destination):
                        with self.assertRaisesRegex(RuntimeError, "DESTINATION_INVALID"):
                            sprint._validate_extra_spec(
                                _spec(root, source, destination), stage, stage.parent
                            )
                overlap = root / "overlap.json"
                overlap.write_text(json.dumps({
                    "format": "sprint-extras/1",
                    "trees": [
                        {"source": str(source.absolute()), "destination": "evidence/sprint/a"},
                        {"source": str(source.absolute()), "destination": "evidence/sprint/a/b"},
                    ],
                }), encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "OVERLAP"):
                    sprint._validate_extra_spec(overlap, stage, stage.parent)

    def test_source_must_be_inside_actual_workspace(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            workspace.mkdir()
            outside = root / "outside"
            outside.mkdir()
            (outside / "result.json").write_text("{}", encoding="utf-8")
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                with self.assertRaisesRegex(RuntimeError, "SOURCE_INVALID"):
                    sprint._validate_extra_spec(_spec(root, outside), stage, stage.parent)

    def test_stage_source_and_source_containing_stage_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            (stage / "result.json").write_text("{}", encoding="utf-8")
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                for source in (stage, workspace / "out"):
                    with self.subTest(source=source):
                        with self.assertRaisesRegex(RuntimeError, "STAGE_LOOP"):
                            sprint._validate_extra_spec(_spec(root, source), stage, stage.parent)

    def test_symlink_file_blocked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            source = workspace / "chosen"
            source.mkdir(parents=True)
            target = workspace / "actual.txt"
            target.write_text("fixture", encoding="utf-8")
            link = source / "linked.txt"
            try:
                link.symlink_to(target)
            except OSError as exc:
                permission_denied = exc.errno in {errno.EACCES, errno.EPERM} or getattr(exc, "winerror", None) == 1314
                if os.name == "nt" and permission_denied:
                    self.skipTest("Windows account lacks symlink permission")
                raise
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                with self.assertRaisesRegex(RuntimeError, "LINK|BLOCKED"):
                    sprint._validate_extra_spec(_spec(root, source), stage, stage.parent)

    def test_linked_source_ancestor_is_checked_before_resolve(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            actual = workspace / "actual" / "chosen"
            actual.mkdir(parents=True)
            (actual / "result.json").write_text("{}", encoding="utf-8")
            linked_parent = workspace / "linked"
            try:
                linked_parent.symlink_to(actual.parent, target_is_directory=True)
            except OSError as exc:
                permission_denied = exc.errno in {errno.EACCES, errno.EPERM} or getattr(exc, "winerror", None) == 1314
                if os.name == "nt" and permission_denied:
                    self.skipTest("Windows account lacks symlink permission")
                raise
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            source = linked_parent / "chosen"
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                with self.assertRaisesRegex(RuntimeError, "LINK"):
                    sprint._validate_extra_spec(_spec(root, source), stage, stage.parent)

    def test_whole_localdata_rejected_but_explicit_child_allowed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            local = workspace / ".localdata"
            chosen = local / "explicit-run"
            chosen.mkdir(parents=True)
            (chosen / "artifact.pdb").write_text("ATOM fixture", encoding="utf-8")
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                with self.assertRaisesRegex(RuntimeError, "WHOLE_LOCALDATA"):
                    sprint._validate_extra_spec(_spec(root, local), stage, stage.parent)
                plans = sprint._validate_extra_spec(_spec(root, chosen), stage, stage.parent)
            self.assertEqual(plans[0][0], chosen.resolve())

    def test_blocked_source_basename_and_secret_filename_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            blocked = workspace / ".venv"
            blocked.mkdir()
            (blocked / "result.txt").write_text("fixture", encoding="utf-8")
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                with self.assertRaisesRegex(RuntimeError, "BLOCKED_COMPONENT"):
                    sprint._validate_extra_spec(_spec(root, blocked), stage, stage.parent)
                source = workspace / "chosen"
                source.mkdir()
                (source / ".env.production").write_text("fixture", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "SUSPECTED_SECRET|BLOCKED"):
                    sprint._validate_extra_spec(_spec(root, source), stage, stage.parent)

    def test_actual_key_scans_binary_before_any_copy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            first = source / "first.bin"
            first.write_bytes(b"safe")
            second = source / "second.bin"
            second.write_bytes(b"prefix-real-secret-key-suffix")
            stage = root / "stage"
            stage.mkdir()
            plan = [(source, sprint.PurePosixPath("evidence/sprint/run"), [_record(first), _record(second)])]
            with self.assertRaisesRegex(RuntimeError, "CONFIGURED_KEY_DETECTED"):
                sprint._copy_extras(stage, plan, b"real-secret-key")
            self.assertFalse((stage / "evidence").exists())

    def test_copy_refuses_overwrite_and_manifest_has_security_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            artifact = source / "candidate.json"
            artifact.write_text('{"artifact_id":"legitimate-keylike-id"}', encoding="utf-8")
            stage = root / "stage"
            stage.mkdir()
            plan = [(source, sprint.PurePosixPath("evidence/sprint/run"), [_record(artifact)])]
            sprint._copy_extras(stage, plan)
            manifest = json.loads((stage / sprint.EXTRAS_MANIFEST_NAME).read_text(encoding="utf-8"))
            copied = stage / "evidence/sprint/run/candidate.json"
            self.assertEqual(manifest["files"][0]["source"], "source")
            self.assertEqual(manifest["files"][0]["bytes"], artifact.stat().st_size)
            self.assertEqual(manifest["files"][0]["sha256"], hashlib.sha256(artifact.read_bytes()).hexdigest())
            self.assertEqual(manifest["files"][0]["sha256"], hashlib.sha256(copied.read_bytes()).hexdigest())
            self.assertNotIn(str(source.resolve()), json.dumps(manifest))
            self.assertFalse(manifest["security"]["actual_key_scan_performed"])
            with self.assertRaises((FileExistsError, RuntimeError)):
                sprint._copy_extras(stage, plan)

    def test_changed_source_after_validation_rejected_without_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workspace = root / "workspace"
            source = workspace / "chosen"
            source.mkdir(parents=True)
            artifact = source / "candidate.bin"
            artifact.write_bytes(b"frozen-input")
            stage = workspace / "out" / "TPD"
            stage.mkdir(parents=True)
            with mock.patch.object(sprint, "WORKSPACE", workspace):
                plans = sprint._validate_extra_spec(_spec(root, source), stage, stage.parent)
            artifact.write_bytes(b"changed-input")
            with self.assertRaisesRegex(RuntimeError, "SOURCE_CHANGED"):
                sprint._copy_extras(stage, plans)
            self.assertFalse((stage / "evidence").exists())
            self.assertFalse((stage / sprint.EXTRAS_MANIFEST_NAME).exists())

    def test_later_destination_conflict_is_found_before_first_copy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            first = source / "a.bin"
            second = source / "b.bin"
            first.write_bytes(b"first")
            second.write_bytes(b"second")
            stage = root / "stage"
            conflict = stage / "evidence/sprint/run/b.bin"
            conflict.parent.mkdir(parents=True)
            conflict.write_bytes(b"existing")
            plan = [(source, sprint.PurePosixPath("evidence/sprint/run"), [
                _record(first), _record(second),
            ])]
            with self.assertRaisesRegex(RuntimeError, "DESTINATION_EXISTS"):
                sprint._copy_extras(stage, plan)
            self.assertFalse((stage / "evidence/sprint/run/a.bin").exists())
            self.assertEqual(conflict.read_bytes(), b"existing")
            self.assertFalse((stage / sprint.EXTRAS_MANIFEST_NAME).exists())

    def test_safe_binary_actual_key_scan_reports_performed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            artifact = source / "candidate.bin"
            artifact.write_bytes(b"binary-safe-payload\x00\xff")
            stage = root / "stage"
            stage.mkdir()
            plan = [(source, sprint.PurePosixPath("evidence/sprint/run"), [_record(artifact)])]
            sprint._copy_extras(stage, plan, b"configured-key-not-present")
            manifest = json.loads((stage / sprint.EXTRAS_MANIFEST_NAME).read_text(encoding="utf-8"))
            self.assertTrue(manifest["security"]["actual_key_scan_performed"])
            self.assertEqual(
                (stage / "evidence/sprint/run/candidate.bin").read_bytes(),
                artifact.read_bytes(),
            )

    def test_resume_rejects_mismatched_stage_without_calling_v5(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / sprint.STAGE_NAME
            stage.mkdir()
            called = []
            with mock.patch.object(sprint.v5, "stage_package", side_effect=lambda *args: called.append(args)), \
                    mock.patch.object(sprint, "_validate_v6_stage", side_effect=RuntimeError("PACKAGE_SOURCE_MISMATCH")):
                with self.assertRaisesRegex(RuntimeError, "SOURCE_MISMATCH"):
                    sprint.stage_package(root, True)
            self.assertFalse(called)

    def test_archive_restores_make_manifest_and_v5_globals(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / sprint.STAGE_NAME
            stage.mkdir()
            for rel in (
                "evidence/final-reports/research-campaign/index.html",
                "evidence/final-reports/demo/TPD_Navigator_시연.mp4",
                "evidence/final-reports/report.html",
            ):
                path = stage / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture")
            original_make = sprint.scientific.make_manifest
            before = sprint._wrapper_globals()
            with mock.patch.object(sprint, "_validate_v6_stage", return_value={}), \
                    mock.patch.object(sprint.v5, "archive_package", side_effect=RuntimeError("archive fixture")):
                with self.assertRaisesRegex(RuntimeError, "archive fixture"):
                    sprint.archive_package(root, None, None, None)
            self.assertIs(sprint.scientific.make_manifest, original_make)
            self.assertEqual(sprint._wrapper_globals(), before)


if __name__ == "__main__":
    unittest.main()

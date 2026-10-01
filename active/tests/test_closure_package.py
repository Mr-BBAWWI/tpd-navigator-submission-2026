"""Synthetic packaging tests; no fixture is a scientific result or approval."""
from __future__ import annotations

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

import build_closure_package as closure
import build_sprint_package as sprint


class ClosurePackageTests(unittest.TestCase):
    def test_identity_and_no_video_policy(self):
        payload = closure._release_payload()
        text = closure._readme().decode("utf-8")
        self.assertEqual(payload["release_version"], "v7")
        self.assertEqual(payload["archive_name"], "TPD_closure_v7.zip")
        self.assertFalse(payload["reports"]["video_required"])
        self.assertIn("demo-runbook.md", text)
        self.assertNotIn("시연영상_열기", text)

    def test_stage_validation_uses_v5_not_v6_video_validator(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary)
            files = {
                closure.RELEASE_METADATA_NAME: b"{}",
                "program/sprint_cpu_portable.py": b"fixture",
                "02_CPU_직접시험.cmd": b"program\\sprint_cpu_portable.py",
                "09_최종보고서_열기.cmd": (b'set "TARGET=%~dp0evidence\\final-reports\\report.html"\n'
                    b'pathlib.Path(sys.argv[1]).resolve().as_uri()+sys.argv[2]\n'
                    b'"%TARGET%" "#report"'),
                "10_연구비교실_열기.cmd": b"campaign",
                "11_수동시연가이드_열기.cmd": (b'set "TARGET=%~dp0evidence\\final-reports\\report.html"\n'
                    b'pathlib.Path(sys.argv[1]).resolve().as_uri()+sys.argv[2]\n'
                    b'"%TARGET%" "#demo"'),
            }
            for name, raw in files.items():
                path = stage / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(raw)
            manifest = {"files": {"program/sprint_cpu_portable.py": "0" * 64}}
            with mock.patch.object(closure.v5, "_validate_v5_stage", return_value=manifest), \
                    mock.patch.object(sprint, "_validate_v6_stage", side_effect=AssertionError("must not run")):
                self.assertIs(closure._validate_v7_stage(stage), manifest)

    def test_v7_reports_require_manual_guide_and_forbid_video(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary)
            for name in (
                "evidence/final-reports/final-report.md",
                "evidence/final-reports/demo-runbook.md",
                "evidence/final-reports/report.html",
                "evidence/final-reports/research-campaign/index.html",
            ):
                path = stage / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic fixture", encoding="utf-8")
            closure._require_final_reports(stage)
            video = stage / "evidence/final-reports/demo/fake.mp4"
            video.parent.mkdir(parents=True)
            video.write_bytes(b"synthetic")
            with self.assertRaisesRegex(RuntimeError, "VIDEO_FORBIDDEN"):
                closure._require_final_reports(stage)
            video.unlink()
            (stage / "evidence/final-reports/demo-runbook.md").unlink()
            with self.assertRaisesRegex(RuntimeError, "REPORT_MISSING"):
                closure._require_final_reports(stage)

    def test_v6_still_requires_video(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary)
            for name in (
                "evidence/final-reports/research-campaign/index.html",
                "evidence/final-reports/report.html",
            ):
                path = stage / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "TPD_Navigator"):
                sprint._require_final_reports(stage)

    def test_stage_globals_restore_on_failure(self):
        before = closure._wrapper_globals()
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(closure.v5, "stage_package", side_effect=RuntimeError("synthetic")):
            with self.assertRaisesRegex(RuntimeError, "synthetic"):
                closure.stage_package(Path(temporary), False)
        self.assertEqual(closure._wrapper_globals(), before)
        self.assertEqual(closure.v6.RELEASE_VERSION, "v6")
        self.assertEqual(closure.v5.RELEASE_VERSION, "v5")

    def test_archive_restores_manifest_hook_on_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = closure.scientific.make_manifest
            before = closure._wrapper_globals()
            with mock.patch.object(closure, "_validate_v7_stage", return_value={}), \
                    mock.patch.object(closure.v5, "archive_package", side_effect=RuntimeError("synthetic")):
                with self.assertRaisesRegex(RuntimeError, "synthetic"):
                    closure.archive_package(root, None, None, None)
            self.assertIs(closure.scientific.make_manifest, original)
            self.assertEqual(closure._wrapper_globals(), before)

    def test_archive_descriptor_fails_before_core_mutation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / closure.STAGE_NAME
            descriptor = stage / sprint.PROTOCOL_IMPORT_DESCRIPTOR
            descriptor.parent.mkdir(parents=True)
            descriptor.write_text("{}", encoding="utf-8")
            called = []
            with mock.patch.object(closure, "_validate_v7_stage", return_value={}), \
                    mock.patch.object(closure.v5, "archive_package", side_effect=lambda *a: called.append(a)):
                with self.assertRaisesRegex(RuntimeError, "DESCRIPTOR_FORBIDDEN"):
                    closure.archive_package(root, None, None, None)
            self.assertEqual(called, [])

    def test_stage_no_overwrite_delegates_unchanged_core(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / closure.STAGE_NAME
            stage.mkdir()
            (stage / "existing").write_text("keep", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                closure.stage_package(root, False)
            self.assertEqual((stage / "existing").read_text(encoding="utf-8"), "keep")


if __name__ == "__main__":
    unittest.main()

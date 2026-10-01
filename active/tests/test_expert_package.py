from __future__ import annotations

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

import build_acceptance_package as acceptance
import build_scientific_package as scientific
import build_expert_package as expert


class ExpertPackageTests(unittest.TestCase):
    def test_schema3_and_v4_defaults_are_intentional(self):
        self.assertEqual(expert.FORMAT, "tpd-scientific-delivery/3")
        self.assertEqual(expert.SEED_FORMAT, "scientific-seed/3")
        self.assertEqual(expert.DEFAULT_OUT.name, "최종실행본_20260930_v4")
        self.assertEqual(expert.STAGE_NAME, "TPD")
        self.assertEqual(expert.ARCHIVE_NAME, "TPD_expert_v4.zip")

    def test_scientific_globals_restored_when_delegate_errors(self):
        original = expert._wrapper_globals()

        def fail(*args, **kwargs):
            self.assertEqual(scientific.DEFAULT_OUT, expert.DEFAULT_OUT)
            self.assertEqual(scientific.STAGE_NAME, "TPD")
            self.assertEqual(scientific.ARCHIVE_NAME, "TPD_expert_v4.zip")
            self.assertIs(scientific._readme, expert._readme)
            raise RuntimeError("boom")

        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            scientific, "stage_package", side_effect=fail
        ):
            with self.assertRaisesRegex(RuntimeError, "boom"):
                expert.stage_package(Path(temporary) / "new-v4", False)
        self.assertEqual(expert._wrapper_globals(), original)

    def test_scientific_globals_restored_for_seed_and_archive_errors(self):
        original = expert._wrapper_globals()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(
                scientific, "export_scientific_seed", side_effect=RuntimeError("seed")
            ):
                with self.assertRaisesRegex(RuntimeError, "seed"):
                    expert.export_scientific_seed(root, root / "data", "job")
            self.assertEqual(expert._wrapper_globals(), original)
            with mock.patch.object(
                scientific, "archive_package", side_effect=RuntimeError("archive")
            ):
                with self.assertRaisesRegex(RuntimeError, "archive"):
                    expert.archive_package(root, None, None)
            self.assertEqual(expert._wrapper_globals(), original)

    def test_evidence_launcher_is_relative_checked_and_quoted(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "한글 경로" / "TPD"
            stage.mkdir(parents=True)
            expert._evidence_launcher(stage)
            text = (stage / "06_검수자료_열기.cmd").read_text(encoding="utf-8")
        target = "%~dp0evidence\\final-reports\\expert-packet\\index.html"
        self.assertIn(f'set "TARGET={target}"', text)
        self.assertIn('if not exist "%TARGET%" (', text)
        self.assertIn('start "" "%TARGET%"', text)
        self.assertNotIn("http://", text)
        self.assertNotIn("https://", text)

    def test_stage_refuses_nonempty_package_before_delegation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old_v3 = root / "실행본_과학검토_v3"
            old_v3.mkdir()
            marker = old_v3 / "keep.txt"
            marker.write_text("unchanged", encoding="utf-8")
            output = root / "new-v4"
            stage = output / expert.STAGE_NAME
            stage.mkdir(parents=True)
            sentinel = stage / "existing.txt"
            sentinel.write_text("do-not-overwrite", encoding="utf-8")

            delegate = mock.Mock(side_effect=AssertionError("must not run"))
            with mock.patch.object(scientific, "stage_package", delegate):
                with self.assertRaisesRegex(
                    RuntimeError, "INCOMPLETE_NONEMPTY_STAGE_EXISTS"
                ):
                    expert.stage_package(output, False)
            delegate.assert_not_called()
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "do-not-overwrite")
            self.assertEqual(marker.read_text(encoding="utf-8"), "unchanged")

    def test_readme_has_plain_korean_scope_and_no_all_pass_claim(self):
        text = expert._readme().decode("utf-8")
        self.assertIn("warhead analog 3개", text)
        self.assertIn("조립된 구조는 36개", text)
        self.assertIn("신규 GPU 계산의 실제 완료·실패 수와 조건은 최종보고서 실측값에서 확인합니다", text)
        self.assertNotIn("신규 GPU 계산 30건은 아직 실행 대기", text)
        self.assertIn("FX5/6HAZ", text)
        self.assertIn("site 19의 최소값은 1", text)
        self.assertIn("site 8, 11, 12는 UNKNOWN", text)
        self.assertIn("전체 통과 또는 자동 승인을 주장하지 않습니다", text)
        self.assertNotIn("FX5onlyscope", text)

    def test_release_metadata_explicitly_denies_all_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary)
            expert._write_release_metadata(stage)
            payload = json.loads(
                (stage / "RELEASE-VERSION-v4.json").read_text(encoding="utf-8")
            )
        self.assertFalse(payload["scientific_approval"])
        self.assertFalse(payload["automatic_approval"])
        self.assertFalse(payload["all_pass_claimed"])
        self.assertEqual(payload["new_gpu_runs"], {"count": 30, "status": "see_final_report"})


if __name__ == "__main__":
    unittest.main()

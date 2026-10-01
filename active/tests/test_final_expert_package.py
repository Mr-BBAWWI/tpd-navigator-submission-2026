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

import build_expert_package as expert_v4
import build_scientific_package as scientific
import build_final_expert_package as final_expert


class FinalExpertPackageTests(unittest.TestCase):
    def test_v5_defaults_retain_v3_package_and_seed_compatibility(self):
        self.assertEqual(final_expert.FORMAT, "tpd-scientific-delivery/3")
        self.assertEqual(final_expert.SEED_FORMAT, "scientific-seed/3")
        self.assertEqual(final_expert.DEFAULT_OUT.name, "최종실행본_20261001_v5")
        self.assertEqual(final_expert.STAGE_NAME, "TPD")
        self.assertEqual(final_expert.ARCHIVE_NAME, "TPD_expert_v5.zip")
        self.assertEqual(final_expert.RELEASE_VERSION, "v5")

    def test_v4_builder_defaults_are_not_modified(self):
        before = (
            expert_v4.DEFAULT_OUT,
            expert_v4.STAGE_NAME,
            expert_v4.ARCHIVE_NAME,
            expert_v4.RELEASE_VERSION,
        )
        self.assertEqual(expert_v4.DEFAULT_OUT.name, "최종실행본_20260930_v4")
        self.assertEqual(expert_v4.ARCHIVE_NAME, "TPD_expert_v4.zip")
        self.assertEqual(
            before,
            (
                expert_v4.DEFAULT_OUT,
                expert_v4.STAGE_NAME,
                expert_v4.ARCHIVE_NAME,
                expert_v4.RELEASE_VERSION,
            ),
        )

    def test_scientific_globals_restore_when_stage_delegate_fails(self):
        original = final_expert._wrapper_globals()

        def fail(*args, **kwargs):
            self.assertEqual(scientific.DEFAULT_OUT, final_expert.DEFAULT_OUT)
            self.assertEqual(scientific.STAGE_NAME, "TPD")
            self.assertEqual(scientific.ARCHIVE_NAME, "TPD_expert_v5.zip")
            self.assertIs(scientific._readme, final_expert._readme)
            raise RuntimeError("delegated-stage-failure")

        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            scientific, "stage_package", side_effect=fail
        ):
            with self.assertRaisesRegex(RuntimeError, "delegated-stage-failure"):
                final_expert.stage_package(Path(temporary) / "v5", False)
        self.assertEqual(final_expert._wrapper_globals(), original)

    def test_scientific_globals_restore_for_seed_and_archive_failures(self):
        original = final_expert._wrapper_globals()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(
                final_expert, "_validate_v5_stage", side_effect=RuntimeError("seed-failure")
            ):
                with self.assertRaisesRegex(RuntimeError, "seed-failure"):
                    final_expert.export_scientific_seed(root, root / "data", "job")
            self.assertEqual(final_expert._wrapper_globals(), original)

            with mock.patch.object(
                final_expert, "_validate_v5_stage", side_effect=RuntimeError("archive-failure")
            ):
                with self.assertRaisesRegex(RuntimeError, "archive-failure"):
                    final_expert.archive_package(root, None, None)
            self.assertEqual(final_expert._wrapper_globals(), original)

    def test_stage_refuses_nonempty_v5_and_preserves_existing_v4_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old_v4 = root / "최종실행본_20260930_v4"
            old_v4.mkdir()
            old_source = old_v4 / "source.txt"
            old_store = old_v4 / "store.txt"
            old_result = old_v4 / "result.txt"
            old_source.write_text("source-v4", encoding="utf-8")
            old_store.write_text("store-v4", encoding="utf-8")
            old_result.write_text("result-v4", encoding="utf-8")

            output = root / "최종실행본_20261001_v5"
            stage = output / final_expert.STAGE_NAME
            stage.mkdir(parents=True)
            sentinel = stage / "existing.txt"
            sentinel.write_text("do-not-overwrite", encoding="utf-8")

            delegate = mock.Mock(side_effect=AssertionError("must not run"))
            with mock.patch.object(scientific, "stage_package", delegate):
                with self.assertRaisesRegex(
                    RuntimeError, "INCOMPLETE_NONEMPTY_STAGE_EXISTS"
                ):
                    final_expert.stage_package(output, False)

            delegate.assert_not_called()
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "do-not-overwrite")
            self.assertEqual(old_source.read_text(encoding="utf-8"), "source-v4")
            self.assertEqual(old_store.read_text(encoding="utf-8"), "store-v4")
            self.assertEqual(old_result.read_text(encoding="utf-8"), "result-v4")

    def test_fresh_stage_adds_v5_metadata_and_launchers_without_replacing_01_to_07(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "v5"
            stage = output / "TPD"

            def delegated(root, resume):
                self.assertEqual(Path(root), output)
                self.assertFalse(resume)
                stage.mkdir(parents=True)
                for name in (
                    "01_프로그램_실행.cmd",
                    "02_CPU_직접시험.cmd",
                    "03_프로그램_종료.cmd",
                    "04_파일검사.cmd",
                    "05_과학검토_열기.cmd",
                    "07_기존_C01_C02_결과_보기.cmd",
                ):
                    (stage / name).write_text("delegated", encoding="utf-8")
                return stage

            with mock.patch.object(
                scientific, "stage_package", side_effect=delegated
            ), mock.patch.object(scientific, "make_manifest") as make_manifest:
                result = final_expert.stage_package(output, False)

            self.assertEqual(result, stage)
            make_manifest.assert_called_once_with(stage)
            for name in (
                "01_프로그램_실행.cmd",
                "02_CPU_직접시험.cmd",
                "03_프로그램_종료.cmd",
                "04_파일검사.cmd",
                "05_과학검토_열기.cmd",
                "07_기존_C01_C02_결과_보기.cmd",
            ):
                self.assertEqual((stage / name).read_text(encoding="utf-8"), "delegated")
            for name in (
                "06_검수자료_열기.cmd",
                "08_상태비교_열기.cmd",
                "09_최종보고서_열기.cmd",
                "RELEASE-VERSION-v5.json",
            ):
                self.assertTrue((stage / name).is_file(), name)

    def test_report_launchers_are_relative_exist_checked_quoted_and_utf8_safe(self):
        expected = {
            "06_검수자료_열기.cmd": "evidence\\final-reports\\expert-packet\\index.html",
            "08_상태비교_열기.cmd": "evidence\\final-reports\\analog-states\\index.html",
            "09_최종보고서_열기.cmd": "evidence\\final-reports\\report.html",
        }
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "한글 경로" / "TPD"
            stage.mkdir(parents=True)
            final_expert._stage_release_files(stage)
            for name, relative in expected.items():
                text = (stage / name).read_text(encoding="utf-8")
                self.assertIn("chcp 65001 >nul", text)
                self.assertIn('cd /d "%~dp0"', text)
                self.assertIn(f'set "TARGET=%~dp0{relative}"', text)
                self.assertIn('if not exist "%TARGET%" (', text)
                self.assertIn('start "" "%TARGET%"', text)
                self.assertNotIn("http://", text)
                self.assertNotIn("https://", text)

    def test_release_metadata_denies_scientific_completion_and_records_followup(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary)
            final_expert._write_release_metadata(stage)
            payload = json.loads(
                (stage / "RELEASE-VERSION-v5.json").read_text(encoding="utf-8")
            )

        self.assertEqual(payload["release_version"], "v5")
        self.assertEqual(payload["package_schema"], "tpd-scientific-delivery/3")
        self.assertEqual(payload["seed_schema"], "scientific-seed/3")
        self.assertEqual(payload["source_followup"], "cases/expert_opinions/20261001")
        self.assertTrue(payload["source_followup_rules_installed"])
        self.assertFalse(payload["scientific_accepted"])
        self.assertFalse(payload["scientific_approval"])
        self.assertFalse(payload["automatic_approval"])
        self.assertFalse(payload["all_pass_claimed"])
        self.assertFalse(payload["secrets_accounts_or_gpu_weights_included"])
        self.assertEqual(
            payload["reports"]["expert_packet"],
            "evidence/final-reports/expert-packet/index.html",
        )

    def test_readme_marks_demo_incomplete_and_has_no_hardcoded_pass_counts(self):
        text = final_expert._readme().decode("utf-8")
        self.assertIn("cases/expert_opinions/20261001", text)
        self.assertIn("DEMO", text)
        self.assertIn("INCOMPLETE", text)
        self.assertIn("과학적 전체 통과", text)
        self.assertIn("통과·실패 건수를 README에 고정하지 않습니다", text)
        self.assertIn("사람의 과학적 판정과 최종 승인은 완료되지 않았", text)
        self.assertNotIn("신규 GPU 계산 30건", text)
        self.assertNotIn("전체 통과했습니다", text)
        self.assertNotIn("과학적으로 승인되었습니다", text)

    def test_resume_requires_matching_v5_metadata_before_delegation(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "v5"
            stage = output / "TPD"
            stage.mkdir(parents=True)
            wrong = final_expert._release_payload()
            wrong["release_version"] = "v4"
            (stage / final_expert.RELEASE_METADATA_NAME).write_text(
                json.dumps(wrong, ensure_ascii=False), encoding="utf-8"
            )
            delegate = mock.Mock(side_effect=AssertionError("must not run"))
            with mock.patch.object(scientific, "stage_package", delegate):
                with self.assertRaisesRegex(
                    RuntimeError, "V5_RELEASE_METADATA_MISMATCH"
                ):
                    final_expert.stage_package(output, True)
            delegate.assert_not_called()

    def test_resume_verifies_manifest_and_current_bound_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "v5"
            stage = output / "TPD"
            stage.mkdir(parents=True)
            (stage / final_expert.RELEASE_METADATA_NAME).write_text(
                json.dumps(final_expert._release_payload(), ensure_ascii=False),
                encoding="utf-8",
            )
            manifest = {"format": final_expert.FORMAT}
            with mock.patch.object(
                scientific, "verify_manifest", return_value=manifest
            ) as verify, mock.patch.object(
                scientific, "_verify_source"
            ) as verify_source, mock.patch.object(
                scientific, "stage_package", return_value=stage
            ) as delegate:
                result = final_expert.stage_package(output, True)

            self.assertEqual(result, stage)
            verify.assert_called_once_with(stage.resolve())
            verify_source.assert_called_once_with(stage.resolve(), manifest, current=True)
            delegate.assert_called_once_with(output, True)

    def test_seed_and_archive_use_custom_hardened_paths(self):
        original_copy_usage = scientific._copy_usage
        original_copy_reports = scientific._copy_report_dir
        original_export = scientific.export_scientific_seed
        original_archive = scientific.archive_package

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(
                final_expert, "_validate_v5_stage", side_effect=RuntimeError("custom-seed-path")
            ), mock.patch.object(scientific, "export_scientific_seed") as export:
                with self.assertRaisesRegex(RuntimeError, "custom-seed-path"):
                    final_expert.export_scientific_seed(root, root / "data", "JOB-1")
                export.assert_not_called()
            with mock.patch.object(
                final_expert, "_validate_v5_stage", side_effect=RuntimeError("custom-archive-path")
            ), mock.patch.object(scientific, "archive_package") as archive:
                with self.assertRaisesRegex(RuntimeError, "custom-archive-path"):
                    final_expert.archive_package(root, root / "reports", root / "usage.json")
                archive.assert_not_called()

        self.assertIs(scientific._copy_usage, original_copy_usage)
        self.assertIs(scientific._copy_report_dir, original_copy_reports)
        self.assertIs(scientific.export_scientific_seed, original_export)
        self.assertIs(scientific.archive_package, original_archive)


if __name__ == "__main__":
    unittest.main()

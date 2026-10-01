import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_sprint_report.py"


class SprintReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.assessment = self.root / "assessment.json"
        criteria = []
        for i in range(14):
            criteria.append({"id": f"gate-{i}", "title": f"합성 기준 {i}",
                             "status": "failed" if i == 0 else "pass",
                             "reason": "<script>alert(1)</script> https://evil.example/x" if i == 0 else "확인",
                             "observed": {"value": i}, "required": "전문가 데이터 필요"})
        self.assessment.write_text(json.dumps({"id": "synthetic-label-only", "scientific_accepted": False,
                                               "criteria": criteria}), encoding="utf-8")
        self.usage = self.root / "usage.json"
        self.usage.write_text(json.dumps({
            "developer_calls": {"calls": 2, "known_usage_calls": 2, "unknown_usage_calls": 0,
                                "input_tokens": 30, "output_tokens": 10, "total_tokens": 40},
            "product_calls": {"calls": 1, "known_usage_calls": 0, "unknown_usage_calls": 1,
                              "input_tokens": 0, "output_tokens": 0, "total_tokens": 0},
            "unknown_usage_failures": 1,
            "nested_taint": {"total_tokens": 999999},
            "note": "https://bad.example"
        }), encoding="utf-8")
        self.log = self.root / "tests.log"
        self.log.write_text("diagnostic output only\n", encoding="utf-8")
        self.dirs = {}
        for name in ("campaign", "dashboard", "parent", "science", "prior", "demo"):
            p = self.root / name
            p.mkdir()
            self.dirs[name] = p
        (self.dirs["campaign"] / "rounds").mkdir()
        (self.dirs["campaign"] / "rounds" / "result.json").write_text(
            json.dumps({"rounds_completed": 2, "ready_for_experiment": True,
                        "nested_taint": {"calls": 900, "tokens": 800}}), encoding="utf-8")
        for round_number in (1, 2):
            round_dir = self.dirs["campaign"] / "rounds" / f"round-{round_number}"
            round_dir.mkdir()
            for role in ("proposer", "adversarial_critic", "experiment_judge"):
                (round_dir / f"{role}.response.metadata.json").write_text(
                    json.dumps({"request_id": f"{round_number}-{role}"}), encoding="utf-8")
        (self.dirs["campaign"] / "token-report.json").write_text(json.dumps({
            "observed_usage": {"input_tokens": 30, "output_tokens": 10, "total_tokens": 40},
            "unknown_usage_calls": 1,
            "nested_taint": {"total_tokens": 777777}
        }), encoding="utf-8")
        evidence = {"parents": [{"parent_id": "P1", "counts_recomputed_from_records": {
                                  "docking_attempts": 234, "docking_completed": 231,
                                  "failed_docking": 3, "selected_analogs": 27,
                                  "assemblies": 324}}],
                    "campaign_counts_recomputed_from_records": {"exploratory_assemblies": 1,
                                                                  "families": 1},
                    "criteria": [{"id": "gate-0", "status": "pass"}],
                    "scientific_accepted": True}
        (self.dirs["campaign"] / "evidence.json").write_text(json.dumps(evidence), encoding="utf-8")
        (self.dirs["dashboard"] / "index.html").write_text("dashboard", encoding="utf-8")
        (self.dirs["parent"] / "manifest.json").write_text("{}", encoding="utf-8")
        (self.dirs["science"] / "result.html").write_text("science", encoding="utf-8")
        for sub in ("expert-packet", "analog-states"):
            d = self.dirs["prior"] / sub
            d.mkdir()
            (d / "index.html").write_text("old", encoding="utf-8")
        (self.dirs["demo"] / "TPD_Navigator_시연.mp4").write_bytes(b"synthetic-video")
        (self.dirs["demo"] / "caption.srt").write_text("synthetic", encoding="utf-8")
        (self.dirs["demo"] / "narration.txt").write_text("synthetic", encoding="utf-8")
        (self.dirs["demo"] / "manifest.json").write_text("{}", encoding="utf-8")
        (self.dirs["prior"] / "report-manifest.json").write_text("{}", encoding="utf-8")
        self.output = self.root / "out"

    def tearDown(self):
        self.tmp.cleanup()

    def command(self, output=None):
        return [sys.executable, str(SCRIPT), "--assessment", str(self.assessment),
                "--campaign-dir", str(self.dirs["campaign"]), "--dashboard-dir", str(self.dirs["dashboard"]),
                "--parent-pack", str(self.dirs["parent"]), "--science-pack", str(self.dirs["science"]),
                "--prior-report-dir", str(self.dirs["prior"]), "--demo-dir", str(self.dirs["demo"]),
                "--usage-json", str(self.usage), "--test-log", str(self.log),
                "--output-dir", str(output or self.output)]

    def build(self):
        result = subprocess.run(self.command(), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads((self.output / "facts.json").read_text(encoding="utf-8"))

    def test_malicious_html_is_escaped_and_external_urls_not_injected(self):
        self.build()
        report = (self.output / "report.html").read_text(encoding="utf-8")
        self.assertNotIn("<script>alert(1)</script>", report)
        self.assertIn("&lt;script&gt;alert(1)&lt;/script&gt;", report)
        self.assertNotIn("https://evil.example", report)
        self.assertNotIn("https://bad.example", report)

    def test_unknown_test_log_does_not_pass(self):
        facts = self.build()
        self.assertEqual(facts["test_result"]["status"], "unverified")
        self.assertFalse(facts["artifact_status"]["software_tests_pass"])

    def test_diagnostic_cannot_override_original_criterion(self):
        facts = self.build()
        self.assertEqual(facts["criteria"][0]["status"], "failed")
        self.assertFalse(facts["artifact_status"]["scientific_all_pass"])

    def test_actual_parent_record_counts_are_used(self):
        facts = self.build()
        self.assertEqual(facts["parents"], [{
            "parent": "P1", "attempt": 234, "completed": 231,
            "failure": 3, "selected": 27, "assembled": 324
        }])

    def test_nested_approval_taint_is_never_accepted(self):
        assessment = json.loads(self.assessment.read_text(encoding="utf-8"))
        assessment["nested"] = {"human_approval": True, "formal_approval": True}
        self.assessment.write_text(json.dumps(assessment), encoding="utf-8")
        facts = self.build()
        self.assertEqual(facts["artifact_status"]["human_approval"], "unverified")

    def test_formal_top_decision_is_the_only_human_approval_source(self):
        assessment = json.loads(self.assessment.read_text(encoding="utf-8"))
        assessment["decision"] = {"formal_confirmation": True, "human_approval": True}
        self.assessment.write_text(json.dumps(assessment), encoding="utf-8")
        facts = self.build()
        self.assertIs(facts["artifact_status"]["human_approval"], True)

    def test_explicit_api_ledgers_and_token_report_are_used(self):
        facts = self.build()
        self.assertEqual(facts["api_usage"]["developer_calls"]["total_tokens"], 40)
        self.assertEqual(facts["api_usage"]["product_calls"]["unknown_usage_calls"], 1)
        self.assertEqual(facts["api_usage"]["unknown_usage_failures"], 1)
        self.assertEqual(facts["ai_campaign"]["actual_call_count"], 6)
        self.assertEqual(facts["ai_campaign"]["observed_usage"]["total_tokens"], 40)
        self.assertTrue((self.output / "raw" / "token-report.json").is_file())

    def test_skipped_tests_are_not_reported_as_passed(self):
        self.log.write_text("Ran 5 tests in 0.01s\n\nOK (skipped=2)\n", encoding="utf-8")
        facts = self.build()
        self.assertEqual(facts["test_result"]["tests_run"], 5)
        self.assertEqual(facts["test_result"]["skipped"], 2)
        self.assertEqual(facts["test_result"]["passed"], 3)
        self.assertEqual(facts["test_result"]["failed"], 0)

    def test_failed_test_count_is_reported(self):
        self.log.write_text("Ran 5 tests in 0.01s\n\nFAILED (failures=1, errors=1, skipped=1)\n",
                            encoding="utf-8")
        facts = self.build()
        self.assertEqual(facts["test_result"]["status"], "failed")
        self.assertEqual(facts["test_result"]["tests_run"], 5)
        self.assertEqual(facts["test_result"]["skipped"], 1)
        self.assertEqual(facts["test_result"]["failed"], 2)
        self.assertIsNone(facts["test_result"]["passed"])

    def test_launcher_steps_and_query_buttons_match_actual_names(self):
        self.build()
        report = (self.output / "report.html").read_text(encoding="utf-8")
        for expected in (
            "04: 파일 해시 확인", "01: 프로그램 시작", "10: 연구 비교",
            "02: CPU 실험", "07: 이전 C01/C02 확인", "03: 종료",
            "CRBN·VHL 후보와 linker 비교", "도킹 실패와 재실험 항목 확인",
            "14개 기준과 AI개선제안 확인"
        ):
            self.assertIn(expected, report)
        self.assertIn("table-layout:fixed", report)
        self.assertIn("overflow-x:auto", report)

    def test_submission_strings_are_never_inferred(self):
        assessment = json.loads(self.assessment.read_text(encoding="utf-8"))
        assessment["nested"] = {"public_url": "https://example.invalid/report",
                                "youtube_url": "https://example.invalid/video",
                                "slides_pdf": "slides.pdf"}
        self.assessment.write_text(json.dumps(assessment), encoding="utf-8")
        facts = self.build()
        self.assertEqual(set(facts["submission_items_provided"].values()), {"unverified/미제공"})

    def test_output_overwrite_is_blocked(self):
        self.build()
        second = subprocess.run(self.command(), text=True, capture_output=True)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("FileExistsError", second.stderr)

    def test_tree_copy_hash_and_manifest(self):
        self.build()
        copied = self.output / "demo" / "TPD_Navigator_시연.mp4"
        self.assertEqual(copied.read_bytes(), b"synthetic-video")
        manifest = json.loads((self.output / "report-manifest.json").read_text(encoding="utf-8"))
        paths = {row["path"] for row in manifest["outputs"]}
        self.assertIn("demo/TPD_Navigator_시연.mp4", paths)
        self.assertNotIn("report-manifest.json", paths)
        self.assertIn("previous-v5/report-manifest.json", paths)
        self.assertTrue((self.output / "expert-packet" / "index.html").is_file())


if __name__ == "__main__":
    unittest.main()

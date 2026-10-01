import hashlib
import importlib.util
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts/render_final_delivery_report.py"
SPEC = importlib.util.spec_from_file_location("renderer", SCRIPT)
renderer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(renderer)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FinalDeliveryReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.evidence = self.root / "evidence"
        self.evidence.mkdir()

        expert = self.evidence / "expert-packet"
        expert.mkdir()
        expert_index = expert / "index.html"
        expert_index.write_text("<h1>expert-packet</h1>", encoding="utf-8")
        self.write_manifest(expert, [expert_index])

        analog = self.evidence / "analog-states"
        analog.mkdir()
        analog_index = analog / "index.html"
        analog_index.write_text("<h1>analog-states</h1>", encoding="utf-8")
        sources = analog / "sources"
        sources.mkdir()
        self.source = sources / "reply.docx"
        self.source.write_bytes(b"minimal reply docx fixture")
        self.write_manifest(analog, [analog_index, self.source])

        raw = self.evidence / "raw"
        raw.mkdir()
        (raw / "policy.json").write_text("{}", encoding="utf-8")
        (raw / "api-usage.json").write_text("{}", encoding="utf-8")
        (raw / "runtime-strict.json").write_text("{}", encoding="utf-8")
        (raw / "tests.log").write_text(
            "Ran 8 tests\nFAILED (failures=1, skipped=2)\n",
            encoding="utf-8",
        )
        receipts = raw / "parent-receipts"
        receipts.mkdir()
        (receipts / "p1.json").write_text("{}", encoding="utf-8")

        criteria = []
        statuses = ["pass", "failed", "pending"]
        for index, criterion_id in enumerate(renderer.IDS):
            criteria.append({
                "id": criterion_id,
                "title": f"title {index}",
                "status": statuses[index % len(statuses)],
                "reason": "reason",
                "observed_text": (
                    "<script>alert(1)</script>" if index == 0 else "observed"
                ),
                "required_text": "required",
            })

        self.facts = {
            "format": "tpd-final-delivery-facts/1",
            "scientific_accepted": False,
            "criteria": criteria,
            "tests": {
                "command": "python -m unittest",
                "tests_run": 8,
                "failures": 1,
                "errors": 0,
                "skipped": 2,
                "exit_code": 1,
                "log_sha256": digest(raw / "tests.log"),
            },
            "software_features": ["reporting", "<escaping>"],
            "measurements": [{
                "title": "strict current live",
                "headers": ["kind", "count"],
                "rows": [["dock", 8], ["verified facts", 36]],
            }],
            "api_usage": {
                "development_known_tokens": 10,
                "product_known_tokens": 20,
                "total_known_tokens": 30,
                "unknown_usage_calls": 1,
                "latest_quota_estimate": 100,
                "quota_note": "estimate <not balance>",
            },
            "freshness": {
                "current_policy": True,
                "current_implementation": True,
                "source_runtime_current": False,
            },
            "source": {
                "sha256": digest(self.source),
                "bytes": self.source.stat().st_size,
                "path": renderer.SOURCE_PATH,
            },
            "delivery_status": "DEMO/INCOMPLETE",
        }

        (raw / "test-summary.json").write_text(
            json.dumps(self.facts["tests"], ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        self.write_assessment()
        self.facts_path = self.root / "input-facts.json"
        self.write_facts()

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self, base, files):
        manifest = {
            "files": [{
                "path": path.relative_to(base).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            } for path in files]
        }
        (base / "manifest.json").write_text(
            json.dumps(manifest),
            encoding="utf-8",
        )

    def write_assessment(self):
        assessment_source = {
            "sha256": self.facts["source"]["sha256"],
            "bytes": self.facts["source"]["bytes"],
            "path": "cases/expert_opinions/20261001/reply.docx",
            "name": "TPD_Navigator_expert_followup_검토회신안_20261001.docx",
            "verified_extracts": {"/body/13": "a" * 64},
        }
        assessment = {
            "format": "scientific-acceptance/20260930.3",
            "scientific_accepted": False,
            "expert_questions": [{
                "id": "parent-selection",
                "parent_id": "SMARCA2-FX5",
                "question": "실제 알려진 부모 ID 중 어떤 후보가 선택됩니까?",
            }, {
                "id": "site-counts",
                "parent_id": "SMARCA2-FX5",
                "atom_maps": [19],
                "question": "MODIFIABLE map [19]가 기준을 충족합니까?",
            }],
            "pending_actions": [{
                "criterion_id": "parent_funnel",
                "action": "Only server-trusted parent funnel rows count.",
            }, {
                "criterion_id": "formal_expert_decision",
                "action": "Formal approval remains pending.",
            }],
            "source_runtime_current": False,
            "current_policy": True,
            "current_implementation": True,
            "criteria": [
                {
                    **item,
                    "observed": {"fixture": True},
                    "additional_assessment_key": f"preserved-{index}",
                }
                for index, item in enumerate(self.facts["criteria"])
            ],
            "expert_followup": {
                "source": assessment_source,
                "additional_followup_key": "preserved",
            },
            "additional_top_level_key": {"preserved": True},
        }
        (self.evidence / "raw/assessment.json").write_text(
            json.dumps(assessment, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    def write_facts(self):
        self.facts_path.write_text(
            json.dumps(self.facts, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    def run_renderer(self, name="out"):
        output = self.root / name
        renderer.main([
            "--facts", str(self.facts_path),
            "--evidence-dir", str(self.evidence),
            "--output", str(output),
        ])
        return output

    def test_renders_copies_escaped_features_api_measurements_and_questions(self):
        original = self.facts_path.read_bytes()
        original_assessment = (self.evidence / "raw/assessment.json").read_bytes()
        output = self.run_renderer()
        self.assertEqual(original, (output / "facts.json").read_bytes())
        self.assertEqual(
            original_assessment,
            (output / "raw/assessment.json").read_bytes(),
        )
        self.assertEqual(
            self.source.read_bytes(),
            (output / renderer.SOURCE_PATH).read_bytes(),
        )
        report = (output / "report.html").read_text(encoding="utf-8")
        report_md = (output / "report.md").read_text(encoding="utf-8")
        self.assertIn("&lt;script&gt;alert(1)&lt;/script&gt;", report)
        self.assertNotIn("<script>alert(1)</script>", report)
        self.assertIn("&lt;escaping&gt;", report)
        self.assertIn("estimate &lt;not balance&gt;", report)
        self.assertIn("실패 1건, 오류 0건", report)
        self.assertIn("TEST SKIPPED", report)
        self.assertNotIn("<button", report)
        self.assertIn("01로 애플리케이션을 실행", report)
        self.assertIn("실제 알려진 부모 ID 중 어떤 후보가 선택됩니까?", report)
        self.assertIn("MODIFIABLE map [19]가 기준을 충족합니까?", report)
        self.assertIn("Only server-trusted parent funnel rows count.", report)
        self.assertIn("source runtime: 현재 기준과 불일치", report)
        self.assertIn("Formal approval remains pending.", report)
        self.assertIn("## API 사용량", report_md)
        self.assertIn("## 측정값", report_md)
        self.assertIn("## 전문가 후속 질문 및 조치", report_md)
        self.assertIn("실제 알려진 부모 ID 중 어떤 후보가 선택됩니까?", report_md)
        self.assertIn("Formal approval remains pending.", report_md)
        self.assertIn("verified facts", report_md)
        self.assertIn("&lt;escaping&gt;", report_md)
        manifest = json.loads(
            (output / "report-manifest.json").read_text(encoding="utf-8")
        )
        self.assertIn("report.html", manifest["files"])
        self.assertIn(renderer.SOURCE_PATH, manifest["files"])
        self.assertNotIn("report-manifest.json", manifest["files"])

    def test_three_scientific_statuses_are_exact_and_skipped_is_test_only(self):
        self.assertEqual(renderer.STATUSES, {"pass", "failed", "pending"})
        counts = renderer.status_counts(self.facts["criteria"])
        self.assertEqual(set(counts), {"pass", "failed", "pending"})
        for invalid in ("fail", "unknown", "skipped", "PASS"):
            self.facts["criteria"][0]["status"] = invalid
            with self.assertRaisesRegex(ValueError, "invalid criterion status"):
                renderer.validate_facts(self.facts, self.evidence)
        self.assertEqual(self.facts["tests"]["skipped"], 2)

    def test_tampered_manifest_artifact_is_refused_upfront(self):
        (self.evidence / "expert-packet/index.html").write_text(
            "tampered", encoding="utf-8"
        )
        output = self.root / "out"
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            renderer.main([
                "--facts", str(self.facts_path),
                "--evidence-dir", str(self.evidence),
                "--output", str(output),
            ])
        self.assertFalse(output.exists())

    def test_manifest_path_escape_and_windows_colon_are_refused(self):
        manifest = self.evidence / "analog-states/manifest.json"
        data = json.loads(manifest.read_text(encoding="utf-8"))
        data["files"][0]["path"] = "../index.html"
        manifest.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unsafe path"):
            renderer.validate_evidence(self.evidence)
        self.assertRaisesRegex(ValueError, "unsafe path", renderer.safe_rel, "C:/x")

    def test_bad_fourteen_ids_are_refused(self):
        self.facts["criteria"][3]["id"] = "invented"
        self.write_facts()
        with self.assertRaisesRegex(ValueError, "IDs or ordering"):
            renderer.validate_facts(self.facts, self.evidence)

    def test_real_assessment_shape_is_accepted_with_original_notes(self):
        assessment = renderer.validate_facts(self.facts, self.evidence)
        self.assertEqual(
            assessment["format"],
            "scientific-acceptance/20260930.3",
        )
        self.assertEqual(
            renderer.original_notes(assessment),
            [
                "실제 알려진 부모 ID 중 어떤 후보가 선택됩니까?",
                "MODIFIABLE map [19]가 기준을 충족합니까?",
                "Only server-trusted parent funnel rows count.",
                "Formal approval remains pending.",
            ],
        )

    def test_assessment_status_tamper_is_refused(self):
        assessment_path = self.evidence / "raw/assessment.json"
        assessment = json.loads(assessment_path.read_text(encoding="utf-8"))
        assessment["criteria"][1]["status"] = "pass"
        assessment_path.write_text(json.dumps(assessment), encoding="utf-8")
        with self.assertRaisesRegex(
            ValueError, "assessment criteria do not exactly match facts"
        ):
            renderer.validate_facts(self.facts, self.evidence)

    def test_assessment_source_tamper_is_refused(self):
        assessment_path = self.evidence / "raw/assessment.json"
        assessment = json.loads(assessment_path.read_text(encoding="utf-8"))
        assessment["expert_followup"]["source"]["sha256"] = "0" * 64
        assessment_path.write_text(json.dumps(assessment), encoding="utf-8")
        with self.assertRaisesRegex(
            ValueError, "assessment source hash or bytes do not exactly match"
        ):
            renderer.validate_facts(self.facts, self.evidence)

    def test_source_and_test_bindings_are_required_before_output_creation(self):
        self.source.write_bytes(b"tampered source")
        output = self.root / "source-out"
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            renderer.main([
                "--facts", str(self.facts_path),
                "--evidence-dir", str(self.evidence),
                "--output", str(output),
            ])
        self.assertFalse(output.exists())

        self.setUp_after_source_tamper_not_supported()

    def setUp_after_source_tamper_not_supported(self):
        pass

    def test_test_log_and_summary_tampering_are_refused(self):
        log = self.evidence / "raw/tests.log"
        original_log = log.read_text(encoding="utf-8")
        log.write_text(original_log + "tampered\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "tests.log checksum mismatch"):
            renderer.validate_facts(self.facts, self.evidence)
        log.write_text(original_log, encoding="utf-8")

        summary = self.evidence / "raw/test-summary.json"
        value = json.loads(summary.read_text(encoding="utf-8"))
        value["skipped"] = 0
        summary.write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "does not exactly match"):
            renderer.validate_facts(self.facts, self.evidence)

    def test_invalid_nonfinite_api_and_known_sum_are_refused(self):
        self.facts["api_usage"]["latest_quota_estimate"] = math.inf
        with self.assertRaisesRegex(ValueError, "non-finite"):
            renderer.validate_facts(self.facts, self.evidence)
        self.facts["api_usage"]["latest_quota_estimate"] = 100
        self.facts["api_usage"]["total_known_tokens"] = 31
        with self.assertRaisesRegex(ValueError, "known-token sum"):
            renderer.validate_facts(self.facts, self.evidence)

    def test_all_pass_inconsistency_and_existing_output_are_refused(self):
        for criterion in self.facts["criteria"]:
            criterion["status"] = "pass"
        self.write_assessment()
        with self.assertRaisesRegex(ValueError, "all-pass"):
            renderer.validate_facts(self.facts, self.evidence)

        self.facts["criteria"][-1]["status"] = "pending"
        self.write_assessment()
        self.write_facts()
        output = self.root / "exists"
        output.mkdir()
        with self.assertRaisesRegex(ValueError, "new and nonexistent"):
            renderer.main([
                "--facts", str(self.facts_path),
                "--evidence-dir", str(self.evidence),
                "--output", str(output),
            ])

    @unittest.skipUnless(hasattr(os, "symlink"), "symlink unsupported")
    def test_symlink_evidence_root_is_refused(self):
        linked = self.root / "linked-evidence"
        try:
            linked.symlink_to(self.evidence, target_is_directory=True)
        except OSError:
            self.skipTest("symlink creation unavailable")
        with self.assertRaisesRegex(ValueError, "symlink or junction root"):
            renderer.validate_evidence(linked)


if __name__ == "__main__":
    unittest.main()

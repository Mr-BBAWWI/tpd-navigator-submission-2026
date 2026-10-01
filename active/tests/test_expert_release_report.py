import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_expert_release_report.py"
SPEC = importlib.util.spec_from_file_location("build_expert_release_report", MODULE_PATH)
reporter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(reporter)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def criteria_fixture():
    statuses = ["pass"] * 3 + ["failed"] * 3 + ["pending"] * 8
    return [
        {
            "id": f"criterion_{index:02d}",
            "title": f"기준 {index}",
            "status": status,
            "required": True,
            "reason": f"source reason {index}",
        }
        for index, status in enumerate(statuses, 1)
    ]


def policy_fixture():
    return {
        "numerical_criteria": {
            "parent_funnel_min": 5,
            "parent_funnel_max": 10,
            "modifiable_sites_min": 1,
            "distinct_graphs_per_site_min": 30,
            "broad_families_min": 6,
            "panel_min": 10,
            "panel_max": 20,
            "calibration_seeds_min": 5,
            "novel_ternary_repeats_min": 5,
        }
    }




def assessment_fixture():
    return {
        "id": "assessment-dba30f0a65004dd28a79472d94551c97",
        "job_id": "job-4ce0b3d57ecb4be5b795504ca22762a0",
        "project_id": "local-research",
        "revision": 3,
        "source_binding": {
            "input_sha256": "9324d846a8a7a563f7c4ebf7f113880ab66cf9e474838531cda703cb3ed5b4c1",
            "job_id": "job-4ce0b3d57ecb4be5b795504ca22762a0",
            "project": "local-research",
            "result_binding_kind": "exact_archived_design_json_bytes",
            "result_sha256": "51c85df73558b5ce87c3999c96fca937d4bab272e8e4b5fb36509bb3f2742369",
        },
        "confirmation_status": "unconfirmed",
        "scientific_accepted": False,
        "criteria": criteria_fixture(),
        "counts": {"pass": 3, "failed": 3, "pending": 8},
        "current_policy": policy_fixture(),
    }


def benchmark_fixture():
    seeds = []
    values_e3 = [5.2205889006556205, 39.19674979017389, 3.220835929717891, 27.807100353950336, 10.621024449350978]
    values_ligand = [1.5420394953522754, 12.38776540052114, 1.1780957439827384, 7.4289082360252525, 2.505346661961516]
    for index, seed in enumerate(reporter.EXPECTED_SEEDS):
        seeds.append({
            "seed": seed,
            "original_status": "failure" if seed == 23 else "success",
            "original_failure_reason": "PREDICTION_CHAIN_MAPPING_AMBIGUOUS" if seed == 23 else None,
            "metrics": {
                "comparison.metrics.e3_CA_RMSD_after_target_alignment_A": values_e3[index],
                "comparison.metrics.ligand_heavy_atom_RMSD_after_target_alignment_A": values_ligand[index],
            },
        })
    return {
        "seeds": seeds,
        "technical_success_count": 5,
        "scientific_model_quality_success_count": None,
        "expert_geometric_success_cutoff": None,
        "distribution_strength": "weak_single_known_case",
        "outlier_deletion_performed": False,
        "statistics": {
            "comparison.metrics.e3_CA_RMSD_after_target_alignment_A": {
                "count": 5, "median": 10.621024449350978, "min": 3.220835929717891,
                "max": 39.19674979017389, "q1": 5.2205889006556205,
                "q3": 27.807100353950336, "iqr": 22.586511453294715,
            },
            "comparison.metrics.ligand_heavy_atom_RMSD_after_target_alignment_A": {
                "count": 5, "median": 2.505346661961516, "min": 1.1780957439827384,
                "max": 12.38776540052114, "q1": 1.5420394953522754,
                "q3": 7.4289082360252525, "iqr": 5.886868740672977,
            },
        },
    }


def packet_fixture(root):
    root.mkdir()
    index = root / "index.html"
    index.write_text("<!doctype html><title>packet</title>", encoding="utf-8")
    representatives = []
    number = 0
    for analog in reporter.EXPECTED_ANALOGS:
        for e3 in ("CRBN", "VHL"):
            number += 1
            representatives.append({
                "candidate_id": f"D-{number:02d}",
                "warhead_analog_id": analog,
                "e3_type": e3,
                "linker_id": "alkyl_c6",
                "candidate_graph_sha256": hashlib.sha256(str(number).encode()).hexdigest(),
            })
    manifest = {
        "format": "expert-evidence-packet/test",
        "representatives": representatives,
        "artifacts": [{"path": "index.html", "sha256": digest(index), "bytes": index.stat().st_size}],
    }
    write_json(root / "manifest.json", manifest)
    return manifest


def batch_fixture(root, failed=False, aggregates=False):
    assessment = assessment_fixture()
    policy_digest = "8" * 64
    number = 0
    for analog in reporter.EXPECTED_ANALOGS:
        for e3 in ("CRBN", "VHL"):
            number += 1
            candidate_id = f"D-{number:012d}"
            graph_sha = hashlib.sha256(f"graph-{number}".encode()).hexdigest()
            plan_digest = hashlib.sha256(f"plan-{number}".encode()).hexdigest()
            run_name = f"{analog}--{e3}"
            plan = {
                "format": "tpd-novel-ternary-plan/1.0",
                "bindings": {
                    "job": {
                        "candidate_graph_sha256": graph_sha,
                        "input_sha256": assessment["source_binding"]["input_sha256"],
                        "job_id": assessment["job_id"],
                        "project": assessment["project_id"],
                        "result_sha256": assessment["source_binding"]["result_sha256"],
                    },
                    "policy": {"digest": policy_digest, "module": "test", "revision": "0"},
                },
                "candidate_graph": {
                    "candidate_id": candidate_id,
                    "warhead_analog_id": analog,
                    "e3_type": e3,
                    "graph_sha256": graph_sha,
                },
                "plan_digest": plan_digest,
                "msa_mode": "single_sequence",
                "seeds": list(reporter.EXPECTED_SEEDS),
            }
            write_json(root / f"{run_name}.plan.json", plan)
            for seed in reporter.EXPECTED_SEEDS:
                is_failure = failed and number == 1 and seed == 41
                receipt = {
                    "actual_computation": True,
                    "bindings": {
                        "job_id": assessment["job_id"],
                        "project": assessment["project_id"],
                        "policy_digest": policy_digest,
                    },
                    "candidate_id": candidate_id,
                    "e3_type": e3,
                    "exit_code": 9 if is_failure else 0,
                    "inspection": {
                        "actual_computation": True,
                        "reference_free": True,
                        "status": "computed_hypothesis",
                    },
                    "plan_digest": plan_digest,
                    "reference_free": True,
                    "scientific_status": "computed_hypothesis",
                    "seed": seed,
                    "status": "failed" if is_failure else "completed",
                }
                write_json(root / run_name / f"seed-{seed}" / "receipt.json", receipt)
            if aggregates:
                write_json(root / run_name / "receipt.json", {"run": run_name, "receipts": 5})
    if aggregates:
        write_json(root / "batch-plan.json", {"runs": 6})


class ExpertReleaseReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.packet = self.root / "packet"
        packet_fixture(self.packet)

    def tearDown(self):
        self.temp.cleanup()

    def test_packet_manifest_accepts_exact_files(self):
        manifest, files = reporter.validate_packet(self.packet)
        self.assertEqual(manifest["format"], "expert-evidence-packet/test")
        self.assertEqual(files, ["index.html", "manifest.json"])

    def test_packet_tamper_is_rejected(self):
        (self.packet / "index.html").write_text("tampered", encoding="utf-8")
        with self.assertRaisesRegex(reporter.ReleaseError, "크기 불일치|SHA-256 불일치"):
            reporter.validate_packet(self.packet)

    def test_unmanifested_secret_file_is_rejected(self):
        (self.packet / ".env").write_text("API_KEY=not-for-release", encoding="utf-8")
        with self.assertRaisesRegex(reporter.ReleaseError, "manifest에 없는 파일"):
            reporter.validate_packet(self.packet)

    def test_status_counts_are_derived_not_fake_all_pass(self):
        assessment = assessment_fixture()
        criteria = reporter.extract_criteria(assessment)
        counts = reporter.status_counts(criteria)
        self.assertEqual(counts, {"pass": 3, "failed": 3, "pending": 8})
        assessment["counts"] = {"pass": 14, "failed": 0, "pending": 0}
        batch = self.root / "batch"
        batch_fixture(batch)
        packet_manifest, packet_files = reporter.validate_packet(self.packet)
        with self.assertRaisesRegex(reporter.ReleaseError, "criteria-derived count"):
            reporter.build_summary(
                assessment,
                benchmark_fixture(),
                packet_manifest,
                packet_files,
                batch,
                {"developer_calls": {}, "product_calls": {}, "unknown_usage_failures": 0},
                policy_fixture(),
                "assessment.current_policy",
                None,
            )

    def test_failed_actual_run_is_preserved(self):
        batch = self.root / "batch"
        batch_fixture(batch, failed=True)
        novel = reporter.summarize_novel_batches(assessment_fixture(), batch)
        self.assertEqual(novel["completed_exit_zero"], 29)
        self.assertEqual(len(novel["failed_or_nonzero_receipts"]), 1)
        failure = novel["failed_or_nonzero_receipts"][0]
        self.assertEqual(failure["seed"], 41)
        self.assertEqual(failure["exit_code"], 9)
        self.assertEqual(novel["computed_hypothesis"], "technical_not_scientific_pass")

    def test_batch_root_preserves_failed_receipt(self):
        batch = self.root / "batch"
        batch_fixture(batch, failed=True)
        result = reporter.validate_batch_root(batch, assessment_fixture())
        self.assertEqual(result["plan_count"], 6)
        self.assertEqual(result["receipt_count"], 30)
        self.assertEqual(result["exit_zero_count"], 29)
        self.assertEqual(len(result["failed_or_nonzero_receipts"]), 1)

    def test_real_shape_without_batch_files_fails(self):
        batch = self.root / "empty-batch"
        batch.mkdir()
        with self.assertRaisesRegex(reporter.ReleaseError, "root \\*\\.plan\\.json 6개"):
            reporter.validate_batch_root(batch, assessment_fixture())

    def test_only_root_plans_and_seed_receipts_count_as_runs(self):
        batch = self.root / "batch-with-aggregates"
        batch_fixture(batch, aggregates=True)
        result = reporter.validate_batch_root(batch, assessment_fixture())
        self.assertEqual(result["plan_count"], 6)
        self.assertEqual(result["receipt_count"], 30)
        self.assertEqual(result["aggregate_records"]["root_json"], ["batch-plan.json"])
        self.assertEqual(len(result["aggregate_records"]["run_receipts"]), 6)

    def test_plan_graph_result_binding_tamper_is_rejected(self):
        batch = self.root / "tampered-batch"
        batch_fixture(batch)
        plan_path = sorted(batch.glob("*.plan.json"))[0]
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        plan["bindings"]["job"]["result_sha256"] = "0" * 64
        write_json(plan_path, plan)
        with self.assertRaisesRegex(reporter.ReleaseError, "graph/result binding 불일치"):
            reporter.validate_batch_root(batch, assessment_fixture())

    def test_packet_accepts_actual_analog_id_fallback_only(self):
        manifest = json.loads((self.packet / "manifest.json").read_text(encoding="utf-8"))
        for row in manifest["representatives"]:
            row["analog_id"] = row.pop("warhead_analog_id")
        summary = reporter.summarize_packet(manifest)
        self.assertEqual(summary["representative_count"], 6)
        manifest["representatives"][0]["warhead_analog_id"] = reporter.EXPECTED_ANALOGS[1]
        with self.assertRaisesRegex(reporter.ReleaseError, "analog 필드가 충돌"):
            reporter.summarize_packet(manifest)

    def test_nested_usage_counters_and_missing_metrics_render_unavailable(self):
        totals = reporter.summarize_usage_totals({
            "developer_calls": {"calls": 167, "total_tokens": 6375167},
            "product_calls": {"calls": 24, "total_tokens": 201066},
            "unknown_usage_failures": 12,
        })
        self.assertEqual(totals["developer_total_tokens"], 6375167)
        self.assertEqual(totals["developer_calls"], 167)
        self.assertEqual(totals["product_total_tokens"], 201066)
        self.assertEqual(totals["product_calls"], 24)

        missing = reporter.summarize_usage_totals({"developer_calls": {}, "product_calls": None})
        for key in ("developer_total_tokens", "developer_calls", "product_total_tokens", "product_calls"):
            self.assertIsNone(missing[key])
        rendered = reporter.render_report_html({
            "status_counts": {"pass": 3, "failed": 3, "pending": 8},
            "validation": {"provided": False, "note": "not provided"},
            "confirmation_status": "unconfirmed",
            "novel": {"completed_exit_zero": 30},
            "usage_totals": missing,
            "assessment_id": "test",
            "assessment_revision": 1,
            "policy_source": "test",
        }, [])
        self.assertGreaterEqual(rendered.count("확인 불가(unavailable)"), 4)

    def test_full_generation_outputs_source_bound_status(self):
        assessment = self.root / "assessment.json"
        benchmark = self.root / "benchmark.json"
        usage = self.root / "usage.json"
        batch = self.root / "batch"
        output = self.root / "release"
        write_json(assessment, assessment_fixture())
        write_json(benchmark, benchmark_fixture())
        write_json(usage, {
            "developer_calls": {"calls": 167, "total_tokens": 6375167, "unknown_usage_calls": 12},
            "product_calls": {"calls": 24, "total_tokens": 201066, "unknown_usage_calls": 0},
            "unknown_usage_failures": 12,
            "source_hashes": [{"path": "C:\\Users\\private-user\\ledger.json", "sha256": "a" * 64}],
            "latest_quota": {"label": "estimated", "nonmonotonic_observations": True},
            "quota_caveat": "Estimated observations only",
        })
        batch_fixture(batch)
        args = reporter.make_parser().parse_args([
            "--assessment", str(assessment),
            "--benchmark", str(benchmark),
            "--packet", str(self.packet),
            "--batch-root", str(batch),
            "--usage", str(usage),
            "--output", str(output),
        ])
        reporter.generate(args)
        self.assertTrue((output / "report.html").is_file())
        self.assertTrue((output / "report.md").is_file())
        self.assertTrue((output / "expert-followup.html").is_file())
        self.assertTrue((output / "expert-packet" / "index.html").is_file())
        measured = json.loads((output / "measured-summary.json").read_text(encoding="utf-8"))
        self.assertEqual(measured["status_counts"], {"pass": 3, "failed": 3, "pending": 8})
        self.assertFalse(measured["scientific_accepted"])
        self.assertEqual(measured["novel"]["completed_exit_zero"], 30)
        self.assertEqual(measured["usage_totals"]["developer_total_tokens"], 6375167)
        self.assertEqual(measured["usage_totals"]["developer_calls"], 167)
        self.assertEqual(measured["usage_totals"]["product_total_tokens"], 201066)
        self.assertEqual(measured["usage_totals"]["product_calls"], 24)
        report_html = (output / "report.html").read_text(encoding="utf-8")
        self.assertIn("6375167", report_html)
        self.assertIn("201066", report_html)
        followup = (output / "expert-followup.html").read_text(encoding="utf-8")
        self.assertIn("CURRENT STATUS", followup)
        self.assertIn("FAILED", followup)
        self.assertEqual(followup.count("<tr>"), 15)  # header 1 + criteria 14
        usage_copy = (output / "api-usage.json").read_text(encoding="utf-8")
        self.assertNotIn("private-user", usage_copy)
        self.assertIn("<private-path-redacted>", usage_copy)

    def test_json_nan_is_rejected(self):
        bad = self.root / "bad.json"
        bad.write_text('{"value": NaN}', encoding="utf-8")
        with self.assertRaises(reporter.ReleaseError):
            reporter.load_json(bad, "bad")

    def test_symlink_packet_file_rejected_when_supported(self):
        target = self.root / "outside.txt"
        target.write_text("outside", encoding="utf-8")
        link = self.packet / "linked.txt"
        try:
            link.symlink_to(target)
        except (OSError, NotImplementedError):
            self.skipTest("symlink creation is unavailable")
        with self.assertRaisesRegex(reporter.ReleaseError, "symlink|reparse"):
            reporter.validate_packet(self.packet)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_campaign_diagnostics", ROOT / "scripts" / "build_campaign_diagnostics.py"
)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


class CampaignDiagnosticsUnitTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, relative: str, data: bytes = b"evidence") -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def test_rejects_path_escape(self):
        outside = Path(self.temporary.name).parent / "campaign-outside.json"
        outside.write_text("{}", encoding="utf-8")
        try:
            with self.assertRaises(module.BuildError):
                module.safe_source(outside, self.root)
        finally:
            outside.unlink(missing_ok=True)

    def test_rejects_symlink_when_supported(self):
        source = self.write("real/report.json", b"{}")
        link = self.root / "linked.json"
        try:
            link.symlink_to(source)
        except (OSError, NotImplementedError):
            self.skipTest("symlinks unavailable")
        with self.assertRaises(module.BuildError):
            module.safe_source(link, self.root)

    def test_copy_hash_integrity_and_expected_hash(self):
        source = self.write("source/report.json", b'{"ok":true}\n')
        stage = self.root / "stage"
        stage.mkdir()
        expected = hashlib.sha256(source.read_bytes()).hexdigest()
        record = module.copy_verified(source, stage, "evidence/report.json", self.root, expected)
        copied = stage / "evidence/report.json"
        self.assertEqual(copied.read_bytes(), source.read_bytes())
        self.assertEqual(record["source_sha256"], expected)
        with self.assertRaises(module.BuildError):
            module.copy_verified(source, stage, "evidence/other.json", self.root, "0" * 64)

    def synthetic_core(self):
        rows = []
        for index in range(40):
            rows.append({
                "analog_id": "A" if index < 20 else "B",
                "pose_index": index % 5,
                "required_contact_evidence": index < 16,
                "failure": None if index != 39 else "synthetic failure",
            })
        pka = []
        for log in range(2):
            for index in range(44):
                pka.append({
                    "chain": "A",
                    "residue": "RES" + str(index),
                    "sequence_number": index,
                    "pKa": float(index) / 10,
                    "log_source": log,
                })
        return {"rows": rows, "protein_logs": pka, "analogs": []}

    def test_core_preserves_false_null_failure_and_pka_denominators(self):
        compact = module.build_core_compact(self.synthetic_core())
        self.assertEqual(compact["required_contact_true_count"], 16)
        self.assertEqual(compact["required_contact_false_count"], 24)
        self.assertEqual(compact["failure_count"], 1)
        self.assertEqual(compact["failure_null_count"], 39)
        self.assertFalse(compact["ligand_pKa_prediction_performed"])
        self.assertFalse(compact["population_prediction_performed"])
        self.assertEqual(compact["protein_pKa_raw_row_count"], 88)
        self.assertEqual(compact["protein_pKa_unique_chain_residue_sequence_pKa_count"], 44)

    def test_pka_88_rows_are_not_accepted_as_88_unique(self):
        source = self.synthetic_core()
        for index, row in enumerate(source["protein_logs"]):
            row["sequence_number"] = index
        with self.assertRaisesRegex(module.BuildError, "UNIQUE_COUNT_NOT_44"):
            module.build_core_compact(source)

    def test_calibration_raw_and_selected_denominators_are_distinct(self):
        receipts = []
        raw_passes = 0
        for seed_number, seed in enumerate(module.SEEDS):
            models = []
            for index in range(5):
                passed = raw_passes < 16
                raw_passes += int(passed)
                models.append({
                    "model_index": index,
                    "metrics": {"pass": passed,
                                "e3_CA_RMSD_after_target_alignment_A": 3.0 + index},
                })
            selected_index = 0 if seed_number < 2 else 4
            if seed_number >= 2:
                models[selected_index]["metrics"]["pass"] = False
            receipts.append({"seed": seed, "selection": {"model_index": selected_index},
                             "models": models})
        # Restore 16 raw successes after making three selected models fail.
        current = sum(model["metrics"]["pass"] for receipt in receipts for model in receipt["models"])
        for receipt in receipts:
            for model in receipt["models"]:
                if current == 16:
                    break
                if model["model_index"] != receipt["selection"]["model_index"] and not model["metrics"]["pass"]:
                    model["metrics"]["pass"] = True
                    current += 1
        plan = {"preserved_original_baseline": {"geometric_successes": 2, "seed_count": 5},
                "cutoffs": {}, "effective_protocol_settings": {}}
        compact = module.build_calibration_compact(
            plan, receipts, lambda metrics: metrics["pass"]
        )
        self.assertEqual(compact["raw_model_count"], 25)
        self.assertEqual(compact["raw_all_four_pass_count"], 16)
        self.assertEqual(compact["selected_pass_count"], 2)

    def test_compact_limit_is_strictly_under_20_kib(self):
        path = self.root / "compact.json"
        size = module.write_compact(path, {"false_value": False, "null_value": None})
        self.assertLess(size, 20 * 1024)
        loaded = json.loads(path.read_text(encoding="utf-8"))
        self.assertIs(loaded["false_value"], False)
        self.assertIsNone(loaded["null_value"])
        with self.assertRaises(module.BuildError):
            module.write_compact(self.root / "large.json", {"x": "a" * 21000})

    def test_manifest_detects_tamper(self):
        stage = self.root / "pack"
        stage.mkdir()
        evidence = stage / "evidence.json"
        evidence.write_text('{"value":1}\n', encoding="utf-8")
        module.write_json_exclusive(stage / "output-manifest.json", module.make_manifest(stage))
        module.verify_manifest(stage)
        evidence.write_text('{"value":2}\n', encoding="utf-8")
        with self.assertRaisesRegex(module.BuildError, "HASH_MISMATCH"):
            module.verify_manifest(stage)

    def test_manifest_rejects_unlisted_file(self):
        stage = self.root / "pack"
        stage.mkdir()
        self.write_manifest_file(stage)
        (stage / "extra.json").write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(module.BuildError, "FILE_SET_MISMATCH"):
            module.verify_manifest(stage)

    def test_nested_output_manifest_is_hashed_and_tamper_detected(self):
        stage = self.root / "nested-pack"
        nested = stage / "evidence"
        nested.mkdir(parents=True)
        (nested / "output-manifest.json").write_text('{"nested":1}\n', encoding="utf-8")
        module.write_json_exclusive(stage / "output-manifest.json", module.make_manifest(stage))
        module.verify_manifest(stage)
        (nested / "output-manifest.json").write_text('{"nested":2}\n', encoding="utf-8")
        with self.assertRaisesRegex(module.BuildError, "HASH_MISMATCH"):
            module.verify_manifest(stage)

    def test_topology_hashed_files_and_exact_set(self):
        report = self.root / "topology"
        report.mkdir()
        summary = report / "summary.json"
        html = report / "report.html"
        summary.write_text("{}", encoding="utf-8")
        html.write_text("report", encoding="utf-8")
        manifest = {
            "hashed_files": {
                "summary.json": hashlib.sha256(summary.read_bytes()).hexdigest(),
                "report.html": hashlib.sha256(html.read_bytes()).hexdigest(),
            },
            "exact_output_file_set": ["output-manifest.json", "report.html", "summary.json"],
        }
        (report / "output-manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        module.verify_artifact_manifest(report, manifest, self.root)
        summary.write_text('{"tampered":true}', encoding="utf-8")
        with self.assertRaisesRegex(module.BuildError, "HASH_MISMATCH"):
            module.verify_artifact_manifest(report, manifest, self.root)

    def test_topology_compact_does_not_require_comparison(self):
        candidates = []
        for candidate_index in range(6):
            candidates.append({
                "candidate_id": "C" + str(candidate_index),
                "e3_type": "CRBN",
                "seeds": [{"seed": seed, "interface_clashes": {
                    "raw_under_2A_pair_count": candidate_index,
                    "vdw_pair_count": candidate_index + 1,
                }} for seed in module.SEEDS],
                "pairwise_contact_jaccard": [{
                    "protein_protein_contact_jaccard": 0.5,
                    "target_warhead_site_jaccard": 0.6,
                    "e3_recruiter_site_jaccard": 0.7,
                } for _ in range(10)],
                "cluster_sensitivity": [{"threshold": value} for value in (0.3, 0.5, 0.7)],
            })
        compact = module.build_topology_compact({
            "candidates": candidates, "raw_seed_result_count": 30, "seeds": module.SEEDS,
        })
        first = compact["candidates"][0]
        self.assertEqual(first["clashes_under2_A_all_five"], [0] * 5)
        self.assertEqual(first["protein_protein_contact_jaccard"],
                         {"min": 0.5, "median": 0.5, "max": 0.5})
        self.assertEqual(first["target_warhead_site_jaccard"],
                         {"min": 0.6, "median": 0.6, "max": 0.6})
        self.assertEqual(first["e3_recruiter_site_jaccard"],
                         {"min": 0.7, "median": 0.7, "max": 0.7})

    def test_novel_indices_come_from_seed_selection_and_schema_is_slim(self):
        quantitative = {
            "finite_selected_ipTM_values": [
                0.8886020183563232, 0.8535048365592957, 0.8920925259590149,
                0.8611037135124207, 0.8600322008132935,
            ],
            "ipTM_IQR": 0.028569817543029785,
            "ipTM_IQR_maximum": 0.2,
            "finite_selected_endpoint_values_A": [
                5.431362069343195, 5.652832642030012, 6.20434067398785,
                5.555010776785226, 6.123707515051645,
            ],
            "endpoint_IQR_A": 0.5686967382664196,
            "endpoint_IQR_maximum_A": 1.5,
            "positive_target_warhead_and_e3_recruiter_contacts": [True] * 5,
            "quantitative_criteria_met": True,
            "automatic_approval": False,
            "qualitative_topology_or_clash_acceptance_used": False,
        }
        candidate = {
            "candidate_id": "N1", "e3_type": "CRBN",
            "old_expert_quantitative_metrics": quantitative,
            "seeds": [{"seed": seed, "selection": {"selected_model_index": index}}
                      for index, seed in enumerate(module.SEEDS)],
            "selected_model_inspections": [{"seed": 23, "inspection": {}}],
            "geometric_diagnostics": {
                "all_five_seeds_including_failures": [{"seed": seed, "interface_clashes": {
                    "raw_under_2A_pair_count": 1, "vdw_pair_count": 2,
                }} for seed in module.SEEDS],
                "pairwise_contact_jaccard": [{
                    "protein_protein_jaccard": 0.1, "target_site_jaccard": 0.2,
                    "e3_site_jaccard": 0.3, "large_array": list(range(100)),
                } for _ in range(10)],
                "cluster_sensitivity": [{"threshold": value} for value in (0.3, 0.5, 0.7)],
            },
        }
        slim = module._slim_novel_candidate(candidate)
        self.assertEqual([row["selected_model_index"] for row in slim["selected_top1_per_seed"]], list(range(5)))
        self.assertEqual(slim["all_five_quantitative_metrics"], quantitative)
        self.assertNotIn("pairwise_contact_jaccard", slim["geometric_diagnostics"])

    def test_copy_novel_uses_exact_receipts_plan_hash_and_protocol_digest(self):
        root = self.root / "novel"
        stage = self.root / "stage"
        root.mkdir()
        stage.mkdir()

        def digest(value):
            encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
            return hashlib.sha256(encoded).hexdigest()

        final_candidates = []
        summary_candidates = []
        for candidate_index, e3_type in enumerate(("CRBN", "VHL")):
            cid = "N" + str(candidate_index + 1)
            candidate_root = root / f"{cid}--{e3_type}"
            plan = self.write(f"novel/{cid}--{e3_type}/new-plan.json", b'{"plan":true}\n')
            final_candidates.append({
                "candidate_id": cid,
                "e3_type": e3_type,
                "new_plan_path": str(plan),
                "new_plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
            })
            receipts = []
            for seed in module.SEEDS:
                seed_root = candidate_root / f"seed-{seed}"
                seed_root.mkdir(parents=True)
                (seed_root / "novel_ternary.yaml").write_text("version: 1\n", encoding="utf-8")
                (seed_root / "settings.json").write_text("{}\n", encoding="utf-8")
                (seed_root / "selection-receipt.json").write_text("{}\n", encoding="utf-8")
                models = []
                for index in range(5):
                    prediction = seed_root / f"source-model-{index}.cif"
                    confidence = seed_root / f"source-confidence-{index}.json"
                    prediction.write_text(f"data_{cid}_{seed}_{index}\n", encoding="utf-8")
                    confidence.write_text(json.dumps({"index": index}) + "\n", encoding="utf-8")
                    models.append({
                        "model_index": index,
                        "raw_prediction_path": str(prediction),
                        "raw_prediction_sha256": hashlib.sha256(prediction.read_bytes()).hexdigest(),
                        "raw_confidence_path": str(confidence),
                        "raw_confidence_sha256": hashlib.sha256(confidence.read_bytes()).hexdigest(),
                    })
                receipt = {"seed": seed, "models": models}
                (seed_root / "receipt.json").write_text(
                    json.dumps(receipt, sort_keys=True), encoding="utf-8")
                receipts.append(receipt)
            summary_candidates.append({
                "candidate_id": cid, "e3_type": e3_type, "seeds": receipts,
            })

        final = {"candidates": final_candidates}
        protocol_digest = digest(final)
        final["protocol_digest"] = protocol_digest
        summary = {"protocol_digest": protocol_digest, "candidates": summary_candidates}
        (root / "final-protocol.json").write_text(json.dumps(final), encoding="utf-8")
        (root / "protocol-summary.json").write_text(json.dumps(summary), encoding="utf-8")

        provenance = {}
        original_import = module.importlib.import_module
        try:
            module.importlib.import_module = lambda name: type(
                "ProtocolModule", (), {"digest": staticmethod(digest)})
            module._copy_novel(root, summary, stage, self.root, provenance)
        finally:
            module.importlib.import_module = original_import

        self.assertIn("novel_final-protocol.json", provenance)
        self.assertIn("novel_protocol-summary.json", provenance)
        copied_plan = stage / "evidence/novel-msa/plans/N1--CRBN/plan.json"
        self.assertEqual(hashlib.sha256(copied_plan.read_bytes()).hexdigest(),
                         final_candidates[0]["new_plan_sha256"])
        copied_receipt = stage / "evidence/novel-msa/raw/N1--CRBN/seed-23/receipt.json"
        self.assertEqual(json.loads(copied_receipt.read_text(encoding="utf-8")),
                         summary_candidates[0]["seeds"][0])
        self.assertFalse(any("msa" in path.name.lower() or "cache" in path.name.lower()
                             for path in stage.rglob("*") if path.is_file()
                             and path.name != "novel_ternary.yaml"))

    def test_invalid_ranking_rule_preserves_null(self):
        evaluations = {
            "rule-" + str(index): ({"selected_pass_count": None, "valid": False, "error": "failed"}
                                    if index == 0 else {"selected_pass_count": index % 6, "valid": True})
            for index in range(7)
        }
        compact = module.build_ranking_compact({"evaluations": evaluations})
        invalid = compact["rules"]["rule-0"]
        self.assertIsNone(invalid["selected_pass_count"])
        self.assertFalse(invalid["valid"])
        self.assertEqual(invalid["error"], "failed")

    def render_fixture(self, include_novel: bool = False):
        compact = {
            "scientific_approved": False,
            "known_calibration": {
                "selected_pass_count": 2, "seed_count": 5,
                "raw_all_four_pass_count": 16, "raw_model_count": 25,
                "per_seed_selected": [{
                    "seed": seed, "selected_model_index": 0,
                    "selected_e3_CA_RMSD_after_target_alignment_A": float(index + 1),
                    "selected_pass": index < 2,
                } for index, seed in enumerate(module.SEEDS)],
            },
            "ranking_development": {"rules": {
                "rule-safe": {"selected_pass_count": 2, "valid": True, "error": None},
            }},
            "core_state": {
                "row_count": 40, "required_contact_true_count": 16,
                "required_contact_false_count": 24, "failure_count": 1,
                "failure_null_count": 39, "protein_pKa_raw_row_count": 88,
                "protein_pKa_unique_chain_residue_sequence_pKa_count": 44,
                "microstate_gate": "pending",
            },
            "topology_single_sequence": {
                "candidate_count": 1, "raw_seed_result_count": 5,
                "seeds": module.SEEDS,
                "candidates": [{
                    "candidate_id": "N1", "e3_type": "CRBN",
                    "priority_diagnostic": False,
                    "protein_protein_contact_jaccard": {"median": 0.010714285714285714},
                    "target_warhead_site_jaccard": {"median": 0.5},
                    "e3_recruiter_site_jaccard": {"median": 0.75},
                    "clashes_under2_A_all_five": [0, 1, 2, 3, 4],
                    "vdw_clashes_all_five": [1, 2, 3, 4, 5],
                }],
            },
            "sources": {
                "summary": {
                    "portable_pack_relative_path": "evidence/summary.json",
                    "active_root_relative_path": "source/summary.json",
                    "source_sha256": "a" * 64,
                },
            },
            "source_provenance": {
                "portable_pack_relative_path": "source-provenance.json",
                "sha256": "b" * 64,
            },
        }
        if include_novel:
            compact["novel_msa"] = {
                "completed_seed_count": 10, "raw_model_count": 50,
                "candidates": [{
                    "candidate_id": "N1", "e3_type": "CRBN",
                    "selected_top1_per_seed": [
                        {"seed": seed, "selected_model_index": index}
                        for index, seed in enumerate(module.SEEDS)
                    ],
                    "all_five_quantitative_metrics": {
                        "finite_selected_ipTM_values": [
                            0.8886020183563232, 0.8535048365592957,
                            0.8920925259590149, 0.8611037135124207,
                            0.8600322008132935,
                        ],
                        "ipTM_IQR": 0.028569817543029785,
                        "ipTM_IQR_maximum": 0.2,
                        "finite_selected_endpoint_values_A": [
                            5.431362069343195, 5.652832642030012,
                            6.20434067398785, 5.555010776785226,
                            6.123707515051645,
                        ],
                        "endpoint_IQR_A": 0.5686967382664196,
                        "endpoint_IQR_maximum_A": 1.5,
                        "positive_target_warhead_and_e3_recruiter_contacts": [True] * 5,
                        "quantitative_criteria_met": True,
                        "automatic_approval": False,
                        "qualitative_topology_or_clash_acceptance_used": False,
                    },
                    "geometric_diagnostics": {
                        "protein_protein_jaccard_stats": {"median": 0.27439},
                        "raw_under_2A_pair_count_all_five": [0, 1, 0, 1, 0],
                        "vdw_pair_count_all_five": [1, 2, 3, 4, 5],
                    },
                }],
            }
        return compact

    def test_render_index_escapes_hostile_dynamic_values_and_paths(self):
        compact = self.render_fixture()
        hostile = '<script>alert("x")</script>'
        compact["topology_single_sequence"]["candidates"][0]["candidate_id"] = hostile
        compact["sources"]["summary"]["portable_pack_relative_path"] = hostile + '.json'
        rendered = module.render_index(compact)
        self.assertNotIn("<script>", rendered)
        self.assertIn("&lt;script&gt;", rendered)
        self.assertNotIn('href="<script>', rendered)
        self.assertNotIn("http://", rendered)
        self.assertNotIn("https://fonts", rendered)
        self.assertNotIn("<script src=", rendered)

    def test_render_index_keeps_calibration_core_and_novel_denominators_separate(self):
        rendered = module.render_index(self.render_fixture(include_novel=True))
        self.assertIn("2 / 5", rendered)
        self.assertIn("16 / 25", rendered)
        self.assertIn("완료 5 / 원시 25 / top-1 선택 5", rendered)
        self.assertIn("50개 원시 모델에서 선택한 10개", rendered)
        self.assertIn("양성", rendered)
        self.assertIn(">16<", rendered)
        self.assertIn("음성", rendered)
        self.assertIn(">24<", rendered)
        self.assertIn(">0.010714<", rendered)
        self.assertIn(">0.27439<", rendered)
        self.assertIn(">0.263676<", rendered)
        self.assertIn("finite_selected_ipTM_values", rendered)
        self.assertIn(">0.861104<", rendered)
        self.assertIn(">0.02857<", rendered)
        self.assertIn("finite_selected_endpoint_values_A", rendered)
        self.assertIn(">5.652833<", rendered)
        self.assertIn(">0.568697<", rendered)
        self.assertIn("quantitative only, not scientific approval", rendered)
        self.assertIn("quantitative_criteria_met</code>: <strong>참</strong>", rendered)
        self.assertIn("automatic_approval</code>: <strong>거짓</strong>", rendered)
        self.assertIn("qualitative_topology_or_clash_acceptance_used</code>: <strong>거짓</strong>", rendered)
        self.assertIn("과학 승인: 거짓", rendered)

    def test_missing_baseline_identity_does_not_invent_old_value_or_delta(self):
        compact = self.render_fixture(include_novel=True)
        candidate = compact["novel_msa"]["candidates"][0]
        candidate["candidate_id"] = "absent"
        self.assertIsNone(module._old_pp_value(candidate, compact["topology_single_sequence"]))
        rendered = module.render_index(compact)
        self.assertIn("<td>미계산</td><td>0.27439</td><td>미계산</td>", rendered)

    def test_render_without_novel_is_explicitly_not_computed(self):
        compact = self.render_fixture(include_novel=False)
        rendered = module.render_index(compact)
        readme = module.render_readme(compact)
        self.assertIn("계산되지 않음", rendered)
        self.assertIn("성공, 실패 또는 승인으로 해석하지 않습니다", rendered)
        self.assertIn("계산되지 않았습니다", readme)
        self.assertNotIn("--novel-protocol-root", readme)

    def test_render_readme_adds_optional_flag_only_when_novel_present(self):
        baseline = module.render_readme(self.render_fixture(False))
        novel = module.render_readme(self.render_fixture(True))
        self.assertNotIn("--novel-protocol-root", baseline)
        self.assertIn("--novel-protocol-root", novel)
        self.assertNotIn(" \\\n  --novel-protocol-root", novel)
        self.assertIn("science-evidence-v2-REPLAY-NEW", novel)
        self.assertIn("novel-msa-r10-d5-rerun-NEW", novel)
        self.assertIn("from sources.scripts.build_campaign_diagnostics import verify_manifest", novel)
        self.assertIn("verify_manifest(Path('.'))", novel)

    def test_generated_reader_files_are_in_manifest(self):
        stage = self.root / "rendered-pack"
        stage.mkdir()
        compact = self.render_fixture()
        module.write_text_exclusive(stage / "index.html", module.render_index(compact))
        module.write_text_exclusive(stage / "README.md", module.render_readme(compact))
        manifest = module.make_manifest(stage)
        paths = {item["path"] for item in manifest["files"]}
        self.assertEqual(paths, {"README.md", "index.html"})

    def write_manifest_file(self, stage: Path):
        module.write_json_exclusive(stage / "output-manifest.json", module.make_manifest(stage))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from scripts import diagnose_calibration_ranking as ranking


class CalibrationRankingChoiceTests(unittest.TestCase):
    def models(self):
        result = []
        for index in range(5):
            result.append({
                "model_index": index,
                "confidence": {
                    "confidence_score": 0.50 + index * 0.01,
                    "iptm": 0.20 + index * 0.10,
                    "protein_iptm": 0.30 + index * 0.05,
                    "ligand_iptm": 0.40 + index * 0.03,
                    "complex_iplddt": 0.60 + index * 0.02,
                    "complex_ipde": 5.0 - index * 0.5,
                },
                "metrics": {
                    "prediction_descriptive_clashes": [
                        {"dummy_clash": clash} for clash in range(10 - index)
                    ],
                    "target_CA_RMSD_A": 1000 + index,
                    "ligand_heavy_atom_RMSD_after_target_alignment_A": 1000 + index,
                    "e3_CA_RMSD_after_target_alignment_A": 1000 + index,
                    "contact_jaccard": 0.0,
                },
            })
        return result

    def test_choice_ignores_injected_reference_values(self):
        models = self.models()
        before = ranking.select_rule(models, "iptm_max")
        for model in models:
            model["metrics"]["target_CA_RMSD_A"] = 0.0 if model["model_index"] == 0 else 9999.0
            model["metrics"]["contact_jaccard"] = 1.0 if model["model_index"] == 0 else 0.0
            model["reference_oracle_index"] = 0
        self.assertEqual(before, 4)
        self.assertEqual(ranking.select_rule(models, "iptm_max"), 4)

    def test_ties_use_confidence_then_lower_index(self):
        models = self.models()
        for model in models:
            model["confidence"]["iptm"] = 0.8
        models[1]["confidence"]["confidence_score"] = 0.99
        models[2]["confidence"]["confidence_score"] = 0.99
        self.assertEqual(ranking.select_rule(models, "iptm_max"), 1)

    def test_baseline_tie_is_deterministic(self):
        models = self.models()
        for model in models:
            model["confidence"]["confidence_score"] = 0.7
        self.assertEqual(ranking.select_rule(models, "confidence_score_max"), 0)

    def test_missing_field_rejected(self):
        models = self.models()
        del models[3]["confidence"]["protein_iptm"]
        with self.assertRaisesRegex(ranking.DiagnosisError, "RULE_FIELD_UNAVAILABLE"):
            ranking.select_rule(models, "protein_iptm_max")

    def test_nonfinite_field_rejected(self):
        models = self.models()
        models[2]["confidence"]["complex_iplddt"] = math.nan
        with self.assertRaisesRegex(ranking.DiagnosisError, "FINITE_NUMERIC_FIELD_REQUIRED"):
            ranking.select_rule(models, "complex_iplddt_max")

    def test_out_of_range_field_rejected(self):
        models = self.models()
        models[0]["confidence"]["iptm"] = 1.01
        with self.assertRaisesRegex(ranking.DiagnosisError, "FIELD_OUT_OF_NATURAL_RANGE"):
            ranking.select_rule(models, "iptm_max")

    def test_lower_ipde_is_correct(self):
        models = self.models()
        models[2]["confidence"]["complex_ipde"] = 0.1
        models[4]["confidence"]["complex_ipde"] = 0.2
        self.assertEqual(ranking.select_rule(models, "complex_ipde_min"), 2)

    def test_clash_rule_uses_only_predicted_clash_field(self):
        models = self.models()
        models[1]["metrics"]["prediction_descriptive_clashes"] = []
        models[4]["metrics"]["target_CA_RMSD_A"] = 0
        self.assertEqual(ranking.select_rule(models, "clash_count_min"), 1)

    def test_clash_rule_distinguishes_empty_from_one_record(self):
        models = self.models()
        models[0]["metrics"]["prediction_descriptive_clashes"] = [{"dummy": True}]
        models[1]["metrics"]["prediction_descriptive_clashes"] = []
        self.assertEqual(ranking.select_rule(models, "clash_count_min"), 1)

    def test_clash_rule_rejects_numeric_malformed_value(self):
        models = self.models()
        models[0]["metrics"]["prediction_descriptive_clashes"] = 0
        with self.assertRaisesRegex(ranking.DiagnosisError, "CLASH_RECORD_LIST_REQUIRED"):
            ranking.select_rule(models, "clash_count_min")

    def test_exactly_five_unique_models_required(self):
        with self.assertRaisesRegex(ranking.DiagnosisError, "ALL_FIVE_MODELS_REQUIRED"):
            ranking.select_rule(self.models()[:4], "confidence_score_max")
        models = self.models()
        models[4]["model_index"] = 3
        with self.assertRaisesRegex(ranking.DiagnosisError, "MODEL_INDEX_SET_INVALID"):
            ranking.select_rule(models, "confidence_score_max")


class CalibrationRankingOutputTests(unittest.TestCase):
    def receipt(self, seed, passing_index=4):
        models = CalibrationRankingChoiceTests().models()
        for model in models:
            passes = model["model_index"] == passing_index
            model["metrics"].update({
                "target_CA_RMSD_A": 1.0 if passes else 2.0,
                "ligand_heavy_atom_RMSD_after_target_alignment_A": 2.0 if passes else 4.0,
                "e3_CA_RMSD_after_target_alignment_A": 5.0 if passes else 11.0,
                "contact_jaccard": 0.5 if passes else 0.2,
            })
        return {"seed": seed, "models": models}

    def test_proposal_keeps_invalid_rule_selection_null(self):
        receipts = [
            self.receipt(seed, passing_index=4 if position < 2 else 0)
            for position, seed in enumerate(ranking.SEEDS)
        ]
        del receipts[0]["models"][0]["confidence"]["ligand_iptm"]
        proposal = ranking.proposals(receipts)
        rule = proposal["rules"]["ligand_iptm_max"]
        self.assertFalse(rule["valid"])
        self.assertIsNone(rule["selection_indices"][str(ranking.SEEDS[0])])
        report = ranking.evaluate(receipts, proposal)
        evaluation = report["evaluations"]["ligand_iptm_max"]
        self.assertIsNone(evaluation["selected_pass_count"])
        self.assertIn("RULE_FIELD_UNAVAILABLE", evaluation["error"])

    def test_diagnosis_has_no_result_approval_or_protocol_pass_flag(self):
        receipts = []
        for position, seed in enumerate(ranking.SEEDS):
            # Confidence top-1 is model 4 and passes only for the first two seeds.
            receipts.append(self.receipt(seed, passing_index=4 if position < 2 else 0))
        proposal = ranking.proposals(receipts)
        report = ranking.evaluate(receipts, proposal)
        encoded = json.dumps(report, sort_keys=True)
        self.assertNotIn("scientifically_approved", encoded)
        self.assertNotIn("protocol_pass", encoded)
        self.assertNotIn("recommended_winner", encoded)
        self.assertEqual(
            report["evaluations"]["confidence_score_max"]["selected_pass_count"], 2
        )
        self.assertIn("No winner is recommended", report["recommendation"])
        self.assertTrue(any("Zero results" in item for item in report["limitations"]))

    def test_write_json_exclusive_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "value.json"
            ranking.write_json_exclusive(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                ranking.write_json_exclusive(path, {"a": 2})


if __name__ == "__main__":
    unittest.main()

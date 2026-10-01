from __future__ import annotations

import math
import shutil
import tempfile
import unittest
from pathlib import Path

from packages.science.expert_followup import (
    EXPECTED_SEEDS,
    ROOT,
    SOURCE_RELATIVE_PATH,
    evaluate_calibration,
    evaluate_novel,
    load_followup,
)


LIMIT_VALUES = {
    "target_CA_RMSD_A": 1.5,
    "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
    "e3_CA_RMSD_after_target_alignment_A": 10.0,
    "contact_jaccard": 0.40,
}


def calibration_row(seed, *, passing=True, metrics=None, technical=True):
    values = dict(LIMIT_VALUES)
    if not passing:
        values["target_CA_RMSD_A"] = 1.5001
    if metrics:
        values.update(metrics)
    return {
        "seed": seed,
        "technical_execution_success": technical,
        "metrics": values,
        "original_status": "failed" if seed == 41 else "completed",
        "original_failure_reason": "retained poor seed" if seed == 41 else None,
    }


def distribution(passing_seeds):
    return {"seeds": [calibration_row(seed, passing=seed in passing_seeds)
                      for seed in EXPECTED_SEEDS]}


def novel_receipt(e3, candidate, seed, *, iptm=.7, endpoint=8.0,
                  target_contacts=2, e3_contacts=3):
    return {
        "candidate_id": candidate,
        "e3_type": e3,
        "seed": seed,
        "status": "completed",
        "execution_success": True,
        "actual_computation": True,
        "reference_free": True,
        "inspection": {
            "model_confidence": {"iptm": iptm},
            "descriptive_metrics": {
                "linker_endpoint_distance_A": endpoint,
                "contacts_by_protein_role_and_ligand_group": {
                    "target": {"warhead": target_contacts},
                    "e3": {"recruiter": e3_contacts},
                },
                "protein_ligand_clashes": {"count": seed % 3},
                "target_e3_heavy_atom_clashes": {"count": seed % 2},
            },
        },
    }


class CalibrationTests(unittest.TestCase):
    def test_exact_boundaries_and_prefixed_keys_pass(self):
        rows = []
        for seed in EXPECTED_SEEDS:
            row = calibration_row(seed)
            if seed == 23:
                row["metrics"] = {"comparison.metrics." + key: value
                                  for key, value in LIMIT_VALUES.items()}
            rows.append(row)
        result = evaluate_calibration({"seeds": rows})
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["passing_seed_count"], 5)
        self.assertEqual(result["seed_verdicts"][0]["evaluations"][0]
                         ["e3_geometry_class"], "screening_acceptable")
        self.assertFalse(result["scientific_accepted"])
        self.assertEqual(result["failure_count"], 0)
        self.assertEqual(result["missing_count"], 0)

    def test_three_of_five_pass_with_bad_seed_41_retained(self):
        result = evaluate_calibration(distribution({23, 61, 79}))
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["passing_seed_count"], 3)
        seed41 = next(item for item in result["seed_verdicts"] if item["seed"] == 41)
        self.assertFalse(seed41["passed"])
        self.assertEqual(seed41["evaluations"][0]["original_failure_reason"],
                         "retained poor seed")

    def test_two_of_five_fails_without_averaging(self):
        result = evaluate_calibration(distribution({23, 61}))
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["passing_seed_count"], 2)
        self.assertEqual(result["failure_count"], 3)
        self.assertEqual(result["missing_count"], 0)

    def test_nonfinite_bool_negative_and_jaccard_domain_fail(self):
        rows = distribution({23, 41, 61, 79, 97})["seeds"]
        rows[0]["metrics"]["target_CA_RMSD_A"] = True
        rows[1]["metrics"]["ligand_heavy_atom_RMSD_after_target_alignment_A"] = math.inf
        rows[2]["metrics"]["e3_CA_RMSD_after_target_alignment_A"] = -0.1
        rows[3]["metrics"]["contact_jaccard"] = 1.01
        rows[4]["metrics"]["contact_jaccard"] = False
        result = evaluate_calibration({"seeds": rows})
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["passing_seed_count"], 0)
        self.assertNotIn(True, result["statistics"]["target_CA_RMSD_A"]
                         ["all_finite_values"])

    def test_prefixed_bare_conflict_rejected(self):
        rows = distribution(set(EXPECTED_SEEDS))["seeds"]
        rows[0]["metrics"]["comparison.metrics.target_CA_RMSD_A"] = 0.1
        result = evaluate_calibration({"seeds": rows})
        first = result["seed_verdicts"][0]["evaluations"][0]
        self.assertEqual(first["metric_conflicts"], ["target_CA_RMSD_A"])
        self.assertFalse(first["passed_all_metrics"])

    def test_duplicate_fails_and_missing_is_pending(self):
        rows = distribution(set(EXPECTED_SEEDS))["seeds"]
        duplicate = evaluate_calibration({"seeds": rows + [dict(rows[0])]})
        self.assertEqual(duplicate["status"], "failed")
        self.assertEqual(duplicate["duplicate_seeds"], [23])
        missing = evaluate_calibration({"seeds": rows[:-1]})
        self.assertEqual(missing["status"], "pending")
        self.assertEqual(missing["missing_seeds"], [97])
        self.assertEqual(missing["missing_count"], 1)

    def test_nonfinite_raw_values_are_json_safe_and_represented(self):
        rows = distribution(set(EXPECTED_SEEDS))["seeds"]
        rows[0]["metrics"]["target_CA_RMSD_A"] = math.nan
        result = evaluate_calibration({"seeds": rows})
        import json
        json.dumps(result, allow_nan=False)
        self.assertEqual(result["all_outcomes"][0]["metrics"]
                         ["target_CA_RMSD_A"]["nonfinite_numeric"], "nan")


class NovelTests(unittest.TestCase):
    def _priority_receipts(self):
        rows = []
        for index, seed in enumerate(EXPECTED_SEEDS):
            rows.append(novel_receipt("CRBN", "D-99b12e64986a", seed,
                                      iptm=.60 + index * .01,
                                      endpoint=8.0 + index * .1))
            rows.append(novel_receipt("VHL", "D-b39273b7a53b", seed,
                                      iptm=.75 + index * .01,
                                      endpoint=9.0 + index * .1))
        for index in range(20):
            rows.append(novel_receipt("CRBN" if index % 2 else "VHL",
                                      "D-nonpriority-%02d" % index,
                                      1000 + index))
        return rows

    def test_low_iqr_nonzero_contacts_still_geometry_pending(self):
        rows = self._priority_receipts()
        candidates = {"candidates": [
            {"candidate_id": "D-99b12e64986a", "e3_type": "CRBN",
             "warhead_analog_id": "W-c2afc5e73c1a"},
            {"candidate_id": "D-b39273b7a53b", "e3_type": "VHL",
             "warhead_analog_id": "W-80f8f4a11b5d"},
        ]}
        result = evaluate_novel(rows, candidates)
        self.assertEqual(result["quantitative_status"], "pass")
        self.assertEqual(result["geometry_review"]["status"], "pending")
        self.assertFalse(result["scientific_accepted"])
        self.assertEqual(len(result["all_raw_receipts"]), 30)
        self.assertFalse(result["cross_e3_raw_score_winner_selected"])
        self.assertEqual(result["priority_groups"]["CRBN"]["receipts"][0]
                         ["protein_ligand_clashes"]["count"], 2)

    def test_failure_and_duplicate_not_accepted(self):
        rows = self._priority_receipts()[:10]
        rows[0]["status"] = "failed"
        result = evaluate_novel(rows)
        self.assertEqual(result["quantitative_status"], "failed")
        rows = self._priority_receipts()[:10]
        rows.append(dict(rows[0]))
        duplicate = evaluate_novel(rows)
        self.assertEqual(duplicate["quantitative_status"], "failed")
        self.assertEqual(duplicate["priority_groups"]["CRBN"]["duplicate_seeds"], [23])

    def test_execution_success_domains_contacts_and_invalid_priority_seed(self):
        rows = self._priority_receipts()[:10]
        rows[0]["execution_success"] = False
        rows[1]["inspection"]["model_confidence"]["iptm"] = 1.01
        rows[2]["inspection"]["descriptive_metrics"]\
            ["contacts_by_protein_role_and_ligand_group"]["target"]["warhead"] = 1.5
        rows[3]["seed"] = True
        result = evaluate_novel(rows)
        self.assertEqual(result["quantitative_status"], "failed")
        self.assertTrue(result["priority_groups"]["VHL"]["invalid_priority_rows"])

    def test_completed_text_cannot_override_failure_and_nonfinite_is_safe(self):
        rows = self._priority_receipts()[:10]
        rows[0]["failure"] = {"reason": "producer failed"}
        rows[1]["inspection"]["model_confidence"]["iptm"] = math.inf
        result = evaluate_novel(rows)
        import json
        json.dumps(result, allow_nan=False)
        self.assertEqual(result["quantitative_status"], "failed")
        self.assertTrue(result["priority_groups"]["CRBN"]["receipts"][0]
                        ["contradictory_failure_status"])

    def test_missing_seed_pending_and_identity_mismatch_failed(self):
        rows = self._priority_receipts()[:8]
        pending = evaluate_novel(rows)
        self.assertEqual(pending["quantitative_status"], "pending")
        candidates = [{"candidate_id": "D-99b12e64986a", "e3_type": "CRBN",
                       "warhead_analog_id": "wrong"},
                      {"candidate_id": "D-b39273b7a53b", "e3_type": "VHL",
                       "warhead_analog_id": "W-80f8f4a11b5d"}]
        failed = evaluate_novel(self._priority_receipts()[:10], candidates)
        self.assertEqual(failed["quantitative_status"], "failed")
        self.assertEqual(failed["candidate_identity_failures"], ["CRBN"])


class SourceBindingTests(unittest.TestCase):
    def test_real_source_loads_without_mutation(self):
        result = load_followup()
        self.assertFalse(result["authority"])
        self.assertFalse(result["scientific_accepted"])
        self.assertEqual(result["source"]["bytes"], 41532)
        self.assertTrue(result["rules"])
        self.assertTrue(all("locator" in rule and "extract_sha256" in rule
                            for rule in result["rules"]))
        parents = next(rule for rule in result["rules"]
                       if rule["id"] == "exploratory_parents")["allowed"]
        self.assertTrue(all(parent.startswith("SMARCA2-") for parent in parents))

    def test_copied_source_tamper_is_rejected(self):
        source = ROOT / SOURCE_RELATIVE_PATH
        original = source.read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SOURCE_RELATIVE_PATH
            target.parent.mkdir(parents=True)
            shutil.copyfile(source, target)
            tampered = bytearray(original)
            tampered[-1] ^= 1
            target.write_bytes(tampered)
            with self.assertRaises(Exception):
                load_followup(root)
        self.assertEqual(source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()

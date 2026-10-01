import copy
import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from packages.contracts import ContractError
from packages.science import expert_followup


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads((
    ROOT / "contracts" / "drafts" / "scientific_acceptance.schema.json"
).read_text(encoding="utf-8"))


class FollowupGeometryReviewTests(unittest.TestCase):
    def _receipt(self, e3, candidate, seed, iptm, endpoint):
        return {
            "candidate_id": candidate,
            "e3_type": e3,
            "seed": seed,
            "actual_computation": True,
            "execution_success": True,
            "reference_free": True,
            "failure": None,
            "inspection": {
                "model_confidence": {"iptm": iptm},
                "descriptive_metrics": {
                    "linker_endpoint_distance_A": endpoint,
                    "contacts_by_protein_role_and_ligand_group": {
                        "target": {"warhead": 3},
                        "e3": {"recruiter": 2},
                    },
                    "protein_ligand_clashes": 1,
                    "target_e3_heavy_atom_clashes": 1,
                },
            },
        }

    def _quantitative(self):
        receipts = []
        candidates = []
        for e3, (candidate, warhead) in expert_followup._PRIORITIES.items():
            candidates.append({
                "candidate_id": candidate,
                "e3_type": e3,
                "warhead_analog_id": warhead,
            })
            for index, seed in enumerate(expert_followup.EXPECTED_SEEDS):
                receipts.append(self._receipt(
                    e3, candidate, seed, .70 + index * .01, 8.0 + index * .1
                ))
        result = expert_followup.evaluate_novel(receipts, candidates)
        self.assertEqual(result["quantitative_status"], "pass")
        return result

    def _ref(self, name):
        return {
            "artifact_id": name,
            "version": 1,
            "sha256": "a" * 64,
            "media_type": "application/json",
            "schema_id": "saved-result/1",
            "provenance": "computed",
        }

    def _decision(self, labels=None):
        labels = labels or ["A", "A", "A", "A", "B"]
        branches = []
        for e3, (candidate, _warhead) in expert_followup._PRIORITIES.items():
            reviews = []
            for seed, label in zip(expert_followup.EXPECTED_SEEDS, labels):
                evidence_ref = self._ref(f"report-{candidate}-{seed}")
                reviews.append({
                    "seed": seed,
                    "topology_label": label,
                    "severe_clash_acceptable": True,
                    "rationale": "Authenticated per-seed expert review.",
                    "evidence_ref": evidence_ref,
                    "evidence_binding": {
                        "evidence_id": f"evidence-{candidate}-{seed}",
                        "evidence_ref": evidence_ref,
                        "evidence_sha256": "b" * 64,
                        "artifact_sha256": "a" * 64,
                        "raw_artifact_sha256": {
                            "receipt": "c" * 64,
                            "prediction": "d" * 64,
                            "plan_ref": "e" * 64,
                            "aggregate_receipt_ref": "f" * 64,
                        },
                        "candidate_graph_sha256": "1" * 64,
                        "plan_digest": "plan-digest",
                        "project_id": "project",
                        "job_id": "job",
                        "input_sha256": "2" * 64,
                        "result_sha256": "3" * 64,
                        "parent_id": "SMARCA2-FX5",
                        "candidate_id": candidate,
                        "e3_type": e3,
                        "seed": seed,
                        "actual_computation": True,
                        "execution_success": True,
                        "reference_free": True,
                    },
                })
            branches.append({
                "candidate_id": candidate,
                "e3_type": e3,
                "clash_definition": "Expert-defined severe-clash recurrence.",
                "rationale": "Authenticated branch review.",
                "seed_reviews": reviews,
            })
        source_ref = self._ref("statement")
        return {
            "_validated_by_platform": True,
            "decision_id": "decision-1",
            "actor": {"id": "actor-1"},
            "scope": {"project_id": "project", "job_id": "job"},
            "source_ref": source_ref,
            "source_binding": {
                "project_id": "project",
                "job_id": "job",
                "result_sha256": "3" * 64,
                "source_ref": source_ref,
                "source_sha256": "a" * 64,
                "followup_source_policy_digest":
                    expert_followup.SOURCE_POLICY_DIGEST,
            },
            "branches": branches,
        }

    def test_platform_gate_is_required(self):
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), self._decision())

    def test_four_of_five_and_quantitative_pass(self):
        result = expert_followup.evaluate_geometry_review(
            self._quantitative(), self._decision(), True)
        self.assertEqual(result["status"], "pass")
        self.assertEqual(
            result["disposition"], "internal_geometry_review_candidate")
        self.assertFalse(result["scientific_accepted"])
        self.assertFalse(result["formal_acceptance"])

    def test_five_of_five_passes(self):
        result = expert_followup.evaluate_geometry_review(
            self._quantitative(), self._decision(["A"] * 5), True)
        self.assertEqual(result["status"], "pass")

    def test_no_four_matching_labels_fails(self):
        result = expert_followup.evaluate_geometry_review(
            self._quantitative(),
            self._decision(["A", "A", "A", "B", "B"]), True)
        self.assertEqual(result["status"], "failed")
        self.assertTrue(any(value.startswith("insufficient_same_topology")
                            for value in result["failures"]))

    def test_explicit_false_clash_review_fails(self):
        decision = self._decision()
        decision["branches"][0]["seed_reviews"][2][
            "severe_clash_acceptable"] = False
        result = expert_followup.evaluate_geometry_review(
            self._quantitative(), decision, True)
        self.assertEqual(result["status"], "failed")
        self.assertTrue(any(value.startswith("severe_clash_rejected")
                            for value in result["failures"]))

    def test_duplicate_seed_is_rejected(self):
        decision = self._decision()
        decision["branches"][0]["seed_reviews"][4]["seed"] = 23
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), decision, True)

    def test_duplicate_branch_is_rejected_not_overwritten(self):
        decision = self._decision()
        decision["branches"][1] = copy.deepcopy(decision["branches"][0])
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), decision, True)

    def test_wrong_candidate_evidence_identity_is_rejected(self):
        decision = self._decision()
        decision["branches"][0]["seed_reviews"][0][
            "evidence_binding"]["candidate_id"] = "D-b39273b7a53b"
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), decision, True)

    def test_wrong_seed_evidence_identity_is_rejected(self):
        decision = self._decision()
        decision["branches"][0]["seed_reviews"][0][
            "evidence_binding"]["seed"] = 41
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), decision, True)

    def test_wrong_source_or_result_identity_is_rejected(self):
        decision = self._decision()
        decision["branches"][0]["seed_reviews"][0][
            "evidence_binding"]["result_sha256"] = "4" * 64
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                self._quantitative(), decision, True)

    def test_topline_pass_without_group_coverage_is_rejected(self):
        quantitative = self._quantitative()
        quantitative["priority_groups"]["CRBN"]["receipts"].pop()
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                quantitative, self._decision(), True)

    def test_nonfinite_or_excess_iqr_is_rejected(self):
        quantitative = self._quantitative()
        quantitative["priority_groups"]["CRBN"]["iptm"]["iqr"] = .21
        with self.assertRaises(ContractError):
            expert_followup.evaluate_geometry_review(
                quantitative, self._decision(), True)

    def test_schema_rejects_duplicate_seed_coverage(self):
        decision = self._decision()
        submission = {
            "kind": "followup_geometry_review",
            "action": "accept",
            "intent": "record authenticated review",
            "reason": "source-bound review",
            "expected_policy_revision": 0,
            "data": {
                "branches": [],
                "source_ref": decision["source_ref"],
            },
        }
        for branch in decision["branches"]:
            clean = copy.deepcopy(branch)
            for review in clean["seed_reviews"]:
                review.pop("evidence_binding")
            submission["data"]["branches"].append(clean)
        wrapper = {"kind": "decision", "payload": submission}
        errors = list(Draft202012Validator(SCHEMA).iter_errors(wrapper))
        self.assertEqual(errors, [])
        submission["data"]["branches"][0]["seed_reviews"][4]["seed"] = 23
        self.assertEqual(
            list(Draft202012Validator(SCHEMA).iter_errors(wrapper)), [])

    def test_records_are_not_mutated(self):
        quantitative = self._quantitative()
        decision = self._decision()
        original_quantitative = copy.deepcopy(quantitative)
        original_decision = copy.deepcopy(decision)
        expert_followup.evaluate_geometry_review(
            quantitative, decision, True)
        self.assertEqual(quantitative, original_quantitative)
        self.assertEqual(decision, original_decision)


if __name__ == "__main__":
    unittest.main()

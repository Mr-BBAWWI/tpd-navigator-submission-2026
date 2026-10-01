import unittest

from packages.science.scientific_assessment import evaluate_design


def criterion(assessment, identifier):
    return {row["id"]: row for row in assessment["criteria"]}[identifier]


def analog(identifier, smiles):
    return {
        "id": identifier,
        "mapped_smiles": smiles,
        "canonical_smiles": smiles,
        "attachment_site_atom_maps": [19],
        "transformation_class": "ring_expansion",
        "cheap_filter": {"valid": True},
        "selected": True,
        "pipeline_status": "qualified",
        "assembly_eligible": True,
        "parent_redocking_supported": True,
        "docking": {
            "status": "completed_with_limits",
            "pose_preserved": True,
            "passing_pose_count": 1,
        },
    }


def base_result():
    return {
        "sites": {"atoms": [{"atom_map": 19, "state": "MODIFIABLE"}]},
        "analogs": [],
        "protac_candidates": [],
        "calibration": {"CRBN": {"seed_receipts": []}},
    }


def trusted_policy(decisions=None, scope=None):
    return {
        "_validated_by_platform": True,
        "numerical_criteria": {"panel_min": 1, "panel_max": 20},
        "decisions": decisions or {},
        "scope": scope,
    }


def bound_report(digest, **values):
    report = {
        "_validated_by_platform": True,
        "job_id": values.pop("job_id", "worker-job"),
        "policy_digest": digest,
    }
    report.update(values)
    return report


class ScientificCoverageRegressionTests(unittest.TestCase):
    def test_one_interaction_report_cannot_cover_two_selected_analogs(self):
        result = base_result()
        result["analogs"] = [analog("A-1", "CC"), analog("A-2", "CCC")]
        policy = trusted_policy()
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        report = bound_report(
            digest,
            analog_id="A-1",
            report_digest="interaction-a1",
            requirements=[{
                "id": "required-contact",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
            }],
            failures=[],
        )

        assessment = evaluate_design(
            result,
            policy=policy,
            supplements={"interactions": [report]},
        )
        row = criterion(assessment, "core_interaction_preservation")

        self.assertEqual(row["status"], "pending")
        self.assertEqual(row["observed"]["missing_analogs"], ["A-2"])
        self.assertIs(assessment["scientific_accepted"], False)

    def test_scoped_report_uses_actual_job_id_without_current_job_id(self):
        result = base_result()
        result["analogs"] = [analog("A-1", "CC")]
        policy = trusted_policy(scope={"job_id": "service-job"})
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        report = bound_report(
            digest,
            job_id="service-job",
            analog_id="A-1",
            requirements=[{
                "id": "required-contact",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
            }],
            failures=[],
        )

        row = criterion(
            evaluate_design(result, policy=policy, supplements={"interactions": [report]}),
            "core_interaction_preservation",
        )

        self.assertEqual(row["status"], "pass")
        self.assertEqual(row["observed"]["covered_analogs"], ["A-1"])

    def test_scoped_report_rejects_conflicting_optional_current_job_id(self):
        result = base_result()
        result["analogs"] = [analog("A-1", "CC")]
        policy = trusted_policy(scope={"job_id": "service-job"})
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        report = bound_report(
            digest,
            job_id="service-job",
            current_job_id="different-job",
            analog_id="A-1",
            requirements=[{
                "id": "required-contact",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
            }],
            failures=[],
        )

        row = criterion(
            evaluate_design(result, policy=policy, supplements={"interactions": [report]}),
            "core_interaction_preservation",
        )

        self.assertEqual(row["status"], "pending")
        self.assertEqual(row["observed"]["missing_analogs"], ["A-1"])
        self.assertEqual(row["observed"]["unbound_reports"], 1)

    def test_nested_actual_hydrogen_receipts_require_review(self):
        result = base_result()
        result["analogs"] = [analog("A-1", "CC")]
        policy = trusted_policy({
            "microstate_expert_decision": {
                "accepted": True,
                "state_choices": {"A-1": 0},
            }
        })
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        report = bound_report(
            digest,
            analog_id="A-1",
            hydrogen_receipts={
                "parent": {"requires_review": False},
                "pose": {"requires_review": True},
            },
            requirements=[{
                "id": "required-contact",
                "required": True,
                "status": "preserved",
                "computed_pass": True,
            }],
            failures=[],
            state_alternatives=[{
                "status": "computed",
                "pH_conditions": {"pH": 7.4},
                "requirements": [{
                    "id": "required-contact",
                    "required": True,
                    "computed_pass": True,
                    "unambiguous": True,
                    "pending_missing_protein_hydrogen": False,
                }],
            }],
        )

        assessment = evaluate_design(
            result,
            policy=policy,
            supplements={"interactions": [report]},
        )
        core = criterion(assessment, "core_interaction_preservation")
        microstate = criterion(assessment, "microstates_h_direction")

        self.assertEqual(core["status"], "pending")
        self.assertTrue(any(
            item.get("reason") == "parent_or_pose_hydrogen_receipt_requires_review"
            for item in core["observed"]["pending"]
        ))
        self.assertEqual(microstate["status"], "pending")
        self.assertTrue(any(
            item.get("reason") == "selected_state_failed_or_minimizer_requires_review"
            for item in microstate["observed"]["pending"]
        ))

    def test_optional_only_interaction_requirements_do_not_vacuously_pass(self):
        result = base_result()
        result["analogs"] = [analog("A-1", "CC")]
        policy = trusted_policy()
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        report = bound_report(
            digest,
            analog_id="A-1",
            requirements=[{
                "id": "optional-contact",
                "required": False,
                "status": "preserved",
                "computed_pass": True,
            }],
            failures=[],
        )

        row = criterion(
            evaluate_design(result, policy=policy, supplements={"interactions": [report]}),
            "core_interaction_preservation",
        )

        self.assertEqual(row["status"], "pending")
        self.assertTrue(any(item["reason"] == "nonempty_required_requirements_required" for item in row["observed"]["pending"]))

    def test_ternary_report_cannot_claim_candidate_under_wrong_e3(self):
        result = base_result()
        result["protac_candidates"] = [
            {"candidate_id": "C-CRBN", "e3_type": "CRBN", "mapped_smiles": "CC"},
            {"candidate_id": "C-VHL", "e3_type": "VHL", "mapped_smiles": "CCC"},
        ]
        policy = trusted_policy({
            "ternary_criterion": {
                "chosen_per_e3": {"CRBN": "C-CRBN", "VHL": "C-VHL"},
                "expected_seeds": [23, 41],
            },
            "ternary_geometry_criterion": {
                "disposition": "accept",
                "metric": "collision_score",
                "threshold": 1.0,
                "comparison": "lte",
            },
        })
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        reports = [
            bound_report(
                digest,
                candidate_id="C-CRBN",
                e3_type="VHL",
                seed=seed,
                actual_computation=True,
                execution_success=True,
                reference_free=True,
                inspection_summary={
                    "compatibility": "compatible",
                    "metrics": {"collision_score": 0.1},
                },
            )
            for seed in (23, 41)
        ]

        row = criterion(
            evaluate_design(result, policy=policy, supplements={"ternary": reports}),
            "novel_ternary_repeats",
        )

        self.assertEqual(row["status"], "pending")
        self.assertTrue(all(not outcome["candidate_e3_match"] for outcome in row["observed"]["all_outcomes"]))

    def test_calibration_requires_metric_for_every_expected_seed(self):
        result = base_result()
        result["calibration"]["CRBN"]["seed_receipts"] = [
            {
                "seed": seed,
                "execution_success": True,
                "inspection": {"comparison": {"metrics": metrics}},
            }
            for seed, metrics in (
                (23, {"rmsd": 1.0}),
                (41, {}),
                (67, {"rmsd": 1.5}),
            )
        ]
        policy = trusted_policy({
            "calibration_criterion": {
                "metric": "rmsd",
                "threshold": 2.0,
                "expected_seeds": [23, 41, 67],
                "scope": "CRBN",
                "comparison": "max_lte",
            }
        })

        row = criterion(evaluate_design(result, policy=policy), "known_crbn_calibration")

        self.assertEqual(row["status"], "pending")
        self.assertEqual(row["observed"]["missing_metric_seeds"], [41])
        self.assertTrue(row["observed"]["known_seed_41_retained"])
        self.assertEqual(
            row["observed"]["metric_distributions"]["rmsd"]["by_seed"],
            {"23": [1.0], "67": [1.5]},
        )

    def test_complete_diagnostic_ternary_pool_without_chosen_geometry_is_pending(self):
        result = base_result()
        result["protac_candidates"] = [
            {"candidate_id": "C-CRBN", "e3_type": "CRBN", "mapped_smiles": "CC"},
            {"candidate_id": "C-VHL", "e3_type": "VHL", "mapped_smiles": "CCC"},
        ]
        policy = trusted_policy()
        digest = evaluate_design(result, policy=policy)["policy_digest"]
        reports = []
        for candidate_id, e3 in (("C-CRBN", "CRBN"), ("C-VHL", "VHL")):
            for seed in (23, 41):
                reports.append(bound_report(
                    digest,
                    candidate_id=candidate_id,
                    e3_type=e3,
                    seed=seed,
                    actual_computation=True,
                    execution_success=True,
                    reference_free=True,
                    inspection_summary={
                        "compatibility": "compatible",
                        "metrics": {"collision_score": 0.1},
                    },
                ))

        assessment = evaluate_design(result, policy=policy, supplements={"ternary": reports})

        self.assertEqual(criterion(assessment, "novel_ternary_repeats")["status"], "pending")
        self.assertEqual(criterion(assessment, "novel_ternary_geometry")["status"], "pending")
        self.assertIs(assessment["scientific_accepted"], False)


if __name__ == "__main__":
    unittest.main()

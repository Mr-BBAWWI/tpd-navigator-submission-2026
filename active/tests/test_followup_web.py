import pathlib
import re
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
JS = (ROOT / "apps/web/scientific-review.js").read_text(encoding="utf-8")
HTML = (ROOT / "apps/web/scientific-review.html").read_text(encoding="utf-8")
CSS = (ROOT / "apps/web/scientific-review.css").read_text(encoding="utf-8")


class FollowupWebContractTests(unittest.TestCase):
    def test_followup_fallback_is_source_only_and_freshness_bound(self):
        self.assertIn('id="expertFollowup"', HTML)
        self.assertIn('Legacy relayed expert opinion', HTML)
        self.assertLess(HTML.index('id="expertFollowup"'), HTML.index('id="expertOpinion"'))
        self.assertIn("state.policy?.followup", JS)
        self.assertIn("current_policy===true", JS)
        self.assertIn("current_implementation===true", JS)
        self.assertIn("source-only policy fallback", JS)
        self.assertIn("이전 assessment의 정량 subcheck를 현재 결과처럼 표시하지 않습니다", JS)
        self.assertIn("extractHash(value)", JS)

    def test_quantitative_tables_use_dynamic_priority_ids_and_compact_clashes(self):
        for token in (
            "target_CA_RMSD_A",
            "ligand_heavy_atom_RMSD_after_target_alignment_A",
            "e3_CA_RMSD_after_target_alignment_A",
            "contact_jaccard",
            "Poor seed 41 is retained",
            "range · median · IQR",
            "group?.candidate_id",
            "Target-warhead contacts",
            "E3-recruiter contacts",
            "Protein-ligand clash",
            "Target-E3 clash",
            "heavyClashValue(x.protein_ligand_clashes)",
            "authenticated_followup_geometry_review",
        ):
            self.assertIn(token, JS)
        self.assertNotIn("D-99b12e64986a", JS)
        self.assertNotIn("D-b39273b7a53b", JS)
        self.assertIn("cell-details", CSS)

    def test_state_set_note_requires_current_selected_source_and_pose_range(self):
        for identifier in (
            "stateSetReviewForm",
            "stateSetAnalog",
            "stateSetPose",
            "statePopulationEvidence",
            "stateSetObservations",
            "stateSetReviewArtifact",
        ):
            self.assertIn(f'id="{identifier}"', HTML)
        self.assertIn('max="4"', HTML)
        handler = re.search(
            r"\$\('stateSetReviewForm'\)\.onsubmit=async e=>\{(.+?)\n\$\('geometryTemplateButton'\)",
            JS,
            re.DOTALL,
        )
        self.assertIsNotNone(handler)
        body = handler.group(1)
        self.assertIn("analogRows().some", body)
        self.assertIn("x?.selected===true", body)
        self.assertIn("poseIndex>4", body)
        self.assertIn("verifiedFollowupSourceSHA()", body)
        self.assertIn("if(!observations)", body)
        self.assertIn("followup_source_sha256:followupSHA", body)
        self.assertNotIn("followup_source_sha256:null", body)
        self.assertIn("selected_state:null", body)
        self.assertIn("approval:false", body)
        self.assertIn("policy_decision:false", body)
        self.assertIn("microstate_gate_effect:'none'", body)
        self.assertNotIn("innerHTML", JS)

    def test_geometry_form_elements_and_handlers_match(self):
        identifiers = (
            "geometryReviewForm",
            "geometryTemplateButton",
            "geometryBranches",
            "geometryReviewRationale",
            "geometryReviewButton",
            "geometryReviewStatus",
        )
        for identifier in identifiers:
            self.assertEqual(HTML.count(f'id="{identifier}"'), 1)
            self.assertIn(f"$('{identifier}')", JS)
        self.assertIn("decision('followup_geometry_review'", JS)
        self.assertIn("{branches,source_ref:sourceRef()}", JS)
        self.assertIn("button.disabled=submitDisabled", JS)
        self.assertIn("template.disabled=loadDisabled", JS)
        self.assertIn("loadDisabled=!auth||ro||!assessment||!fresh", JS)
        self.assertIn("submitDisabled=loadDisabled||!runtimeCurrent", JS)

    def test_geometry_template_uses_registered_dynamic_schema_without_posting(self):
        loader = re.search(
            r"async function loadGeometryTemplate\(\)\{(.+?)\}function validateGeometryBranches",
            JS,
            re.DOTALL,
        )
        self.assertIsNotNone(loader)
        body = loader.group(1)
        self.assertIn("state.assessment?.registered_evidence", JS)
        self.assertIn("priority_groups", JS)
        self.assertIn("groups[e3]?.candidate_id", JS)
        self.assertIn("/api/lab/artifacts/", body)
        self.assertIn("/preview", body)
        self.assertIn("seed_reviews", body)
        self.assertIn("topology_label:''", body)
        self.assertIn("severe_clash_acceptable:null", body)
        self.assertIn("clash_definition:''", body)
        self.assertIn("evidence_ref:found.get", body)
        self.assertNotIn("post(", body)
        self.assertNotRegex(body, r"D-[0-9a-f]{12}")

    def test_geometry_validation_enforces_two_by_five_typed_reviews(self):
        validator = re.search(
            r"function validateGeometryBranches\(value\)\{(.+?)return value\}",
            JS,
            re.DOTALL,
        )
        self.assertIsNotNone(validator)
        body = validator.group(1)
        for token in (
            "value.length!==2",
            "expected.length!==5",
            "branch.seed_reviews.length!==5",
            "seenBranches.has",
            "seenSeeds.has",
            "typeof review.topology_label!=='string'",
            "typeof review.severe_clash_acceptable!=='boolean'",
            "allowedRefs.get(artifactKey(review.evidence_ref))",
            "branch.candidate_id!==groups[branch.e3_type]?.candidate_id",
        ):
            self.assertIn(token, body)

    def test_archived_runtime_blocks_policy_and_compute_but_keeps_template_read_only(self):
        notice = "이전 실행 버전의 저장 결과입니다. 현재 코드로 새 설계를 실행할 수 있습니다. 이 기록의 원본과 판정은 보존합니다."
        self.assertIn(notice, HTML)
        self.assertIn('id="runtimeArchiveNotice"', HTML)
        self.assertIn('href="/design"', HTML)
        self.assertIn("function runtimeCurrentForPolicy()", JS)
        self.assertIn("requireRuntimeCurrentForPolicy();if(!state.assessment?.id)", JS)
        self.assertIn("policyWritable=policyReady&&runtimeCurrent", JS)
        self.assertIn("$('computeButton').disabled=!auth||ro||review||!policyWritable", JS)
        self.assertIn("template.disabled=loadDisabled", JS)
        self.assertIn("button.disabled=submitDisabled", JS)
        self.assertIn("읽기 전용 템플릿", JS)

    def test_branch_cards_use_bound_priority_groups_before_legacy_fallback(self):
        branch = re.search(
            r"function renderBranches\(\)\{(.+?)\}\s*function renderQuestions",
            JS,
            re.DOTALL,
        )
        self.assertIsNotNone(branch)
        body = branch.group(1)
        self.assertIn("currentAssessmentFresh()?state.assessment?.expert_followup?.quantitative_subchecks?.novel:null", body)
        self.assertIn("groups[name]", body)
        self.assertIn("group.candidate_id", body)
        self.assertIn("5-seed coverage", body)
        self.assertIn("group.iptm?.iqr", body)
        self.assertIn("group.linker_endpoint_distance_A?.iqr", body)
        self.assertIn("authenticatedGeometryStatus(novel)", body)
        self.assertIn("branchDiagnostics(name)", body)

    def test_followup_source_hashes_are_collapsed_before_metrics(self):
        self.assertIn("검증된 source · 검토일", JS)
        self.assertIn("검증된 source hash·locator 세부 정보", JS)
        self.assertIn("sourceDetails=el('details'", JS)
        self.assertLess(JS.index("sourceDetails=el('details'"), JS.index("Core interaction map · /body/6"))
        self.assertIn("Math.round((value+Number.EPSILON)*1e6)/1e6", JS)

    def test_no_automatic_geometry_decision_and_packaged_state_instructions(self):
        self.assertIn("버튼을 직접 누른 경우에만", HTML)
        self.assertIn("자동 저장하지 않습니다", JS)
        self.assertIn("패키지 release evidence", HTML)
        self.assertNotIn("구현되지 않은 URL", HTML)
        self.assertIn("같은 topology_label이 최소 4/5 seed", HTML)
        self.assertIn("모든 seed", HTML)
        self.assertIn("quantitative IQR/contact 기준도 통과", HTML)
        self.assertIn("formal approval은 생성되지 않습니다", HTML)
        self.assertIn(".table-scroll", CSS)
        self.assertIn("overflow-x:auto", CSS)


if __name__ == "__main__":
    unittest.main()

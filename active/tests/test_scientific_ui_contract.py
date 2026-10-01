import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
JS_PATH = ROOT / "apps" / "web" / "scientific-review.js"
HTML_PATH = ROOT / "apps" / "web" / "scientific-review.html"
SCHEMA_PATH = ROOT / "contracts" / "drafts" / "scientific_acceptance.schema.json"


class ScientificReviewUiContractTests(unittest.TestCase):
    def setUp(self):
        self.javascript = JS_PATH.read_text(encoding="utf-8")
        self.html = HTML_PATH.read_text(encoding="utf-8")
        self.schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))

    def test_ternary_form_uses_exact_schema_names_and_shapes(self):
        ternary = self.schema["$defs"]["ternaryDecision"]["allOf"][1]["properties"]["data"]
        self.assertEqual(
            ternary["required"],
            ["chosen_candidates", "expected_seeds", "geometry"],
        )
        handler = re.search(
            r"\$\('ternaryForm'\)\.onsubmit=.*?\n\$\('synthesisForm'\)",
            self.javascript,
            re.DOTALL,
        )
        self.assertIsNotNone(handler)
        source = handler.group(0)
        self.assertIn("chosen_candidates:{CRBN:[", source)
        self.assertIn("VHL:[", source)
        self.assertIn("expected_seeds:seeds", source)
        self.assertIn("accepted:true", source)
        self.assertNotIn("chosen_per_e3", source)
        self.assertNotIn("expected_unique_seeds", source)
        self.assertNotIn("disposition:", source)

    def test_contact_options_exactly_match_requirement_schema_enum(self):
        requirement_kinds = self.schema["$defs"]["requirement"]["properties"]["kind"]["enum"]
        add_contact = re.search(
            r"function addContact\(\).*?\$\('contactForm'\)",
            self.javascript,
            re.DOTALL,
        )
        self.assertIsNotNone(add_contact)
        option_pairs = re.findall(
            r"\['([^']*)','([^']+)'\]",
            add_contact.group(0),
        )
        self.assertEqual(option_pairs[0], ("", "선택하세요"))
        offered = dict(option_pairs[1:])
        self.assertEqual(list(offered), requirement_kinds)
        self.assertEqual(
            offered,
            {
                "directional_hbond": "방향성 수소 결합",
                "salt_bridge": "염다리",
                "aromatic": "방향족 상호작용",
                "hydrophobic": "소수성 접촉",
                "proximity": "근접 접촉",
            },
        )

    def test_microstate_editor_accumulates_and_requires_full_selected_coverage(self):
        self.assertIn("microChoices:{}", self.javascript)
        self.assertIn("expectedMicroReports()", self.javascript)
        self.assertIn("state.microChoices[reportId]=", self.javascript)
        self.assertIn("누락된 선택", self.javascript)
        self.assertIn("id=\"microChoices\"", self.html)

    def test_statement_and_artifact_views_are_text_only(self):
        self.assertIn("id=\"statementSource\"", self.html)
        self.assertIn("id=\"artifactFacts\"", self.html)
        self.assertIn("id=\"freshnessFacts\"", self.html)
        self.assertIn("textContent", self.javascript)
        self.assertNotIn("innerHTML", self.javascript)
        self.assertIn("locator:$('statementLocator').value.trim()", self.javascript)

    def test_protein_hydrogen_review_uses_registered_assessment_reference(self):
        self.assertIn('id="proteinReviewForm"', self.html)
        self.assertIn("decision('protein_hydrogen_review'", self.javascript)
        self.assertIn("evidence_ref:evidence.ref", self.javascript)
        self.assertIn("source_ref:sourceRef()", self.javascript)
        self.assertIn("accepted:accepted==='accept'", self.javascript)
        self.assertIn("requires_review", self.javascript)
        self.assertIn("original source flags", self.html)
        self.assertNotIn("evidence_ref:{", self.javascript)

    def test_site_policy_separates_actual_and_explicit_desired_state(self):
        self.assertIn('id="siteState" readonly', self.html)
        self.assertIn('id="siteDesiredState"', self.html)
        self.assertIn('<option value="">선택하지 않음</option>', self.html)
        self.assertIn("stateValue=$('siteDesiredState').value", self.javascript)
        self.assertIn("changed=stateValue!==actualState", self.javascript)
        self.assertIn("original_state:actualState", self.javascript)
        self.assertIn("expert_classification:classification", self.javascript)
        self.assertIn("id=\"expertClassification\"", self.html)
        self.assertIn("새 실험실 사실을 생성하지 않습니다", self.html)
        self.assertIn("release_source_ref=sourceRef()", self.javascript)

    def test_running_job_does_not_request_scientific_policy_or_show_fake_ready_state(self):
        self.assertIn("if(job?.state!=='completed')", self.javascript)
        self.assertIn("Promise.allSettled", self.javascript)
        self.assertIn("policyReady=completed&&!!state.assessment&&!!state.policy", self.javascript)
        self.assertIn('id="assessmentLoadState"', self.html)

    def test_branch_renderer_prefers_current_priority_data_and_keeps_fallback(self):
        self.assertIn("function renderBranches()", self.javascript)
        self.assertIn("currentAssessmentFresh()", self.javascript)
        self.assertIn(
            "expert_followup?.quantitative_subchecks?.novel",
            self.javascript,
        )
        self.assertIn("novel?.priority_groups", self.javascript)
        self.assertIn("['CRBN','VHL'].some(name=>groups[name])", self.javascript)
        self.assertIn("현재 정책·구현에 결합됨", self.javascript)
        self.assertIn(
            "5-seed coverage ${covered}/${expected.length||5}",
            self.javascript,
        )
        self.assertIn(
            "group.status||novel.quantitative_status||'pending'",
            self.javascript,
        )
        self.assertIn("const rows=state.assessment?.criteria||[]", self.javascript)
        self.assertIn("실제 후보 ${count}개", self.javascript)
        self.assertIn("연결된 criterion ${xs.length}개", self.javascript)
        self.assertIn("분기 criterion 미연결", self.javascript)

    def test_policy_rerun_uses_workbench_envelope_and_freshness_gates(self):
        self.assertIn("post('/api/lab/jobs',body)", self.javascript)
        self.assertIn("request_key:crypto.randomUUID()", self.javascript)
        self.assertIn("scientific_policy_id:state.assessment.id", self.javascript)
        self.assertIn("parent_id:state.job.result.parent_scope.actual_design_parent_id", self.javascript)
        self.assertIn("exploratory:false,dock:true", self.javascript)
        self.assertIn("panel_size:12", self.javascript)
        self.assertIn("state.health?.llm_enabled===true", self.javascript)
        self.assertIn("a.source_inputs_current!==true", self.javascript)
        self.assertIn("a.current_implementation!==true", self.javascript)
        self.assertIn("active_policy_decision_ids", self.javascript)
        self.assertIn("['queued','running']", self.javascript)
        self.assertIn("새 결과는 자동 승인되지 않습니다.", self.html)


if __name__ == "__main__":
    unittest.main()

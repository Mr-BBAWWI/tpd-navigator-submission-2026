"""Safety/correctness checks, NOT scientific quality scores for LLMs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import httpx

from packages.agents.bundle import ACTIVE, build_bundle, json_digest
from packages.agents.provider import AgentError, Completion, DaconProvider, load_key
from packages.agents.runtime import ReviewRuntime, validate_claims, validate_response, output_schema


def small_bundle():
    sources = {"paper:F1": {"kind": "article_caption", "text": "PROTAC 1 DC50 is 300 nM."},
               "b:geometry": {"kind": "b_artifact", "content": {"SASA_A2": 21.0}}}
    facts = [{"id": "C01:DC50", "candidate_id": "C01", "domain": "literature",
              "kind": "literature_observation", "value": {"value": 300, "unit": "nM", "review_status": "pending"}, "source_ids": ["paper:F1"]},
             {"id": "START:geometry", "candidate_id": "START", "domain": "molecule",
              "kind": "computed_descriptive_geometry", "value": {"SASA_A2": 21.0}, "source_ids": ["b:geometry"]}]
    b = {"facts": facts, "sources": sources, "coverage": "Synthetic test only", "not_completed": ["human review"]}
    b["input_digest"] = json_digest(b)
    return b


def claims_for(facts, sources, prefix="test"):
    return [{"id": prefix + ":" + f["id"], "candidate_id": f["candidate_id"], "evidence_id": f["id"],
             "kind": f["kind"], "interpretation": "테스트용 해석이며 효능·승인을 뜻하지 않습니다.",
             "citations": [{"source_id": sid, "quote": sources[sid].get("text")} for sid in f["source_ids"]]}
            for f in facts]


class ScriptedProvider:
    mode = "synthetic_test_fixture_no_llm"
    def __init__(self, revision=False, unresolved=False):
        self.critic_calls = 0
        self.revision, self.unresolved = revision, unresolved

    def complete(self, *, instructions, context, **kwargs):
        result = {"input_digest": context["input_digest"]}
        if "Choose the order" in instructions:
            result.update(order=["literature", "molecule"], focus={"literature": "Attribution", "molecule": "Geometry"})
        elif "Independently compare" in instructions:
            self.critic_calls += 1
            problem = self.unresolved or (self.revision and self.critic_calls == 1)
            result.update(verdict="needs_revision" if problem else "no_issues_found",
                          findings=[{"claim_id": context["claims"][0]["id"], "evidence_ids": [context["claims"][0]["evidence_id"]], "problem": "Revise test wording"}] if problem else [])
        elif "Order all validated" in instructions:
            result["claim_order"] = [c["id"] for c in context["claims"]]
        elif not context["read_sources"]:
            result.update(action="read_sources", source_ids=[s["id"] for s in context["source_catalog"]])
        else:
            result.update(action="submit", claims=claims_for(context["facts"], context["read_sources"]), open_questions=[])
        return Completion(json.dumps(result, ensure_ascii=False), {"usage": {"total_tokens": 100}})


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.b = small_bundle()
        self.claims = claims_for(self.b["facts"], self.b["sources"])

    def check(self):
        validate_claims(self.claims, self.b["facts"], self.b["sources"], set(self.b["sources"]))

    def test_cross_candidate_attribution_rejected(self):
        self.claims[0]["candidate_id"] = "C02"
        with self.assertRaisesRegex(AgentError, "CANDIDATE"): self.check()

    def test_computation_cannot_be_experimental(self):
        self.claims[1]["kind"] = "literature_observation"
        with self.assertRaisesRegex(AgentError, "KIND"): self.check()

    def test_invented_quote_rejected(self):
        self.claims[0]["citations"][0]["quote"] = "PROTAC 2 DC50 is 300 nM."
        with self.assertRaisesRegex(AgentError, "QUOTE"): self.check()

    def test_repeated_source_cannot_fill_citation_slots(self):
        self.claims[0]['citations'].append(dict(self.claims[0]['citations'][0]))
        with self.assertRaisesRegex(AgentError, 'DUPLICATE_CITATION_SOURCE'):
            self.check()

    def test_unknown_evidence_rejected(self):
        self.claims[0]["evidence_id"] = "missing"
        with self.assertRaisesRegex(AgentError, "EVIDENCE"): self.check()

    def test_unread_source_rejected(self):
        with self.assertRaisesRegex(AgentError, "UNREAD"):
            validate_claims(self.claims, self.b["facts"], self.b["sources"], set())

    def test_nested_text_does_not_make_structured_record_quoteable(self):
        self.b['sources']['b:geometry']['content']['text'] = 'Nested original sentence.'
        self.check()
        self.claims[1]['citations'][0]['quote'] = 'Nested original sentence.'
        with self.assertRaisesRegex(AgentError, 'STRUCTURED_CITATION'):
            self.check()

    def test_structured_and_original_sources_are_both_required(self):
        self.b['sources']['a:record'] = {'kind':'extraction', 'content':{'value':300}}
        self.b['facts'][0]['source_ids'].append('a:record')
        with self.assertRaisesRegex(AgentError, 'REQUIRED_SOURCE_NOT_CITED'):
            self.check()
        self.claims[0]['citations'].append({'source_id':'a:record','quote':None})
        self.check()

    def test_llm_numeric_field_not_allowed(self):
        self.claims[0]["value"] = 99
        payload = {"input_digest": self.b["input_digest"], "action": "submit", "claims": self.claims, "open_questions": []}
        with self.assertRaisesRegex(AgentError, "SCHEMA"): validate_response(payload, "assessment", self.b["input_digest"])

    def test_stale_response_rejected(self):
        payload = {"input_digest": "0" * 64, "action": "read_sources", "source_ids": ["paper:F1"]}
        with self.assertRaisesRegex(AgentError, "STALE"): validate_response(payload, "assessment", self.b["input_digest"])

    def test_read_action_forbidden_in_submission_phase(self):
        payload = {"input_digest": self.b["input_digest"], "action": "read_sources", "source_ids": ["paper:F1"]}
        with self.assertRaisesRegex(AgentError, "SCHEMA"):
            validate_response(payload, "submission", self.b["input_digest"])

    def test_submission_prompt_schema_excludes_read_action(self):
        schema = output_schema("submission")
        self.assertEqual(schema["properties"]["action"]["const"], "submit")
        self.assertNotIn("read_sources", json.dumps(schema))
        self.assertNotIn("$ref", json.dumps(schema))

    def test_semantic_errors_require_critic_or_human(self):
        # Deliberately demonstrates the boundary; code cannot prove entailment of prose.
        self.claims[1]["interpretation"] = "이 값은 분해 효능을 증명한다."
        self.check()


class RuntimeTests(unittest.TestCase):
    def run_case(self, provider=None, **kwargs):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        runtime = ReviewRuntime(provider or ScriptedProvider(), Path(folder.name) / "run", **kwargs)
        return runtime, small_bundle()

    def test_hierarchy_and_fixture_label(self):
        rt, b = self.run_case()
        result = rt.run(b)
        self.assertEqual(result["status"], "draft_pending_human_review")
        self.assertEqual(result["calls"], 7)
        self.assertEqual(result["execution_mode"], "synthetic_test_fixture_no_llm")
        self.assertIsNone(result["official_approval"])

    def test_revision_has_one_bounded_round(self):
        rt, b = self.run_case(ScriptedProvider(revision=True))
        result = rt.run(b)
        self.assertEqual(result["status"], "draft_pending_human_review")
        self.assertEqual(result["calls"], 10)
        self.assertEqual(result["revision_rounds"], 1)

    def test_unresolved_critique_never_publishes_report(self):
        rt, b = self.run_case(ScriptedProvider(unresolved=True))
        result = rt.run(b)
        self.assertEqual(result["error"], "UNRESOLVED_CRITIC_FINDINGS")
        self.assertFalse((rt.output / "REVIEW.md").exists())

    def test_budget_blocks_before_next_request(self):
        rt, b = self.run_case(max_calls=1)
        result = rt.run(b)
        self.assertEqual(result["error"], "LOCAL_REVIEW_BUDGET_EXHAUSTED")
        self.assertEqual(result["calls"], 1)

    def test_input_changes_block_completion(self):
        rt, b = self.run_case()
        result = rt.run(b, revalidate=lambda: "0" * 64)
        self.assertEqual(result["error"], "STALE_BUNDLE_AT_COMPLETION")
        self.assertFalse((rt.output / "REVIEW.md").exists())

    def test_single_and_ablation_mark_unreviewed(self):
        for architecture, calls in [("single", 2), ("no_critic", 6)]:
            rt, b = self.run_case()
            result = rt.run(b, architecture=architecture)
            self.assertEqual(result["status"], "unreviewed_comparison_draft")
            self.assertEqual(result["calls"], calls)

    def test_altered_bundle_rejected_before_api(self):
        rt, b = self.run_case()
        b["facts"][0]["value"]["value"] = 400
        result = rt.run(b)
        self.assertEqual(result["error"], "BUNDLE_DIGEST_MISMATCH")
        self.assertEqual(result["calls"], 0)

    def test_critic_source_id_is_not_a_fact_id(self):
        class ConfusedCritic(ScriptedProvider):
            def complete(self, **kwargs):
                if "Independently compare" in kwargs['instructions']:
                    context = kwargs['context']
                    data = {"input_digest": context['input_digest'], "verdict": "needs_revision",
                            "findings": [{"claim_id": context['claims'][0]['id'],
                                          "evidence_ids": ["paper:F1"], "problem": "test"}]}
                    return Completion(json.dumps(data), {"usage": {"total_tokens": 100}})
                return super().complete(**kwargs)
        rt, b = self.run_case(ConfusedCritic())
        result = rt.run(b)
        self.assertEqual(result["error"], "CRITIC_UNKNOWN_REFERENCE")
        self.assertFalse((rt.output / "REVIEW.md").exists())


class ProviderTests(unittest.TestCase):
    def provider(self, handler):
        provider = DaconProvider("SYNTHETIC-SECRET-TEST", httpx.MockTransport(handler))
        self.addCleanup(provider.close)
        return provider

    def invoke(self, provider):
        return provider.complete(model="gpt-5.6-sol", instructions="JSON only", context={}, max_output_tokens=256)

    def test_auth_and_redaction_and_usage(self):
        def handler(request):
            self.assertEqual(request.headers["api-key"], "SYNTHETIC-SECRET-TEST")
            self.assertEqual(request.url.path, "/hackathon/openai/v1/responses")
            self.assertFalse(json.loads(request.content)["store"])
            return httpx.Response(200, json={"status": "completed", "output": [{"type": "message", "role": "assistant", "content": [{"type": "output_text", "text": "SYNTHETIC-SECRET-TEST"}]}], "usage": {"total_tokens": 42}}, headers={"x-team-tokens-consumed": "42"})
        result = self.invoke(self.provider(handler))
        self.assertEqual(result.text, "[REDACTED]")
        self.assertEqual(result.metadata["usage"]["total_tokens"], 42)

    def test_http_failures_do_not_echo_body_or_retry(self):
        for status in (401, 403, 429, 500, 307):
            calls = []
            def handler(request):
                calls.append(request)
                return httpx.Response(status, text="SYNTHETIC-SECRET-TEST")
            with self.assertRaises(AgentError) as ctx: self.invoke(self.provider(handler))
            self.assertNotIn("SYNTHETIC-SECRET", str(ctx.exception))
            self.assertEqual(len(calls), 1)

    def test_refusal_and_incomplete(self):
        for data in ({"status": "incomplete"}, {"status": "completed", "output": [{"type": "message", "role": "assistant", "content": [{"type": "refusal"}]}]}):
            with self.assertRaises(AgentError): self.invoke(self.provider(lambda request: httpx.Response(200, json=data)))

    @patch.dict("os.environ", {}, clear=True)
    def test_markdown_key_and_placeholder(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "key.md"
            path.write_text("```dacon-api-key\nPASTE_YOUR_DACON_API_KEY_HERE\n```", encoding="utf-8")
            with self.assertRaises(AgentError): load_key(path)
            path.write_text("```dacon-api-key\nSYNTHETIC-SECRET-TEST\n```", encoding="utf-8")
            self.assertEqual(load_key(path), "SYNTHETIC-SECRET-TEST")




if __name__ == "__main__":
    unittest.main()

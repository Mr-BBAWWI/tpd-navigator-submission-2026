"""Synthetic source/LLM contract and concurrency tests, never scientific validation."""
import copy
import json
import tempfile
import unittest
import uuid
from concurrent.futures import ThreadPoolExecutor

from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.agents.provider import Completion, AgentError
from packages.contracts import missing, ContractError
from packages.platform.analysis import AnalysisService
from packages.platform.orchestrator import Orchestrator
from packages.platform.store import Store
from packages.review_contracts import ReviewReader
from packages.literature.assessment import check_literals, document_label
from test_i2 import FakeTransport


class ScriptedLiterature:
    mode = "synthetic_test_fixture_no_llm"

    def __init__(self, mutation=None, before_return=None):
        self.calls = 0
        self.mutation, self.before_return = mutation, before_return

    def close(self):
        pass

    def complete(self, *, context, **kwargs):
        self.calls += 1
        if "claims" in context:
            result = {"input_digest": context["input_digest"], "reviewed_claim_ids": [c["claim_id"] for c in context["claims"]],
                      "findings": [], "limitations": ["Synthetic critique only"]}
        else:
            record = {"basis_kind": "author_interpretation", "statement": "가상 문헌은 실험 주장을 하지 않는다.",
                      "segment_id": context["sources"][0]["segment_id"], "quote": context["sources"][0]["text"],
                      "compound_label": None, "reported_target": None,
                      "assay_context": {"assay_id": missing(), "assay_type": missing(), "cell_line": missing(),
                          "tissue": missing(), "duration": {"kind": "missing", **missing()}, "conditions": []},
                      "observations": [], "limitations": ["Synthetic only"]}
            result = {"input_digest": context["input_digest"], "records": [record], "limitations": []}
        if self.mutation:
            self.mutation(result)
        if self.before_return:
            self.before_return()
        return Completion(json.dumps(result), {"usage": {"total_tokens": 100}, "returned_model": "scripted-fixture"})


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(self.temp.name)
        request = {"query": "TESTGENE", "query_kind": "gene_symbol", "objective": "", "cell_line": "", "pdf_artifact_ids": [], "max_documents": 6}
        self.run, _ = self.store.create_run(PROJECT, uuid.uuid4().hex, request)
        Orchestrator(self.store, PROJECT, FakeTransport).execute(self.run["id"])
        self.run = self.store.update(PROJECT, self.run["id"], {"data_mode": "test_fixture"})
        self.port = self.store.scope(PROJECT)
        self.bundle = self.port.json(self.run["bundle_ref"])
        self.ids = [s["segment_id"] for s in self.bundle["segments"][:1]]
        self.provider = ScriptedLiterature()
        self.service = AnalysisService(self.store, PROJECT, lambda: self.provider)

    def start(self, key="synthetic-request"):
        return self.service.create(self.run["id"], self.run["revision"], self.ids, key)[0]

    def execute(self):
        job = self.start()
        self.service.execute(job["id"])
        return self.service.get(job["id"])

    def test_complete_graph_read_history_and_no_approval(self):
        result = self.execute()
        self.assertEqual(result["state"], "review_ready", result)
        view = self.service.view(result["id"])
        self.assertEqual(view["annotated_bundle"]["coverage"]["llm_read_segment_count"], 1)
        self.assertEqual(view["annotated_bundle"]["evidence_records"][0]["review"]["state"], "unreviewed")
        self.assertEqual(view["assessment"]["proposals"], [])
        self.assertEqual(self.port.json(self.run["bundle_ref"]), self.bundle)
        self.assertFalse(self.service.g1_status(self.run["id"])["ready"])
        self.assertEqual(result["total_tokens"], 200)
        for ref in (result["assessment_ref"], result["review_ref"], result["code_check_ref"]):
            self.assertTrue(ReviewReader(self.port.read).verify(ref)["contract_valid"])

    def test_pdb_identifier_is_not_assay_identifier(self):
        record={'compound_label':None,'reported_target':None,'assay_context':{'assay_id':{'state':'known','value':'6HAX','reason':None}},'observations':[]}
        with self.assertRaisesRegex(AgentError,'STRUCTURE_ID_IS_NOT_ASSAY_ID'):
            check_literals(record,'PDB ID: 6HAX')

    def test_document_number_alias_requires_explicit_unique_source_pair(self):
        seg={'s1':{'locator':{'document_id':'doc1'}},'s2':{'locator':{'document_id':'doc2'}}}
        records=[{'compound_label':'Drug Z','segment_id':'s1'}]
        sources=[{'segment_id':'s1','text':'Drug Z ( 3 ) was studied.'}]
        name,note=document_label('Drug Z ( 3 )','doc1',records,seg,sources)
        self.assertEqual(name,'Drug Z');self.assertIsNotNone(note)
        sources[0]['text']+=' Drug Z ( 4 ) also appears.'
        self.assertEqual(document_label('Drug Z ( 3 )','doc1',records,seg,sources)[0],'Drug Z ( 3 )')
        self.assertEqual(document_label('Drug Z ( 3 )','doc2',records,seg,sources)[0],'Drug Z ( 3 )')

    def test_research_context_change_invalidates_saved_analysis(self):
        result=self.execute()
        query=self.port.json(self.run['target_query_ref'])
        query['research_context']['objective']=missing('Synthetic changed context','not_provided')
        ref=self.port.put_json(query,'TargetQuery','human_annotation')
        self.store.update(PROJECT,self.run['id'],{'target_query_ref':ref})
        self.assertFalse(self.service.view(result['id'])['current_input'])

    def test_same_request_and_duplicate_workers_call_only_once(self):
        job = self.start()
        again, fresh = self.service.create(self.run["id"], 1, self.ids, "synthetic-request")
        self.assertFalse(fresh)
        self.assertEqual(job["id"], again["id"])
        with ThreadPoolExecutor(2) as pool:
            list(pool.map(self.service.execute, [job["id"], job["id"]]))
        self.assertEqual(self.provider.calls, 2)

    def test_same_document_label_shared_without_claiming_structure_identity(self):
        def mutation(v):
            if "records" in v:
                v["records"][0]["compound_label"] = "Synthetic"
                v["records"].append(copy.deepcopy(v["records"][0]))
        self.provider.mutation = mutation
        result = self.execute()
        self.assertEqual(result["state"], "review_ready", result)
        bundle = self.service.view(result["id"])["annotated_bundle"]
        self.assertEqual(len(bundle["compounds"]), 1)
        self.assertEqual(bundle["evidence_records"][0]["compound_ids"], bundle["evidence_records"][1]["compound_ids"])
        self.assertEqual(bundle["compounds"][0]["structure_state"], "unknown")

    def test_concurrent_distinct_requests_admit_only_one(self):
        def create(i):
            try:
                return self.start("request-" + str(i))["id"]
            except ValueError:
                return None
        with ThreadPoolExecutor(2) as pool:
            results = list(pool.map(create, [1, 2]))
        self.assertEqual(sum(r is not None for r in results), 1)

    def test_wrong_quote_keeps_rejected_answer(self):
        self.provider.mutation = lambda v: v["records"][0].update(quote="invented observation")
        result = self.execute()
        self.assertEqual(result["state"], "failed")
        self.assertEqual(self.provider.calls, 1)
        self.assertEqual(len(result["trace"]), 1)
        self.assertIn("invented observation", self.port.json(result["trace"][0]["response_ref"])["text"])

    def test_unread_segment_rejected(self):
        self.provider.mutation = lambda v: v["records"][0].update(segment_id=self.bundle["segments"][1]["segment_id"])
        self.assertEqual(self.execute()["state"], "failed")

    def test_unquoted_compound_and_self_approval_rejected(self):
        for mutation in (lambda v: v["records"][0].update(compound_label="Not in source"),
                         lambda v: v.update(approved=True)):
            self.provider = ScriptedLiterature(mutation)
            job = self.start(uuid.uuid4().hex)
            self.service.execute(job["id"])
            self.assertEqual(self.service.get(job["id"])["state"], "failed")

    def test_blocking_critic_preserved_and_held(self):
        def mutation(v):
            if "findings" in v:
                v["findings"] = [{"severity": "blocking", "message": "Synthetic unsupported statement", "claim_ids": v["reviewed_claim_ids"]}]
        self.provider.mutation = mutation
        result = self.execute()
        self.assertEqual(result["state"], "held")
        self.assertEqual(self.service.view(result["id"])["review"]["findings"][0]["severity"], "blocking")

    def test_critic_must_review_all_claims(self):
        def mutation(v):
            if "reviewed_claim_ids" in v:
                v["reviewed_claim_ids"] = []
        self.provider.mutation = mutation
        result = self.execute()
        self.assertEqual(result["state"], "failed")
        self.assertEqual(len(result["trace"]), 2)

    def test_input_changed_between_calls_prevents_critic(self):
        self.provider.before_return = lambda: self.store.update(PROJECT, self.run["id"], {"revision": 2})
        result = self.execute()
        self.assertEqual(result["state"], "stale")
        self.assertEqual(self.provider.calls, 1)

    def test_changed_bundle_same_revision_invalidates_saved_result(self):
        result = self.execute()
        new_ref = self.port.put_json(copy.deepcopy(self.bundle), "EvidenceBundle")
        self.store.update(PROJECT, self.run["id"], {"bundle_ref": new_ref})
        self.assertFalse(self.service.view(result["id"])["current_input"])
        self.assertIn("CURRENT_LITERATURE_REVIEW_NOT_READY", self.service.g1_status(self.run["id"])["blocking_reasons"])

    def test_tamper_detected_on_replay_without_provider_call(self):
        result = self.execute()
        ref = self.bundle["segments"][0]["text_artifact_ref"]
        (self.store.root / "blobs" / ref["artifact_id"]).write_bytes(b"tampered")
        with self.assertRaises(ContractError):
            self.service.view(result["id"])
        self.assertEqual(self.provider.calls, 2)

    def test_restart_does_not_retry_and_preserves_partial_trace(self):
        job = self.start()
        self.service.update(job["id"], {"state": "running", "pending_prompt_ref": self.run["bundle_ref"]})
        self.service.interrupt_pending()
        self.service.execute(job["id"])
        result = self.service.get(job["id"])
        self.assertEqual(result["state"], "interrupted")
        self.assertEqual(result["usage_status"], "unreconciled")
        self.assertEqual(self.provider.calls, 0)

    def test_second_app_startup_does_not_interrupt_live_worker(self):
        self.provider.before_return = lambda: AnalysisService(self.store, PROJECT).interrupt_pending()
        result = self.execute()
        self.assertEqual(result["state"], "review_ready", result)
        self.assertEqual(self.provider.calls, 2)

    def test_fixture_provider_cannot_masquerade_as_real(self):
        self.run = self.store.update(PROJECT, self.run["id"], {"data_mode": "real"})
        self.assertEqual(self.execute()["error_code"], "FIXTURE_PROVIDER_IN_REAL_ANALYSIS")
        self.assertEqual(self.provider.calls, 0)

    def test_api_is_opt_in_and_gets_never_invoke_model(self):
        app = create_app(self.temp.name, enable_worker=False)
        client = TestClient(app)
        self.assertFalse(client.get("/api/health").json()["llm_enabled"])
        body = {"revision": 1, "segment_ids": self.ids, "request_key": "synthetic-request"}
        self.assertEqual(client.post(f"/api/runs/{self.run['id']}/analyses", json=body, headers={"x-tpd-local": "1"}).status_code, 409)
        job = self.execute()
        for _ in range(2):
            self.assertEqual(client.get(f"/api/analyses/{job['id']}").status_code, 200)
            self.assertEqual(client.get(f"/api/runs/{self.run['id']}").status_code, 200)
        self.assertEqual(self.provider.calls, 2)
        other = AnalysisService(self.store, "other-project", lambda: self.provider)
        with self.assertRaises(KeyError):
            other.get(job["id"])

    def test_api_rejects_forged_authority_and_requires_local_header(self):
        app = create_app(self.temp.name, enable_worker=False, provider_factory=lambda: self.provider)
        client = TestClient(app)
        url = f"/api/runs/{self.run['id']}/analyses"
        body = {"revision": 1, "segment_ids": self.ids, "request_key": "synthetic-request"}
        self.assertEqual(client.post(url, json=body).status_code, 403)
        self.assertEqual(client.post(url, json={**body, "actor_name": "fake"}, headers={"x-tpd-local": "1"}).status_code, 422)
        self.assertEqual(client.post(url, json=body, headers={"x-tpd-local": "1"}).status_code, 202)
        self.assertEqual(self.provider.calls, 0)

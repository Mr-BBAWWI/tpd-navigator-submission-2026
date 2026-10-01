"""Scripted harness faults, not model accuracy measurements."""
import json
import tempfile
import unittest
from pathlib import Path

from packages.agents.provider import Completion,AgentError
from packages.agents.review_harness import ClaimReviewHarness


class Provider:
    mode='synthetic_test_fixture_no_llm'
    def __init__(self,verdicts=('supported','supported'),mutation=None):
        self.prompts=[];self.verdicts=verdicts;self.mutation=mutation
    def complete(self,context,**kwargs):
        self.prompts.append(context)
        value={'input_digest':context['input_digest'],'verdict':self.verdicts[len(self.prompts)-1],
            'evidence_ids':['f1'],'rationale':'Synthetic explanation','missing_fields':[]}
        if self.mutation:self.mutation(value)
        return Completion(json.dumps(value),{'usage':{'total_tokens':100},'returned_model':'fixture'})


class HarnessTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.input={'proposition':'Synthetic claim','evidence':[{'id':'f1','value':'Synthetic evidence'}],'context_note':''}
    def make(self,provider,**kw):
        return ClaimReviewHarness(provider,Path(self.temp.name)/'run',**kw)
    def test_blind_pair_never_passes_first_answer_or_gold(self):
        p=Provider();r=self.make(p).run('blind_pair',self.input)
        self.assertEqual(r['status'],'reviewed_pending_human');self.assertEqual(p.prompts[0],p.prompts[1])
        self.assertNotIn('expected',p.prompts[0]);self.assertFalse(r['human_approved']);self.assertFalse(r['dispatch_authorized'])
        self.assertEqual(r['calls'],2)
    def test_disagreement_held_without_retry(self):
        p=Provider(('supported','contradicted'));r=self.make(p).run('blind_pair',self.input)
        self.assertEqual(r['error_code'],'HARNESS_REVIEW_DISAGREEMENT');self.assertIsNone(r['verdict']);self.assertEqual(r['calls'],2)
    def test_stale_answer_kept_but_not_used(self):
        p=Provider(mutation=lambda v:v.update(input_digest='0'*64));r=self.make(p).run('direct',self.input)
        self.assertEqual(r['error_code'],'HARNESS_STALE_OUTPUT');self.assertEqual(r['calls'],1)
        self.assertTrue((Path(self.temp.name)/'run/0-answer.json').is_file())
    def test_unread_evidence_rejected(self):
        p=Provider(mutation=lambda v:v.update(evidence_ids=['secret-unread']));r=self.make(p).run('direct',self.input)
        self.assertEqual(r['error_code'],'HARNESS_UNDELIVERED_EVIDENCE')
    def test_model_self_approval_rejected(self):
        p=Provider(mutation=lambda v:v.update(human_approved=True));r=self.make(p).run('direct',self.input)
        self.assertEqual(r['error_code'],'HARNESS_OUTPUT_SCHEMA')
    def test_budget_stops_before_provider_call(self):
        p=Provider();r=self.make(p,max_tokens=1).run('direct',self.input)
        self.assertEqual(r['calls'],0);self.assertFalse(p.prompts)
    def test_gold_cannot_enter_runtime_input(self):
        p=Provider()
        with self.assertRaisesRegex(AgentError,'INPUT_FIELDS'):self.make(p).run('direct',{**self.input,'expected':'supported'})
        self.assertFalse(p.prompts)
    def test_source_manifest_has_no_duplicate_ids(self):
        p=Provider();self.input['evidence']*=2
        with self.assertRaisesRegex(AgentError,'EVIDENCE_INVALID'):self.make(p).run('direct',self.input)

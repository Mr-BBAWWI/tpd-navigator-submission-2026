"""Versioned, read-only candidate-claim harness shared by local evaluation runners.

The supervisor is code. Blind reviewers receive the same sealed evidence, never
one another's verdict or the evaluator's answer. Disagreement stops publication;
there is no majority vote, approval, filesystem/tool execution, or automatic retry.
"""
import json
import time
from pathlib import Path

from jsonschema import Draft202012Validator

from packages.contracts import encoded, parse_json
from packages.platform.dossiers import digest
from .provider import AgentError

ROOT = Path(__file__).resolve().parents[2]
PROFILE_PATH = ROOT/'cases/harness/review_profiles.json'
SCHEMA_PATH = ROOT/'contracts/drafts/claim_review.schema.json'


class ClaimReviewHarness:
    def __init__(self, provider, output, model='gpt-5.6-sol', max_tokens=40000):
        self.provider, self.output, self.model = provider, Path(output), model
        self.output.mkdir(parents=True,exist_ok=False)
        self.profiles = parse_json(PROFILE_PATH.read_bytes())
        self.schema = parse_json(SCHEMA_PATH.read_bytes())
        self.max_tokens, self.tokens, self.calls = max_tokens, 0, 0
        self.trace = []
        self.unreconciled = False

    def save(self, name, value):
        with (self.output/name).open('xb') as handle:
            handle.write(encoded(value))

    def call(self, role, instruction, context):
        prompt = {**context,'output_schema':self.schema}
        # UTF-8 byte count plus output allowance is a conservative admission estimate,
        # not provider billing. Every completion's actual usage is separately retained.
        reserve = len(encoded(prompt))+len(instruction.encode())+2500+2048
        if self.calls>=2 or self.tokens+reserve>self.max_tokens:
            raise AgentError('HARNESS_BUDGET_LIMIT')
        index = self.calls
        self.save(f'{index}-prompt.json',{'role':role,'instructions':instruction,'context':prompt,
            'profile_version':self.profiles['version'],'model':self.model,'reservation_estimate':reserve})
        self.calls += 1
        self.unreconciled = True
        result = self.provider.complete(model=self.model,instructions=instruction,context=prompt,max_output_tokens=2500)
        self.save(f'{index}-answer.json',{'role':role,'text':result.text,'metadata':result.metadata})
        n = result.metadata.get('usage',{}).get('total_tokens')
        if type(n) is not int or n<=0:
            raise AgentError('HARNESS_USAGE_UNKNOWN')
        self.tokens += n
        self.unreconciled = False
        self.trace.append({'role':role,'prompt':f'{index}-prompt.json','answer':f'{index}-answer.json',
            'tokens':n,'returned_model':result.metadata.get('returned_model')})
        try:
            value = parse_json(result.text.encode())
            if list(Draft202012Validator(self.schema).iter_errors(value)):
                raise ValueError()
        except (ValueError,TypeError):
            raise AgentError('HARNESS_OUTPUT_SCHEMA') from None
        if value['input_digest'] != context['input_digest']:
            raise AgentError('HARNESS_STALE_OUTPUT')
        available = {f['id'] for f in context['evidence']}
        if not set(value['evidence_ids']) <= available:
            raise AgentError('HARNESS_UNDELIVERED_EVIDENCE')
        if self.tokens>self.max_tokens:
            raise AgentError('HARNESS_BUDGET_EXHAUSTED')
        return value

    def run(self, profile, input_data):
        if profile not in self.profiles['profiles']:
            raise AgentError('HARNESS_PROFILE_UNKNOWN')
        if set(input_data) != {'proposition','evidence','context_note'}:
            raise AgentError('HARNESS_INPUT_FIELDS')  # Evaluator labels cannot enter prompts.
        ids = [f['id'] for f in input_data['evidence']]
        if not ids or len(ids)!=len(set(ids)) or len(encoded(input_data))>20000:
            raise AgentError('HARNESS_EVIDENCE_INVALID')
        context = {**input_data,'input_digest':digest(encoded(input_data))}
        self.save('input.json',context)
        spec = self.profiles['profiles'][profile]
        result = {'version':self.profiles['version'],'profile':profile,'input_digest':context['input_digest'],
            'profile_sha256':digest(PROFILE_PATH.read_bytes()),'schema_sha256':digest(SCHEMA_PATH.read_bytes()),
            'mode':self.provider.mode,'model':self.model,'status':'held','verdict':None,'reviews':[],
            'human_approved':False,'dispatch_authorized':False}
        started=time.monotonic()
        try:
            for role in spec['roles']:
                instruction=self.profiles['common_instruction']+'\n'+self.profiles['role_instructions'][role]
                # Each role starts with the original evidence only. No reasoning transcript relay.
                result['reviews'].append(self.call(role,instruction,context))
            votes={r['verdict'] for r in result['reviews']}
            if len(votes)==1:
                result.update(status='reviewed_pending_human',verdict=votes.pop())
            else:
                result['error_code']='HARNESS_REVIEW_DISAGREEMENT'
        except AgentError as exc:
            result['error_code']=str(exc)
        except Exception:
            result['error_code']='HARNESS_EXECUTION_FAILED'
        result.update(calls=self.calls,total_tokens=self.tokens,usage_status='unreconciled' if self.unreconciled else 'observed',
            elapsed_seconds=round(time.monotonic()-started,3),trace=self.trace)
        self.save('result.json',result)
        return result

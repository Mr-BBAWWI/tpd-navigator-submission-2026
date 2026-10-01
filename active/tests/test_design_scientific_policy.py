"""Synthetic recorded-policy and call-journal tests; fixtures are not human approval."""
import json
from pathlib import Path
import tempfile
import unittest

from rdkit import Chem

from packages.agents.provider import AgentError, Completion
from packages.platform.store import Store
from packages.platform.workbench import CallJournal, WorkbenchService, LIMITS
from packages.science.analog_generation import generate
from packages.science.dual_e3 import apply_scientific_policy, validate_scientific_policy


class ScientificPolicyTests(unittest.TestCase):
    @staticmethod
    def artifact_ref(identifier='decision-1'):
        return {
            'artifact_id': identifier,
            'version': 1,
            'sha256': 'a' * 64,
            'media_type': 'application/json',
            'schema_id': 'scientific-decision/v1',
            'provenance': 'source',
        }

    def policy(self, parent_id='PARENT', site_policy=None):
        return {
            'scope': {
                'project_id': 'synthetic-project',
                'job_id': 'recorded-assessment-job',
                'parent_id': parent_id,
                'assessment_id': 'assessment-1',
                'policy_revision': 3,
                'policy_digest': 'sha256:synthetic-policy',
            },
            'site_policy': site_policy or [],
            'parent_funnel': [{'id': 'PARENT', 'state': 'reviewed_fixture'}],
        }

    def test_scope_mismatch_is_rejected(self):
        policy = self.policy()
        with self.assertRaisesRegex(ValueError, 'PARENT_SCOPE_MISMATCH'):
            validate_scientific_policy(
                policy,
                project_id='synthetic-project',
                policy_id='assessment-1',
                parent_id='OTHER',
            )
        with self.assertRaisesRegex(ValueError, 'PROJECT_SCOPE_MISMATCH'):
            validate_scientific_policy(policy, project_id='other-project')

    def test_protected_site_requires_explicit_release_and_justification(self):
        sites = [{
            'atom_map': 1,
            'state': 'PROTECTED',
            'protected': True,
            'evidence': {'SAR': {}},
        }]
        source_ref = self.artifact_ref()
        declaration = {
            'atom_map': 1,
            'state': 'MODIFIABLE',
            'allowed_transforms': ['LH_OH'],
            'rationale': '',
            'source_ref': source_ref,
            'allow_release_protected': True,
        }
        retained = apply_scientific_policy(
            sites,
            self.policy(site_policy=[declaration]),
            [1],
        )[0]
        self.assertEqual(retained['original_state'], 'PROTECTED')
        self.assertEqual(retained['state'], 'PROTECTED')
        declaration.update(rationale='recorded expert rationale')
        released = apply_scientific_policy(
            sites,
            self.policy(site_policy=[declaration]),
            [1],
        )[0]
        self.assertEqual(released['state'], 'MODIFIABLE')
        self.assertEqual(
            released['evidence']['scientific_policy']['source_ref'],
            source_ref,
        )

    def test_exact_rule_and_class_permissions_are_enforced(self):
        mol = Chem.MolFromSmiles('[CH3:1]')
        exact = {
            'atom_map': 1,
            'state': 'MODIFIABLE',
            'allowed_transforms': ['LH_OH'],
            'rationale': 'synthetic exact-rule review',
            'source_ref': self.artifact_ref('exact'),
            'allow_release_protected': False,
        }
        sites = apply_scientific_policy(
            [{
                'atom_map': 1,
                'state': 'UNKNOWN',
                'protected': False,
                'evidence': {'SAR': {}},
            }],
            self.policy(site_policy=[exact]),
        )
        result = generate(mol, sites, False)
        self.assertTrue(any(
            row['rule_id'] == 'LH_OH' for row in result['analogs']
        ))
        self.assertFalse(any(
            row['rule_id'] == 'LH_NH2' for row in result['analogs']
        ))

        class_rule = dict(
            exact,
            allowed_transforms=['linker_handle_introduction'],
        )
        class_sites = apply_scientific_policy(
            [{
                'atom_map': 1,
                'state': 'UNKNOWN',
                'protected': False,
                'evidence': {'SAR': {}},
            }],
            self.policy(site_policy=[class_rule]),
        )
        class_result = generate(mol, class_sites, False)
        self.assertTrue(any(
            row['rule_id'] == 'LH_NH2' for row in class_result['analogs']
        ))

    def test_unknown_outside_policy_stays_exploratory_and_pending(self):
        mol = Chem.MolFromSmiles('[CH3:1]')
        sites = apply_scientific_policy(
            [{
                'atom_map': 1,
                'state': 'UNKNOWN',
                'protected': False,
                'evidence': {'SAR': {}},
            }],
            self.policy(),
        )
        self.assertTrue(sites[0]['scientific_pending'])
        self.assertEqual(generate(mol, sites, False)['analogs'], [])
        exploratory = generate(mol, sites, True)['analogs']
        self.assertTrue(exploratory)
        self.assertTrue(all(
            row['scientific_pending'] for row in exploratory
        ))


class CallJournalPolicyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.service = WorkbenchService(
            Store(Path(self.temp.name)),
            'synthetic-project',
        )
        self.job_id = 'job-synthetic-journal'
        row = {
            'id': self.job_id,
            'project_id': 'synthetic-project',
            'run_id': 'run',
            'result_id': 'design:SMARCA2',
            'operation': 'review',
            'title': 'synthetic',
            'parameters': {},
            'input_digest': 'input',
            'binding': {},
            'runtime': {},
            'state': 'running',
            'stage': 'running',
            'tool_inputs': [],
            'created_at': 'synthetic',
            'updated_at': 'synthetic',
            'retry_of': None,
            'attempt': 1,
            'calls': 0,
            'total_tokens': 0,
            'usage_status': 'not_started',
            'limits': dict(LIMITS),
            'outputs': [],
            'error_code': None,
            'authority': {},
            'result': None,
        }
        with self.service.store.db() as db:
            db.execute(
                'INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)',
                (
                    self.job_id,
                    'synthetic-project',
                    'run',
                    'result',
                    'request-key',
                    'fingerprint',
                    'running',
                    json.dumps(row),
                ),
            )

        def active(identifier):
            if self.service.get(identifier)['state'] != 'running':
                raise AgentError('JOB_CANCELLED_OR_NOT_RUNNING')

        self.service.check_active = active

    def viewed(self, identifier):
        row = self.service.get(identifier)
        with self.service.store.db() as db:
            row['calls_detail'] = [
                json.loads(record[0])
                for record in db.execute(
                    'SELECT body FROM workbench_calls '
                    'WHERE job_id=? ORDER BY ordinal',
                    (identifier,),
                )
            ]
        for call in row['calls_detail']:
            for key in ('prompt_ref', 'answer_ref'):
                if call.get(key):
                    self.service.port.read(call[key])
        return row

    @staticmethod
    def error(code, metadata=None):
        error = AgentError(code)
        if metadata is not None:
            error.metadata = metadata
        return error

    def test_failure_known_total_records_sanitized_blob_and_total(self):
        error = self.error(
            'API_PROVIDER_ERROR',
            {
                'usage': {'total_tokens': 17},
                'provider_status': 'rejected',
            },
        )

        class Provider:
            mode = 'synthetic'

            def complete(self, **request):
                raise error

        with self.assertRaises(AgentError):
            CallJournal(
                self.service,
                self.job_id,
                Provider(),
            ).complete(model='synthetic')

        viewed = self.viewed(self.job_id)
        call = viewed['calls_detail'][0]
        self.assertEqual(viewed['total_tokens'], 17)
        self.assertEqual(call['usage_status'], 'observed')
        self.assertIsNotNone(call['answer_ref'])
        blob = self.service.port.read(call['answer_ref'])
        answer = json.loads(blob)
        self.assertEqual(answer['text'], '')
        self.assertNotIn('message', answer['metadata'])
        self.assertNotIn('key', answer['metadata'])

    def test_incomplete_usage_is_partial_but_network_unknown_unreconciled(self):
        partial = self.error(
            'API_PROVIDER_ERROR',
            {'usage': {'input_tokens': 9}},
        )

        class Partial:
            mode = 'synthetic'

            def complete(self, **request):
                raise partial

        with self.assertRaises(AgentError):
            CallJournal(
                self.service,
                self.job_id,
                Partial(),
            ).complete(model='synthetic')

        viewed = self.viewed(self.job_id)
        self.assertEqual(viewed['usage_status'], 'observed_partial')
        self.assertEqual(viewed['total_tokens'], 0)

        second_id = 'job-network-unknown'
        row = self.service.get(self.job_id)
        row.update(id=second_id, calls=0, usage_status='not_started')
        with self.service.store.db() as db:
            db.execute(
                'INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)',
                (
                    second_id,
                    'synthetic-project',
                    'run',
                    'result',
                    'request-key-2',
                    'fingerprint-2',
                    'running',
                    json.dumps(row),
                ),
            )

        unknown = self.error('API_TRANSPORT_ERROR_USAGE_UNKNOWN')

        class Unknown:
            mode = 'synthetic'

            def complete(self, **request):
                raise unknown

        with self.assertRaises(AgentError):
            CallJournal(
                self.service,
                second_id,
                Unknown(),
            ).complete(model='synthetic')

        viewed = self.viewed(second_id)
        self.assertEqual(viewed['usage_status'], 'unreconciled')
        self.assertIsNotNone(viewed['calls_detail'][0]['answer_ref'])

    def test_cancellation_after_response_preserves_usage_without_retry(self):
        service = self.service
        calls = []

        class Provider:
            mode = 'synthetic'

            def complete(self, **request):
                calls.append(request)
                service.cancel(service_job)
                return Completion(
                    'synthetic response',
                    {'usage': {'total_tokens': 23}},
                )

        service_job = self.job_id
        with self.assertRaisesRegex(AgentError, 'CANCELLED'):
            CallJournal(
                service,
                service_job,
                Provider(),
            ).complete(model='synthetic')

        viewed = self.viewed(service_job)
        self.assertEqual(len(calls), 1)
        self.assertEqual(viewed['total_tokens'], 23)
        self.assertEqual(viewed['usage_status'], 'observed')
        self.assertIsNotNone(viewed['calls_detail'][0]['answer_ref'])


if __name__ == '__main__':
    unittest.main()

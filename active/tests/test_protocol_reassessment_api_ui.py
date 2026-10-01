"""Protocol reassessment API and UI integration contracts."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fastapi.testclient import TestClient

from apps.api.main import create_app

ROOT = Path(__file__).resolve().parents[1]
HTML = (ROOT / 'apps/web/scientific-review.html').read_text(encoding='utf-8')
JS = (ROOT / 'apps/web/scientific-review.js').read_text(encoding='utf-8')


class ProtocolReassessmentApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def client(self, **kwargs):
        return TestClient(create_app(
            data_root=self.tmp.name, enable_worker=False, **kwargs))

    def authenticated_post(self, client, body):
        identities = client.app.state.report_reviews.identities
        with patch.object(
            identities, 'session', return_value={'id': 'reviewer-server'}
        ):
            return client.post(
                '/api/lab/jobs/job-synthetic/scientific-assessments',
                headers={'X-TPD-Local': '1'}, json=body,
            )

    def test_explicit_protocol_ids_and_digest_are_forwarded_with_server_actor(self):
        service_result = ({'id': 'assessment-new'}, True)
        with self.client() as client, patch.object(
            client.app.state.scientific_acceptance,
            'create', return_value=service_result,
        ) as create:
            response = self.authenticated_post(client, {
                'expected_result_sha256': 'a' * 64,
                'protocol_ids': ['known-protocol', 'novel-protocol'],
            })
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['assessment_id'], 'assessment-new')
        create.assert_called_once_with(
            'job-synthetic', 'a' * 64,
            reviewer_id='reviewer-server',
            protocol_ids=['known-protocol', 'novel-protocol'],
        )

    def test_omitted_and_empty_protocol_ids_preserve_legacy_call(self):
        for body in ({}, {'protocol_ids': []}):
            with self.subTest(body=body), self.client() as client, patch.object(
                client.app.state.scientific_acceptance,
                'create', return_value=({'id': 'assessment-legacy'}, True),
            ) as create:
                response = self.authenticated_post(client, body)
                self.assertEqual(response.status_code, 200)
                create.assert_called_once_with(
                    'job-synthetic', None, reviewer_id='reviewer-server')

    def test_protocol_ids_body_is_strict_bounded_and_unique(self):
        invalid = [
            {'protocol_ids': 'known-protocol'},
            {'protocol_ids': [1]},
            {'protocol_ids': ['']},
            {'protocol_ids': ['   ']},
            {'protocol_ids': ['x' * 257]},
            {'protocol_ids': ['same', 'same']},
            {'protocol_ids': [str(index) for index in range(9)]},
            {'protocol_ids': ['ok'], 'reviewer_id': 'browser-claim'},
        ]
        with self.client() as client, patch.object(
            client.app.state.scientific_acceptance, 'create'
        ) as create:
            for body in invalid:
                with self.subTest(body=body):
                    response = self.authenticated_post(client, body)
                    self.assertEqual(response.status_code, 422)
            create.assert_not_called()

    def test_authentication_and_read_only_boundary_remain_enforced(self):
        path = '/api/lab/jobs/job-synthetic/scientific-assessments'
        with self.client() as client:
            response = client.post(
                path, headers={'X-TPD-Local': '1'},
                json={'protocol_ids': ['known-protocol']},
            )
            self.assertEqual(response.status_code, 401)
        with self.client(read_only=True) as client, patch.object(
            client.app.state.scientific_acceptance, 'create'
        ) as create:
            response = client.post(
                path, headers={'X-TPD-Local': '1'},
                json={'protocol_ids': ['known-protocol']},
            )
            self.assertEqual(response.status_code, 403)
            create.assert_not_called()


class ProtocolReassessmentUiTests(unittest.TestCase):
    def test_explicit_selection_action_and_old_assessment_fallback_exist(self):
        self.assertEqual(HTML.count('id="protocolReassessmentButton"'), 1)
        self.assertEqual(HTML.count('id="protocolEvaluation"'), 1)
        self.assertIn('자동 선택, 최신값 선택, 최고값 선택은 하지 않습니다', HTML)
        self.assertIn("input.type='checkbox'", JS)
        self.assertIn('state.selectedProtocolIds.add(row.protocol_id)', JS)
        self.assertIn('같은 종류(', JS)
        self.assertIn('선택한 계산으로 새 평가', HTML)
        self.assertIn('이전 형식의 평가', JS)

    def test_reassessment_posts_only_existing_source_bound_ids_and_selects_new_revision(self):
        self.assertIn('const protocol_ids=[...state.selectedProtocolIds]', JS)
        self.assertIn('protocol_ids},120000', JS)
        self.assertIn('state.selectedAssessmentId=newId', JS)
        self.assertIn('await reloadAssessments(newId)', JS)
        self.assertIn('state.assessment?.id!==newId', JS)
        self.assertIn('state.selectedProtocolIds.clear()', JS)
        self.assertIn('state.protocolReassessmentBusy', JS)
        self.assertNotIn('setTimeout(createProtocolReassessment', JS)
        self.assertNotIn('setInterval(createProtocolReassessment', JS)

    def test_server_evaluation_status_authority_sources_and_categories_are_rendered(self):
        self.assertIn('state.assessment?.protocol_evaluation', JS)
        self.assertIn('Array.isArray(closure.records)', JS)
        self.assertIn('record?.evaluation||{}', JS)
        self.assertIn('badge(item.effective_status)', JS)
        self.assertIn('badge(check?.status)', JS)
        self.assertIn('record.record_ref', JS)
        self.assertIn('record.archive_ref', JS)
        self.assertIn('protocolCandidateIqrTable(check?.candidates)', JS)
        self.assertIn('scope_review_not_implemented:', JS)
        self.assertIn('추가 개발 필요 상태', JS)
        self.assertIn('Scientific authority: false', JS)
        self.assertIn('Source runtime current: false', JS)
        self.assertIn("if(state.assessment?.source_runtime_current===false)head.append(badge('Archived'))", JS)
        self.assertNotIn("badge('archived diagnostic unchanged')", JS)
        for label in ('추가 개발 필요', '계산 결과 미달', '근거 자료 부족', '전문가 검수 필요'):
            self.assertIn(label, JS)
        self.assertIn("PENDING_ACTION_LABELS[x.category]||x.category", JS)
        self.assertNotIn('innerHTML', JS)


if __name__ == '__main__':
    unittest.main()

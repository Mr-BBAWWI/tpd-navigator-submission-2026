import tempfile
import unittest
from unittest.mock import patch

from fastapi.testclient import TestClient

from apps.api.main import create_app
from packages.platform.scientific_acceptance import ScientificAcceptanceService


class ScientificApiBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def client(self, **kwargs):
        return TestClient(
            create_app(
                data_root=self.tmp.name,
                enable_worker=False,
                **kwargs,
            )
        )

    def test_real_service_class_is_installed(self):
        with self.client() as client:
            self.assertIsInstance(
                client.app.state.scientific_acceptance,
                ScientificAcceptanceService,
            )

    def test_existing_routes_and_review_page_remain_available(self):
        with self.client() as client:
            self.assertEqual(client.get('/api/health').status_code, 200)
            self.assertEqual(client.get('/design').status_code, 200)
            page = client.get('/scientific-review/job-missing')
            self.assertEqual(page.status_code, 200)
            self.assertIn('과학적 준비도', page.text)

    def test_anonymous_post_is_401_before_lookup(self):
        with self.client() as client:
            response = client.post(
                '/api/lab/jobs/job-missing/scientific-assessments',
                headers={'X-TPD-Local': '1'},
                json={},
            )
            self.assertEqual(response.status_code, 401)

    def test_strict_actor_payload_is_rejected(self):
        with self.client() as client:
            response = client.post(
                '/api/lab/jobs/job-missing/scientific-assessments',
                headers={'X-TPD-Local': '1'},
                json={'reviewer_id': 'browser-claim'},
            )
            self.assertEqual(response.status_code, 422)

    def test_read_only_blocks_every_scientific_post(self):
        paths = [
            ('/api/lab/jobs/x/scientific-assessments', {}),
            ('/api/scientific-statements', {
                'text': 'x',
                'locator': 'doi:x',
            }),
            ('/api/scientific-assessments/x/decisions', {}),
            ('/api/scientific-assessments/x/compute', {
                'candidate_ids': None,
                'analog_ids': None,
            }),
            ('/api/scientific-assessments/x/export', {}),
            ('/api/scientific-decisions/x/withdraw', {}),
        ]
        with self.client(read_only=True) as client:
            for path, body in paths:
                with self.subTest(path=path):
                    response = client.post(
                        path,
                        headers={'X-TPD-Local': '1'},
                        json=body,
                    )
                    self.assertEqual(response.status_code, 403)

    def test_review_only_compute_is_403_and_review_mutation_reaches_auth(self):
        with self.client(review_only=True) as client:
            compute = client.post(
                '/api/scientific-assessments/x/compute',
                headers={'X-TPD-Local': '1'},
                json={'candidate_ids': None, 'analog_ids': None},
            )
            self.assertEqual(compute.status_code, 403)
            statement = client.post(
                '/api/scientific-statements',
                headers={'X-TPD-Local': '1'},
                json={'text': 'x', 'locator': 'doi:x'},
            )
            self.assertEqual(statement.status_code, 401)

    def test_compute_resolves_assessment_job_id(self):
        with self.client() as client, patch.object(
            client.app.state.scientific_acceptance,
            'view',
            return_value={'job_id': 'job-real'},
        ), patch.object(
            client.app.state.scientific_acceptance,
            'compute',
            return_value={'job_id': 'job-real'},
        ) as compute:
            with patch.object(
                client.app.state.report_reviews.identities,
                'session',
                return_value={'id': 'reviewer-1'},
            ):
                response = client.post(
                    '/api/scientific-assessments/assessment-1/compute',
                    headers={'X-TPD-Local': '1'},
                    json={'candidate_ids': None, 'analog_ids': None},
                )
            self.assertEqual(response.status_code, 200)
            self.assertEqual(compute.call_args.args[0], 'job-real')

    def test_nonlocal_host_is_rejected(self):
        with TestClient(
            create_app(data_root=self.tmp.name, enable_worker=False),
            base_url='http://public.example',
        ) as client:
            self.assertEqual(client.get('/api/health').status_code, 403)


if __name__ == '__main__':
    unittest.main()

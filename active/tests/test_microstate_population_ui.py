"""Microstate population UI and real API boundary checks."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fastapi.testclient import TestClient

from apps.api.main import create_app

ROOT = Path(__file__).resolve().parents[1]
HTML = (ROOT / "apps/web/scientific-review.html").read_text(encoding="utf-8")
JS = (ROOT / "apps/web/scientific-review.js").read_text(encoding="utf-8")


class MicrostatePopulationUiTests(unittest.TestCase):
    def test_controls_and_validation_are_wired(self):
        identifiers = (
            "populationProject", "populationJob", "populationResult",
            "populationSourceFile", "populationSourceText",
            "populationUploadButton", "populationUploadStatus",
            "populationSource", "populationLocator", "populationReport",
            "populationReview", "populationInterpretation",
            "addPopulationRow", "cancelPopulationEdit", "populationRows",
            "populationWarning",
        )
        for identifier in identifiers:
            self.assertEqual(HTML.count(f'id="{identifier}"'), 1)
            self.assertIn(f"$('{identifier}')", JS)
        self.assertIn("file.size>1024*1024", JS)
        self.assertIn("population_evidence:state.populationRows", JS)
        self.assertIn("같은 analog의 population 행을 중복 추가할 수 없습니다.", JS)
        self.assertIn("pH<0||pH>14", JS)
        self.assertIn("state.populationEditingIndex=index", JS)
        self.assertNotIn("innerHTML", JS)

    def test_placeholders_do_not_fabricate_measurements(self):
        self.assertIn('"selected_state_index": null', HTML)
        self.assertIn('"pH": null', HTML)
        self.assertIn('"state_population": null', HTML)
        self.assertIn("그대로는 유효하지 않습니다", HTML)
        self.assertIn("graph·ID·수량을 추측하지 않습니다", HTML)
        self.assertIn("인간 전문가 판단이나 승인을 자동 생성하지 않습니다", HTML)


class MicrostatePopulationApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def client(self, **kwargs):
        return TestClient(create_app(
            data_root=self.tmp.name, enable_worker=False, **kwargs))

    def test_anonymous_get_and_post_are_401(self):
        path = '/api/lab/jobs/job-missing/microstate-population-sources'
        with self.client() as client:
            self.assertEqual(client.get(path).status_code, 401)
            response = client.post(
                path, headers={'X-TPD-Local': '1'}, json={'text': '{}'})
            self.assertEqual(response.status_code, 401)

    def test_read_only_post_is_403(self):
        path = '/api/lab/jobs/job-missing/microstate-population-sources'
        with self.client(read_only=True) as client:
            response = client.post(
                path, headers={'X-TPD-Local': '1'}, json={'text': '{}'})
            self.assertEqual(response.status_code, 403)

    def test_missing_local_header_is_403(self):
        path = '/api/lab/jobs/job-missing/microstate-population-sources'
        with self.client() as client:
            response = client.post(path, json={'text': '{}'})
            self.assertEqual(response.status_code, 403)

    def test_invalid_pydantic_record_is_422_on_real_route(self):
        path = '/api/lab/jobs/job-missing/microstate-population-sources'
        with self.client() as client, patch.object(
            client.app.state.report_reviews.identities,
            'session', return_value={'id': 'reviewer-1'},
        ):
            response = client.post(
                path, headers={'X-TPD-Local': '1'}, json={'text': 7})
            self.assertEqual(response.status_code, 422)


if __name__ == '__main__':
    unittest.main()

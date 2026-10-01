"""Read-only protocol diagnostics API and scientific-review UI contracts."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from fastapi.testclient import TestClient

from apps.api.main import PROJECT, create_app
from packages.contracts import ContractError
from packages.platform.protocol_evidence import ProtocolEvidenceService
from packages.science.protocol_evidence import ProtocolEvidenceError

ROOT = Path(__file__).resolve().parents[1]
HTML = (ROOT / 'apps/web/scientific-review.html').read_text(encoding='utf-8')
JS = (ROOT / 'apps/web/scientific-review.js').read_text(encoding='utf-8')


def synthetic(identifier='pdiag-' + 'a' * 32, job_id='job-synthetic'):
    return {
        'id': identifier,
        'project_id': PROJECT,
        'job_id': job_id,
        'protocol_id': 'known-calibration',
        'kind': 'calibration',
        'created_at': '2026-10-01T00:00:00Z',
        'status': 'computed_diagnostic_unreviewed',
        'source_binding': {'result_sha256': '1' * 64},
        'assessment_link': {'assessment_id': 'assessment-synthetic'},
        'summary': {'observed': 2, 'secret_key': 'remove-me'},
        'measurement_verification': {
            'stage': {'graph': 'bound', 'local_path': '/private/raw'}
        },
        'pack_manifest_sha256': '2' * 64,
        'compact_sha256': '3' * 64,
        'archive_ref': {
            'artifact_id': 'artifact-synthetic', 'version': 1,
            'sha256': '4' * 64,
        },
        'scientific_approved': False,
        'formal_acceptance': False,
        'gates_affected': False,
        'operator_kind': 'explicit_local_cli',
    }


class ProtocolDiagnosticsApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def client(self, **kwargs):
        return TestClient(create_app(
            data_root=self.tmp.name, enable_worker=False, **kwargs))

    def test_real_service_class_is_installed(self):
        with self.client() as client:
            self.assertIsInstance(
                client.app.state.protocol_evidence, ProtocolEvidenceService)

    def test_list_and_view_delegate_anonymously_and_redact_recursively(self):
        row = synthetic()
        with self.client() as client:
            service = client.app.state.protocol_evidence
            with patch.object(service, 'list', return_value=[row]) as listed, \
                    patch.object(service, 'view', return_value=row) as viewed:
                listing = client.get(
                    '/api/lab/jobs/job-synthetic/protocol-diagnostics')
                detail = client.get('/api/protocol-diagnostics/' + row['id'])
        self.assertEqual(listing.status_code, 200)
        self.assertEqual(detail.status_code, 200)
        listed.assert_called_once_with('job-synthetic')
        viewed.assert_called_once_with(row['id'])
        self.assertNotIn('secret_key', listing.text)
        self.assertNotIn('local_path', detail.text)
        self.assertNotIn('/private/raw', detail.text)

    def test_gets_are_allowed_in_read_only_and_review_only_modes(self):
        for mode in ({'read_only': True}, {'review_only': True}):
            with self.subTest(mode=mode), self.client(**mode) as client, patch.object(
                client.app.state.protocol_evidence, 'list', return_value=[]
            ):
                response = client.get(
                    '/api/lab/jobs/job-synthetic/protocol-diagnostics')
                self.assertEqual(response.status_code, 200)

    def test_missing_job_and_identifier_use_generic_404(self):
        with self.client() as client:
            service = client.app.state.protocol_evidence
            for path, method in (
                ('/api/lab/jobs/job-missing/protocol-diagnostics', 'list'),
                ('/api/protocol-diagnostics/pdiag-' + 'f' * 32, 'view'),
            ):
                with self.subTest(path=path), patch.object(
                    service, method, side_effect=KeyError('foreign-or-missing')):
                    response = client.get(path)
                    self.assertEqual(response.status_code, 404)
                    self.assertEqual(response.json(), {
                        'detail': '이 프로젝트에서 자료를 찾지 못했습니다.'})

    def test_real_missing_job_and_identifier_use_generic_404(self):
        with self.client() as client:
            for path in (
                '/api/lab/jobs/job-missing/protocol-diagnostics',
                '/api/protocol-diagnostics/pdiag-' + 'f' * 32,
            ):
                with self.subTest(path=path):
                    response = client.get(path)
                    self.assertEqual(response.status_code, 404)
                    self.assertEqual(response.json(), {
                        'detail': '이 프로젝트에서 자료를 찾지 못했습니다.'})

    def test_service_is_constructed_with_current_project_scope(self):
        fake = MagicMock()
        fake.list.return_value = []
        with patch('apps.api.main.ProtocolEvidenceService', return_value=fake) as cls:
            with self.client() as client:
                self.assertIs(client.app.state.protocol_evidence, fake)
                self.assertEqual(client.get(
                    '/api/lab/jobs/job-synthetic/protocol-diagnostics').status_code, 200)
        self.assertEqual(cls.call_args.args[1], PROJECT)

    def test_zip_artifact_download_preserves_bytes_filename_and_project_scope(self):
        artifact_id = 'artifact-synthetic'
        reference = {'artifact_id': artifact_id, 'media_type': 'application/zip'}
        archive = b'PK\x03\x04synthetic-zip-fixture'
        with self.client() as client:
            store = client.app.state.store
            with patch.object(
                store, 'reference', return_value=reference,
            ) as referenced, patch.object(
                store, 'read', return_value=archive,
            ) as read:
                response = client.get(f'/api/artifacts/{artifact_id}/download')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.headers['content-type'], 'application/octet-stream')
        self.assertEqual(
            response.headers['content-disposition'],
            'attachment; filename="artifact-synthetic.zip"',
        )
        self.assertEqual(response.content, archive)
        referenced.assert_called_once_with(PROJECT, artifact_id)
        read.assert_called_once_with(PROJECT, reference)

    def test_contract_failures_are_409(self):
        for error in (
            ContractError('malformed contract record'),
            ProtocolEvidenceError('malformed protocol evidence record'),
        ):
            with self.subTest(error=type(error).__name__), self.client() as client, patch.object(
                client.app.state.protocol_evidence, 'list', side_effect=error,
            ):
                response = client.get(
                    '/api/lab/jobs/job-synthetic/protocol-diagnostics')
                self.assertEqual(response.status_code, 409)

    def test_no_mutation_route_or_service_call(self):
        path = '/api/lab/jobs/job-synthetic/protocol-diagnostics'
        with self.client() as client:
            service = client.app.state.protocol_evidence
            with patch.object(service, 'list') as listed, patch.object(
                    service, 'view') as viewed:
                response = client.post(
                    path, headers={'X-TPD-Local': '1'}, json={})
                self.assertEqual(response.status_code, 405)
                listed.assert_not_called()
                viewed.assert_not_called()
        with self.client(read_only=True) as client:
            response = client.post(
                path, headers={'X-TPD-Local': '1'}, json={})
            self.assertEqual(response.status_code, 403)

    def test_nonlocal_host_is_rejected(self):
        with TestClient(create_app(
            data_root=self.tmp.name, enable_worker=False),
            base_url='http://public.example',
        ) as client:
            response = client.get(
                '/api/lab/jobs/job-synthetic/protocol-diagnostics')
            self.assertEqual(response.status_code, 403)


class ProtocolDiagnosticsUiContractTests(unittest.TestCase):
    def test_section_and_controls_are_unique_and_wired(self):
        for identifier in (
            'protocolDiagnostics', 'protocolDiagnosticsState',
            'protocolDiagnosticsRetry',
        ):
            self.assertEqual(HTML.count(f'id="{identifier}"'), 1)
            self.assertIn(f"$('{identifier}')", JS)
        self.assertIn('프로토콜별 추가 계산', HTML)
        self.assertIn('미검수 · 과학 승인 아님', HTML)
        self.assertIn('기존 registry 기록은 미검수 상태로 그대로 보존됩니다. 명시적으로 선택해 새 평가에 사용하면 서버가 검증한 입력과 서버 criterion 상태가 새 revision에 반영되지만 과학적 승인이나 인간 정책 결정은 자동으로 생기지 않습니다.', HTML)
        self.assertIn('Raw model, 선택된 seed 또는 계산, 후보 수와 분모는 각 서버 record의 실제 summary 값이 있을 때만 표시합니다. 화면은 14개 gate의 통과 여부를 계산량에서 추론하지 않습니다.', HTML)

    def test_fail_closed_flags_scope_and_artifact_only_download(self):
        self.assertIn("PROTOCOL_DIAGNOSTIC_STATUS='computed_diagnostic_unreviewed'", JS)
        self.assertIn("/^pdiag-[a-f0-9]{32}$/.test(row.id)", JS)
        self.assertIn('row.job_id!==jobId', JS)
        self.assertIn("Array.isArray(row.source_binding)", JS)
        self.assertIn("Array.isArray(row.summary)", JS)
        self.assertIn("Array.isArray(row.measurement_verification)", JS)
        self.assertIn('row.scientific_approved!==false', JS)
        self.assertIn('row.formal_acceptance!==false', JS)
        self.assertIn('row.gates_affected!==false', JS)
        self.assertIn('artifactRef(full.archive_ref)', JS)
        self.assertIn("link('검증된 raw pack ZIP 다운로드',full.archive_ref)", JS)
        self.assertIn("protocolDetails('전체 source binding',full.source_binding)", JS)
        self.assertIn("protocolDetails('전체 assessment link',full.assessment_link)", JS)
        self.assertIn("protocolDetails('전체 진단 레코드',full)", JS)
        self.assertIn("pre.textContent='미관측'", JS)
        self.assertNotIn('innerHTML', JS)

    def test_optional_loading_is_generation_guarded_and_manual_retry_only(self):
        self.assertIn('protocolGeneration', JS)
        self.assertIn('generation!==state.protocolGeneration', JS)
        self.assertIn("$('protocolDiagnosticsRetry').onclick=loadProtocolDiagnostics", JS)
        self.assertIn('loadProtocolDiagnostics();if(job?.state', JS)
        self.assertNotIn('setInterval(loadProtocolDiagnostics', JS)
        self.assertNotIn('setTimeout(loadProtocolDiagnostics', JS)


if __name__ == '__main__':
    unittest.main()

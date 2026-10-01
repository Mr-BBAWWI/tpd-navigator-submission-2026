"""Fresh synthetic fixtures in temporary storage; never scientific/replay evidence."""
import copy
import io
import tempfile
import unittest
import uuid
import ast
from pathlib import Path

from fastapi.testclient import TestClient
from lxml import etree
from pypdf import PdfWriter
from pypdf.generic import DictionaryObject, NameObject, DecodedStreamObject

from apps.api.main import create_app, PROJECT
from packages.contracts import ArtifactReader, ContractError, encoded, validate_exchange
from packages.platform.store import Store
from packages.platform.orchestrator import Orchestrator
from packages.transport import Fetch, FetchError
from packages.literature.extract import extract_xml, extract_pdf, strip_markup


XML = b'''<article><front><article-meta><title-group><article-title>Synthetic source</article-title></title-group></article-meta></front><body><sec><title>Results</title><p id="p1">Synthetic text; no experimental claim.</p><table-wrap id="t1"><caption><p>Test table</p></caption><table><tr><td>Label</td><td>Value</td></tr></table></table-wrap><fig id="f1"><caption><p>Test caption</p></caption></fig></sec></body></article>'''


def pdf_bytes():
    writer = PdfWriter()
    page = writer.add_blank_page(width=300, height=300)
    font = DictionaryObject({NameObject("/Type"): NameObject("/Font"), NameObject("/Subtype"): NameObject("/Type1"), NameObject("/BaseFont"): NameObject("/Helvetica")})
    page[NameObject("/Resources")] = DictionaryObject({NameObject("/Font"): DictionaryObject({NameObject("/F1"): writer._add_object(font)})})
    content = DecodedStreamObject()
    content.set_data(b"BT /F1 12 Tf 20 250 Td (Synthetic attachment test only) Tj ET")
    page[NameObject("/Contents")] = writer._add_object(content)
    output = io.BytesIO()
    writer.write(output)
    return output.getvalue()


def target_row(accession="Q16539", membrane=False):
    return {"primaryAccession": accession, "organism": {"taxonId": 9606},
            "proteinDescription": {"recommendedName": {"fullName": {"value": "Synthetic target"}}},
            "genes": [{"geneName": {"value": "TESTGENE"}}], "features": [{"type": "Transmembrane"}] if membrane else []}


class FakeTransport:
    def __init__(self, rows=None, failed_search=False, failed_xml=False, failed_target=False):
        self.rows = [target_row()] if rows is None else rows
        self.failed_search, self.failed_xml, self.failed_target = failed_search, failed_xml, failed_target
        self.calls = []

    def get(self, url, params=None, max_bytes=2_000_000):
        self.calls.append((url, params))
        if "uniprot" in url:
            if self.failed_target:
                raise FetchError("NETWORK_ERROR", "Synthetic network failure")
            body = encoded({"results": self.rows})
            headers = {"x-total-results": str(len(self.rows))}
        elif url.endswith("/search"):
            if self.failed_search:
                raise FetchError("NETWORK_ERROR", "Synthetic network failure")
            rows = [{"source": "MED", "id": "1", "pmcid": "PMC1", "title": "Synthetic fulltext", "isOpenAccess": "Y", "abstractText": "Fallback abstract"},
                    {"source": "MED", "id": "2", "title": "Synthetic abstract", "isOpenAccess": "N", "abstractText": "Synthetic abstract only"},
                    {"source": "MED", "id": "3", "title": "Synthetic metadata", "isOpenAccess": "N"}]
            body = encoded({"hitCount": 3, "resultList": {"result": rows[:params["pageSize"]]}})
            headers = {}
        else:
            if self.failed_xml:
                raise FetchError("HTTP_404", "Synthetic missing fulltext")
            body, headers = XML, {}
        return Fetch(body, url, 200, headers)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(self.temp.name)

    def run_pipeline(self, transport=None, uploads=None):
        self.transport = transport or FakeTransport()
        self.runner = Orchestrator(self.store, PROJECT, lambda: self.transport)
        request = {"query": "TESTGENE", "query_kind": "gene_symbol", "objective": "", "cell_line": "", "pdf_artifact_ids": uploads or [], "max_documents": 6}
        run, _ = self.store.create_run(PROJECT, uuid.uuid4().hex, request)
        self.runner.execute(run["id"])
        return self.runner.view(run["id"])

    def assert_contracts(self, run):
        port = self.store.scope(PROJECT)
        for job in run["jobs"]:
            verdict = validate_exchange(port.json(job["job_ref"]), ArtifactReader(port.read), port.json(job["result_ref"]))
            self.assertTrue(verdict["contract_valid"])
            self.assertFalse(verdict["dispatch_authorized"])

    def test_pipeline_source_levels_and_xml_locators(self):
        run = self.run_pipeline()
        self.assertEqual(run["state"], "partial", run["issues"])
        bundle = run["bundle"]
        self.assertEqual([d["access_level"] for d in bundle["documents"]], ["fulltext", "abstract", "metadata"])
        self.assertEqual(bundle["coverage"]["llm_read_segment_count"], 0)
        self.assertEqual(bundle["evidence_records"], [])
        root = etree.fromstring(XML)
        for segment in bundle["segments"]:
            loc = segment["locator"]
            if loc["kind"] == "jats_xml":
                node = root.xpath(loc["xpath"])
                self.assertEqual(len(node), 1)
                expected = " ".join(" ".join(node[0].itertext()).split())
                self.assertEqual(self.store.read(PROJECT, segment["text_artifact_ref"]).decode(), expected)
        self.assert_contracts(run)

    def test_ambiguous_requires_current_human_selection(self):
        run = self.run_pipeline(FakeTransport(rows=[target_row(), target_row("Q06187")]))
        self.assertEqual(run["state"], "awaiting_selection", run["issues"])
        self.assertEqual(len(self.transport.calls), 1)
        with self.assertRaises(ValueError):
            self.runner.choose(run["id"], "NOT_SHOWN", 1)
        self.runner.choose(run["id"], "Q06187", 1)
        with self.assertRaises(ValueError):
            self.runner.choose(run["id"], "Q16539", 1)
        self.runner.execute(run["id"])
        selected = self.runner.view(run["id"])
        self.assertEqual(selected["target"]["selected_candidate_id"], "Q06187")
        self.assertEqual(selected["revision"], 2)
        self.assertEqual(len([c for c in self.transport.calls if "uniprot" in c[0]]), 1)
        self.assert_contracts(selected)

    def test_unsupported_membrane_does_not_search(self):
        run = self.run_pipeline(FakeTransport(rows=[target_row(membrane=True)]))
        self.assertEqual(run["state"], "unsupported", run["issues"])
        self.assertEqual(len(self.transport.calls), 1)

    def test_not_found_does_not_search(self):
        run = self.run_pipeline(FakeTransport(rows=[]))
        self.assertEqual(run["state"], "not_found", run["issues"])
        self.assertEqual(len(self.transport.calls), 1)

    def test_target_failure_has_saved_failed_result(self):
        run = self.run_pipeline(FakeTransport(failed_target=True))
        self.assertEqual(run["state"], "failed")
        self.assertEqual(len(run["jobs"]), 1)
        self.assertEqual(run["jobs"][0]["status"], "failed")
        self.assert_contracts(run)

    def test_fulltext_failure_keeps_abstract_and_reason(self):
        run = self.run_pipeline(FakeTransport(failed_xml=True))
        self.assertEqual(run["state"], "partial", run["issues"])
        document = run["bundle"]["documents"][0]
        self.assertEqual(document["access_level"], "abstract")
        self.assertIn("HTTP_404", [x["code"] for x in document["issues"]])
        self.assert_contracts(run)

    def test_search_failure_preserves_uploaded_pdf(self):
        reference = self.store.scope(PROJECT).put_raw(pdf_bytes(), "application/pdf")
        run = self.run_pipeline(FakeTransport(failed_search=True), [reference["artifact_id"]])
        self.assertEqual(run["state"], "partial", run["issues"])
        bundle = run["bundle"]
        self.assertEqual(bundle["searches"][0]["status"], "failed")
        self.assertEqual(bundle["documents"][0]["access_level"], "fulltext")
        self.assertEqual(bundle["segments"][0]["locator"]["page"], 1)
        self.assert_contracts(run)

    def test_bad_pdf_does_not_claim_read_fulltext(self):
        reference = self.store.scope(PROJECT).put_raw(b"%PDF-1.7\nbroken", "application/pdf")
        run = self.run_pipeline(uploads=[reference["artifact_id"]])
        self.assertEqual(run["state"], "partial", run["issues"])
        document = run["bundle"]["documents"][0]
        self.assertEqual(document["extraction"]["status"], "failed")
        self.assertEqual(document["access_level"], "unverified")
        self.assertIsNone(document["extraction"]["snapshot_ref"])
        self.assert_contracts(run)

    def test_saved_view_never_calls_provider(self):
        run = self.run_pipeline()
        count = len(self.transport.calls)
        other = Orchestrator(Store(self.temp.name), PROJECT, lambda: (_ for _ in ()).throw(AssertionError("Network called")))
        self.assertEqual(other.view(run["id"])["bundle"], run["bundle"])
        self.assertEqual(len(self.transport.calls), count)

    def test_duplicate_worker_does_not_overwrite_finished_run(self):
        run = self.run_pipeline()
        count = len(self.transport.calls)
        self.runner.execute(run["id"])
        self.assertEqual(self.runner.view(run["id"]), run)
        self.assertEqual(len(self.transport.calls), count)

    def test_rejected_cross_project_output_is_audited(self):
        run = self.run_pipeline()
        port = self.store.scope(PROJECT)
        job = port.json(run["jobs"][-1]["job_ref"])
        request_ref = port.json(job["input_manifest"])["payload_ref"]
        wrong = copy.deepcopy(run["bundle"])
        wrong["project_id"] = "wrong-project"
        with self.assertRaises(FetchError):
            self.runner.record_job(run, "collect_literature", request_ref, job["input_artifacts"], job["policy_ref"], lambda: wrong)
        result = port.json(self.store.get_run(PROJECT, run["id"])["jobs"][-1]["result_ref"])
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["output_artifacts"], [])
        self.assertEqual(result["issues"][0]["code"], "OUTPUT_CONTRACT_REJECTED")

    def test_store_scope_integrity_and_idempotency(self):
        reference = self.store.scope("p1").put_raw(b"abc", "text/plain")
        with self.assertRaises(KeyError):
            self.store.read("p2", reference)
        forged = {**reference, "provenance": "computed"}
        with self.assertRaises(ContractError):
            self.store.read("p1", forged)
        (Path(self.temp.name) / "blobs" / reference["artifact_id"]).write_bytes(b"corrupted")
        with self.assertRaises(ContractError):
            self.store.read("p1", reference)
        run, fresh = self.store.create_run("p1", "same-key", {"query": "one"})
        same, fresh2 = self.store.create_run("p1", "same-key", {"query": "one"})
        self.assertTrue(fresh)
        self.assertFalse(fresh2)
        self.assertEqual(run["id"], same["id"])
        with self.assertRaises(ValueError):
            self.store.create_run("p1", "same-key", {"query": "two"})
        self.store.interrupt_pending("p1")
        self.assertEqual(self.store.get_run("p1", run["id"])["state"], "interrupted")


class ExtractionTests(unittest.TestCase):
    def test_escaped_title_markup_preserves_text_and_ampersand(self):
        self.assertEqual(strip_markup("Title &lt;b&gt;X&lt;/b&gt; &amp; Y"), "Title X & Y")

    def test_pdf_page_and_text(self):
        chunks, warnings = extract_pdf(pdf_bytes(), 4)
        self.assertEqual(chunks[0]["page"], 1)
        self.assertIn("Synthetic attachment", chunks[0]["text"])
        self.assertIn("PDF_LAYOUT_AND_IMAGES_NOT_INTERPRETED", warnings)

    def test_xml_limits_and_entities(self):
        chunks, warnings = extract_xml(XML, 1)
        self.assertEqual(len(chunks), 1)
        self.assertIn("SEGMENT_LIMIT", warnings)
        with self.assertRaises(ValueError):
            extract_xml(b'<!DOCTYPE article [<!ENTITY x "expansion">]><article><p>&x;</p></article>', 4)

    def test_new_runtime_imports_keep_module_boundary(self):
        root = Path(__file__).resolve().parents[1]
        forbidden = {"Study", "study", "tpd_navigator", "degradeAI"}
        for directory in [root / "apps", root / "packages"]:
            for path in directory.rglob("*.py"):
                for item in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
                    names = [n.name for n in item.names] if isinstance(item, ast.Import) else [item.module or ""] if isinstance(item, ast.ImportFrom) else []
                    for name in names:
                        self.assertNotIn(name.split(".")[0], forbidden, str(path))
                        if "literature" in path.parts or "science" in path.parts:
                            self.assertFalse(name.startswith(("packages.platform", "apps")), str(path))


class ApiTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.app = create_app(self.temp.name, FakeTransport, enable_worker=False)
        self.client = TestClient(self.app)
        self.client.__enter__()
        self.addCleanup(self.client.__exit__, None, None, None)
        self.headers = {"X-TPD-Local": "1"}

    def test_upload_submission_saved_view_and_input_boundaries(self):
        response = self.client.post("/api/uploads", files={"file": ("test.pdf", pdf_bytes(), "application/pdf")}, headers=self.headers)
        self.assertEqual(response.status_code, 200, response.text)
        artifact_id = response.json()["artifact"]["artifact_id"]
        body = {"query": "TESTGENE", "query_kind": "gene_symbol", "pdf_artifact_ids": [artifact_id], "request_key": "test-idempotency"}
        first = self.client.post("/api/runs", json=body, headers=self.headers)
        self.assertEqual(first.status_code, 202, first.text)
        identifier = first.json()["run_id"]
        second = self.client.post("/api/runs", json=body, headers=self.headers)
        self.assertEqual(second.json(), {"run_id": identifier, "created": False})
        changed = self.client.post("/api/runs", json={**body, "query": "different"}, headers=self.headers)
        self.assertEqual(changed.status_code, 409)
        self.app.state.orchestrator.execute(identifier)
        result = self.client.get("/api/runs/" + identifier).json()
        self.assertEqual(result["state"], "partial", result["issues"])
        segment_id = result["bundle"]["segments"][0]["text_artifact_ref"]["artifact_id"]
        self.assertIn("Synthetic attachment", self.client.get("/api/artifacts/" + segment_id + "/text").json()["text"])
        download = self.client.get("/api/artifacts/" + artifact_id + "/download")
        self.assertEqual(download.content, pdf_bytes())
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.client.get("/static/app.js").status_code, 200)

    def test_mutation_origin_upload_type_scope_and_schema_rejection(self):
        self.assertEqual(self.client.post("/api/runs", json={}).status_code, 403)
        self.assertEqual(self.client.post("/api/runs", json={}, headers={**self.headers, "Origin": "https://example.org"}).status_code, 403)
        self.assertEqual(self.client.post("/api/runs", json={}, headers=self.headers).status_code, 422)
        self.assertEqual(self.client.post("/api/uploads", files={"file": ("fake.pdf", b"not pdf")}, headers=self.headers).status_code, 422)
        self.assertEqual(self.client.post("/api/uploads", files={"file": ("huge.pdf", b"%PDF-" + b"x" * 6_100_000)}, headers=self.headers).status_code, 413)
        self.assertEqual(self.client.get("/api/artifacts/absent/download").status_code, 404)
        self.assertFalse(self.client.get("/api/health").json()["public_deployment_ready"])


if __name__ == "__main__":
    unittest.main()

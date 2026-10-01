"""Authenticated structured human-source statements remain strict and non-approving."""
from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from packages.contracts import ContractError
from packages.platform.review_identity import ReviewAuthError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.science import synthesis_review


class ScientificStructuredSourceTests(unittest.TestCase):
    project = "scientific-structured-source-test"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.store = Store(Path(self.temp.name))
        self.service = ScientificAcceptanceService(self.store, self.project)
        key = self.service.auth.register("reviewer", "Reviewer")
        self.token = self.service.auth.login("reviewer", key)

    def _document(self, nested_ref=None):
        route = {
            "candidate_id": "fixture-candidate",
            "canonical_smiles": "CC",
            "source": {
                "doi": "10.fixture/example",
                "title": "Synthetic fixture; not a scientific experiment",
            },
            "locator": "fixture route record",
            "steps": [{
                "reagents": ["fixture reagent"],
                "conditions": "fixture conditions",
                "purification": "fixture purification",
                "characterization": "fixture characterization",
            }],
        }
        if nested_ref is not None:
            route["supporting_artifact"] = nested_ref
        return {"routes": [route]}

    def test_authenticated_structured_statement_resolves_and_validates_exact_fixture(self):
        registered = self.service.port.put_raw(b"registered fixture", "text/plain", "source")
        text = json.dumps(self._document(registered), ensure_ascii=False)
        made = self.service.create_statement(text, "/document/routes/0", self.token)
        self.assertEqual(made["text"], text)
        self.assertEqual(made["document"], self._document(registered))
        self.assertEqual(made["evidence_provenance"], "authenticated human-supplied")
        self.assertFalse(made["scientific_approved"])
        self.assertFalse(made["route_approved"])
        self.assertFalse(made["synthesized"])
        archived = self.service.port.json(made["ref"])
        route = self.service._resolve_locator(archived, "/document/routes/0")
        self.assertEqual(route["supporting_artifact"], registered)
        report = synthesis_review._route_records([route], "fixture-candidate", "CC")
        self.assertTrue(report["documented_route_for_graph"])
        self.assertFalse(report["approved"])
        foreign = copy.deepcopy(route)
        foreign["canonical_smiles"] = "CCC"
        rejected = synthesis_review._route_records([foreign], "fixture-candidate", "CC")
        self.assertFalse(rejected["documented_route_for_graph"])

    def test_strict_json_and_locator_rejections(self):
        bad = [
            ('{"routes":[],"routes":[]}', "SCIENTIFIC_STRUCTURED_STATEMENT_DUPLICATE_KEY"),
            ('{"value":NaN}', "SCIENTIFIC_STRUCTURED_STATEMENT_NONFINITE"),
            ('{"value":1e1000}', "SCIENTIFIC_STRUCTURED_STATEMENT_NONFINITE"),
            ('[]', "SCIENTIFIC_STRUCTURED_STATEMENT_OBJECT"),
            ('{"routes":[]}', "SCIENTIFIC_SOURCE_LOCATOR"),
        ]
        for text, code in bad:
            with self.subTest(code=code), self.assertRaisesRegex(ContractError, code):
                self.service.create_statement(text, "/document/routes/0", self.token)

    def test_structured_size_precedes_schema_and_json_validation(self):
        oversized_invalid = '{"route":"' + ("x" * 50000)
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_STRUCTURED_STATEMENT_SIZE"):
            self.service.create_statement(
                oversized_invalid,
                "/document/route",
                self.token,
            )
        oversized_utf8 = json.dumps({"route": "가" * 20000}, ensure_ascii=False)
        self.assertLessEqual(len(oversized_utf8), 50000)
        self.assertGreater(len(oversized_utf8.encode("utf-8")), 50000)
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_STRUCTURED_STATEMENT_SIZE"):
            self.service.create_statement(
                oversized_utf8,
                "/document/route",
                self.token,
            )

    def test_ordinary_doi_text_remains_unstructured(self):
        made = self.service.create_statement(
            "Synthetic bibliographic fixture only.",
            "doi:10.fixture/ordinary-text",
            self.token,
        )
        self.assertEqual(made["locator"], "doi:10.fixture/ordinary-text")
        self.assertNotIn("document", made)
        self.assertNotIn("structured_evidence", made)

    def test_unregistered_tampered_cross_project_and_other_author_refs_rejected(self):
        fake = {
            "artifact_id": "missing", "version": 1, "sha256": "0" * 64,
            "media_type": "text/plain", "schema_id": "urn:tpd-navigator:raw:1",
            "provenance": "source",
        }
        with self.assertRaises(ContractError):
            self.service.create_statement(json.dumps(self._document(fake)), "/document/routes/0", self.token)

        ref = self.service.port.put_raw(b"bound", "text/plain", "source")
        tampered = copy.deepcopy(ref)
        tampered["sha256"] = "f" * 64
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_ARTIFACT_METADATA"):
            self.service.create_statement(json.dumps(self._document(tampered)), "/document/routes/0", self.token)

        other = ScientificAcceptanceService(self.store, "other-project")
        other_key = other.auth.register("other", "Other")
        other_token = other.auth.login("other", other_key)
        with self.assertRaisesRegex(ContractError, "SCIENTIFIC_ARTIFACT_PROJECT"):
            other.create_statement(json.dumps(self._document(ref)), "/document/routes/0", other_token)

        author_key = self.service.auth.register("author", "Author")
        author_token = self.service.auth.login("author", author_key)
        author_statement = self.service.create_statement("author fixture", "fixture:author", author_token)
        with self.assertRaises(ContractError):
            self.service.create_statement(
                json.dumps(self._document(author_statement["ref"])),
                "/document/routes/0", self.token,
            )

    def test_unauthenticated_and_failed_write_rollback(self):
        text = json.dumps(self._document())
        with self.assertRaises(ReviewAuthError):
            self.service.create_statement(text, "/document/routes/0", "invalid-session")
        before = {p.name for p in (self.store.root / "blobs").iterdir()}
        original = self.service.saved._write

        def write_then_fail(db, paths, raw, media="application/json", provenance="computed"):
            original(db, paths, raw, media, provenance)
            raise RuntimeError("structured rollback")

        with mock.patch.object(self.service.saved, "_write", side_effect=write_then_fail):
            with self.assertRaisesRegex(RuntimeError, "structured rollback"):
                self.service.create_statement(text, "/document/routes/0", self.token)
        after = {p.name for p in (self.store.root / "blobs").iterdir()}
        self.assertEqual(after, before)
        with self.store.db() as db:
            count = db.execute(
                "SELECT COUNT(*) FROM scientific_statements WHERE project=?",
                (self.project,),
            ).fetchone()[0]
        self.assertEqual(count, 0)


if __name__ == "__main__":
    unittest.main()

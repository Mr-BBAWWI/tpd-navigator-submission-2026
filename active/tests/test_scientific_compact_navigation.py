"""Compact scientific assessment navigation remains scoped and non-authoritative."""
from __future__ import annotations

import sqlite3
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest import mock

from apps.api.main import create_app
from packages.platform.scientific_acceptance import ScientificAcceptanceService


class _Store:
    def __init__(self):
        self.connection = sqlite3.connect(":memory:")
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("""CREATE TABLE scientific_assessments(
            id TEXT PRIMARY KEY,
            project TEXT NOT NULL,
            job_id TEXT NOT NULL,
            revision INTEGER NOT NULL,
            policy_revision INTEGER NOT NULL,
            policy_digest TEXT NOT NULL,
            fingerprint TEXT NOT NULL,
            ref TEXT NOT NULL)""")

    @contextmanager
    def db(self):
        yield self.connection


class ScientificCompactNavigationTests(unittest.TestCase):
    def setUp(self):
        self.store = _Store()
        self.addCleanup(self.store.connection.close)
        self.service = ScientificAcceptanceService.__new__(ScientificAcceptanceService)
        self.service.store = self.store
        self.service.project = "project-a"
        rows = [
            ("assessment-a1", "project-a", "job-1", 1),
            ("assessment-a3", "project-a", "job-1", 3),
            ("assessment-a2-other-job", "project-a", "job-2", 2),
            ("assessment-foreign", "project-b", "job-1", 99),
        ]
        self.store.connection.executemany(
            """INSERT INTO scientific_assessments(
                id,project,job_id,revision,policy_revision,policy_digest,
                fingerprint,ref) VALUES(?,?,?,?,0,'digest','fingerprint','{}')""",
            rows,
        )

    def test_compact_list_is_project_scoped_ordered_and_does_not_expand_views(self):
        self.service._job = mock.Mock(return_value={"id": "job-1"})
        self.service._verified_result = mock.Mock(
            side_effect=AssertionError("compact navigation must not verify full result artifacts")
        )
        self.service._assessment = mock.Mock(
            side_effect=AssertionError("compact navigation must not read assessment artifacts")
        )
        self.service._effective = mock.Mock(
            side_effect=AssertionError("compact navigation must not build effective views")
        )

        rows = self.service.list("job-1", compact=True)

        self.assertEqual(rows, [
            {"id": "assessment-a3", "revision": 3, "summary_only": True},
            {"id": "assessment-a1", "revision": 1, "summary_only": True},
        ])
        self.service._job.assert_called_once()
        self.service._verified_result.assert_not_called()
        self.service._assessment.assert_not_called()
        self.service._effective.assert_not_called()

    def test_default_list_still_returns_full_effective_views(self):
        self.service._verified_result = mock.Mock(return_value=None)
        self.service._assessment = mock.Mock(
            side_effect=lambda _db, identifier: {"id": identifier, "criteria": [identifier]}
        )
        self.service._effective = mock.Mock(
            side_effect=lambda _db, value: {**value, "validated_full_view": True}
        )

        rows = self.service.list("job-1")

        self.assertEqual(
            [row["id"] for row in rows],
            ["assessment-a3", "assessment-a1"],
        )
        self.assertTrue(all(row["validated_full_view"] for row in rows))
        self.service._verified_result.assert_called_once()
        self.assertEqual(self.service._assessment.call_count, 2)
        self.assertEqual(self.service._effective.call_count, 2)

    def test_individual_view_still_uses_assessment_and_effective_validation_path(self):
        self.service._assessment = mock.Mock(return_value={"id": "assessment-a3"})
        self.service._effective = mock.Mock(return_value={"id": "assessment-a3", "criteria": []})

        value = self.service.view("assessment-a3")

        self.assertEqual(value["id"], "assessment-a3")
        self.service._assessment.assert_called_once()
        self.service._effective.assert_called_once()

    def test_api_compact_is_opt_in_and_default_call_remains_full(self):
        with tempfile.TemporaryDirectory() as directory:
            app = create_app(Path(directory), enable_worker=False)
            route = next(
                item for item in app.routes
                if getattr(item, "path", None) ==
                "/api/lab/jobs/{identifier}/scientific-assessments"
                and "GET" in getattr(item, "methods", set())
            )
            service = app.state.scientific_acceptance
            with mock.patch.object(service, "list", return_value=[]) as listing:
                self.assertEqual(route.endpoint("job-1", compact=True), [])
                listing.assert_called_once_with("job-1", compact=True)
                listing.reset_mock()
                self.assertEqual(route.endpoint("job-1"), [])
                listing.assert_called_once_with("job-1")

    def test_frontend_requests_compact_navigation_then_one_bounded_full_view(self):
        javascript = (
            Path(__file__).resolve().parents[1]
            / "apps" / "web" / "scientific-review.js"
        ).read_text(encoding="utf-8")
        self.assertIn("/scientific-assessments?compact=true", javascript)
        self.assertIn("ASSESSMENT_VIEW_TIMEOUT=30000", javascript)
        self.assertIn(
            "api(`/api/scientific-assessments/${encodeURIComponent(id)}`,{},ASSESSMENT_VIEW_TIMEOUT)",
            javascript,
        )
        self.assertNotIn("scientific-assessments?compact=false", javascript)


if __name__ == "__main__":
    unittest.main()

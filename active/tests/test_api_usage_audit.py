import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

import httpx

from packages.agents.provider import AgentError, DaconProvider
from scripts.summarize_api_usage import summarize


class ProviderUsageTests(unittest.TestCase):
    def call(self, body, headers=None, status=200, key="test-secret-key"):
        transport = httpx.MockTransport(
            lambda request: httpx.Response(status, json=body, headers=headers or {})
        )
        provider = DaconProvider(key, transport=transport)
        self.addCleanup(provider.close)
        return provider.complete(
            model="gpt-5.6-sol", instructions="i", context={}, max_output_tokens=10
        )

    def test_incomplete_preserves_code_and_observed_usage(self):
        body = {
            "id": "resp_1", "model": "gpt-5.6-sol", "status": "incomplete",
            "usage": {"input_tokens": 4, "output_tokens": 2, "total_tokens": 6},
            "output": [{"type": "reasoning", "text": "do-not-copy"}],
        }
        with self.assertRaises(AgentError) as raised:
            self.call(body)
        self.assertEqual(raised.exception.code, "API_INCOMPLETE_USAGE_UNRECONCILED")
        self.assertEqual(raised.exception.metadata["usage"]["total_tokens"], 6)
        self.assertNotIn("do-not-copy", json.dumps(raised.exception.metadata))

    def test_http_code_is_compatible_and_usage_is_attached(self):
        body = {
            "id": "resp_http", "status": "failed",
            "usage": {"total_tokens": 9}, "request": {"authorization": "test-secret-key"},
        }
        with self.assertRaises(AgentError) as raised:
            self.call(body, status=429)
        self.assertEqual(raised.exception.code, "API_HTTP_429_NO_RETRY_USAGE_UNRECONCILED")
        self.assertEqual(raised.exception.metadata["usage"], {"total_tokens": 9})
        self.assertNotIn("test-secret-key", json.dumps(raised.exception.metadata))

    def test_key_is_redacted_from_output_and_metadata(self):
        key = "EXACT-SECRET-123"
        body = {
            "id": "resp_" + key, "model": key, "status": "completed",
            "usage": {"total_tokens": 3},
            "output": [{"type": "message", "role": "assistant", "content": [
                {"type": "output_text", "text": "before " + key + " after"}
            ]}],
        }
        result = self.call(body, key=key)
        self.assertEqual(result.text, "before [REDACTED] after")
        self.assertNotIn(key, json.dumps(result.metadata))

    def test_empty_refusal_and_invalid_usage_attach_partial_observation(self):
        cases = [
            ([], "API_EMPTY_OUTPUT"),
            ([{"type": "message", "role": "assistant", "content": [
                {"type": "refusal", "text": "no"}
            ]}], "API_REFUSAL"),
        ]
        for output, code in cases:
            with self.subTest(code=code), self.assertRaises(AgentError) as raised:
                self.call({
                    "id": "resp", "status": "completed",
                    "usage": {"input_tokens": 4, "output_tokens": True, "total_tokens": -1},
                    "output": output,
                })
            self.assertEqual(raised.exception.code, code)
            self.assertEqual(raised.exception.metadata["usage"], {"input_tokens": 4})
            self.assertFalse(raised.exception.metadata["usage_known"])

    def test_total_only_is_accepted(self):
        result = self.call({
            "id": "resp_total", "status": "completed", "usage": {"total_tokens": 7},
            "output": [{"type": "message", "role": "assistant", "content": [
                {"type": "output_text", "text": "ok"}
            ]}],
        })
        self.assertEqual(result.metadata["usage"], {"total_tokens": 7})


class SummaryTests(unittest.TestCase):
    def make_store(self, root, calls, project="project-a"):
        database = root / "index.sqlite3"
        connection = sqlite3.connect(database)
        connection.executescript("""
            CREATE TABLE workbench_jobs(id TEXT PRIMARY KEY, body TEXT NOT NULL);
            CREATE TABLE workbench_calls(job_id TEXT, ordinal INTEGER, body TEXT);
            CREATE TABLE artifacts(id TEXT PRIMARY KEY, project TEXT, metadata TEXT);
        """)
        connection.execute(
            "INSERT INTO workbench_jobs VALUES(?,?)",
            ("job-1", json.dumps({"id": "job-1", "project_id": project})),
        )
        for ordinal, call in enumerate(calls, 1):
            connection.execute(
                "INSERT INTO workbench_calls VALUES(?,?,?)",
                ("job-1", ordinal, json.dumps(call)),
            )
        connection.commit()
        connection.close()
        return database

    def register_answer(self, root, reference, project, payload):
        blob = json.dumps(payload).encode()
        reference = dict(reference, sha256=hashlib.sha256(blob).hexdigest())
        (root / "blobs").mkdir(exist_ok=True)
        (root / "blobs" / reference["artifact_id"]).write_bytes(blob)
        connection = sqlite3.connect(root / "index.sqlite3")
        connection.execute(
            "INSERT INTO artifacts VALUES(?,?,?)",
            (reference["artifact_id"], project, json.dumps(reference)),
        )
        connection.execute(
            "UPDATE workbench_calls SET body=? WHERE ordinal=1",
            (json.dumps({
                "requested_at": "2025-01-01T00:00:00Z",
                "answer_ref": reference, "usage_status": "observed",
                "usage": {"total_tokens": 99},
            }),),
        )
        connection.commit()
        connection.close()
        return reference

    def test_registered_project_hash_blob_metadata_is_authoritative(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_store(root, [{}])
            self.register_answer(
                root, {"artifact_id": "answer-1", "version": 1}, "project-a",
                {"text": "never summarize this", "metadata": {
                    "response_id": "resp-product", "usage": {"total_tokens": 8}
                }},
            )
            result = summarize([], [root])
            self.assertEqual(result["product_calls"]["calls"], 1)
            self.assertEqual(result["product_calls"]["total_tokens"], 8)
            self.assertNotIn("never summarize this", json.dumps(result))

    def test_tamper_cross_project_and_traversal_are_corrupted_refs(self):
        for mode in ("tamper", "cross-project", "traversal"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self.make_store(root, [{}])
                reference = self.register_answer(
                    root, {"artifact_id": "answer-1", "version": 1}, "project-a",
                    {"text": "secret", "metadata": {"usage": {"total_tokens": 8}}},
                )
                connection = sqlite3.connect(root / "index.sqlite3")
                if mode == "tamper":
                    (root / "blobs" / "answer-1").write_text("tampered", encoding="utf-8")
                elif mode == "cross-project":
                    connection.execute("UPDATE artifacts SET project='project-b'")
                else:
                    bad = dict(reference, artifact_id="../answer-1")
                    connection.execute(
                        "UPDATE workbench_calls SET body=? WHERE ordinal=1",
                        (json.dumps({"answer_ref": bad, "usage": {"total_tokens": 5}}),),
                    )
                connection.commit()
                connection.close()
                result = summarize([], [root])
                self.assertEqual(result["source_errors"].get("corruptedref"), 1)
                self.assertEqual(result["product_calls"]["source_error_calls"], 1)

    def test_partial_total_only_invalid_bool_and_conflict(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "a.metadata.json").write_text(json.dumps({
                "response_id": "partial", "usage": {"input_tokens": 2, "output_tokens": True}
            }), encoding="utf-8")
            (root / "b.metadata.json").write_text(json.dumps({
                "response_id": "total", "usage": {"total_tokens": 10}
            }), encoding="utf-8")
            (root / "c.metadata.json").write_text(json.dumps({
                "response_id": "conflict", "usage": {"total_tokens": 4}
            }), encoding="utf-8")
            (root / "d.metadata.json").write_text(json.dumps({
                "response_id": "conflict", "usage": {"total_tokens": 5}
            }), encoding="utf-8")
            result = summarize([root], [])
            calls = result["developer_calls"]
            self.assertEqual(calls["calls"], 3)
            self.assertEqual(calls["conflict_calls"], 1)
            self.assertEqual(calls["input_tokens"], 2)
            self.assertEqual(calls["total_tokens"], 10)
            self.assertEqual(calls["partial_usage_calls"], 2)

    def test_non_api_parent_metadata_and_unrelated_failure_are_ignored(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "parent.metadata.json").write_text(json.dumps({"parent_id": "P-1", "status": "failed"}), encoding="utf-8")
            (root / "unrelated.failure.json").write_text(json.dumps({"message": "ordinary pipeline failure"}), encoding="utf-8")
            result = summarize([root], [])
            self.assertEqual(result["developer_calls"]["calls"], 0)
            self.assertEqual(result["unknown_usage_failures"], 0)
            self.assertEqual(result["source_hashes"], [])

    def test_malformed_api_receipt_is_a_parse_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "broken.metadata.json"
            path.write_text('{"response_id":', encoding="utf-8")
            result = summarize([root], [])
            self.assertEqual(result["developer_calls"]["calls"], 0)
            self.assertEqual(result["parse_errors"][0]["path"], str(path.resolve()))

    def test_chronological_nonmonotonic_quota_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "later.metadata.json").write_text(json.dumps({"response_id": "later", "observed_at": "2025-01-02T00:00:00Z", "quota_headers": {"x-team-remaining-tokens": 20}}), encoding="utf-8")
            (root / "earlier.metadata.json").write_text(json.dumps({"response_id": "earlier", "observed_at": "2025-01-01T00:00:00Z", "quota_headers": {"x-team-remaining-tokens": 10}}), encoding="utf-8")
            result = summarize([root], [])
            self.assertTrue(result["latest_quota"]["nonmonotonic_observations"])
            self.assertEqual(result["latest_quota"]["values"]["x-team-remaining-tokens"], 20)

    def test_unknown_failure_reconciles_and_overlapping_dirs_dedupe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            nested = root / "nested"
            nested.mkdir()
            (nested / "one.failure.json").write_text(json.dumps({
                "api_metadata": {"response_id": "same", "usage": {}}
            }), encoding="utf-8")
            (nested / "two.metadata.json").write_text(json.dumps({
                "response_id": "same", "usage": {"total_tokens": 6}
            }), encoding="utf-8")
            result = summarize([root, nested], [])
            self.assertEqual(result["unknown_usage_failures"], 0)
            self.assertEqual(result["developer_calls"]["calls"], 1)
            paths = [item["path"] for item in result["source_hashes"]]
            self.assertEqual(len(paths), len(set(paths)))


if __name__ == "__main__":
    unittest.main()

"""Provider metadata retention and token-accounting regression tests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from packages.agents.research_campaign import run_campaign
from scripts.summarize_api_usage import _extract, summarize
from tests.test_research_campaign import (
    FakeProvider, critic, evidence, judge, proposals,
)


class MetadataProvider(FakeProvider):
    def __init__(self, values):
        super().__init__(values)
        self.index = 0

    def complete(self, **kwargs):
        response = super().complete(**kwargs)
        self.index += 1
        response.metadata = {
            "requested_model": "gpt-5.6-sol",
            "returned_model": "gpt-5.6-sol",
            "observed_at": f"2026-10-01T00:00:0{self.index}Z",
            "response_id": f"response-{self.index}",
            "status": "completed",
            "http_status": 200,
            "usage": {
                "input_tokens": 2,
                "output_tokens": 3,
                "total_tokens": 5,
            },
            "usage_known": True,
            "quota_headers": {},
            "quota_is_estimate": True,
            "elapsed_seconds": 0.125,
            "Authorization": "Bearer fake-do-not-retain",
            "api_key": "fake-do-not-retain",
            "private_extra": {"internal": "fake-do-not-retain"},
        }
        return response


class CampaignMetadataTests(unittest.TestCase):
    def test_safe_allowlist_and_flat_nested_usage_without_duplication(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "campaign"
            provider = MetadataProvider([proposals(), critic(), judge()])
            result = run_campaign(evidence(), output, provider=provider, max_rounds=1)
            self.assertEqual(result["status"], "complete")

            token_report = json.loads((output / "token-report.json").read_text(encoding="utf-8"))
            self.assertEqual(token_report["observed_usage"]["input_tokens"], 6)
            self.assertEqual(token_report["observed_usage"]["output_tokens"], 9)
            self.assertEqual(token_report["observed_usage"]["total_tokens"], 15)

            metadata_files = sorted(output.glob("rounds/**/*.metadata.json"))
            self.assertEqual(len(metadata_files), 3)
            response_ids = set()
            for path in metadata_files:
                payload = json.loads(path.read_text(encoding="utf-8"))
                self.assertLessEqual(set(payload), {
                    "format", "requested_model", "returned_model", "response_id",
                    "status", "http_status", "observed_at", "usage",
                    "usage_known", "quota_headers", "quota_is_estimate",
                    "elapsed_seconds", "request_sha256", "output_sha256",
                    "applied_automatically",
                })
                self.assertNotIn("Authorization", payload)
                self.assertNotIn("api_key", payload)
                self.assertNotIn("private_extra", payload)
                self.assertEqual(payload["requested_model"], "gpt-5.6-sol")
                self.assertEqual(payload["returned_model"], "gpt-5.6-sol")
                response_id = payload["response_id"]
                response_ids.add(response_id)
                response_index = response_id.removeprefix("response-")
                self.assertEqual(
                    payload["observed_at"],
                    f"2026-10-01T00:00:0{response_index}Z",
                )
                self.assertEqual(payload["usage"], {
                    "input_tokens": 2,
                    "output_tokens": 3,
                    "total_tokens": 5,
                })
                for key in ("request_sha256", "output_sha256"):
                    digest = payload[key]
                    self.assertEqual(len(digest), 64)
                    self.assertTrue(all(character in "0123456789abcdef" for character in digest))
                self.assertIs(payload["applied_automatically"], False)

            self.assertEqual(response_ids, {
                "response-1", "response-2", "response-3",
            })

    def test_raw_nested_v1_and_flat_v2_deduplicate_without_journal_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            audit = Path(temporary)
            journal = audit / "journal.jsonl"
            journal.write_text('{"event":"immutable"}\n', encoding="utf-8")
            before = hashlib.sha256(journal.read_bytes()).hexdigest()
            fallback = "2026-10-01T00:00:00Z"

            for index in range(1, 4):
                usage = {
                    "input_tokens": 2,
                    "output_tokens": 3,
                    "total_tokens": 5,
                }
                provider_metadata = {
                    "format": "provider-response-metadata/1",
                    "requested_model": "gpt-5.6-sol",
                    "returned_model": "gpt-5.6-sol",
                    "response_id": f"raw-response-{index}",
                    "observed_at": f"2026-10-01T00:00:0{index}Z",
                    "usage": usage,
                }
                receipt = {
                    "format": "provider-response-metadata/2",
                    "response_id": f"raw-response-{index}",
                    "observed_at": f"2026-10-01T00:00:0{index}Z",
                    "usage": usage,
                    "provider_metadata": provider_metadata,
                }
                self.assertEqual(_extract(receipt, fallback)["usage"], usage)
                (audit / f"receipt-{index}.metadata.json").write_text(
                    json.dumps(receipt), encoding="utf-8",
                )

            totals = summarize([audit], [])["developer_calls"]
            self.assertEqual(totals["calls"], 3)
            self.assertEqual(totals["input_tokens"], 6)
            self.assertEqual(totals["output_tokens"], 9)
            self.assertEqual(totals["total_tokens"], 15)
            self.assertEqual(hashlib.sha256(journal.read_bytes()).hexdigest(), before)


if __name__ == "__main__":
    unittest.main()

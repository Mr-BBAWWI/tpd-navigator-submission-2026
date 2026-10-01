import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from packages.agents.research_campaign import (
    build_context, load_diagnostic_evidence, run_campaign,
)


def evidence():
    criteria = [{"id": f"c{i}", "status": "failed" if i == 0 else "pending"} for i in range(14)]
    candidates = []
    for parent in ("P1", "P2"):
        for e3 in ("CRBN", "VHL"):
            for i in range(6):
                candidates.append({"candidate_key": f"{parent}/{e3}/D{i}", "candidate_id": f"D{i}",
                    "parent_id": parent, "job_id": "J", "evidence_source_id": "E1",
                    "canonical_smiles": f"C{i}", "e3_type": e3,
                    "selected_analog_docking_metrics": {"best_core_rmsd_A": i + .25}})
    return {"format": "research-campaign-evidence/1", "digest": "abc",
        "scientific_accepted": False,
        "assessment_snapshot": {"criteria14": criteria, "original_scientific_accepted": False},
        "evidence_sources": [{"source_id": "E1", "closure_verified": True}],
        "parents": [{"parent_id": x, "job_id": "J", "evidence_source_id": "E1",
                     "mode": "strict", "counts_recomputed_from_records": {}} for x in ("P1", "P2")],
        "candidates": candidates,
        "analogs": [{"analog_id": f"W{i}", "selected": True, "metric": i} for i in range(35)]}


def proposals(ids=("PX1", "PX2", "PX3"), unknown=False):
    return {"role": "proposer", "proposals": [{"proposal_id": pid,
        "criterion_ids": ["c0"], "candidate_keys": ["BAD" if unknown else "P1/CRBN/D0"],
        "parent_ids": ["P1"], "evidence_source_ids": ["E1"], "hypothesis": "가설",
        "action": "측정 실험을 설계", "measurable_success_criteria": ["값을 사전 기준과 비교"],
        "rationale": "E1 근거", "limits": "새 측정값 없음"} for pid in ids]}


def critic(ids=("PX1", "PX2", "PX3"), status="supported"):
    return {"role": "adversarial_critic", "reviews": [{"proposal_id": x,
        "evidence_source_ids": ["E1"], "findings": ["근거 확인"],
        "counterarguments": ["한계 존재"], "status": status} for x in ids]}


def judge(ids=("PX1", "PX2", "PX3"), verdict="ready_for_experiment"):
    return {"role": "experiment_judge", "decisions": [{"proposal_id": x,
        "verdict": verdict, "reason": "실험 준비도 판정", "conditions": []} for x in ids]}


class FakeProvider:
    mode = "mocked_test_only"
    def __init__(self, values): self.values, self.calls, self.requests = list(values), 0, []
    def complete(self, **kwargs):
        self.calls += 1
        self.requests.append(kwargs)
        self.assert_request(kwargs)
        value = self.values.pop(0)
        if isinstance(value, Exception): raise value
        return SimpleNamespace(text=value if isinstance(value, str) else json.dumps(value),
            metadata={"usage": {"input_tokens": 2, "output_tokens": 3, "total_tokens": 5}})

    @staticmethod
    def assert_request(request):
        assert request["model"] == "gpt-5.6-sol"
        context = request["context"]
        contract = context["response_contract"]
        assert contract["role"] == context["role"]
        assert contract["json_schema"]["additionalProperties"] is False
        assert "expected_ids" in contract
        if context["role"] != "proposer" or context["round"] > 1:
            assert "current_proposals" in context


class CampaignTests(unittest.TestCase):
    def run_temp(self, provider, rounds=1):
        root = Path(tempfile.mkdtemp()) / "out"
        return root, run_campaign(evidence(), root, provider, rounds)

    def test_full_flow_three_roles(self):
        root, result = self.run_temp(FakeProvider([proposals(), critic(), judge()]))
        self.assertEqual(result["status"], "complete")
        self.assertFalse(result["scientific_final_approval"])
        self.assertEqual(result["original_criterion_statuses"]["c0"], "failed")
        self.assertEqual(json.loads((root / "token-report.json").read_text())["observed_usage"]["total_tokens"], 15)

    def test_malformed_retained_with_usage(self):
        root, result = self.run_temp(FakeProvider(["not-json"]))
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["reason"], "MALFORMED_JSON")
        self.assertEqual((root / "rounds/round-01/proposer.response.json").read_text(), "not-json")
        report = json.loads((root / "token-report.json").read_text())
        self.assertEqual(report["observed_usage"]["total_tokens"], 5)
        self.assertEqual(report["unknown_usage_calls"], 0)
        note = json.loads((root / "rounds/round-01/proposer.failure-note.json").read_text())
        self.assertEqual(note["protocol"], "research-campaign/4")
        self.assertEqual(note["reason"], "MALFORMED_JSON")
        self.assertTrue(note["response_received"])
        self.assertEqual(note["raw_stored_bytes"], len(b"not-json"))
        self.assertEqual(note["partial_round_successful_roles"], [])

    def test_request_hash_matches_pretty_saved_request_bytes(self):
        root, _ = self.run_temp(FakeProvider(["not-json"]))
        rd = root / "rounds/round-01"
        request_path = rd / "proposer.request.json"
        request_bytes = request_path.read_bytes()
        request = json.loads(request_bytes)
        metadata = json.loads((rd / "proposer.response.metadata.json").read_text())
        self.assertEqual(metadata["request_sha256"], hashlib.sha256(request_bytes).hexdigest())
        canonical = json.dumps(
            request, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("utf-8")
        self.assertNotEqual(metadata["request_sha256"], hashlib.sha256(canonical).hexdigest())
        self.assertEqual(request["artifact_version"], "prompt/4")
        self.assertEqual(request["context"]["prompt_version"], "research-campaign-prompts/4")

    def test_oversized_unicode_response_preserved_and_counted_once(self):
        raw = '{"payload":"' + ("한" * 40_000) + '"}'
        raw_bytes = raw.encode("utf-8")
        self.assertGreater(len(raw_bytes), 100_000)
        root, result = self.run_temp(FakeProvider([raw]))
        rd = root / "rounds/round-01"
        response_path = rd / "proposer.response.json"
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["reason"], "RESPONSE_TOO_LARGE")
        self.assertEqual(response_path.read_bytes(), raw_bytes)
        metadata = json.loads((rd / "proposer.response.metadata.json").read_text())
        self.assertEqual(metadata["output_sha256"], hashlib.sha256(raw_bytes).hexdigest())
        note = json.loads((rd / "proposer.failure-note.json").read_text())
        self.assertEqual(note["reason"], "RESPONSE_TOO_LARGE")
        self.assertTrue(note["response_received"])
        self.assertEqual(note["raw_stored_bytes"], len(raw_bytes))
        self.assertEqual(note["structured_parse_limit"], 100_000)
        report = json.loads((root / "token-report.json").read_text())
        self.assertEqual(report["observed_usage"]["total_tokens"], 5)
        self.assertEqual(report["unknown_usage_calls"], 0)

    def test_provider_exception_writes_empty_raw_and_unknown_usage_note(self):
        root, result = self.run_temp(FakeProvider([RuntimeError("private provider detail")]))
        rd = root / "rounds/round-01"
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["reason"], "PROVIDER_ERROR")
        self.assertEqual((rd / "proposer.response.json").read_bytes(), b"")
        note = json.loads((rd / "proposer.failure-note.json").read_text())
        self.assertEqual(note["reason"], "PROVIDER_ERROR")
        self.assertFalse(note["response_received"])
        self.assertEqual(note["raw_stored_bytes"], 0)
        self.assertNotIn("private provider detail", json.dumps(note))
        report = json.loads((root / "token-report.json").read_text())
        self.assertEqual(report["observed_usage"]["total_tokens"], 0)
        self.assertEqual(report["unknown_usage_calls"], 1)

    def test_unknown_candidate_rejected(self):
        _, result = self.run_temp(FakeProvider([proposals(unknown=True)]))
        self.assertEqual(result["reason"], "UNKNOWN_CANDIDATE_KEYS")

    def test_judge_cannot_override_unsupported(self):
        _, result = self.run_temp(FakeProvider([proposals(), critic(status="unsupported"), judge()]))
        self.assertEqual(result["status"], "failed")

    def test_revision_round_bounded(self):
        values = [proposals(), critic(status="needs_revision"), judge(verdict="revise"),
                  proposals(), critic(), judge()]
        provider = FakeProvider(values)
        _, result = self.run_temp(provider, 2)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(provider.calls, 6)

    def test_single_proposal_revision_supported(self):
        first_judge = judge(verdict="ready_for_experiment")
        first_judge["decisions"][0]["verdict"] = "revise"
        first_critic = critic()
        first_critic["reviews"][0]["status"] = "needs_revision"
        values = [proposals(), first_critic, first_judge,
                  proposals(("PX1",)), critic(("PX1",)), judge(("PX1",))]
        provider = FakeProvider(values)
        _, result = self.run_temp(provider, 2)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(provider.calls, 6)

    def test_output_overwrite_blocked(self):
        root = Path(tempfile.mkdtemp()) / "output-v3"
        root.mkdir()
        sentinel = root / "sentinel.txt"
        sentinel.write_bytes(b"immutable-v3-output")
        with self.assertRaises(FileExistsError):
            run_campaign(evidence(), root, FakeProvider([proposals()]))
        self.assertEqual(sentinel.read_bytes(), b"immutable-v3-output")
        self.assertEqual(list(root.iterdir()), [sentinel])

    def test_offline_no_api(self):
        root, result = self.run_temp(None)
        self.assertEqual(result["mode"], "comparison_only")
        self.assertEqual(result["rounds_completed"], 0)
        self.assertTrue((root / "journal.jsonl").exists())

    def test_cli_help_from_isolated_directory(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "run_research_campaign.py"
        completed = subprocess.run(
            [sys.executable, str(script), "--help"], cwd=tempfile.mkdtemp(),
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--rounds", completed.stdout)
        self.assertIn("--diagnostic-json", completed.stdout)

    def test_sample_group_coverage(self):
        context = build_context(evidence())
        self.assertEqual(len(context["candidate_sample_coverage"]), 4)
        self.assertTrue(all(x["full_count"] == 6 and x["sample_count"] == 4
                            and not x["sample_is_ranked_best"]
                            for x in context["candidate_sample_coverage"]))
        self.assertEqual(len(context["selected_analog_sample"]), 30)

    def test_candidate_attachment_nested_fields_preserved(self):
        supplied = evidence()
        supplied["candidates"][0]["attachment_metadata"] = {
            "warhead": {"atom_index": 7}, "linker": {"bond": {"kind": "single"}}
        }
        context = build_context(supplied)
        row = next(x for x in context["candidate_sample"]
                   if x["candidate_key"] == "P1/CRBN/D0")
        self.assertEqual(row["attachment_metadata"]["warhead"]["atom_index"], 7)
        self.assertEqual(row["attachment_metadata"]["linker"]["bond"]["kind"], "single")

    def test_diagnostic_bytes_source_binding_and_copy(self):
        directory = Path(tempfile.mkdtemp())
        source = directory / "calibration recent.json"
        raw = b'{"calibration":{"run":"N3"},"topology":{"ok":true}}\n'
        source.write_bytes(raw)
        prepared = load_diagnostic_evidence([source])
        expected = "diag-" + __import__("hashlib").sha256(raw).hexdigest()
        self.assertEqual(prepared[0]["source_id"], expected)
        output = directory / "out"
        result = run_campaign(evidence(), output, supplemental_evidence=prepared)
        copied = next((output / "diagnostics").iterdir())
        self.assertEqual(copied.read_bytes(), raw)
        context = json.loads((output / "context.json").read_text(encoding="utf-8"))
        self.assertEqual(context["supplemental_diagnostics"][0]["source_id"], expected)
        self.assertIn(expected, {x["source_id"] for x in context["evidence_sources"]})
        self.assertEqual(result["original_criterion_statuses"]["c0"], "failed")
        self.assertEqual(len(result["original_criterion_statuses"]), 14)

    def test_diagnostic_malformed_rejected_before_output(self):
        directory = Path(tempfile.mkdtemp())
        source = directory / "bad.json"
        source.write_text('{"x":1,"x":2}', encoding="utf-8")
        output = directory / "out"
        with self.assertRaisesRegex(ValueError, "DIAGNOSTIC_DUPLICATE_JSON_KEY"):
            run_campaign(evidence(), output, supplemental_evidence=[source])
        self.assertFalse(output.exists())

    def test_diagnostic_tamper_rejected_before_output(self):
        directory = Path(tempfile.mkdtemp())
        source = directory / "diagnostic.json"
        source.write_text('{"core":{"status":"recent"}}', encoding="utf-8")
        prepared = load_diagnostic_evidence([source])
        prepared[0]["raw_bytes"] += b" "
        output = directory / "out"
        with self.assertRaisesRegex(ValueError, "DIAGNOSTIC_TAMPERED"):
            run_campaign(evidence(), output, supplemental_evidence=prepared)
        self.assertFalse(output.exists())

    def test_diagnostic_nonobject_and_nan_rejected(self):
        directory = Path(tempfile.mkdtemp())
        for name, raw, code in (
            ("array.json", "[]", "DIAGNOSTIC_NOT_OBJECT"),
            ("nan.json", '{"value":NaN}', "DIAGNOSTIC_NONFINITE_JSON"),
        ):
            source = directory / name
            source.write_text(raw, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, code):
                load_diagnostic_evidence([source])


if __name__ == "__main__":
    unittest.main()

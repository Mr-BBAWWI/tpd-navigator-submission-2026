"""Synthetic I3-1 artifacts linked to existing new-implementation I1 fixtures."""
import copy
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ACTIVE = ROOT.parents[2]
sys.path.insert(0, str(ACTIVE))
from packages.contracts import encoded, parse_json
from validation import URI, validate


class Fixture:
    def __init__(self):
        self.index, self.raw, self.refs = [], {}, {}
        base = ACTIVE / "contracts/payloads/v0.1.0/examples"
        for row in parse_json((base / "artifact-index.json").read_bytes())["artifacts"]:
            reference = row["reference"]
            self.raw[reference["artifact_id"]] = (reference, (base / row["path"]).read_bytes())
            self.refs[reference["artifact_id"].removeprefix("fixture-")] = reference
            self.index.append({"reference": reference, "path": "../../../payloads/v0.1.0/examples/" + row["path"]})
        self.counter = 0

    def put(self, name, value, kind=None):
        if kind: validate(kind, value)
        self.counter += 1
        raw = encoded(value)
        reference = {"artifact_id": "i3-fixture-" + name + "-" + str(self.counter), "version": 1,
                     "sha256": hashlib.sha256(raw).hexdigest(), "media_type": "application/json",
                     "schema_id": URI + "#/$defs/" + kind if kind else "urn:tpd-navigator:raw:test-fixture", "provenance": "test_fixture"}
        self.raw[reference["artifact_id"]] = (reference, raw)
        self.refs[name] = reference
        self.index.append({"reference": reference, "path": "artifacts/" + reference["artifact_id"] + ".json"})
        return reference

    def read(self, reference):
        registered, data = self.raw[reference["artifact_id"]]
        if registered != reference: raise ValueError("FIXTURE_REGISTERED_REFERENCE_MISMATCH")
        return data

    def value(self, name): return parse_json(self.read(self.refs[name]))

    def add(self, name, kind, **fields):
        return self.put(name, dict(payload_type=kind, payload_version="0.1.0", project_id="fixture-project",
                                   run_id="fixture-run", data_mode="test_fixture", **fields), kind)


def build():
    f = Fixture()
    scope = f.put("scope", {"fixture_only": True, "meaning": "No scientific setting; contract exercise"})
    fields = [{"field_id": "fixture_limit", "mutability": "reviewable", "value": {"kind": "integer", "value": 1, "unit": "count", "endpoint": None},
               "scope_ref": scope, "effect": "fixed_universe_filter", "required_for_execution": True, "rationale": "Synthetic technical limit"},
              {"field_id": "fixture_integrity", "mutability": "locked", "value": {"kind": "boolean", "value": True, "unit": None, "endpoint": None},
               "scope_ref": scope, "effect": "fixed_universe_filter", "required_for_execution": True, "rationale": "Synthetic locked field"}]
    target = f.refs["target-resolution"]
    policy = f.add("policy", "PolicySnapshot", target_resolution_ref=target, fields=fields, limitations=["Synthetic fixture only"])
    evidence = f.refs["evidence-bundle"]
    bundle = f.value("evidence-bundle")
    segment = bundle["segments"][0]
    segment_id = segment["segment_id"]
    request = f.add("analysis-request", "LiteratureAssessmentRequest", input_revision=1, evidence_bundle_ref=evidence, target_resolution_ref=target,
        base_policy_ref=policy, allowed_field_ids=["fixture_limit"], requested_segment_ids=[segment_id], research_context_ref=f.refs["target-query"],
        limits={"max_delivered_segments": 1, "max_claims": 3, "max_proposals": 3})
    response = f.put("model-response", {"fixture_only": True, "note": "No LLM was called"})
    receipt = f.add("receipt", "DeliveryReceipt", request_ref=request, role="LiteratureAnalyst", status="completed",
        delivered_segment_ids=[segment_id], response_ref=response, model="fixture-no-model", prompt_version="fixture-no-prompt", issues=[])
    source_text = f.read(segment["text_artifact_ref"]).decode()
    claim = {"claim_id": "claim-fixture", "statement": "가상 원문의 인용 연결을 검사하는 문장입니다.", "basis_kind": "researcher_hypothesis",
             "evidence_ids": [], "citations": [{"segment_id": segment_id, "start": 0, "end": min(12, len(source_text)), "quote": source_text[:12]}], "limitations": ["Scientific meaning not evaluated"]}
    proposal = {"proposal_id": "proposal-fixture", "field_id": "fixture_limit", "before": copy.deepcopy(fields[0]["value"]),
                "proposed": {"kind": "integer", "value": 2, "unit": "count", "endpoint": None}, "scope_ref": scope,
                "basis_kind": "researcher_hypothesis", "supporting_claim_ids": ["claim-fixture"], "rationale": "Contract fixture", "expected_tradeoff": "Synthetic candidate count change"}
    assessment = f.add("assessment", "LiteratureAssessment", request_ref=request, receipt_ref=receipt, annotated_evidence_ref=evidence,
        claims=[claim], summary=[{"text": "가상 요약", "claim_ids": ["claim-fixture"]}], conflicts=[], proposals=[proposal],
        unassessed_segment_ids=[s["segment_id"] for s in bundle["segments"] if s["segment_id"] != segment_id], limitations=["Synthetic only"], issues=[])
    pre_request = f.add("proposal-review-request", "ReviewRequest", subject_kind="LiteratureAssessment", subject_ref=assessment, required_claim_ids=["claim-fixture"])
    pre_review = f.add("proposal-review", "ReviewRecord", request_ref=pre_request, role="Critic", status="completed", reviewed_claim_ids=["claim-fixture"],
        response_ref=response, findings=[], limitations=["No Critic was run"])
    new_fields = copy.deepcopy(fields)
    new_fields[0]["value"] = copy.deepcopy(proposal["proposed"])
    proposed = f.add("proposed-policy", "PolicySnapshot", target_resolution_ref=target, fields=new_fields, limitations=["Synthetic fixture only"])
    molecule = f.put("molecule", {"fixture_only": True, "not_a_molecular_structure": True})
    mapping = f.put("mapping", {"fixture_only": True, "not_a_validated_atom_map": True})
    calculation = f.put("calculation", {"fixture_only": True, "not_a_scientific_result": True})
    context = f.put("evaluation-context", {"engine": "fixture-no-calculation", "version": "1"})
    universe = f.add("universe", "CandidateUniverse", target_resolution_ref=target,
        candidates=[{"candidate_id": name, "structure_ref": molecule} for name in ["candidate-a", "candidate-b"]], generation_context_ref=context, limitations=["Synthetic identities"])
    preview_request = f.add("preview-request", "PreviewRequest", assessment_ref=assessment, selected_proposal_ids=["proposal-fixture"], base_policy_ref=policy,
        proposed_policy_ref=proposed, universe_ref=universe, evaluation_context_ref=context, proposal_review_ref=pre_review)
    outcomes = [{"candidate_id": name, "disposition": "eligible", "mapping_status": "resolved", "mapping_ref": mapping, "chemical_status": "valid",
                 "evidence_kind": "computed_hypothesis", "assessment_refs": [calculation], "reasons": ["Synthetic outcome; no chemistry executed"], "limitations": ["Fixture only"]} for name in ["candidate-a", "candidate-b"]]
    baseline = copy.deepcopy(outcomes)
    baseline[1]["disposition"] = "excluded"
    before = f.add("baseline", "Phase1Result", universe_ref=universe, policy_ref=policy, assessment_ref=assessment, evaluation_context_ref=context,
        status="complete", candidates=baseline, issues=[])
    after = f.add("proposed-result", "Phase1Result", universe_ref=universe, policy_ref=proposed, assessment_ref=assessment, evaluation_context_ref=context,
        status="complete", candidates=outcomes, issues=[])
    preview = f.add("preview", "CpuPreview", request_ref=preview_request, baseline_result_ref=before, proposed_result_ref=after,
        added_candidate_ids=["candidate-b"], removed_candidate_ids=[], held_candidate_ids=[], limitations=["No CPU chemistry was executed"])
    review_request = f.add("review-request", "ReviewRequest", subject_kind="CpuPreview", subject_ref=preview, required_claim_ids=["claim-fixture"])
    critic = f.add("critic", "ReviewRecord", request_ref=review_request, role="Critic", status="completed", reviewed_claim_ids=["claim-fixture"],
        response_ref=response, findings=[], limitations=["Synthetic review"])
    from validation import REQUIRED_CHECKS
    checks = f.add("code-check", "CodeCheckReport", request_ref=review_request, validator_version="fixture-not-runtime", status="complete",
        checks=[{"check_id": n, "outcome": "pass"} for n in sorted(REQUIRED_CHECKS)], findings=[])
    quote = f.put("quote", {"fixture_only": True, "not_a_provider_quote": True})
    f.add("packet", "GatePacket", gate="G1", input_revision=1, preview_ref=preview, critic_ref=critic, code_check_ref=checks,
        selected_candidate_ids=["candidate-a", "candidate-b"], cost_scope={"execution_mode": "development_live", "quote_ref": quote, "max_gpu_seconds": 0, "max_cost_krw": 0}, purpose="phase2_candidate_design")
    return f


def save(f):
    out = ROOT / "examples"
    (out / "artifacts").mkdir(parents=True, exist_ok=True)
    for item in f.index:
        if item["path"].startswith("artifacts/"):
            (out / item["path"]).write_bytes(f.read(item["reference"]))
    (out / "artifact-index.json").write_bytes(encoded({"data_mode": "test_fixture", "artifacts": f.index, "named_refs": f.refs}))


if __name__ == "__main__":
    fixture = build()
    save(fixture)
    print(f"Saved {fixture.counter} new fixture artifacts; I1 sources remain external, unchanged.")

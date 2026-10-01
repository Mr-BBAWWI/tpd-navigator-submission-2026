"""Executable I3-1 contract cases. All input/results are synthetic, no model/science run."""
import copy
import json
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2]))
from packages.contracts import parse_json, encoded, ContractError
from build_schema import SCHEMA, ROOTS
from build_fixtures import build, save
from validation import ReviewReader, validate, validate_decision, URI
from jsonschema import Draft202012Validator


def revise(f, name, mutation=None, root="packet"):
    """Build a new immutable graph, preserving every original fixture byte."""
    changes = name if isinstance(name, dict) else {name: mutation}
    targets = {f.refs[n]["artifact_id"]: action for n, action in changes.items()}
    memo = {}
    def visit(reference):
        identity = reference["artifact_id"]
        if identity in memo: return memo[identity]
        if not reference["schema_id"].startswith(URI): return reference
        value = parse_json(f.read(reference))
        if identity in targets: targets[identity](value)
        def transform(item):
            if isinstance(item, dict):
                if "artifact_id" in item and "sha256" in item: return visit(item)
                return {k: transform(v) for k, v in item.items()}
            if isinstance(item, list): return [transform(v) for v in item]
            return item
        revised = transform(value)
        new = f.put("revision-" + identity, revised, revised["payload_type"])
        memo[identity] = new
        return new
    return visit(f.refs[root])


def decision(f, reference=None, **changes):
    reference = reference or f.refs["packet"]
    reader = ReviewReader(f.read)
    ready = reader.readiness(reference)
    body = {"packet_ref": reference, "input_revision": 1, "decision": "approve", "acknowledged_warning_ids": ready["warning_ids"],
            "reason": "Synthetic approval contract exercise only", "idempotency_key": "fixture-decision-1"}
    body.update(changes)
    return body


ACTOR = {"kind": "human", "id": "fixture-researcher", "name": "Fixture Researcher", "project_id": "fixture-project"}
TIME = "2026-09-11T00:00:00Z"


def main():
    Draft202012Validator.check_schema(SCHEMA)
    require_schema = parse_json((ROOT / "workflow.schema.json").read_bytes())
    assert require_schema == SCHEMA
    positives, negatives = [], []
    f = build()
    save(f)

    unchanged = build()
    no_change = revise(unchanged, {
        "assessment": lambda x: x.update(proposals=[]),
        "preview-request": lambda x: x.update(selected_proposal_ids=[], proposal_review_ref=None),
        "proposed-policy": lambda x: x["fields"][0]["value"].update(value=1),
        "proposed-result": lambda x: x["candidates"][1].update(disposition="excluded"),
        "preview": lambda x: x.update(added_candidate_ids=[]),
        "packet": lambda x: x.update(selected_candidate_ids=["candidate-a"]),
    })
    assert ReviewReader(unchanged.read).readiness(no_change)["approval_contract_ready"]
    positives.append("zero_proposals_keep_baseline")
    reader = ReviewReader(f.read)
    for name, reference in f.refs.items():
        if reference["schema_id"].startswith(URI):
            reader.verify(reference)
            positives.append("graph:" + name)
    state = reader.readiness(f.refs["packet"])
    assert state["approval_contract_ready"] and not state["dispatch_authorized"]
    record = validate_decision(decision(f), reader, current_packet_ref=f.refs["packet"], current_revision=1, authenticated_actor=ACTOR, decided_at=TIME)
    assert record["data_mode"] == "test_fixture" and not record["dispatch_authorized"]
    approval = f.put("approval-record", record, "ApprovalRecord")
    reader.verify(approval)
    positives.extend(["gate_readiness", "synthetic_approval_record"])
    save(f)

    def bad(name, execute, expected):
        try: execute()
        except (ContractError, ValueError) as exc:
            assert expected in str(exc), (name, expected, str(exc))
            negatives.append({"case": name, "rejected_with": expected})
        else: raise AssertionError("Invalid case accepted: " + name)

    def mutation(name, target, change, expected, root="packet"):
        def run():
            test = build()
            changed = revise(test, target, change, root)
            ReviewReader(test.read).verify(changed)
        bad(name, run, expected)

    mutation("invented quote", "assessment", lambda x: x["claims"][0]["citations"][0].update(quote="invented"), "CITATION_TEXT_MISMATCH")
    mutation("wrong text offset", "assessment", lambda x: x["claims"][0]["citations"][0].update(end=999999), "INVALID_CITATION_SPAN")
    mutation("unread citation", "receipt", lambda x: x.update(delivered_segment_ids=[]), "CITATION_NOT_DELIVERED")
    mutation("unknown summary claim", "assessment", lambda x: x["summary"][0].update(claim_ids=["absent"]), "UNKNOWN_SUMMARY_CLAIM")
    mutation("invented experiment", "assessment", lambda x: x["claims"][0].update(basis_kind="experiment"), "CLAIM_BASIS_WITHOUT_OBSERVATION")
    mutation("false coverage", "assessment", lambda x: x.update(unassessed_segment_ids=["absent"]), "ANALYSIS_COVERAGE_MISMATCH")
    mutation("output limit", "analysis-request", lambda x: x["limits"].update(max_claims=0), "ANALYSIS_OUTPUT_LIMIT")
    mutation("refused model promoted", "receipt", lambda x: x.update(status="refused", issues=[{"code":"REFUSED","category":"invalid_output","message":"Fixture refusal","affected_artifact_ids":[],"retry_hint":"do_not_retry"}]), "INCOMPLETE_MODEL_OUTPUT_PROMOTED")
    mutation("locked field allowlist", "analysis-request", lambda x: x.update(allowed_field_ids=["fixture_integrity"]), "LOCKED_OR_UNKNOWN_ALLOWED_FIELD")
    mutation("unknown policy field", "assessment", lambda x: x["proposals"][0].update(field_id="absent"), "UNAUTHORIZED_POLICY_FIELD")
    mutation("stale before value", "assessment", lambda x: x["proposals"][0]["before"].update(value=999), "PROPOSAL_BASE_OR_SCOPE_MISMATCH")
    mutation("changed endpoint domain", "assessment", lambda x: x["proposals"][0].update(proposed={"kind":"number","value":2,"unit":"nM","endpoint":"IC50"}), "POLICY_VALUE_DOMAIN_CHANGED")
    mutation("hidden hard change", "proposed-policy", lambda x: x["fields"][1]["value"].update(value=False), "UNDECLARED_POLICY_CHANGE")
    mutation("candidate generation disguised as filter", "policy", lambda x: x["fields"][0].update(effect="requires_regeneration"), "CANDIDATE_REGENERATION_REQUIRED")
    mutation("missing pre-preview Critic", "preview-request", lambda x: x.update(proposal_review_ref=None), "PROPOSAL_REVIEW_REQUIRED")
    mutation("invented delta", "preview", lambda x: x.update(added_candidate_ids=[]), "PREVIEW_DELTA_MISMATCH")
    mutation("candidate omitted", "proposed-result", lambda x: x["candidates"].pop(), "CANDIDATE_COVERAGE_MISMATCH")
    mutation("unresolved mapping eligible", "proposed-result", lambda x: x["candidates"][0].update(mapping_status="ambiguous"), "INVALID_CANDIDATE_ELIGIBLE")
    mutation("missing map artifact", "proposed-result", lambda x: x["candidates"][0].update(mapping_ref=None), "RESOLVED_MAPPING_WITHOUT_ARTIFACT")
    mutation("invalid chemistry eligible", "proposed-result", lambda x: x["candidates"][0].update(chemical_status="invalid"), "INVALID_CANDIDATE_ELIGIBLE")
    mutation("Critic hides unreviewed claim", "critic", lambda x: x.update(reviewed_claim_ids=[]), "INCOMPLETE_CRITIC_COVERAGE")
    mutation("omitted required review", "review-request", lambda x: x.update(required_claim_ids=[]), "REVIEW_CLAIMS_OMITTED")
    mutation("unknown selected candidate", "packet", lambda x: x.update(selected_candidate_ids=["absent"]), "SELECTED_CANDIDATE_UNKNOWN")
    mutation("public GPU scope", "packet", lambda x: x["cost_scope"].update(execution_mode="public_explore", max_gpu_seconds=1), "PUBLIC_GPU_FORBIDDEN")
    mutation("wrong packet revision", "packet", lambda x: x.update(input_revision=2), "PACKET_REVISION_MISMATCH")
    mutation("cross-project graph", "packet", lambda x: x.update(project_id="another-project"), "GRAPH_CONTEXT_MISMATCH")
    mutation("fixture promoted to real", "packet", lambda x: x.update(data_mode="real"), "FIXTURE_IN_REAL_GRAPH")
    mutation("unknown payload field", "assessment", lambda x: x.update(auto_approve=True), "SCHEMA")
    mutation("duplicate proposal ID", "assessment", lambda x: x["proposals"].append({**x["proposals"][0], "rationale":"different text"}), "DUPLICATE_PROPOSAL_ID")

    def gated(name, target, change, reason):
        test=build()
        changed=revise(test,target,change)
        reading=ReviewReader(test.read)
        readiness=reading.readiness(changed)
        assert readiness["review_ready"] and reason in readiness["blocking_reasons"]
        bad(name, lambda: validate_decision(decision(test,changed), reading, current_packet_ref=changed, current_revision=1, authenticated_actor=ACTOR, decided_at=TIME), "GATE_APPROVAL_BLOCKED")
        rejected=validate_decision(decision(test,changed,decision="reject"),reading,current_packet_ref=changed,current_revision=1,authenticated_actor=ACTOR,decided_at=TIME)
        assert rejected["decision"]=="reject"
        positives.append("reject_packet:"+name)
    gated("empty selection", "packet", lambda x:x.update(selected_candidate_ids=[]), "NO_VALID_SELECTED_CANDIDATES")
    gated("incomplete Critic", "critic", lambda x:x.update(status="incomplete"), "CRITIC_NOT_COMPLETED")
    gated("missing code check", "code-check", lambda x:x["checks"].pop(), "REQUIRED_CODE_CHECK_NOT_PASSED")
    gated("failed code check", "code-check", lambda x:x["checks"][0].update(outcome="fail"), "REQUIRED_CODE_CHECK_NOT_PASSED")
    finding={"finding_id":"block-fixture","severity":"blocking","message":"Synthetic blocker","subject_refs":[f.refs["assessment"]],"claim_ids":["claim-fixture"]}
    gated("blocking finding", "critic", lambda x:x["findings"].append(finding), "BLOCKING_FINDING")
    base=decision(f)
    def decide(body, **kwargs):
        return validate_decision(body, ReviewReader(f.read), current_packet_ref=f.refs["packet"], current_revision=kwargs.pop("revision",1), authenticated_actor=kwargs.pop("actor",ACTOR), decided_at=TIME, **kwargs)
    bad("stale input revision",lambda:decide(base,revision=2),"STALE_GATE_DECISION")
    bad("different displayed packet",lambda:decide({**base,"packet_ref":f.refs["preview"]}),"STALE_GATE_DECISION")
    bad("LLM cannot approve",lambda:decide(base,actor={**ACTOR,"kind":"agent"}),"HUMAN_ACTOR_REQUIRED")
    bad("client provided actor",lambda:decide({**base,"actor_name":"Claimed user"}),"SCHEMA")
    bad("other project actor",lambda:decide(base,actor={**ACTOR,"project_id":"other"}),"HUMAN_ACTOR_REQUIRED")
    bad("unacknowledged uncertainty",lambda:decide({**base,"acknowledged_warning_ids":[]}),"WARNING_NOT_ACKNOWLEDGED")
    bad("made up warning",lambda:decide({**base,"acknowledged_warning_ids":["absent"]}),"UNKNOWN_WARNING_ACKNOWLEDGEMENT")
    bad("tampered artifact bytes",lambda:ReviewReader(lambda ref:f.read(ref)+b" ").verify(f.refs["packet"]),"ARTIFACT_HASH_MISMATCH")
    bad("nonfinite number",lambda:validate("ParameterValue",{"kind":"number","value":float('nan'),"unit":"nM","endpoint":"Kd"}),"NON_FINITE_NUMBER")
    missing = build()
    unset = {"kind":"unset","value":None,"unit":None,"endpoint":None,"reason":"Fixture required setting missing"}
    missing_ref = revise(missing,{"policy":lambda x:x["fields"][1].update(value=unset), "proposed-policy":lambda x:x["fields"][1].update(value=unset)})
    missing_reader = ReviewReader(missing.read)
    assert "REQUIRED_POLICY_UNSET" in missing_reader.readiness(missing_ref)["blocking_reasons"]
    bad("required policy not configured",lambda:validate_decision(decision(missing,missing_ref),missing_reader,current_packet_ref=missing_ref,current_revision=1,authenticated_actor=ACTOR,decided_at=TIME),"GATE_APPROVAL_BLOCKED")
    pre_warning = build()
    warning={"finding_id":"pre-warning","severity":"warning","message":"Synthetic prior uncertainty","subject_refs":[pre_warning.refs["assessment"]],"claim_ids":["claim-fixture"]}
    warning_ref=revise(pre_warning,"proposal-review",lambda x:x["findings"].append(warning))
    warning_reader=ReviewReader(pre_warning.read)
    assert "pre-warning" in warning_reader.readiness(warning_ref)["warning_ids"]
    positives.append("pre_preview_warning_preserved")
    body=decision(pre_warning,warning_ref)
    body["acknowledged_warning_ids"].remove("pre-warning")
    bad("prior uncertainty not acknowledged",lambda:validate_decision(body,warning_reader,current_packet_ref=warning_ref,current_revision=1,authenticated_actor=ACTOR,decided_at=TIME),"WARNING_NOT_ACKNOWLEDGED")
    disk_index=parse_json((ROOT/'examples/artifact-index.json').read_bytes())
    for item in disk_index['artifacts']:
        raw=(ROOT/'examples'/item['path']).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==item['reference']['sha256']
    report={"scope":"I3-1 offline contracts only; no real LLM, chemistry, approval persistence, or dispatch", "definitions":len(SCHEMA['$defs']),"entrypoints":len(ROOTS),
            "positive_checks":len(positives),"negative_cases":len(negatives),"artifact_hashes_checked":len(disk_index['artifacts']),"positive":positives,"negative":negatives}
    (ROOT/'verification.json').write_bytes(encoded(report))
    print(f"PASS I3-1: {len(positives)} positive, {len(negatives)} negative, {len(disk_index['artifacts'])} artifact hashes; no live execution.")


if __name__ == "__main__": main()

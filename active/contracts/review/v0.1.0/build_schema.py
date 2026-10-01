"""I3-1 offline review contracts; no runtime capability is enabled."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
URI = "urn:tpd-navigator:review:0.1.0"
I1 = "urn:tpd-navigator:payloads:0.1.0"
BOUNDARY = "urn:tpd-navigator:boundary:0.1.0"
TEXT = {"type": "string", "minLength": 1}
ID = {"type": "string", "pattern": "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$"}
INT = {"type": "integer", "minimum": 0}
POS = {"type": "integer", "minimum": 1}
NULL = {"type": "null"}
BOOL = {"type": "boolean"}
D = {}


def ref(n): return {"$ref": "#/$defs/" + n}
def enum(*v): return {"enum": list(v)}
def const(v): return {"const": v}
def nullable(v): return {"anyOf": [v, NULL]}
def arr(v, minimum=0): return {"type": "array", "items": v, "minItems": minimum, "uniqueItems": True}
def obj(**p): return {"type": "object", "properties": p, "required": list(p), "additionalProperties": False}
def union(*v): return {"oneOf": list(v)}
def payload(n, **p):
    return obj(payload_type=const(n), payload_version=const("0.1.0"), project_id=ID, run_id=ID,
               data_mode=enum("real", "test_fixture"), **p)


D["ArtifactRef"] = {"$ref": BOUNDARY + "#/$defs/ArtifactRef"}
D["Issue"] = {"$ref": BOUNDARY + "#/$defs/Issue"}
A = ref("ArtifactRef")
D["ParameterValue"] = union(
    obj(kind=const("number"), value={"type": "number"}, unit=TEXT, endpoint=nullable(enum("Kd", "Ki", "IC50", "DC50", "Dmax", "other"))),
    obj(kind=const("integer"), value=INT, unit=const("count"), endpoint=NULL),
    obj(kind=const("boolean"), value=BOOL, unit=NULL, endpoint=NULL),
    obj(kind=const("choice"), value=TEXT, unit=NULL, endpoint=NULL),
    obj(kind=const("unset"), value=NULL, unit=NULL, endpoint=NULL, reason=TEXT))
D["PolicyField"] = obj(field_id=ID, mutability=enum("locked", "reviewable"), value=ref("ParameterValue"),
    scope_ref=A, effect=enum("fixed_universe_filter", "requires_regeneration"), required_for_execution=BOOL, rationale=TEXT)
D["PolicySnapshot"] = payload("PolicySnapshot", target_resolution_ref=A, fields=arr(ref("PolicyField")), limitations=arr(TEXT))
D["Citation"] = obj(segment_id=ID, start=INT, end=POS, quote=TEXT)
D["Claim"] = obj(claim_id=ID, statement=TEXT,
    basis_kind=enum("experiment", "computation", "author_interpretation", "researcher_hypothesis"),
    evidence_ids=arr(ID), citations=arr(ref("Citation"), 1), limitations=arr(TEXT))
D["PolicyProposal"] = obj(proposal_id=ID, field_id=ID, before=ref("ParameterValue"), proposed=ref("ParameterValue"),
    scope_ref=A, basis_kind=enum("literature_informed", "researcher_hypothesis"), supporting_claim_ids=arr(ID, 1),
    rationale=TEXT, expected_tradeoff=TEXT)
D["LiteratureAssessmentRequest"] = payload("LiteratureAssessmentRequest", input_revision=POS,
    evidence_bundle_ref=A, target_resolution_ref=A, base_policy_ref=A, allowed_field_ids=arr(ID),
    requested_segment_ids=arr(ID), research_context_ref=A,
    limits=obj(max_delivered_segments=INT, max_claims=INT, max_proposals=INT))
D["DeliveryReceipt"] = payload("DeliveryReceipt", request_ref=A, role=const("LiteratureAnalyst"),
    status=enum("completed", "incomplete", "failed", "refused"), delivered_segment_ids=arr(ID),
    response_ref=nullable(A), model=TEXT, prompt_version=TEXT, issues=arr(ref("Issue")))
D["LiteratureAssessment"] = payload("LiteratureAssessment", request_ref=A, receipt_ref=A, annotated_evidence_ref=A,
    claims=arr(ref("Claim")), summary=arr(obj(text=TEXT, claim_ids=arr(ID, 1))),
    conflicts=arr(obj(claim_ids=arr(ID, 2), reason=TEXT)), proposals=arr(ref("PolicyProposal")),
    unassessed_segment_ids=arr(ID), limitations=arr(TEXT), issues=arr(ref("Issue")))
D["CandidateIdentity"] = obj(candidate_id=ID, structure_ref=A)
D["CandidateUniverse"] = payload("CandidateUniverse", target_resolution_ref=A, candidates=arr(ref("CandidateIdentity")),
    generation_context_ref=A, limitations=arr(TEXT))
D["CandidateOutcome"] = obj(candidate_id=ID, disposition=enum("eligible", "excluded", "hold"),
    mapping_status=enum("resolved", "ambiguous", "unavailable", "failed"), mapping_ref=nullable(A),
    chemical_status=enum("valid", "invalid", "unassessed"),
    evidence_kind=enum("direct_precedent", "local_sar", "experimental_structure", "computed_hypothesis", "cross_target_reference", "unknown"),
    assessment_refs=arr(A, 1), reasons=arr(TEXT, 1), limitations=arr(TEXT))
D["Phase1Result"] = payload("Phase1Result", universe_ref=A, policy_ref=A, assessment_ref=A, evaluation_context_ref=A,
    status=enum("complete", "partial", "failed"), candidates=arr(ref("CandidateOutcome")), issues=arr(ref("Issue")))
D["PreviewRequest"] = payload("PreviewRequest", assessment_ref=A, selected_proposal_ids=arr(ID),
    base_policy_ref=A, proposed_policy_ref=A, universe_ref=A, evaluation_context_ref=A, proposal_review_ref=nullable(A))
D["CpuPreview"] = payload("CpuPreview", request_ref=A, baseline_result_ref=A, proposed_result_ref=A,
    added_candidate_ids=arr(ID), removed_candidate_ids=arr(ID), held_candidate_ids=arr(ID),
    limitations=arr(TEXT))
D["ReviewRequest"] = payload("ReviewRequest", subject_kind=enum("LiteratureAssessment", "CpuPreview"), subject_ref=A, required_claim_ids=arr(ID))
D["Finding"] = obj(finding_id=ID, severity=enum("blocking", "warning", "info"), message=TEXT,
    subject_refs=arr(A, 1), claim_ids=arr(ID))
D["ReviewRecord"] = payload("ReviewRecord", request_ref=A, role=const("Critic"),
    status=enum("completed", "incomplete", "failed", "refused"), reviewed_claim_ids=arr(ID),
    response_ref=nullable(A), findings=arr(ref("Finding")), limitations=arr(TEXT))
D["CodeCheckReport"] = payload("CodeCheckReport", request_ref=A, validator_version=TEXT,
    status=enum("complete", "incomplete", "failed"), checks=arr(obj(check_id=ID, outcome=enum("pass", "fail", "not_run"))),
    findings=arr(ref("Finding")))
D["CostScope"] = obj(execution_mode=enum("development_live", "public_explore"),
    quote_ref=A, max_gpu_seconds=INT, max_cost_krw=INT)
D["GatePacket"] = payload("GatePacket", gate=const("G1"), input_revision=POS, preview_ref=A, critic_ref=A, code_check_ref=A,
    selected_candidate_ids=arr(ID), cost_scope=ref("CostScope"), purpose=const("phase2_candidate_design"))
D["GateDecisionRequest"] = obj(packet_ref=A, input_revision=POS, decision=enum("approve", "reject"),
    acknowledged_warning_ids=arr(ID), reason=TEXT, idempotency_key=TEXT)
D["ApprovalRecord"] = payload("ApprovalRecord", packet_ref=A, input_digest={"type": "string", "pattern": "^[a-f0-9]{64}$"},
    input_revision=POS, decision=enum("approve", "reject"), acknowledged_warning_ids=arr(ID), reason=TEXT,
    actor_id=ID, actor_name=TEXT, decided_at={"type": "string", "format": "date-time"},
    idempotency_key=TEXT, dispatch_authorized=const(False))

ROOTS = [n for n, s in D.items() if "payload_type" in s.get("properties", {})] + ["GateDecisionRequest"]
SCHEMA = {"$schema": "https://json-schema.org/draft/2020-12/schema", "$id": URI,
          "title": "I3-1 offline review and G1 contracts", "$defs": D, "oneOf": [ref(n) for n in ROOTS]}

if __name__ == "__main__":
    (ROOT / "workflow.schema.json").write_text(json.dumps(SCHEMA, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Generated {len(D)} definitions, {len(ROOTS)} entrypoints")

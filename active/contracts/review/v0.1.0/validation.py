"""Offline referential/approval checks. Never dispatches, authenticates, or writes state."""
import math
from pathlib import Path
from jsonschema import Draft202012Validator, FormatChecker
from referencing import Resource
from packages.contracts import module as i1
from packages.contracts import ContractError, parse_json

ROOT = Path(__file__).resolve().parent
URI = "urn:tpd-navigator:review:0.1.0"
SCHEMA = parse_json((ROOT / "workflow.schema.json").read_bytes())
REGISTRY = i1.REGISTRY.with_resource(URI, Resource.from_contents(SCHEMA))
REQUIRED_CHECKS = {"artifact_integrity", "source_linkage", "policy_diff", "preview_consistency", "candidate_validity", "cost_scope", "execution_scope"}
require = i1.require


def validate(kind, value):
    require(kind in SCHEMA["$defs"], "UNKNOWN_REVIEW_KIND")
    validator = Draft202012Validator({"$ref": URI + "#/$defs/" + kind}, registry=REGISTRY, format_checker=FormatChecker())
    errors = list(validator.iter_errors(value))
    require(not errors, "SCHEMA:" + kind + (":" + errors[0].message if errors else ""))
    for item in i1.walk(value):
        for number in item.values():
            require(not isinstance(number, float) or math.isfinite(number), "NON_FINITE_NUMBER")
    for field, identifier in [("fields", "field_id"), ("claims", "claim_id"), ("proposals", "proposal_id"),
                              ("candidates", "candidate_id"), ("findings", "finding_id"), ("checks", "check_id")]:
        if field in value:
            i1.unique(value[field], identifier)
    return value


class ReviewReader(i1.ArtifactReader):
    def __init__(self, read_bytes):
        super().__init__(read_bytes)
        self.cache = {}

    def load(self, reference, kind):
        if kind not in SCHEMA["$defs"]:
            return self.typed(reference, kind)
        require(reference["schema_id"] == URI + "#/$defs/" + kind, "REVIEW_SCHEMA_MISMATCH")
        require(reference["media_type"] == "application/json", "REVIEW_MEDIA_MISMATCH")
        return validate(kind, parse_json(self.raw(reference)))

    def verify(self, root_ref):
        root = parse_json(self.raw(root_ref))
        require(root_ref["schema_id"] == URI + "#/$defs/" + root["payload_type"], "ROOT_SCHEMA_MISMATCH")
        pending, seen = [root_ref], {}
        while pending:
            reference = pending.pop()
            key = i1.key(reference)
            if key in seen:
                require(seen[key] == reference, "CONFLICTING_REFERENCE")
                continue
            seen[key] = reference
            require(root["data_mode"] != "real" or reference["provenance"] != "test_fixture", "FIXTURE_IN_REAL_GRAPH")
            raw = self.raw(reference)
            schema_id = reference["schema_id"]
            if schema_id.startswith((URI + "#/$defs/", i1.DOMAIN_URI + "#/$defs/")):
                kind = schema_id.split("#/$defs/")[1]
                value = self.load(reference, kind)
                for field in ("project_id", "run_id", "data_mode"):
                    if field in value:
                        require(value[field] == root[field], "GRAPH_CONTEXT_MISMATCH:" + field)
                pending.extend(i1.refs(value))
                if schema_id.startswith(URI):
                    self.check(reference, kind)
            else:
                require(schema_id.startswith("urn:tpd-navigator:raw:"), "UNSUPPORTED_GRAPH_SCHEMA")
        return {"contract_valid": True, "artifact_count": len(seen), "dispatch_authorized": False}

    def check(self, reference, kind):
        key = (reference["artifact_id"], reference["sha256"])
        if key in self.cache:
            return self.cache[key]
        value = self.load(reference, kind)
        method = getattr(self, "check_" + kind, None)
        if method:
            method(value)
        self.cache[key] = value
        return value

    def check_PolicySnapshot(self, value):
        target = self.typed(value["target_resolution_ref"], "TargetResolution")
        require(target["resolution_status"] == "resolved", "POLICY_TARGET_UNRESOLVED")

    def check_LiteratureAssessmentRequest(self, value):
        bundle = self.typed(value["evidence_bundle_ref"], "EvidenceBundle")
        policy = self.check(value["base_policy_ref"], "PolicySnapshot")
        require(bundle["target_resolution_ref"] == value["target_resolution_ref"] == policy["target_resolution_ref"], "ANALYSIS_TARGET_MISMATCH")
        target = self.typed(value["target_resolution_ref"], "TargetResolution")
        require(value["research_context_ref"] == target["query_ref"], "RESEARCH_CONTEXT_MISMATCH")
        allowed = {f["field_id"] for f in policy["fields"] if f["mutability"] == "reviewable"}
        require(set(value["allowed_field_ids"]) <= allowed, "LOCKED_OR_UNKNOWN_ALLOWED_FIELD")
        require(set(value["requested_segment_ids"]) <= {s["segment_id"] for s in bundle["segments"]}, "UNKNOWN_REQUESTED_SEGMENT")

    def check_DeliveryReceipt(self, value):
        request = self.check(value["request_ref"], "LiteratureAssessmentRequest")
        require(set(value["delivered_segment_ids"]) <= set(request["requested_segment_ids"]), "UNREQUESTED_DELIVERY")
        require(len(value["delivered_segment_ids"]) <= request["limits"]["max_delivered_segments"], "DELIVERY_LIMIT")
        if value["status"] == "completed":
            require(value["response_ref"] is not None, "MISSING_MODEL_RESPONSE")
        else:
            require(bool(value["issues"]), "INCOMPLETE_RECEIPT_WITHOUT_REASON")

    def check_LiteratureAssessment(self, value):
        request = self.check(value["request_ref"], "LiteratureAssessmentRequest")
        receipt = self.check(value["receipt_ref"], "DeliveryReceipt")
        require(receipt["request_ref"] == value["request_ref"], "ASSESSMENT_RECEIPT_MISMATCH")
        require(receipt["status"] == "completed" or not (value["claims"] or value["proposals"]), "INCOMPLETE_MODEL_OUTPUT_PROMOTED")
        original = self.typed(request["evidence_bundle_ref"], "EvidenceBundle")
        annotated = self.typed(value["annotated_evidence_ref"], "EvidenceBundle")
        for field in ("project_id", "request_ref", "target_resolution_ref", "documents", "searches"):
            require(original[field] == annotated[field], "SOURCE_BUNDLE_REWRITTEN:" + field)
        segments = {s["segment_id"]: s for s in original["segments"]}
        delivered = set(receipt["delivered_segment_ids"])
        require({s["segment_id"] for s in annotated["segments"]} == set(segments), "SOURCE_SEGMENTS_CHANGED")
        for segment in annotated["segments"]:
            old = segments[segment["segment_id"]]
            require({k: v for k, v in segment.items() if k != "read_events"} == {k: v for k, v in old.items() if k != "read_events"}, "SOURCE_SEGMENT_REWRITTEN")
            require(all(e in segment["read_events"] for e in old["read_events"]), "READ_HISTORY_REMOVED")
            for event in segment["read_events"]:
                if event not in old["read_events"]:
                    require(segment["segment_id"] in delivered and event["reader_kind"] == "llm" and event["activity_ref"] == receipt["response_ref"], "FABRICATED_READ_EVENT")
        for field, id_field in [("compounds", "compound_id"), ("evidence_records", "evidence_id")]:
            after = {item[id_field]: item for item in annotated[field]}
            require(all(after.get(item[id_field]) == item for item in original[field]), "EXISTING_EVIDENCE_REWRITTEN")
        old_ids = {r["evidence_id"] for r in original["evidence_records"]}
        for record in annotated["evidence_records"]:
            if record["evidence_id"] not in old_ids:
                require(record["review"]["state"] == "unreviewed", "ANALYST_SELF_APPROVAL")
                require(record["producer"]["kind"] == "llm" and record["producer"]["activity_ref"] == receipt["response_ref"], "EVIDENCE_PRODUCER_MISMATCH")
                require(all(any(s["segment_id"] in delivered and s["locator"] == loc for s in original["segments"])
                            for loc in record["source_locators"]), "NEW_EVIDENCE_SOURCE_NOT_DELIVERED")
        records = {r["evidence_id"] for r in annotated["evidence_records"]}
        claims = {c["claim_id"] for c in value["claims"]}
        for claim in value["claims"]:
            require(set(claim["evidence_ids"]) <= records, "UNKNOWN_EVIDENCE_ID")
            if claim["basis_kind"] in {"experiment", "computation"}:
                matched = [r for r in annotated["evidence_records"] if r["evidence_id"] in claim["evidence_ids"]]
                require(bool(matched) and all(r["basis_kind"] == claim["basis_kind"] for r in matched), "CLAIM_BASIS_WITHOUT_OBSERVATION")
            for citation in claim["citations"]:
                require(citation["segment_id"] in delivered, "CITATION_NOT_DELIVERED")
                text = self.raw(segments[citation["segment_id"]]["text_artifact_ref"]).decode("utf-8")
                require(0 <= citation["start"] < citation["end"] <= len(text), "INVALID_CITATION_SPAN")
                require(text[citation["start"]:citation["end"]] == citation["quote"], "CITATION_TEXT_MISMATCH")
        for item in value["summary"] + value["conflicts"]:
            require(set(item["claim_ids"]) <= claims, "UNKNOWN_SUMMARY_CLAIM")
        require(set(value["unassessed_segment_ids"]) == set(segments) - delivered, "ANALYSIS_COVERAGE_MISMATCH")
        require(len(value["claims"]) <= request["limits"]["max_claims"] and len(value["proposals"]) <= request["limits"]["max_proposals"], "ANALYSIS_OUTPUT_LIMIT")
        policy = self.check(request["base_policy_ref"], "PolicySnapshot")
        fields = {f["field_id"]: f for f in policy["fields"]}
        for proposal in value["proposals"]:
            require(proposal["field_id"] in request["allowed_field_ids"], "UNAUTHORIZED_POLICY_FIELD")
            field = fields[proposal["field_id"]]
            require(proposal["before"] == field["value"] and proposal["scope_ref"] == field["scope_ref"], "PROPOSAL_BASE_OR_SCOPE_MISMATCH")
            require(set(proposal["supporting_claim_ids"]) <= claims, "PROPOSAL_CLAIM_MISSING")
            require(proposal["proposed"]["kind"] != "unset", "PROPOSED_VALUE_UNSET")
            if field["value"]["kind"] != "unset":
                for k in ("kind", "unit", "endpoint"):
                    require(proposal["proposed"][k] == field["value"][k], "POLICY_VALUE_DOMAIN_CHANGED:" + k)

    def check_CandidateUniverse(self, value):
        self.typed(value["target_resolution_ref"], "TargetResolution")

    def check_Phase1Result(self, value):
        universe = self.check(value["universe_ref"], "CandidateUniverse")
        policy = self.check(value["policy_ref"], "PolicySnapshot")
        assessment = self.check(value["assessment_ref"], "LiteratureAssessment")
        request = self.check(assessment["request_ref"], "LiteratureAssessmentRequest")
        require(universe["target_resolution_ref"] == policy["target_resolution_ref"] == request["target_resolution_ref"], "CANDIDATE_TARGET_MISMATCH")
        require({c["candidate_id"] for c in universe["candidates"]} == {c["candidate_id"] for c in value["candidates"]}, "CANDIDATE_COVERAGE_MISMATCH")
        for candidate in value["candidates"]:
            require(candidate["mapping_status"] != "resolved" or candidate["mapping_ref"] is not None, "RESOLVED_MAPPING_WITHOUT_ARTIFACT")
            require(candidate["disposition"] != "eligible" or (candidate["chemical_status"] == "valid" and candidate["mapping_status"] == "resolved"), "INVALID_CANDIDATE_ELIGIBLE")
        require(value["status"] == "complete" or bool(value["issues"]), "PARTIAL_CANDIDATES_WITHOUT_REASON")

    def check_PreviewRequest(self, value):
        assessment = self.check(value["assessment_ref"], "LiteratureAssessment")
        analysis_request = self.check(assessment["request_ref"], "LiteratureAssessmentRequest")
        require(value["base_policy_ref"] == analysis_request["base_policy_ref"], "PREVIEW_BASE_POLICY_MISMATCH")
        before = self.check(value["base_policy_ref"], "PolicySnapshot")
        after = self.check(value["proposed_policy_ref"], "PolicySnapshot")
        proposals = {p["proposal_id"]: p for p in assessment["proposals"]}
        require(set(value["selected_proposal_ids"]) <= set(proposals), "UNKNOWN_SELECTED_PROPOSAL")
        selected = [proposals[n] for n in value["selected_proposal_ids"]]
        if selected or value["proposal_review_ref"]:
            require(value["proposal_review_ref"] is not None, "PROPOSAL_REVIEW_REQUIRED")
            review = self.check(value["proposal_review_ref"], "ReviewRecord")
            subject = self.check(review["request_ref"], "ReviewRequest")
            require(subject["subject_kind"] == "LiteratureAssessment" and subject["subject_ref"] == value["assessment_ref"], "PROPOSAL_REVIEW_TARGET_MISMATCH")
            require(review["status"] == "completed" and not any(f["severity"] == "blocking" for f in review["findings"]), "PROPOSAL_REVIEW_BLOCKED")
        require(len({p["field_id"] for p in selected}) == len(selected), "CONFLICTING_FIELD_PROPOSALS")
        expected = {f["field_id"]: dict(f) for f in before["fields"]}
        for proposal in selected:
            field = expected[proposal["field_id"]]
            require(field["mutability"] == "reviewable", "LOCKED_POLICY_CHANGED")
            require(field["effect"] == "fixed_universe_filter", "CANDIDATE_REGENERATION_REQUIRED")
            field["value"] = proposal["proposed"]
        require({k: v for k, v in before.items() if k != "fields"} == {k: v for k, v in after.items() if k != "fields"}, "UNDECLARED_POLICY_CHANGE")
        require({f["field_id"]: f for f in after["fields"]} == expected, "UNDECLARED_POLICY_CHANGE")

    def check_CpuPreview(self, value):
        request = self.check(value["request_ref"], "PreviewRequest")
        before = self.check(value["baseline_result_ref"], "Phase1Result")
        after = self.check(value["proposed_result_ref"], "Phase1Result")
        require(before["policy_ref"] == request["base_policy_ref"] and after["policy_ref"] == request["proposed_policy_ref"], "PREVIEW_RESULT_POLICY_MISMATCH")
        for key in ("universe_ref", "assessment_ref", "evaluation_context_ref"):
            require(before[key] == after[key] == request[key], "PREVIEW_CONTEXT_MISMATCH:" + key)
        if self.load(before["policy_ref"], "PolicySnapshot") == self.load(after["policy_ref"], "PolicySnapshot"):
            require(before["candidates"] == after["candidates"], "UNCHANGED_POLICY_DIFFERENT_RESULT")
        a = {c["candidate_id"] for c in before["candidates"] if c["disposition"] == "eligible"}
        b = {c["candidate_id"] for c in after["candidates"] if c["disposition"] == "eligible"}
        hold = {c["candidate_id"] for c in after["candidates"] if c["disposition"] == "hold"}
        require(set(value["added_candidate_ids"]) == b-a and set(value["removed_candidate_ids"]) == a-b and set(value["held_candidate_ids"]) == hold, "PREVIEW_DELTA_MISMATCH")

    def check_ReviewRequest(self, value):
        if value["subject_kind"] == "CpuPreview":
            preview = self.check(value["subject_ref"], "CpuPreview")
            request = self.check(preview["request_ref"], "PreviewRequest")
            assessment = self.check(request["assessment_ref"], "LiteratureAssessment")
        else:
            assessment = self.check(value["subject_ref"], "LiteratureAssessment")
        require(set(value["required_claim_ids"]) == {c["claim_id"] for c in assessment["claims"]}, "REVIEW_CLAIMS_OMITTED")

    def check_ReviewRecord(self, value):
        request = self.check(value["request_ref"], "ReviewRequest")
        require(set(value["reviewed_claim_ids"]) <= set(request["required_claim_ids"]), "UNKNOWN_REVIEWED_CLAIM")
        for finding in value["findings"]:
            require(set(finding["claim_ids"]) <= set(request["required_claim_ids"]), "UNKNOWN_FINDING_CLAIM")
        if value["status"] == "completed":
            require(value["response_ref"] is not None and set(value["reviewed_claim_ids"]) == set(request["required_claim_ids"]), "INCOMPLETE_CRITIC_COVERAGE")

    def check_CodeCheckReport(self, value):
        self.check(value["request_ref"], "ReviewRequest")

    def check_GatePacket(self, value):
        preview = self.check(value["preview_ref"], "CpuPreview")
        request = self.check(preview["request_ref"], "PreviewRequest")
        assessment = self.check(request["assessment_ref"], "LiteratureAssessment")
        analysis_request = self.check(assessment["request_ref"], "LiteratureAssessmentRequest")
        require(value["input_revision"] == analysis_request["input_revision"], "PACKET_REVISION_MISMATCH")
        critic = self.check(value["critic_ref"], "ReviewRecord")
        code = self.check(value["code_check_ref"], "CodeCheckReport")
        require(critic["request_ref"] == code["request_ref"], "CODE_CRITIC_TARGET_MISMATCH")
        review_request = self.check(critic["request_ref"], "ReviewRequest")
        require(review_request["subject_kind"] == "CpuPreview" and review_request["subject_ref"] == value["preview_ref"], "REVIEW_PREVIEW_MISMATCH")
        after = self.check(preview["proposed_result_ref"], "Phase1Result")
        require(set(value["selected_candidate_ids"]) <= {c["candidate_id"] for c in after["candidates"]}, "SELECTED_CANDIDATE_UNKNOWN")
        ids = [f["finding_id"] for f in critic["findings"] + code["findings"]]
        require(len(ids) == len(set(ids)), "DUPLICATE_FINDING_ID")
        require(value["cost_scope"]["execution_mode"] != "public_explore" or value["cost_scope"]["max_gpu_seconds"] == 0, "PUBLIC_GPU_FORBIDDEN")

    def readiness(self, packet_ref):
        self.verify(packet_ref)
        packet = self.check(packet_ref, "GatePacket")
        critic = self.check(packet["critic_ref"], "ReviewRecord")
        code = self.check(packet["code_check_ref"], "CodeCheckReport")
        preview = self.check(packet["preview_ref"], "CpuPreview")
        before = self.check(preview["baseline_result_ref"], "Phase1Result")
        after = self.check(preview["proposed_result_ref"], "Phase1Result")
        reasons = []
        if critic["status"] != "completed": reasons.append("CRITIC_NOT_COMPLETED")
        if code["status"] != "complete": reasons.append("CODE_CHECK_NOT_COMPLETED")
        passed = {c["check_id"] for c in code["checks"] if c["outcome"] == "pass"}
        if not REQUIRED_CHECKS <= passed or any(c["outcome"] != "pass" for c in code["checks"]): reasons.append("REQUIRED_CODE_CHECK_NOT_PASSED")
        if before["status"] != "complete" or after["status"] != "complete": reasons.append("PREVIEW_NOT_COMPLETE")
        findings = critic["findings"] + code["findings"]
        preview_request = self.check(preview["request_ref"], "PreviewRequest")
        if preview_request["proposal_review_ref"]:
            findings += self.check(preview_request["proposal_review_ref"], "ReviewRecord")["findings"]
        identifiers = [f["finding_id"] for f in findings]
        require(len(identifiers) == len(set(identifiers)), "DUPLICATE_FINDING_ID")
        if any(f["severity"] == "blocking" for f in findings): reasons.append("BLOCKING_FINDING")
        selected = set(packet["selected_candidate_ids"])
        eligible = {c["candidate_id"] for c in after["candidates"] if c["disposition"] == "eligible"}
        if not selected or not selected <= eligible: reasons.append("NO_VALID_SELECTED_CANDIDATES")
        policy = self.check(preview_request["proposed_policy_ref"], "PolicySnapshot")
        if any(f["required_for_execution"] and f["value"]["kind"] == "unset" for f in policy["fields"]): reasons.append("REQUIRED_POLICY_UNSET")
        warnings = [f["finding_id"] for f in findings if f["severity"] == "warning"]
        warnings += ["candidate:" + c["candidate_id"] + ":uncertainty" for c in after["candidates"]
                     if c["candidate_id"] in selected and (c["limitations"] or c["evidence_kind"] in {"computed_hypothesis", "cross_target_reference", "unknown"})]
        return {"review_ready": True, "approval_contract_ready": not reasons, "blocking_reasons": reasons,
                "warning_ids": warnings, "dispatch_authorized": False}

    def check_ApprovalRecord(self, value):
        packet = self.check(value["packet_ref"], "GatePacket")
        require(value["input_digest"] == value["packet_ref"]["sha256"] and value["input_revision"] == packet["input_revision"], "APPROVAL_BINDING_MISMATCH")
        if value["decision"] == "approve":
            ready = self.readiness(value["packet_ref"])
            require(ready["approval_contract_ready"] and set(value["acknowledged_warning_ids"]) == set(ready["warning_ids"]), "INVALID_STORED_APPROVAL")


def validate_decision(request, reader, *, current_packet_ref, current_revision, authenticated_actor, decided_at):
    """Return an offline record candidate. M2 must authenticate and persist with CAS later."""
    validate("GateDecisionRequest", request)
    require(request["packet_ref"] == current_packet_ref and request["input_revision"] == current_revision, "STALE_GATE_DECISION")
    packet = reader.check(current_packet_ref, "GatePacket")
    require(packet["input_revision"] == current_revision, "STALE_PACKET_REVISION")
    require(authenticated_actor.get("kind") == "human" and authenticated_actor.get("id") and str(authenticated_actor.get("name", "")).strip()
            and authenticated_actor.get("project_id") == packet["project_id"], "HUMAN_ACTOR_REQUIRED")
    ready = reader.readiness(current_packet_ref)
    warnings = set(request["acknowledged_warning_ids"])
    require(warnings <= set(ready["warning_ids"]), "UNKNOWN_WARNING_ACKNOWLEDGEMENT")
    if request["decision"] == "approve":
        require(ready["approval_contract_ready"], "GATE_APPROVAL_BLOCKED")
        require(warnings == set(ready["warning_ids"]), "WARNING_NOT_ACKNOWLEDGED")
    record = {k: packet[k] for k in ("payload_version", "project_id", "run_id", "data_mode")}
    record.update(payload_type="ApprovalRecord", packet_ref=current_packet_ref, input_digest=current_packet_ref["sha256"],
        input_revision=current_revision, decision=request["decision"], acknowledged_warning_ids=request["acknowledged_warning_ids"],
        reason=request["reason"], idempotency_key=request["idempotency_key"], actor_id=authenticated_actor["id"],
        actor_name=authenticated_actor["name"], decided_at=decided_at, dispatch_authorized=False)
    return validate("ApprovalRecord", record)

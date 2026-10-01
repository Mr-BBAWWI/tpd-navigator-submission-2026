"""Offline I1 contract validation. No dispatch, authorization, or science inference."""
import hashlib
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker
from referencing import Registry, Resource

ROOT = Path(__file__).resolve().parent
DOMAIN_URI = "urn:tpd-navigator:payloads:0.1.0"
BOUNDARY_URI = "urn:tpd-navigator:boundary:0.1.0"
ARTIFACT_KEYS = {"artifact_id", "version", "sha256", "media_type", "schema_id", "provenance"}
OP_PAYLOADS = {"resolve_target": ("TargetQuery", "TargetResolution"),
               "collect_literature": ("LiteratureCollectionRequest", "EvidenceBundle")}


class ContractError(ValueError):
    pass


def require(condition, code):
    if not condition:
        raise ContractError(code)


def parse_json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "DUPLICATE_JSON_KEY")
            result[key] = value
        return result

    def constant(_):
        raise ContractError("NON_FINITE_JSON_NUMBER")

    def finite_float(text):
        value = float(text)
        require(math.isfinite(value), "NON_FINITE_JSON_NUMBER")
        return value

    try:
        return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant, parse_float=finite_float)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("INVALID_JSON") from exc


DOMAIN = parse_json((ROOT / "domain.schema.json").read_bytes())
BOUNDARY = parse_json((ROOT.parents[1] / "v0.1.0/module_boundary.schema.json").read_bytes())
REGISTRY = Registry().with_resources([
    (DOMAIN_URI, Resource.from_contents(DOMAIN)),
    (BOUNDARY_URI, Resource.from_contents(BOUNDARY)),
])


def walk(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk(child)


def refs(value):
    return [item for item in walk(value) if set(item) == ARTIFACT_KEYS]


def key(reference):
    return reference["artifact_id"], reference["version"]


def unique(items, field):
    values = [item[field] for item in items]
    require(len(values) == len(set(values)), "DUPLICATE_" + field.upper())


def validate_payload(kind, value):
    uri = BOUNDARY_URI if kind in {"JobSpec", "ModuleResult"} else DOMAIN_URI
    require(kind in (BOUNDARY if uri == BOUNDARY_URI else DOMAIN)["$defs"], "UNKNOWN_PAYLOAD_KIND")
    validator = Draft202012Validator({"$ref": uri + "#/$defs/" + kind}, registry=REGISTRY, format_checker=FormatChecker())
    errors = sorted(validator.iter_errors(value), key=lambda e: str(list(e.path)))
    if errors:
        error = errors[0]
        raise ContractError(f"SCHEMA:{kind}:{list(error.path)}:{error.message}")
    for item in walk(value):
        require(all(not isinstance(v, float) or math.isfinite(v) for v in item.values()), "NON_FINITE_JSON_NUMBER")
        if item.get("kind") == "interval":
            require(item["lower"] <= item["upper"], "REVERSED_INTERVAL")
            if item["lower"] == item["upper"]:
                require(item["lower_inclusive"] and item["upper_inclusive"], "EMPTY_INTERVAL")
        if item.get("units") == "normalized" and "x0" in item:
            require(item["x0"] < item["x1"] and item["y0"] < item["y1"], "REVERSED_BBOX")
    if kind == "TargetResolution":
        unique(value["candidates"], "candidate_id")
        selected = [c for c in value["candidates"] if c["candidate_id"] == value["selected_candidate_id"]]
        if value["resolution_status"] == "resolved":
            require(len(selected) == 1, "SELECTED_TARGET_NOT_FOUND")
            if value["selection"]["kind"] == "unique_match":
                require(len(value["candidates"]) == 1, "MULTIPLE_TARGETS_AUTO_SELECTED")
            if value["scope_status"] == "within_current_scope":
                require(selected[0]["identity"]["taxon_id"] == 9606, "NON_HUMAN_SCOPE")
    if kind == "LiteratureCollectionRequest":
        unique(value["questions"], "question_id")
        require(all(r["media_type"] == "application/pdf" for r in value["attached_source_refs"]), "UNSUPPORTED_ATTACHMENT_MEDIA")
    if kind == "InputManifest":
        require(value["payload_ref"] in value["input_artifacts"], "PAYLOAD_NOT_DECLARED")
        keys = [key(r) for r in value["input_artifacts"]]
        require(len(keys) == len(set(keys)), "CONFLICTING_INPUT_REFERENCE")
    if kind == "EvidenceBundle":
        validate_bundle(value)
    return value


def validate_bundle(bundle):
    for collection, field in [("documents", "document_id"), ("segments", "segment_id"),
                               ("compounds", "compound_id"), ("evidence_records", "evidence_id"), ("searches", "question_id")]:
        unique(bundle[collection], field)
    docs = {d["document_id"]: d for d in bundle["documents"]}
    compounds = {c["compound_id"]: c for c in bundle["compounds"]}
    records = {r["evidence_id"]: r for r in bundle["evidence_records"]}
    mime = {"jats_xml": {"application/xml", "text/xml", "application/jats+xml"},
            "pdf": {"application/pdf"}, "abstract_text": {"text/plain"}, "metadata": {"application/json"}}
    for doc in docs.values():
        require(doc["project_id"] == bundle["project_id"], "CROSS_PROJECT_DOCUMENT")
        if doc["source_artifact_ref"]:
            require(doc["source_artifact_ref"]["media_type"] in mime[doc["format"]], "SOURCE_MEDIA_MISMATCH")
        for relation in doc["relations"]:
            other = relation["related_document_id"]
            require(other is None or (other in docs and other != doc["document_id"]), "BROKEN_DOCUMENT_RELATION")
    for item in walk(bundle):
        if "source_artifact_ref" not in item or "extraction_snapshot_ref" not in item:
            continue
        require(item["document_id"] in docs, "LOCATOR_DOCUMENT_NOT_FOUND")
        doc = docs[item["document_id"]]
        require(doc["acquisition_status"] == "retrieved", "LOCATOR_SOURCE_UNAVAILABLE")
        require(doc["extraction"]["status"] in {"complete", "partial"}, "LOCATOR_EXTRACTION_UNAVAILABLE")
        require(item["document_revision"] == doc["revision"], "STALE_DOCUMENT_REVISION")
        require(item["source_artifact_ref"] == doc["source_artifact_ref"], "STALE_SOURCE_REFERENCE")
        require(item["extraction_snapshot_ref"] == doc["extraction"]["snapshot_ref"], "STALE_EXTRACTION_REFERENCE")
        require(item["extractor_version"] == doc["extraction"]["extractor_version"], "STALE_EXTRACTOR")
        require(item["kind"] == doc["format"], "LOCATOR_FORMAT_MISMATCH")
    for compound in compounds.values():
        require(compound["document_id"] in docs, "COMPOUND_DOCUMENT_NOT_FOUND")
        require(all(l["document_id"] == compound["document_id"] for l in compound["structure_locators"]), "COMPOUND_LOCATOR_MISMATCH")
    for record in records.values():
        require(all(c in compounds for c in record["compound_ids"]), "COMPOUND_NOT_FOUND")
        source_docs = {l["document_id"] for l in record["source_locators"]}
        require(all(compounds[c]["document_id"] in source_docs for c in record["compound_ids"]), "COMPOUND_SOURCE_MISMATCH")
        if record["attachment"]:
            require(record["attachment"]["compound_id"] in record["compound_ids"], "ATTACHMENT_COMPOUND_MISMATCH")
        for link in record["comparison_links"]:
            require(link["other_evidence_id"] in records and link["other_evidence_id"] != record["evidence_id"], "BROKEN_COMPARISON_LINK")
        if record["basis_kind"] in {"experiment", "computation"}:
            require(any(l["kind"] != "metadata" for l in record["source_locators"]), "METADATA_IS_NOT_EXPERIMENTAL_EVIDENCE")
        if record["basis_kind"] == "experiment":
            require(record["experiment_origin"]["state"] != "not_applicable", "EXPERIMENT_ORIGIN_REQUIRED")
        else:
            require(record["experiment_origin"]["state"] == "not_applicable", "NON_EXPERIMENT_ORIGIN")
        if record["producer"]["kind"] == "llm":
            for loc in record["source_locators"]:
                require(any(s["locator"] == loc and any(e["reader_kind"] == "llm" for e in s["read_events"])
                            for s in bundle["segments"]), "LLM_SOURCE_NOT_MARKED_READ")
    for search in bundle["searches"]:
        require(all(d in docs for d in search["retrieved_document_ids"]), "SEARCH_DOCUMENT_NOT_FOUND")
        require(all(docs[d]["acquisition_status"] == "retrieved" for d in search["retrieved_document_ids"]), "SEARCH_RETRIEVAL_STATUS_MISMATCH")
    expected = {
        "stored_document_count": len(docs),
        "fulltext_acquired_count": sum(d["acquisition_status"] == "retrieved" and d["access_level"] == "fulltext" for d in docs.values()),
        "extracted_segment_count": len(bundle["segments"]),
        "llm_read_segment_count": sum(any(e["reader_kind"] == "llm" for e in s["read_events"]) for s in bundle["segments"]),
        "evidence_record_count": len(records),
        "reviewed_record_count": sum(r["review"]["state"] != "unreviewed" for r in records.values()),
        "unavailable_document_count": sum(d["acquisition_status"] != "retrieved" for d in docs.values()),
    }
    require(all(bundle["coverage"][k] == v for k, v in expected.items()), "COVERAGE_MISMATCH")
    experiments = [r["experiment_origin"] for r in records.values() if r["basis_kind"] == "experiment"]
    reported = bundle["coverage"]["independent_experiment_count"]
    if any(e["state"] == "unknown" for e in experiments):
        require(reported["state"] == "unknown", "UNKNOWN_EXPERIMENT_COUNT_INVENTED")
    else:
        require(reported["state"] == "known" and reported["value"] == len({e["origin_id"] for e in experiments}), "EXPERIMENT_COUNT_MISMATCH")


class ArtifactReader:
    """Adapter port: callback(ref) -> bytes. Caller supplies access control."""
    def __init__(self, read_bytes):
        self.read_bytes = read_bytes

    def raw(self, reference):
        validate_payload("ArtifactRef", reference)
        data = self.read_bytes(reference)
        require(isinstance(data, bytes), "ARTIFACT_BYTES_REQUIRED")
        require(hashlib.sha256(data).hexdigest() == reference["sha256"], "ARTIFACT_HASH_MISMATCH")
        return data

    def typed(self, reference, kind):
        require(reference["schema_id"] == DOMAIN_URI + "#/$defs/" + kind, "ARTIFACT_SCHEMA_MISMATCH")
        require(reference["media_type"] == "application/json", "PAYLOAD_MEDIA_MISMATCH")
        return validate_payload(kind, parse_json(self.raw(reference)))

    def graph(self, roots, data_mode):
        seen = {}
        pending = list(roots)
        while pending:
            reference = pending.pop()
            identity = key(reference)
            if identity in seen:
                require(seen[identity] == reference, "CONFLICTING_ARTIFACT_REFERENCE")
                continue
            seen[identity] = reference
            require(data_mode != "real" or reference["provenance"] != "test_fixture", "FIXTURE_IN_REAL_GRAPH")
            raw = self.raw(reference)
            schema_id = reference["schema_id"]
            if schema_id.startswith("urn:tpd-navigator:payloads:"):
                require(schema_id.startswith(DOMAIN_URI + "#/$defs/"), "UNSUPPORTED_PAYLOAD_VERSION")
                kind = schema_id.split("#/$defs/", 1)[1]
                value = validate_payload(kind, parse_json(raw))
                pending.extend(refs(value))
        return seen


def validate_exchange(job, reader, result=None):
    """Validate two I1 operations against stored bytes; never dispatch work."""
    validate_payload("JobSpec", job)
    require(job["operation_id"] in OP_PAYLOADS, "I1_OPERATION_NOT_IMPLEMENTED")
    manifest = reader.typed(job["input_manifest"], "InputManifest")
    require(job["input_digest"] == job["input_manifest"]["sha256"], "INPUT_DIGEST_MISMATCH")
    for field in ("project_id", "run_id", "job_id", "input_revision", "operation_id", "execution_mode", "data_mode",
                  "input_artifacts", "parameters_ref", "policy_ref"):
        require(manifest[field] == job[field], "MANIFEST_MISMATCH:" + field)
    request_kind, output_kind = OP_PAYLOADS[job["operation_id"]]
    request = reader.typed(manifest["payload_ref"], request_kind)
    declared = job["input_artifacts"] + [job["parameters_ref"], job["policy_ref"]]
    if job["approval_ref"]:
        declared.append(job["approval_ref"])
    require(all(r in declared for r in refs(request)), "UNDECLARED_PAYLOAD_INPUT")
    reader.graph([job["input_manifest"]] + ([job["approval_ref"]] if job["approval_ref"] else []), job["data_mode"])
    if job["operation_id"] == "collect_literature":
        require(request["collection_policy_ref"] == job["policy_ref"], "COLLECTION_POLICY_MISMATCH")
        target = reader.typed(request["target_resolution_ref"], "TargetResolution")
        require(target["resolution_status"] == "resolved", "UNRESOLVED_SEARCH_TARGET")
        require(target["scope_status"] != "outside_current_scope", "UNSUPPORTED_SEARCH_TARGET")
        selected = next(c for c in target["candidates"] if c["candidate_id"] == target["selected_candidate_id"])
        require(selected["identity"]["taxon_id"] == 9606, "NON_HUMAN_SEARCH_TARGET")
        collection_policy = reader.typed(job["policy_ref"], "CollectionPolicy")
        require(len(request["attached_source_refs"]) <= collection_policy["max_documents"], "ATTACHMENT_DOCUMENT_LIMIT")
        require(sum(len(reader.raw(r)) for r in request["attached_source_refs"]) <= collection_policy["max_source_bytes"], "ATTACHMENT_BYTE_LIMIT")
    if result is not None:
        validate_payload("ModuleResult", result)
        for field in ("contract_version", "project_id", "run_id", "job_id", "attempt", "input_revision", "input_digest", "operation_id", "data_mode"):
            require(result[field] == job[field], "RESULT_CORRELATION_MISMATCH:" + field)
        require(result["producer"]["module"] == ("science" if job["operation_id"] == "resolve_target" else "literature"), "RESULT_PRODUCER_MISMATCH")
        if job["data_mode"] == "test_fixture":
            require(all(r["provenance"] == "test_fixture" for r in result["output_artifacts"]), "FIXTURE_OUTPUT_NOT_MARKED")
        reader.graph(result["output_artifacts"], job["data_mode"])
        outputs = [r for r in result["output_artifacts"] if r["schema_id"] == DOMAIN_URI + "#/$defs/" + output_kind]
        if result["status"] in {"succeeded", "partial"}:
            require(len(outputs) == 1, "PRIMARY_OUTPUT_REQUIRED")
        for reference in outputs:
            output = reader.typed(reference, output_kind)
            request_field = "query_ref" if output_kind == "TargetResolution" else "request_ref"
            require(output[request_field] == manifest["payload_ref"], "OUTPUT_REQUEST_MISMATCH")
            if output_kind == "TargetResolution" and output["resolution_status"] == "resolved":
                selected = next(c for c in output["candidates"] if c["candidate_id"] == output["selected_candidate_id"])
                require(selected["identity"]["taxon_id"] == request["taxon_id"], "RESOLVED_TAXON_MISMATCH")
            if output_kind == "EvidenceBundle":
                require(output["project_id"] == job["project_id"], "BUNDLE_PROJECT_MISMATCH")
                require(output["target_resolution_ref"] == request["target_resolution_ref"], "BUNDLE_TARGET_MISMATCH")
                questions = {q["question_id"]: q for q in request["questions"]}
                require({s["question_id"] for s in output["searches"]} == set(questions), "SEARCH_COVERAGE_MISMATCH")
                require(all(s["query"] == questions[s["question_id"]]["query"] for s in output["searches"]), "SEARCH_QUERY_MISMATCH")
                require(len(output["documents"]) <= collection_policy["max_records"], "RECORD_LIMIT_EXCEEDED")
                require(sum(d["acquisition_status"] == "retrieved" for d in output["documents"]) <= collection_policy["max_documents"], "DOCUMENT_LIMIT_EXCEEDED")
                require(len(output["segments"]) <= collection_policy["max_segments"], "SEGMENT_LIMIT_EXCEEDED")
                source_refs = {key(d["source_artifact_ref"]): d["source_artifact_ref"] for d in output["documents"] if d["source_artifact_ref"]}
                require(sum(len(reader.raw(r)) for r in source_refs.values()) <= collection_policy["max_source_bytes"], "SOURCE_BYTE_LIMIT_EXCEEDED")
                uploaded = [d["source_artifact_ref"] for d in output["documents"] if any(p["kind"] == "upload" for p in d["acquisition_paths"])]
                require(all(r in uploaded for r in request["attached_source_refs"]) and all(r in request["attached_source_refs"] for r in uploaded), "ATTACHMENT_COVERAGE_MISMATCH")
    return {"contract_valid": True, "operation_id": job["operation_id"], "dispatch_authorized": False}

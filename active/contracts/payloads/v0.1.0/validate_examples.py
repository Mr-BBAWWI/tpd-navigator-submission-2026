"""Executable I1 contract tests against synthetic artifacts; no external calls."""
import copy
import hashlib
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from build_schemas import SCHEMA
from validation import (ArtifactReader, BOUNDARY, DOMAIN, DOMAIN_URI, ContractError,
                        key, parse_json, validate_exchange, validate_payload, walk)

ROOT = Path(__file__).resolve().parent
EXAMPLES = ROOT / "examples"


class FixtureStore:
    def __init__(self):
        index = parse_json((EXAMPLES / "artifact-index.json").read_bytes())
        self.data = {}
        self.references = {}
        for entry in index["artifacts"]:
            path = (EXAMPLES / entry["path"]).resolve()
            if not path.is_relative_to(EXAMPLES.resolve()):
                raise ContractError("FIXTURE_PATH_OUTSIDE_STORE")
            reference = entry["reference"]
            if key(reference) in self.data:
                raise ContractError("DUPLICATE_FIXTURE_ID")
            self.data[key(reference)] = path.read_bytes()
            self.references[reference["artifact_id"]] = reference

    def read(self, reference):
        try:
            return self.data[key(reference)]
        except KeyError as exc:
            raise ContractError("ARTIFACT_NOT_FOUND") from exc

    def payload(self, name):
        return parse_json(self.read(self.references["fixture-" + name]))


def main():
    Draft202012Validator.check_schema(DOMAIN)
    Draft202012Validator.check_schema(BOUNDARY)
    assert DOMAIN == SCHEMA, "Generated schema is stale; run build_schemas.py"
    store = FixtureStore()
    reader = ArtifactReader(store.read)
    positive = []
    negative = []

    def accept(name, action):
        action()
        positive.append(name)

    def reject(name, action, expected):
        try:
            action()
        except ContractError as exc:
            assert str(exc).startswith(expected), f"{name}: unexpected failure {exc}"
            negative.append({"case": name, "expected_error": expected})
        else:
            raise AssertionError(f"Invalid example accepted: {name}")

    def mutation(name, kind, base, change, expected="SCHEMA:"):
        value = copy.deepcopy(base)
        change(value)
        reject(name, lambda: validate_payload(kind, value), expected)

    for reference in store.references.values():
        accept("hash:" + reference["artifact_id"], lambda r=reference: reader.raw(r))
        if reference["schema_id"].startswith(DOMAIN_URI + "#/$defs/"):
            kind = reference["schema_id"].split("#/$defs/", 1)[1]
            accept("payload:" + reference["artifact_id"], lambda r=reference, k=kind: reader.typed(r, k))
    exchanges = {}
    for name in ("target", "literature"):
        job = parse_json((EXAMPLES / (name + ".job.json")).read_bytes())
        result = parse_json((EXAMPLES / (name + ".result.json")).read_bytes())
        accept("exchange:" + name, lambda j=job, r=result: validate_exchange(j, reader, r))
        exchanges[name] = job, result

    target = store.payload("target-resolution")
    ambiguous = copy.deepcopy(target)
    candidate = copy.deepcopy(ambiguous["candidates"][0])
    candidate["candidate_id"] = "candidate-beta"
    candidate["identity"]["uniprot_accession"] = "FIXTURE_BETA"
    ambiguous["candidates"].append(candidate)
    ambiguous.update(resolution_status="ambiguous", selected_candidate_id=None, selection=None,
                     scope_status="undetermined", scope_reasons=["Multiple synthetic matches."])
    accept("ambiguous target is a valid result", lambda: validate_payload("TargetResolution", ambiguous))
    not_found = copy.deepcopy(ambiguous)
    not_found.update(resolution_status="not_found", candidates=[])
    accept("no target is a valid result", lambda: validate_payload("TargetResolution", not_found))
    nonhuman = copy.deepcopy(target)
    nonhuman["candidates"][0]["identity"]["taxon_id"] = 10090
    nonhuman.update(scope_status="outside_current_scope", scope_reasons=["Human-only product scope."])
    accept("unsupported target identity preserved", lambda: validate_payload("TargetResolution", nonhuman))

    bundle = store.payload("evidence-bundle")
    locator = bundle["segments"][0]["locator"]
    base_locator = {k: v for k, v in locator.items() if k in {"document_id", "document_revision", "source_artifact_ref", "extraction_snapshot_ref", "extractor_version"}}
    pdf_locator = {**base_locator, "kind": "pdf", "page": 1, "segment_id": "pdf-segment", "bbox": {"x0": 0.1, "y0": 0.1, "x1": 0.9, "y1": 0.2, "units": "normalized", "origin": "top_left"}, "table_id": None, "figure_id": None}
    abstract_locator = {**base_locator, "kind": "abstract_text", "record_id": "fixture-abstract", "field": "abstract", "segment_id": "abstract-segment"}
    metadata_locator = {**base_locator, "kind": "metadata", "record_id": "fixture-metadata", "field": "title"}
    for label, value in [("pdf", pdf_locator), ("abstract", abstract_locator), ("metadata", metadata_locator)]:
        accept("locator shape:" + label, lambda v=value: validate_payload("SourceLocator", v))
    interval = {"kind": "interval", "lower": 10, "upper": 20, "lower_inclusive": True, "upper_inclusive": False,
                "unit": {"state": "known", "value": "nM", "reason": None}, "raw_text": "[10,20) nM"}
    accept("interval preserved", lambda: validate_payload("Measurement", interval))
    qualitative = {"kind": "qualitative", "value": "not_detected", "raw_text": "Not detected in this assay."}
    accept("not detected is an observation", lambda: validate_payload("Measurement", qualitative))
    missing = {"kind": "missing", "state": "not_measured", "value": None, "reason": "Assay was not performed."}
    accept("not measured is not zero", lambda: validate_payload("Measurement", missing))
    llm_bundle = copy.deepcopy(bundle)
    llm_bundle["evidence_records"][0]["producer"]["kind"] = "llm"
    llm_bundle["segments"][0]["read_events"] = [{"reader_kind": "llm", "reader_version": "fixture-only", "activity_ref": store.references["fixture-extraction-activity"]}]
    llm_bundle["coverage"]["llm_read_segment_count"] = 1
    accept("LLM producer does not change experimental basis", lambda: validate_payload("EvidenceBundle", llm_bundle))
    other_target = copy.deepcopy(bundle)
    other_target["evidence_records"][0]["target"]["identity"]["uniprot_accession"] = "FIXTURE_OTHER_TARGET"
    accept("other target evidence is preserved", lambda: validate_payload("EvidenceBundle", other_target))
    no_records = copy.deepcopy(bundle)
    no_records.update(evidence_records=[], compounds=[])
    no_records["coverage"].update(evidence_record_count=0, independent_experiment_count={"state": "known", "value": 0, "reason": None})
    accept("collection can return sources before analysis", lambda: validate_payload("EvidenceBundle", no_records))
    empty = copy.deepcopy(no_records)
    empty.update(documents=[], segments=[], issues=[])
    empty["searches"][0].update(total_hits={"state": "known", "value": 0, "reason": None}, retrieved_document_ids=[])
    empty["coverage"].update(stored_document_count=0, fulltext_acquired_count=0, extracted_segment_count=0)
    accept("zero search hits is not a failed search", lambda: validate_payload("EvidenceBundle", empty))
    duplicate_experiment = copy.deepcopy(bundle)
    duplicate_experiment["evidence_records"][0]["experiment_origin"].update(state="resolved", origin_id="original-experiment-1")
    second = copy.deepcopy(duplicate_experiment["evidence_records"][0])
    second["evidence_id"] = "evidence-same-experiment"
    duplicate_experiment["evidence_records"].append(second)
    duplicate_experiment["coverage"].update(evidence_record_count=2, independent_experiment_count={"state": "known", "value": 1, "reason": None})
    accept("duplicate evidence counts one original experiment", lambda: validate_payload("EvidenceBundle", duplicate_experiment))

    mutation("ambiguous cannot select automatically", "TargetResolution", ambiguous, lambda x: x.update(selected_candidate_id="candidate-alpha"))
    mutation("resolved selection must exist", "TargetResolution", target, lambda x: x.update(selected_candidate_id="missing"), "SELECTED_TARGET_NOT_FOUND")
    mutation("unique match cannot hide second candidate", "TargetResolution", target, lambda x: x["candidates"].append(candidate), "MULTIPLE_TARGETS_AUTO_SELECTED")
    mutation("nonhuman cannot be in scope", "TargetResolution", target, lambda x: x["candidates"][0]["identity"].update(taxon_id=10090), "NON_HUMAN_SCOPE")
    mutation("unknown with invented value", "TextValue", {"state": "unknown", "value": None, "reason": "Missing"}, lambda x: x.update(value="invented"))
    mutation("missing measurement cannot become zero", "Measurement", missing, lambda x: x.update(value=0))
    mutation("negative interval ordering", "Measurement", interval, lambda x: x.update(lower=30), "REVERSED_INTERVAL")
    mutation("unit must be explicit even if unknown", "Measurement", interval, lambda x: x.pop("unit"))
    mutation("PDF pages start at one", "SourceLocator", pdf_locator, lambda x: x.update(page=0))
    mutation("XML locator cannot borrow PDF fields", "SourceLocator", locator, lambda x: x.update(page=1))
    mutation("bbox ordering", "SourceLocator", pdf_locator, lambda x: x["bbox"].update(x0=0.95), "REVERSED_BBOX")
    mutation("stale extraction snapshot", "EvidenceBundle", bundle, lambda x: x["segments"][0]["locator"]["extraction_snapshot_ref"].update(version=2), "STALE_EXTRACTION_REFERENCE")
    mutation("stale document version", "EvidenceBundle", bundle, lambda x: x["segments"][0]["locator"].update(document_revision=2), "STALE_DOCUMENT_REVISION")
    mutation("unregistered document locator", "EvidenceBundle", bundle, lambda x: x["segments"][0]["locator"].update(document_id="missing"), "LOCATOR_DOCUMENT_NOT_FOUND")
    mutation("cross project document", "EvidenceBundle", bundle, lambda x: x["documents"][0].update(project_id="other-project"), "CROSS_PROJECT_DOCUMENT")
    mutation("duplicate document IDs", "EvidenceBundle", bundle, lambda x: x["documents"].append(copy.deepcopy(x["documents"][0])), "DUPLICATE_DOCUMENT_ID")
    mutation("coverage must describe returned data", "EvidenceBundle", bundle, lambda x: x["coverage"].update(llm_read_segment_count=1), "COVERAGE_MISMATCH")
    mutation("failed extraction has no complete snapshot", "EvidenceBundle", bundle, lambda x: x["documents"][1]["extraction"].update(status="complete"))
    mutation("LLM evidence needs a read trace", "EvidenceBundle", bundle, lambda x: x["evidence_records"][0]["producer"].update(kind="llm"), "LLM_SOURCE_NOT_MARKED_READ")
    mutation("unknown experiment independence cannot be counted", "EvidenceBundle", bundle, lambda x: x["coverage"].update(independent_experiment_count={"state": "known", "value": 3, "reason": None}), "UNKNOWN_EXPERIMENT_COUNT_INVENTED")
    mutation("broken comparison reference", "EvidenceBundle", bundle, lambda x: x["evidence_records"][0]["comparison_links"].append({"other_evidence_id": "missing", "role": "parent", "basis": "reported_comparison"}), "BROKEN_COMPARISON_LINK")
    mutation("unresolved structure cannot pretend to have a file", "EvidenceBundle", bundle, lambda x: x["compounds"][0].update(structure_ref=store.references["fixture-article"]))
    mutation("literature cannot invent mapped atom", "EvidenceBundle", bundle, lambda x: x["evidence_records"][0]["attachment"].update(atom_map_id=5))
    mutation("active policy cannot be written in evidence", "EvidenceBundle", bundle, lambda x: x.update(active_policy={"cutoff": 10}))
    mutation("unapproved provider is not silently added", "CollectionPolicy", store.payload("collection-policy"), lambda x: x.update(provider="any_web_url"))
    mutation("zero search budget is invalid", "CollectionPolicy", store.payload("collection-policy"), lambda x: x.update(max_documents=0))
    mutation("wrong payload version", "EvidenceBundle", bundle, lambda x: x.update(payload_version="0.2.0"))
    mutation("completed search cannot say not run", "EvidenceBundle", bundle, lambda x: x["searches"][0].update(termination="not_run"))
    mutation("manifest primary payload must be declared", "InputManifest", store.payload("literature-manifest"), lambda x: x.update(input_artifacts=[]))
    reject("duplicate JSON keys", lambda: parse_json(b'{"x":1,"x":2}'), "DUPLICATE_JSON_KEY")
    reject("NaN", lambda: parse_json(b'{"value":NaN}'), "NON_FINITE_JSON_NUMBER")
    reject("floating point overflow", lambda: parse_json(b'{"value":1e999}'), "NON_FINITE_JSON_NUMBER")
    reject("fixture graph forbidden in real run", lambda: reader.graph([store.references["fixture-evidence-bundle"]], "real"), "FIXTURE_IN_REAL_GRAPH")
    bad_hash = copy.deepcopy(store.references["fixture-article"])
    bad_hash["sha256"] = "0" * 64
    reject("artifact byte hash", lambda: reader.raw(bad_hash), "ARTIFACT_HASH_MISMATCH")
    job, response = exchanges["literature"]
    bad_job = copy.deepcopy(job)
    bad_job["input_digest"] = "0" * 64
    reject("manifest digest", lambda: validate_exchange(bad_job, reader), "INPUT_DIGEST_MISMATCH")
    bad_job = copy.deepcopy(job)
    bad_job["input_revision"] = 2
    reject("manifest revision", lambda: validate_exchange(bad_job, reader), "MANIFEST_MISMATCH:input_revision")
    bad_result = copy.deepcopy(response)
    bad_result["attempt"] = 2
    reject("late result attempt", lambda: validate_exchange(job, reader, bad_result), "RESULT_CORRELATION_MISMATCH:attempt")
    bad_result = copy.deepcopy(response)
    bad_result["producer"]["module"] = "science"
    reject("wrong module producer", lambda: validate_exchange(job, reader, bad_result), "RESULT_PRODUCER_MISMATCH")
    bad_result = copy.deepcopy(response)
    bad_result["output_artifacts"] = [store.references["fixture-article"]]
    reject("success needs typed output", lambda: validate_exchange(job, reader, bad_result), "PRIMARY_OUTPUT_REQUIRED")
    metadata_evidence = copy.deepcopy(bundle)
    metadata_evidence["documents"][0].update(format="metadata", access_level="metadata")
    metadata_evidence["documents"][0]["source_artifact_ref"]["media_type"] = "application/json"
    for item in walk(metadata_evidence):
        if "source_artifact_ref" in item and "extraction_snapshot_ref" in item:
            item["source_artifact_ref"]["media_type"] = "application/json"
            for field in ("block_id", "xpath", "table_id", "figure_id"):
                item.pop(field, None)
            item.update(kind="metadata", record_id="fixture-metadata", field="title")
    metadata_evidence["coverage"]["fulltext_acquired_count"] = 0
    reject("metadata cannot support experimental record", lambda: validate_payload("EvidenceBundle", metadata_evidence), "METADATA_IS_NOT_EXPERIMENTAL_EVIDENCE")
    changed_store = FixtureStore()
    changed_ref = copy.deepcopy(changed_store.references["fixture-target-resolution"])
    raw = (json.dumps(nonhuman, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    changed_store.data[key(changed_ref)] = raw
    changed_ref["sha256"] = hashlib.sha256(raw).hexdigest()
    target_job, target_result = exchanges["target"]
    changed_result = copy.deepcopy(target_result)
    changed_result["output_artifacts"] = [changed_ref]
    reject("resolved species must match requested species", lambda: validate_exchange(target_job, ArtifactReader(changed_store.read), changed_result), "RESOLVED_TAXON_MISMATCH")
    report = {"scope": "offline_schema_and_semantic_contract_tests_only", "schema_definitions": len(DOMAIN["$defs"]),
              "positive_checks": len(positive), "negative_cases": len(negative), "positive_names": positive, "negative_details": negative,
              "real_sources_or_models_executed": False, "dispatch_authorized": False}
    (ROOT / "verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(DOMAIN['$defs'])} schema definitions; {len(positive)} positive checks; {len(negative)} negative cases; 2 byte-hashed exchanges.")
    print("No real literature extraction, target lookup, molecular mapping, API, GPU, or authorization was exercised.")


if __name__ == "__main__":
    main()

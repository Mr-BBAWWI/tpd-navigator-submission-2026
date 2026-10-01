"""Source of the I1 JSON Schema. Generation only; no network or scientific policy."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCHEMA_ID = "urn:tpd-navigator:payloads:0.1.0"
BOUNDARY_ID = "urn:tpd-navigator:boundary:0.1.0"
TEXT = {"type": "string", "minLength": 1}
ID = {"type": "string", "pattern": "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$"}
INT = {"type": "integer", "minimum": 0}
POS = {"type": "integer", "minimum": 1}
NUM = {"type": "number"}
NULL = {"type": "null"}
BOOL = {"type": "boolean"}
D = {}


def ref(name):
    return {"$ref": "#/$defs/" + name}


def enum(*values):
    return {"enum": list(values)}


def const(value):
    return {"const": value}


def nullable(value):
    return {"anyOf": [value, NULL]}


def arr(value, minimum=0, unique=False):
    return {"type": "array", "items": value, "minItems": minimum, "uniqueItems": unique}


def obj(**props):
    return {"type": "object", "properties": props, "required": list(props), "additionalProperties": False}


def alternatives(*schemas):
    return {"oneOf": list(schemas)}


def when(field, values, **restrictions):
    return {"if": {"properties": {field: enum(*values)}, "required": [field]},
            "then": {"properties": restrictions}}


def payload(name, **props):
    return obj(payload_type=const(name), payload_version=const("0.1.0"), **props)


D["ArtifactRef"] = {"$ref": BOUNDARY_ID + "#/$defs/ArtifactRef"}
D["Issue"] = {"$ref": BOUNDARY_ID + "#/$defs/Issue"}
A = ref("ArtifactRef")
D["TextValue"] = alternatives(
    obj(state=const("known"), value=TEXT, reason=NULL),
    obj(state=enum("unknown", "not_provided", "not_reported", "not_applicable"), value=NULL, reason=TEXT),
)
T = ref("TextValue")
D["Count"] = alternatives(obj(state=const("known"), value=INT, reason=NULL),
                           obj(state=const("unknown"), value=NULL, reason=TEXT))
D["Identifier"] = obj(namespace=enum("doi", "pmid", "pmcid", "uniprot", "hgnc", "other"), value=TEXT)
D["TargetIdentity"] = obj(uniprot_accession=TEXT, taxon_id=POS, label=TEXT, gene_symbol=T,
                          isoform=T, mutation=T, construct=T)
D["ResearchContext"] = obj(objective=T, cell_line=T, tissue=T)
D["TargetQuery"] = payload("TargetQuery", query=TEXT, query_kind=enum("name", "gene_symbol", "uniprot_accession"),
                           taxon_id=POS, isoform=T, mutation=T, construct=T, research_context=ref("ResearchContext"))
D["TargetCandidate"] = obj(candidate_id=ID, identity=ref("TargetIdentity"), source_refs=arr(A, 1), match_reason=TEXT)
D["TargetSelection"] = alternatives(
    obj(kind=const("unique_match"), annotation_ref=NULL, reason=TEXT),
    obj(kind=const("human_selection"), annotation_ref=A, reason=TEXT),
)
D["TargetResolution"] = payload("TargetResolution", query_ref=A,
    resolution_status=enum("resolved", "ambiguous", "not_found"), candidates=arr(ref("TargetCandidate")),
    selected_candidate_id=nullable(ID), selection=nullable(ref("TargetSelection")),
    scope_status=enum("within_current_scope", "outside_current_scope", "undetermined"),
    scope_reasons=arr(TEXT), issues=arr(ref("Issue")))
D["TargetResolution"]["allOf"] = [
    when("resolution_status", ["resolved"], candidates={"minItems": 1}, selected_candidate_id=ID, selection=ref("TargetSelection")),
    when("resolution_status", ["ambiguous"], candidates={"minItems": 2}, selected_candidate_id=NULL, selection=NULL, scope_status=const("undetermined")),
    when("resolution_status", ["not_found"], candidates={"maxItems": 0}, selected_candidate_id=NULL, selection=NULL, scope_status=const("undetermined")),
    when("scope_status", ["outside_current_scope", "undetermined"], scope_reasons={"minItems": 1}),
]
D["SearchQuestion"] = obj(question_id=ID, question=TEXT, query=TEXT,
    purpose=enum("target_background", "warhead", "attachment_sar", "e3_context", "assay_context"))
D["CollectionPolicy"] = payload("CollectionPolicy", provider=const("europe_pmc"),
    max_records=POS, max_documents=POS, max_segments=POS, max_source_bytes=POS)
D["OperationParameters"] = payload("OperationParameters", profile_id=ID)
D["LiteratureCollectionRequest"] = payload("LiteratureCollectionRequest", target_resolution_ref=A,
    research_context=ref("ResearchContext"), questions=arr(ref("SearchQuestion"), 1),
    attached_source_refs=arr(A, unique=True), collection_policy_ref=A)

operations = json.loads((ROOT.parents[1] / "v0.1.0/operations.json").read_text(encoding="utf-8"))["operations"]
D["InputManifest"] = payload("InputManifest", boundary_version=const("0.1.0"), project_id=ID, run_id=ID,
    job_id=ID, input_revision=POS, operation_id=enum(*(o["operation_id"] for o in operations)),
    execution_mode=enum("development_live", "public_explore"), data_mode=enum("real", "test_fixture"),
    payload_ref=A, input_artifacts=arr(A, 1, True), parameters_ref=A, policy_ref=A)

D["DocumentRelation"] = obj(kind=enum("duplicate_of", "corrects", "retracts", "version_of"),
    related_identifier=ref("Identifier"), related_document_id=nullable(ID))
D["AcquisitionPath"] = obj(kind=enum("search", "upload", "database"), provider=TEXT,
    acquired_at={"type": "string", "format": "date-time"}, external_uri=nullable({"type": "string", "pattern": "^https?://[^\\s]+$"}))
D["Extraction"] = obj(status=enum("not_started", "complete", "partial", "failed", "not_applicable"),
    extractor_version=nullable(TEXT), snapshot_ref=nullable(A), issues=arr(ref("Issue")))
D["Extraction"]["allOf"] = [
    when("status", ["complete", "partial"], extractor_version=TEXT, snapshot_ref=A),
    when("status", ["partial", "failed"], issues={"minItems": 1}),
    when("status", ["not_started", "not_applicable"], extractor_version=NULL, snapshot_ref=NULL),
    when("status", ["failed"], extractor_version=TEXT, snapshot_ref=NULL),
]
D["ReplayRights"] = obj(status=enum("undetermined", "allowed", "restricted"), basis_ref=nullable(A), reason=TEXT)
D["ReplayRights"]["allOf"] = [when("status", ["allowed"], basis_ref=A)]
D["SourceDocument"] = obj(document_id=ID, revision=POS, project_id=ID, title=T,
    identifiers=arr(ref("Identifier")), acquisition_paths=arr(ref("AcquisitionPath"), 1),
    acquisition_status=enum("retrieved", "unavailable", "failed"),
    access_level=enum("fulltext", "abstract", "metadata", "unverified"),
    format=enum("jats_xml", "pdf", "abstract_text", "metadata"), source_artifact_ref=nullable(A),
    extraction=ref("Extraction"), relations=arr(ref("DocumentRelation")), replay_rights=ref("ReplayRights"), issues=arr(ref("Issue")))
D["SourceDocument"]["allOf"] = [
    when("acquisition_status", ["retrieved"], source_artifact_ref=A),
    when("acquisition_status", ["unavailable", "failed"], source_artifact_ref=NULL, issues={"minItems": 1},
         extraction={"properties": {"status": const("not_applicable")}}),
    when("format", ["jats_xml", "pdf"], access_level=enum("fulltext", "unverified")),
    when("format", ["abstract_text"], access_level=const("abstract")),
    when("format", ["metadata"], access_level=const("metadata")),
]
locator_common = dict(document_id=ID, document_revision=POS, source_artifact_ref=A,
                      extraction_snapshot_ref=A, extractor_version=TEXT)
D["BoundingBox"] = obj(x0={"type": "number", "minimum": 0, "maximum": 1},
    y0={"type": "number", "minimum": 0, "maximum": 1}, x1={"type": "number", "minimum": 0, "maximum": 1},
    y1={"type": "number", "minimum": 0, "maximum": 1}, units=const("normalized"), origin=const("top_left"))
D["SourceLocator"] = alternatives(
    obj(**locator_common, kind=const("jats_xml"), block_id=ID, xpath={"type": "string", "pattern": "^/"}, table_id=nullable(TEXT), figure_id=nullable(TEXT)),
    obj(**locator_common, kind=const("pdf"), page=POS, segment_id=ID, bbox=nullable(ref("BoundingBox")), table_id=nullable(TEXT), figure_id=nullable(TEXT)),
    obj(**locator_common, kind=const("abstract_text"), record_id=TEXT, field=const("abstract"), segment_id=ID),
    obj(**locator_common, kind=const("metadata"), record_id=TEXT, field=TEXT),
)
D["ReadEvent"] = obj(reader_kind=enum("code", "llm", "human"), reader_version=TEXT, activity_ref=A)
D["SourceSegment"] = obj(segment_id=ID, locator=ref("SourceLocator"), text_artifact_ref=A,
    content_kind=enum("text", "table", "table_caption", "figure_caption", "figure_reference", "abstract", "metadata"),
    read_events=arr(ref("ReadEvent")))
D["Measurement"] = alternatives(
    obj(kind=const("scalar"), value=NUM, relation=enum("=", "<", "<=", ">", ">=", "~"), unit=T, raw_text=TEXT),
    obj(kind=const("interval"), lower=NUM, upper=NUM, lower_inclusive=BOOL, upper_inclusive=BOOL, unit=T, raw_text=TEXT),
    obj(kind=const("qualitative"), value=TEXT, raw_text=TEXT),
    obj(kind=const("missing"), state=enum("unknown", "not_reported", "not_measured", "not_applicable"), value=NULL, reason=TEXT),
)
D["Observation"] = obj(endpoint=enum("Kd", "Ki", "IC50", "DC50", "Dmax", "binding", "degradation", "other"),
    endpoint_label=TEXT, measurement=ref("Measurement"), replicate_description=T, uncertainty_description=T)
D["AssayContext"] = obj(assay_id=T, assay_type=T, cell_line=T, tissue=T, duration=ref("Measurement"),
    conditions=arr(obj(name=TEXT, reported_value=TEXT)))
D["ReportedTarget"] = obj(reported_label=T, identity=nullable(ref("TargetIdentity")),
    resolution_state=enum("resolved", "unresolved", "not_reported"), reason=nullable(TEXT))
D["ReportedTarget"]["allOf"] = [
    when("resolution_state", ["resolved"], identity=ref("TargetIdentity"), reason=NULL),
    when("resolution_state", ["unresolved", "not_reported"], identity=NULL, reason=TEXT),
]
D["SourceCompound"] = obj(compound_id=ID, document_id=ID, source_label=TEXT,
    structure_state=enum("provided_unverified", "unavailable", "unknown"), structure_ref=nullable(A),
    structure_locators=arr(ref("SourceLocator")), reason=nullable(TEXT))
D["SourceCompound"]["allOf"] = [
    when("structure_state", ["provided_unverified"], structure_ref=A, structure_locators={"minItems": 1}, reason=NULL),
    when("structure_state", ["unavailable", "unknown"], structure_ref=NULL, reason=TEXT),
]
D["AttachmentDescription"] = obj(compound_id=ID, side=enum("warhead", "e3_ligand", "unknown"),
    source_position_label=T, position_description=T, source_locators=arr(ref("SourceLocator"), 1))
D["EvidenceProducer"] = obj(kind=enum("adapter", "llm", "human"), version=TEXT, activity_ref=A)
D["EvidenceReview"] = obj(state=enum("unreviewed", "code_checked", "human_reviewed", "rejected"),
    activity_ref=nullable(A), reason=nullable(TEXT))
D["EvidenceReview"]["allOf"] = [
    when("state", ["unreviewed"], activity_ref=NULL),
    when("state", ["code_checked", "human_reviewed", "rejected"], activity_ref=A, reason=TEXT),
]
D["ExperimentOrigin"] = alternatives(
    obj(state=const("resolved"), origin_id=ID, basis_locators=arr(ref("SourceLocator"), 1), reason=TEXT),
    obj(state=enum("unknown", "not_applicable"), origin_id=NULL, basis_locators=arr(ref("SourceLocator")), reason=TEXT),
)
D["EvidenceRecord"] = obj(evidence_id=ID,
    basis_kind=enum("experiment", "computation", "author_interpretation", "researcher_hypothesis"),
    statement=TEXT, source_locators=arr(ref("SourceLocator"), 1), producer=ref("EvidenceProducer"), review=ref("EvidenceReview"),
    target=ref("ReportedTarget"), compound_ids=arr(ID, unique=True), assay_context=ref("AssayContext"),
    observations=arr(ref("Observation")), attachment=nullable(ref("AttachmentDescription")),
    comparison_links=arr(obj(other_evidence_id=ID, role=enum("parent", "modified"), basis=enum("reported_comparison", "proposed_comparison"))),
    experiment_origin=ref("ExperimentOrigin"), limitations=arr(TEXT))
D["SearchExecution"] = obj(question_id=ID, query=TEXT, provider=const("europe_pmc"),
    status=enum("completed", "failed", "not_run"), total_hits=ref("Count"),
    search_result_ref=nullable(A), retrieved_document_ids=arr(ID, unique=True),
    termination=enum("source_exhausted", "record_limit", "document_limit", "segment_limit", "byte_limit", "time_limit", "network_error", "provider_error", "not_run"),
    issues=arr(ref("Issue")))
D["SearchExecution"]["allOf"] = [
    when("status", ["completed"], search_result_ref=A, termination=enum("source_exhausted", "record_limit", "document_limit", "segment_limit", "byte_limit", "time_limit")),
    when("status", ["failed"], issues={"minItems": 1}, termination=enum("network_error", "provider_error", "time_limit")),
    when("status", ["not_run"], search_result_ref=NULL, termination=const("not_run"), retrieved_document_ids={"maxItems": 0}, total_hits={"properties": {"state": const("unknown")}}, issues={"minItems": 1}),
]
D["Coverage"] = obj(stored_document_count=INT, fulltext_acquired_count=INT, extracted_segment_count=INT,
    llm_read_segment_count=INT, evidence_record_count=INT, reviewed_record_count=INT,
    unavailable_document_count=INT, independent_experiment_count=ref("Count"))
D["EvidenceBundle"] = payload("EvidenceBundle", project_id=ID, request_ref=A, target_resolution_ref=A,
    documents=arr(ref("SourceDocument")), segments=arr(ref("SourceSegment")), compounds=arr(ref("SourceCompound")),
    evidence_records=arr(ref("EvidenceRecord")), searches=arr(ref("SearchExecution"), 1),
    coverage=ref("Coverage"), issues=arr(ref("Issue")))

ENTRYPOINTS = ["InputManifest", "TargetQuery", "TargetResolution", "LiteratureCollectionRequest", "EvidenceBundle", "CollectionPolicy", "OperationParameters"]
SCHEMA = {"$schema": "https://json-schema.org/draft/2020-12/schema", "$id": SCHEMA_ID,
          "title": "TPD Navigator I1 payloads v0.1.0", "oneOf": [ref(n) for n in ENTRYPOINTS], "$defs": D}


def main():
    (ROOT / "domain.schema.json").write_text(json.dumps(SCHEMA, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Generated {len(D)} definitions, {len(ENTRYPOINTS)} payload entrypoints.")


if __name__ == "__main__":
    main()

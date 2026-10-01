"""Create synthetic, byte-hashed contract examples. Never fetch real literature."""
import hashlib
import json
from pathlib import Path

from build_schemas import SCHEMA_ID

ROOT = Path(__file__).resolve().parent / "examples"
INDEX = []


def known(value):
    return {"state": "known", "value": value, "reason": None}


def unknown(reason="Not reported in this synthetic example.", state="not_reported"):
    return {"state": state, "value": None, "reason": reason}


def missing():
    return {"kind": "missing", "state": "not_reported", "value": None, "reason": "Synthetic example has no value."}


def issue(code):
    return {"code": code, "category": "invalid_output", "message": "Synthetic extraction failure; no real PDF was processed.",
            "affected_artifact_ids": ["fixture-upload"], "retry_hint": "requires_changed_input"}


def add(name, value, kind=None, media="application/json"):
    if kind:
        value = {"payload_type": kind, "payload_version": "0.1.0", **value}
    raw = value if isinstance(value, bytes) else (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    suffix = {"application/json": ".json", "application/xml": ".xml", "text/plain": ".txt", "application/pdf": ".pdf"}[media]
    path = Path("artifacts") / (name + suffix)
    (ROOT / path).parent.mkdir(parents=True, exist_ok=True)
    (ROOT / path).write_bytes(raw)
    reference = {"artifact_id": "fixture-" + name, "version": 1, "sha256": hashlib.sha256(raw).hexdigest(),
                 "media_type": media, "schema_id": SCHEMA_ID + "#/$defs/" + kind if kind else "urn:tpd-navigator:raw:test-fixture",
                 "provenance": "test_fixture"}
    INDEX.append({"reference": reference, "path": path.as_posix()})
    return reference


def write(name, value):
    (ROOT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def job(name, operation, request_ref, inputs, policy, parameters):
    shared = {"project_id": "fixture-project", "run_id": "fixture-run", "job_id": "fixture-job-" + name,
              "input_revision": 1, "operation_id": operation, "execution_mode": "development_live", "data_mode": "test_fixture",
              "input_artifacts": inputs, "parameters_ref": parameters, "policy_ref": policy, "approval_ref": None}
    manifest_shared = {k: v for k, v in shared.items() if k != "approval_ref"}
    manifest = add(name + "-manifest", {"boundary_version": "0.1.0", "payload_ref": request_ref, **manifest_shared}, "InputManifest")
    request = {"contract_version": "0.1.0", **shared, "attempt": 1, "input_manifest": manifest,
               "input_digest": manifest["sha256"], "reservation_id": None,
               "limits": {"wall_time_seconds": 30, "max_llm_requests": 0, "max_llm_tokens": 0, "max_gpu_seconds": 0, "max_output_bytes": 1048576}}
    write(name + ".job.json", request)
    return request


def result(name, request, output, status, issues):
    keys = ["contract_version", "project_id", "run_id", "job_id", "attempt", "input_revision", "input_digest", "operation_id", "data_mode"]
    response = {**{k: request[k] for k in keys}, "status": status, "output_artifacts": [output], "issues": issues,
                "actual_usage": {"wall_time_ms": None, "llm_requests": None, "llm_tokens": None, "gpu_seconds": None, "usage_complete": False},
                "producer": {"module": "science" if name == "target" else "literature", "code_version": "synthetic-fixture-only", "tool_versions_ref": None}}
    write(name + ".result.json", response)


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    parameters = add("parameters", {"profile_id": "fixture-only-no-runtime-profile"}, "OperationParameters")
    target_policy = add("target-policy", {"test_fixture": True, "description": "Scope configuration placeholder, not an approved operational policy."})
    context = {"objective": unknown(state="not_provided"), "cell_line": unknown(state="not_provided"), "tissue": unknown(state="not_provided")}
    query = add("target-query", {"query": "FIXTURE_ALPHA (not a real biological target)", "query_kind": "name", "taxon_id": 9606,
                "isoform": unknown(), "mutation": unknown(), "construct": unknown(), "research_context": context}, "TargetQuery")
    identity = {"uniprot_accession": "FIXTURE_ALPHA", "taxon_id": 9606, "label": "Synthetic target Alpha",
                "gene_symbol": known("FIXTURE_ALPHA"), "isoform": unknown(), "mutation": unknown(), "construct": unknown()}
    target_source = add("target-source", {"test_fixture": True, "identity": identity})
    resolution = add("target-resolution", {"query_ref": query, "resolution_status": "resolved",
        "candidates": [{"candidate_id": "candidate-alpha", "identity": identity, "source_refs": [target_source], "match_reason": "Synthetic unique match."}],
        "selected_candidate_id": "candidate-alpha", "selection": {"kind": "unique_match", "annotation_ref": None, "reason": "One synthetic candidate."},
        "scope_status": "within_current_scope", "scope_reasons": [], "issues": []}, "TargetResolution")
    target_job = job("target", "resolve_target", query, [query], target_policy, parameters)
    result("target", target_job, resolution, "succeeded", [])

    policy = add("collection-policy", {"provider": "europe_pmc", "max_records": 100, "max_documents": 3,
                 "max_segments": 20, "max_source_bytes": 10000}, "CollectionPolicy")
    upload = add("upload", b"%PDF-1.4\n% INTENTIONALLY TRUNCATED SYNTHETIC PDF for the extraction-failed contract example.\n", media="application/pdf")
    question = {"question_id": "q-attachment", "question": "What is actually known about attachment?", "query": "FIXTURE_ALPHA attachment", "purpose": "attachment_sar"}
    collection = add("collection-request", {"target_resolution_ref": resolution, "research_context": context,
        "questions": [question], "attached_source_refs": [upload], "collection_policy_ref": policy}, "LiteratureCollectionRequest")
    collection_job = job("literature", "collect_literature", collection, [collection, resolution, upload], policy, parameters)
    source = add("article", b'<article><body><p id="b001">SYNTHETIC: compound C1 has IC50 &gt;1000 nM. R1 is mentioned; its chemical structure is unavailable.</p></body></article>\n', media="application/xml")
    segment_text = add("segment-text", b"SYNTHETIC: compound C1 has IC50 >1000 nM. R1 is mentioned; its chemical structure is unavailable.\n", media="text/plain")
    extraction = add("extraction", {"test_fixture": True, "blocks": [{"block_id": "b001", "text_ref": segment_text}]})
    activity = add("extraction-activity", {"test_fixture": True, "description": "Fabricated adapter trace; no real extraction was run."})
    search = add("search-response", {"test_fixture": True, "hits": 1, "identifiers": []})
    rights = {"status": "undetermined", "basis_ref": None, "reason": "Fixture only; not a public replay source."}
    article_doc = {"document_id": "document-article", "revision": 1, "project_id": "fixture-project", "title": known("Synthetic attachment example"),
        "identifiers": [], "acquisition_paths": [{"kind": "search", "provider": "europe_pmc", "acquired_at": "2026-09-11T00:00:00Z", "external_uri": None}],
        "acquisition_status": "retrieved", "access_level": "fulltext", "format": "jats_xml", "source_artifact_ref": source,
        "extraction": {"status": "complete", "extractor_version": "fixture-extractor-1", "snapshot_ref": extraction, "issues": []},
        "relations": [], "replay_rights": rights, "issues": []}
    failed = issue("PDF_EXTRACTION_FAILED")
    upload_doc = {"document_id": "document-upload", "revision": 1, "project_id": "fixture-project", "title": known("Intentionally corrupt PDF fixture"),
        "identifiers": [], "acquisition_paths": [{"kind": "upload", "provider": "researcher", "acquired_at": "2026-09-11T00:00:00Z", "external_uri": None}],
        "acquisition_status": "retrieved", "access_level": "unverified", "format": "pdf", "source_artifact_ref": upload,
        "extraction": {"status": "failed", "extractor_version": "fixture-extractor-1", "snapshot_ref": None, "issues": [failed]},
        "relations": [], "replay_rights": rights, "issues": []}
    locator = {"document_id": "document-article", "document_revision": 1, "source_artifact_ref": source,
        "extraction_snapshot_ref": extraction, "extractor_version": "fixture-extractor-1", "kind": "jats_xml",
        "block_id": "b001", "xpath": "/article/body/p[@id='b001']", "table_id": None, "figure_id": None}
    evidence = {"evidence_id": "evidence-c1", "basis_kind": "experiment", "statement": "Synthetic reported IC50; no mapped atom or degradation result.",
        "source_locators": [locator], "producer": {"kind": "adapter", "version": "fixture-extractor-1", "activity_ref": activity},
        "review": {"state": "unreviewed", "activity_ref": None, "reason": None},
        "target": {"reported_label": known("Synthetic target Alpha"), "identity": identity, "resolution_state": "resolved", "reason": None},
        "compound_ids": ["compound-c1"], "assay_context": {"assay_id": unknown(), "assay_type": unknown(), "cell_line": unknown(), "tissue": unknown(), "duration": missing(), "conditions": []},
        "observations": [{"endpoint": "IC50", "endpoint_label": "IC50", "measurement": {"kind": "scalar", "value": 1000, "relation": ">", "unit": known("nM"), "raw_text": ">1000 nM"}, "replicate_description": unknown(), "uncertainty_description": unknown()}],
        "attachment": {"compound_id": "compound-c1", "side": "warhead", "source_position_label": known("R1"), "position_description": unknown(), "source_locators": [locator]},
        "comparison_links": [], "experiment_origin": {"state": "unknown", "origin_id": None, "basis_locators": [locator], "reason": "Not independently verified in this fixture."},
        "limitations": ["No chemical structure; no atom mapping; no degradation observation; entirely synthetic."]}
    bundle = add("evidence-bundle", {"project_id": "fixture-project", "request_ref": collection, "target_resolution_ref": resolution,
        "documents": [article_doc, upload_doc], "segments": [{"segment_id": "segment-b001", "locator": locator, "text_artifact_ref": segment_text, "content_kind": "text", "read_events": []}],
        "compounds": [{"compound_id": "compound-c1", "document_id": "document-article", "source_label": "C1", "structure_state": "unavailable", "structure_ref": None, "structure_locators": [], "reason": "Not available in synthetic source."}],
        "evidence_records": [evidence], "searches": [{"question_id": question["question_id"], "query": question["query"], "provider": "europe_pmc", "status": "completed",
            "total_hits": {"state": "known", "value": 1, "reason": None}, "search_result_ref": search, "retrieved_document_ids": ["document-article"], "termination": "source_exhausted", "issues": []}],
        "coverage": {"stored_document_count": 2, "fulltext_acquired_count": 1, "extracted_segment_count": 1, "llm_read_segment_count": 0,
            "evidence_record_count": 1, "reviewed_record_count": 0, "unavailable_document_count": 0,
            "independent_experiment_count": {"state": "unknown", "value": None, "reason": "Experiment origin unresolved."}}, "issues": [failed]}, "EvidenceBundle")
    result("literature", collection_job, bundle, "partial", [failed])
    write("artifact-index.json", {"data_mode": "test_fixture", "artifacts": INDEX})
    print(f"Generated {len(INDEX)} byte-hashed synthetic artifacts and 2 job/result pairs.")


if __name__ == "__main__":
    main()

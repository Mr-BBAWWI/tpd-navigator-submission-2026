"""Collect one DOI for a resolved target through real official APIs, optionally analyze figures.

This is an internal integration entry point; it never creates candidate structures or approvals.
"""
import argparse
import json
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from packages.contracts import payload, missing
from packages.platform.store import Store
from packages.platform.orchestrator import Orchestrator
from packages.platform.analysis import AnalysisService
from packages.science.target import resolve_target
from packages.literature.collector import Collector
from packages.agents.provider import DaconProvider, load_key
from packages.transport import PublicTransport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--doi", required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--analyze-figure", action="append", default=[])
    args = parser.parse_args()
    if args.analyze_figure:
        load_key()
    store = Store(args.data_dir)
    port = store.scope("local-research")
    transport = PublicTransport()
    orchestrator = Orchestrator(store, port.project, lambda: transport)
    run, _ = store.create_run(port.project, uuid.uuid4().hex, {"query": args.uniprot, "query_kind": "uniprot_accession",
        "objective": "선정 문헌에서 원문 화합물·관측·부착점 근거 확인", "cell_line": "", "pdf_artifact_ids": [], "max_documents": 2})
    context = {"objective": missing("아래 선정 문헌 근거 검토", "not_provided"), "cell_line": missing(), "tissue": missing()}
    query_ref = port.put_json(payload("TargetQuery", query=args.uniprot, query_kind="uniprot_accession", taxon_id=9606,
        isoform=missing(), mutation=missing(), construct=missing(), research_context=context), "TargetQuery", "human_annotation")
    policy_ref = port.put_json({"profile": "local-human-target-v1", "taxon_id": 9606, "exclude_annotated_transmembrane": True})
    target_ref, target = orchestrator.record_job(run, "resolve_target", query_ref, [query_ref], policy_ref,
                                               lambda: resolve_target(query_ref, port, transport))
    if target["resolution_status"] != "resolved" or target["scope_status"] == "outside_current_scope":
        raise SystemExit("Target requires selection or is outside current scope; no analysis performed")
    policy_ref = port.put_json(payload("CollectionPolicy", provider="europe_pmc", max_records=2, max_documents=2,
                                      max_segments=180, max_source_bytes=12_000_000), "CollectionPolicy")
    request_ref = port.put_json(payload("LiteratureCollectionRequest", target_resolution_ref=target_ref, research_context=context,
        questions=[{"question_id": "q-selected-doi", "question": "선정 문헌의 관측과 부착점 근거는 무엇인가?",
                    "query": 'DOI:"' + args.doi.replace('"', '') + '"', "purpose": "attachment_sar"}],
        attached_source_refs=[], collection_policy_ref=policy_ref), "LiteratureCollectionRequest")
    bundle_ref, bundle = orchestrator.record_job(run, "collect_literature", request_ref, [request_ref, target_ref], policy_ref,
                                               lambda: Collector(port, transport).collect(request_ref))
    store.update(port.project, run["id"], {"state": "partial" if bundle["issues"] else "completed", "bundle_ref": bundle_ref,
        "target_ref": target_ref, "target_query_ref": query_ref, "issues": bundle["issues"], "message": "선정 DOI 문헌 수집 완료. AI 분석은 별도 상태에서 확인합니다."})
    print(json.dumps({"run_id": run["id"], "documents": len(bundle["documents"]), "segments": len(bundle["segments"])}), flush=True)
    if args.analyze_figure:
        selected = [s["segment_id"] for s in bundle["segments"] if s["locator"].get("figure_id") in args.analyze_figure]
        service = AnalysisService(store, port.project, lambda: DaconProvider(load_key()))
        job, _ = service.create(run["id"], 1, selected, uuid.uuid4().hex)
        service.execute(job["id"])
        result = service.get(job["id"])
        print(json.dumps({k: result.get(k) for k in ("id", "state", "calls", "total_tokens", "error_code")}), flush=True)
        return 0 if result["state"] == "review_ready" else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

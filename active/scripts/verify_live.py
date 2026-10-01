"""Explicit real public-source smoke check. No fixtures, LLM, or GPU."""
import argparse
import json
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lxml import etree
from packages.contracts import ArtifactReader, validate_exchange, now, encoded
from packages.platform.store import Store
from packages.platform.orchestrator import Orchestrator
from packages.transport import PublicTransport


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accession", nargs="+", help="Explicit verification cases, not product defaults")
    parser.add_argument("--data-dir", default=str(ROOT / ".localdata"))
    parser.add_argument("--output", default=str(ROOT / "verification/i2_live.json"))
    args = parser.parse_args()
    store = Store(args.data_dir)
    project = "local-research"
    results = []
    for accession in args.accession:
        transport = PublicTransport()
        runner = Orchestrator(store, project, lambda: transport)
        request = {"query": accession, "query_kind": "uniprot_accession", "objective": "", "cell_line": "", "pdf_artifact_ids": [], "max_documents": 3}
        run, _ = store.create_run(project, "verify-" + uuid.uuid4().hex, request)
        runner.execute(run["id"])
        value = runner.view(run["id"])
        port = store.scope(project)
        for job in value["jobs"]:
            validate_exchange(port.json(job["job_ref"]), ArtifactReader(port.read), port.json(job["result_ref"]))
        roots, locations = {}, 0
        bundle = value["bundle"]
        if bundle:
            for segment in bundle["segments"]:
                locator = segment["locator"]
                if locator["kind"] == "jats_xml":
                    source = locator["source_artifact_ref"]
                    if source["artifact_id"] not in roots:
                        roots[source["artifact_id"]] = etree.fromstring(port.read(source), etree.XMLParser(resolve_entities=False, no_network=True))
                    node = roots[source["artifact_id"]].xpath(locator["xpath"])
                    assert len(node) == 1
                    expected = " ".join(" ".join(node[0].itertext()).split())[:12000]
                    assert port.read(segment["text_artifact_ref"]).decode() == expected
                    locations += 1
        before = transport.calls
        offline = Orchestrator(Store(args.data_dir), project, lambda: (_ for _ in ()).throw(AssertionError("Unexpected network")))
        assert offline.view(run["id"])["bundle"] == bundle
        assert transport.calls == before
        result = {"accession": accession, "run_id": run["id"], "state": value["state"], "http_calls": transport.calls,
                  "received_bytes": transport.bytes_received, "jobs_verified": len(value["jobs"]), "xml_locations_verified": locations,
                  "saved_view_without_network": True, "target_ref": value["target_ref"], "bundle_ref": value["bundle_ref"],
                  "coverage": bundle["coverage"] if bundle else None,
                  "searches": [{"status": s["status"], "hits": s["total_hits"], "termination": s["termination"]} for s in bundle["searches"]] if bundle else [],
                  "issues": sorted(set(i["code"] for i in value["issues"]))}
        results.append(result)
        print(json.dumps(result, ensure_ascii=True), flush=True)
    report = {"verified_at": now(), "data_mode": "real", "scope": "public source integration; not scientific efficacy or generalization validation", "results": results}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(encoded(report))
    assert all(r["state"] in {"completed", "partial"} and r["coverage"]["stored_document_count"] > 0 for r in results), "Live integration incomplete; inspect report"


if __name__ == "__main__":
    main()

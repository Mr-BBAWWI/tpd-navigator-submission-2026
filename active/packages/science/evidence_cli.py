"""Offline candidate evidence export and saved Boltz run collection; never inference."""
import argparse
from pathlib import Path

from .boltz_collection import build_collection, collect_run, empty_collection, load_archive, summarize_result
from .candidate_evidence import build_evidence, markdown, validate_evidence
from .evidence_common import file_index, write_new_directory, safe_relative, schema_check, verify_seal
from .handoff import Snapshot, artifact_ref, build_material, check, encoded, parse, sha


def export_evidence(snapshot, destination, *, receipt=None, runs=(), exchange=None, allow_synthetic=False):
    check(not destination.exists(), "OUTPUT_EXISTS")
    material = build_material(snapshot)
    if receipt is not None:
        collection, files = build_collection(snapshot, receipt, runs, allow_synthetic=allow_synthetic)
    else:
        check(not runs, "RECEIPT_REQUIRED_FOR_RUNS")
        collection, files = empty_collection(material), {}
    evidence = build_evidence(material, collection, exchange, allow_synthetic=allow_synthetic)
    validate_evidence(evidence, material, collection, exchange, allow_synthetic=allow_synthetic)
    files.update({"material.json": encoded(material), "collection.json": encoded(collection),
                  "evidence.json": encoded(evidence), "README.md": markdown(evidence),
                  "cpu/handoff.json": snapshot.handoff_bytes})
    files.update({"cpu/" + name: data for name, data in snapshot.files.items()})
    if receipt is not None:
        files["receipt.json"] = encoded(receipt)
    if exchange is not None:
        files["literature-exchange.json"] = encoded(exchange)
    files["manifest.json"] = encoded({"format": "tpd-evidence-export/0.1.0-draft", "files": file_index(files),
        "evidence_ref": artifact_ref(files["evidence.json"], schema="urn:tpd-navigator:b-candidate-evidence:0.1.0-draft"),
        "registered_in_M2": False, "public_release_ready": False})
    write_new_directory(destination, files)
    return evidence


def read_export(directory, *, allow_synthetic=False):
    """Verify the packaged source files and exact derived views without RDKit/Boltz."""
    manifest = parse((directory / "manifest.json").read_bytes())
    check(manifest["format"] == "tpd-evidence-export/0.1.0-draft", "EXPORT_FORMAT")
    check(manifest["registered_in_M2"] is False and manifest["public_release_ready"] is False, "EXPORT_AUTHORITY")
    files = {}
    for name, entry in manifest["files"].items():
        data = safe_relative(directory, name).read_bytes()
        check(sha(data) == entry["sha256"] and len(data) == entry["bytes"], "EXPORT_FILE_HASH")
        files[name] = data
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()}
    check(actual == set(files) | {"manifest.json"}, "EXPORT_UNINDEXED_FILES")
    check(manifest["evidence_ref"] == artifact_ref(files["evidence.json"], schema="urn:tpd-navigator:b-candidate-evidence:0.1.0-draft"), "EXPORT_EVIDENCE_REF")
    material, collection, evidence = [parse(files[n + ".json"]) for n in ("material", "collection", "evidence")]
    exchange = parse(files["literature-exchange.json"]) if "literature-exchange.json" in files else None
    validate_evidence(evidence, material, collection, exchange, allow_synthetic=allow_synthetic)
    receipt = parse(files["receipt.json"]) if "receipt.json" in files else None
    check((receipt is None) == (collection["receipt_digest"] is None), "EXPORT_RECEIPT_REQUIRED")
    if receipt is not None:
        schema_check(receipt, "b_boltz_parser.schema.json")
        verify_seal(receipt)
        check(receipt["digest"] == collection["receipt_digest"] and receipt["material_digest"] == material["material_digest"],
              "EXPORT_RECEIPT_BINDING")
    # Every displayed sample retains a hash-bound full comparison and raw run archive.
    for run in collection["runs"]:
        archive_path = safe_relative(directory, run["archive_file"]["path"])
        _, info, _, _ = load_archive(archive_path.parent, receipt)
        check(info == run["reported_execution"], "EXPORT_RUN_PROJECTION")
        references = [run["archive_file"]] + [m["comparison_file"] for m in run["models"] if m["comparison_file"]]
        for ref in references:
            check(ref["path"] in files and sha(files[ref["path"]]) == ref["sha256"], "EXPORT_RESULT_REF")
        for model in run["models"]:
            if model["status"] == "parsed":
                result = parse(files[model["comparison_file"]["path"]])
                source = result["source_manifest"]
                check(source["compound_id"] == run["compound_id"] and source["model_rank"] == model["rank"]
                      and source["runtime_reported"]["run_id"] == run["run_id"], "EXPORT_SAMPLE_BINDING")
                check(summarize_result(result) == model["summary"], "EXPORT_SAMPLE_PROJECTION")
    check(files["README.md"] == markdown(evidence), "EXPORT_REPORT_PROJECTION")
    return evidence


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("bundle", "collect"))
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--run", type=Path, action="append", default=[])
    parser.add_argument("--literature-exchange", type=Path)
    parser.add_argument("--predictions", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--run-info", type=Path)
    args = parser.parse_args(argv)
    try:
        snapshot = Snapshot(args.snapshot, args.raw)
        receipt = parse(args.receipt.read_bytes()) if args.receipt else None
        if args.action == "collect":
            check(receipt is not None and args.predictions is not None and args.input is not None and args.run_info is not None,
                  "COLLECT_RECEIPT_PREDICTIONS_INPUT_RUN_INFO_REQUIRED")
            from .boltz_parser import validate_receipt
            validate_receipt(receipt, snapshot)
            info = parse(args.run_info.read_bytes())
            check(info.get("origin") != "synthetic_test", "SYNTHETIC_RUN_NOT_ALLOWED")
            result = collect_run(args.predictions, args.input, receipt, info, args.out)
            print("Saved run archive: " + result["run_id"] + "; reported execution only, not an approval.")
        else:
            exchange = parse(args.literature_exchange.read_bytes()) if args.literature_exchange else None
            evidence = export_evidence(snapshot, args.out, receipt=receipt, runs=args.run, exchange=exchange)
            print("Candidate evidence exported: " + evidence["digest"] + "; human review pending.")
        return 0
    except (ValueError, OSError, ImportError, KeyError, TypeError, RuntimeError) as error:
        print(f"EVIDENCE_FAILED: {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Prepare official parser receipts or inspect saved results, without model inference."""
import argparse
import json
from pathlib import Path

from .boltz_parser import prepare_receipt
from .boltz_results import read_result
from .handoff import Snapshot, encoded, parse, local_file, sha


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "compare"))
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--mols", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.out.exists():
            raise ValueError("OUTPUT_ALREADY_EXISTS")
        snapshot = Snapshot(args.snapshot, args.raw)
        if args.action == "prepare":
            if args.mols is None:
                raise ValueError("MOLS_REQUIRED")
            receipt, _ = prepare_receipt(snapshot, args.mols)
            files = {"receipt.json": encoded(receipt)}
            summary = {"parser_validation": "passed", "prediction_status": "not_run"}
        else:
            if args.receipt is None or args.result is None:
                raise ValueError("RECEIPT_AND_RESULT_REQUIRED")
            receipt_bytes = args.receipt.read_bytes()
            receipt = parse(receipt_bytes)
            result = read_result(args.result, receipt, snapshot)
            files = {"comparison.json": encoded(result), "receipt.json": receipt_bytes,
                     "raw/manifest.json": local_file(args.result, "manifest.json").read_bytes()}
            for info in result["source_manifest"]["files"].values():
                if info is not None:
                    data = local_file(args.result, info["name"]).read_bytes()
                    if sha(data) != info["sha256"]:
                        raise ValueError("RESULT_CHANGED_DURING_READ")
                    files["raw/" + info["name"]] = data
            summary = {"compound_id": result["compound_id"], "origin": result["origin"],
                       "human_review": "pending", "efficacy_claim": "not_established"}
        files["files.json"] = encoded({"files": {n: sha(data) for n, data in files.items()}, "summary": summary})
        args.out.mkdir(parents=True, exist_ok=False)
        for name, data in files.items():
            path = args.out / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as stream:
                stream.write(data)
        print(json.dumps(summary))
        return 0
    except (ValueError, OSError, ImportError, KeyError, StopIteration) as error:
        print(f"BOLTZ_IO_FAILED: {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Build a source-hash-verified known-case benchmark distribution."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.benchmark_distribution import build_benchmark_distribution


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--base", required=True, type=Path)
    result.add_argument("--receipt", required=True, action="append", type=Path)
    result.add_argument("--output", required=True, type=Path)
    return result


def _load_receipt(path: Path) -> dict:
    if not path.is_file() or path.is_symlink():
        raise ValueError("RECEIPT_REGULAR_FILE_REQUIRED")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("RECEIPT_OBJECT_REQUIRED")
    value = dict(value)
    value["provenance"] = {"path": str(path.resolve()), "sha256": _sha256(path)}
    return value


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if not args.base.is_file() or args.base.is_symlink():
        raise ValueError("BASE_REGULAR_FILE_REQUIRED")
    receipts = [_load_receipt(path) for path in args.receipt]
    result = build_benchmark_distribution(Path(args.base), receipts)
    result["sources"] = {
        "base": {"path": str(args.base.resolve()), "sha256": _sha256(args.base)},
        "receipts": [item["provenance"] for item in receipts],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

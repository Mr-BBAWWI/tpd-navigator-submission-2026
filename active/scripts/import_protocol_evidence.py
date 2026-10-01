#!/usr/bin/env python3
"""Explicit local operator import and read CLI for protocol diagnostics."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.contracts import ContractError
from packages.platform.protocol_evidence import ProtocolEvidenceService
from packages.platform.store import Store


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--store", type=Path, required=True)
    result.add_argument("--project", required=True)
    commands = result.add_subparsers(dest="command", required=True)
    register = commands.add_parser("import")
    register.add_argument("--job-id", required=True)
    register.add_argument("--pack", type=Path, required=True)
    register.add_argument("--expected-manifest-sha256", required=True)
    register.add_argument("--assessment-id")
    listing = commands.add_parser("list")
    listing.add_argument("--job-id", required=True)
    view = commands.add_parser("view")
    view.add_argument("--id", required=True)
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        root = args.store.resolve()
        if not (root / "index.sqlite3").is_file():
            raise ContractError("PROTOCOL_EVIDENCE_EXISTING_STORE_REQUIRED")
        service = ProtocolEvidenceService(Store(root), args.project)
        if args.command == "import":
            value = service.register_pack(
                args.job_id, args.pack,
                expected_manifest_sha256=args.expected_manifest_sha256,
                expected_assessment_id=args.assessment_id,
            )
        elif args.command == "list":
            value = service.list(args.job_id)
        else:
            value = service.view(args.id)
        print(json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False))
        return 0
    except (ContractError, KeyError, ValueError, OSError) as error:
        print(json.dumps({"error": str(error)}, ensure_ascii=False, sort_keys=True,
                         allow_nan=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

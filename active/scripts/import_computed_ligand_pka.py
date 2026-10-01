#!/usr/bin/env python3
"""Import a local computed ligand-pKa evidence pack."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.contracts import ContractError  # noqa: E402
from packages.platform.computed_ligand_pka import import_computed_ligand_pka  # noqa: E402
from packages.platform.store import Store  # noqa: E402


def _is_link(path: Path) -> bool:
    junction = getattr(path, "is_junction", None)
    return path.is_symlink() or bool(junction and junction())


def _linked_chain(path: Path) -> bool:
    current = path
    while True:
        if _is_link(current):
            return True
        if current.parent == current:
            return False
        current = current.parent


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--pack", type=Path, required=True)
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.receipt is not None:
            if not args.receipt.parent.is_dir() \
                    or os.path.lexists(args.receipt) \
                    or _linked_chain(args.receipt.parent):
                raise OSError("receipt path must be new with an existing parent")
        if not (args.store / "index.sqlite3").is_file():
            raise OSError("store/index.sqlite3 does not exist")
        store = Store(args.store)
        result = import_computed_ligand_pka(
            store, args.project, args.job_id, args.pack,
            expected_manifest_sha256=args.expected_manifest_sha256)
        text = json.dumps(result, sort_keys=True, ensure_ascii=False) + "\n"
        if args.receipt is not None:
            with args.receipt.open("x", encoding="utf-8") as handle:
                handle.write(text)
        sys.stdout.write(text)
        return 0
    except (ContractError, KeyError, OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

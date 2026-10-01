#!/usr/bin/env python3
"""Run the predeclared native 9D12 parent expansion v7."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.native_parent_expansion import run_expansion


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the real native Vina 20-case expansion, or an explicit no-docking audit."
    )
    parser.add_argument("--native-probe", required=True, type=Path)
    parser.add_argument("--native-cif", required=True, type=Path)
    parser.add_argument("--source-export", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="Run rule applicability and panel graph audit only; NO DOCKING is executed.",
    )
    args = parser.parse_args()
    result = run_expansion(
        native_probe=args.native_probe,
        native_cif=args.native_cif,
        source_export=args.source_export,
        output=args.output,
        audit_only=args.audit_only,
    )
    print(result["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

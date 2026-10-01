#!/usr/bin/env python3
"""Build a verified, read-only nine-parent funnel dossier."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.parent_funnel import build_parent_funnel


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a source-bound nine-parent computed diagnostic; no selection or approval is performed."
    )
    parser.add_argument("--export-root", required=True, help="Directory containing exactly nine immediate parent export directories")
    parser.add_argument("--expert-source", required=True, help="Expert reply .docx to fingerprint and parse")
    parser.add_argument("--output", required=True, help="Fresh output directory; it must not already exist")
    parser.add_argument(
        "--reference-root",
        default=str(ROOT / "cases" / "reference_parents"),
        help="Hash-verified reference parent catalog (default: cases/reference_parents)",
    )
    args = parser.parse_args()
    result = build_parent_funnel(
        export_root=args.export_root,
        expert_source=args.expert_source,
        output=args.output,
        reference_root=args.reference_root,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

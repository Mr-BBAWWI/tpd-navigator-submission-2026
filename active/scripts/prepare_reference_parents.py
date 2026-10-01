#!/usr/bin/env python3
"""Offline builder for selectable, 6HAZ-aligned literature reference parents."""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.reference_parents import build_reference_parents


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build a new hash-verified reference-parent folder without network access."
    )
    parser.add_argument("--catalog-dir", required=True, help="Original collector-v2 output directory")
    parser.add_argument("--reference-cif", required=True, help="Fixed 6HAZ mmCIF receptor source")
    parser.add_argument("--output", required=True, help="New output folder; must not already exist")
    args = parser.parse_args(argv)
    index = build_reference_parents(args.catalog_dir, args.reference_cif, args.output)
    print(json.dumps({
        "output": args.output,
        "ready_parent_count": len(index["parents"]),
        "not_ready_count": len(index["not_ready"]),
        "default_parent_id": index["default_parent_id"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

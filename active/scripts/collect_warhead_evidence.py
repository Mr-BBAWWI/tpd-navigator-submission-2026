#!/usr/bin/env python3
"""Collect an immutable RCSB-backed warhead evidence catalog."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.warhead_evidence import collect_evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", help="Registered target name or accession (SMARCA2/P51531)")
    parser.add_argument("output", type=Path, help="Output directory")
    parser.add_argument("--curated", type=Path,
                        help="Optional provenance-complete local curated evidence JSON")
    parser.add_argument("--max-entries", type=int, default=32)
    parser.add_argument("--max-final", type=int, default=10)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--domain", choices=("bromodomain", "all_domains"),
                        default="bromodomain",
                        help="Default uses curated 6HAZ construct overlap; all_domains partitions constructs")
    parser.add_argument("--cache-root", type=Path,
                        help="Optional existing immutable RCSB cache; referenced bytes are copied into output")
    parser.add_argument("--fresh", action="store_true",
                        help="Fetch again instead of reusing URL-indexed immutable snapshots")
    args = parser.parse_args()

    catalog = collect_evidence(
        args.target,
        args.output,
        max_entries=args.max_entries,
        max_final=args.max_final,
        timeout=args.timeout,
        curated_path=args.curated,
        reuse_cache=not args.fresh,
        domain=args.domain,
        cache_root=args.cache_root,
    )
    print(json.dumps({
        "status": catalog["status"],
        "catalog": str(args.output / "catalog.json"),
        "cards": len(catalog.get("cards", [])),
        "records": len(catalog.get("all_records", [])),
        "exclusions": len(catalog.get("exclusions", [])),
        "errors": len(catalog.get("errors", [])),
    }, sort_keys=True))
    return 0 if catalog["status"] in {"complete", "partial", "unsupported"} else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Read-only verification CLI for a parent-funnel diagnostic bundle."""
from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence, Set
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.parent_funnel import verify_funnel_bundle


def _plain(value: Any) -> Any:
    """Recursively convert immutable verifier values to JSON-compatible values."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, Set) and not isinstance(value, (str, bytes, bytearray)):
        return [_plain(item) for item in sorted(value, key=repr)]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_plain(item) for item in value]
    return value


def _parent_projection(parent: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "parent_id": _plain(parent.get("parent_id")),
        "parent_redocking": _plain(parent.get("parent_redocking")),
        "actual_counts": _plain(parent.get("actual_counts", {})),
        "evidence_tiers": _plain(parent.get("evidence_tiers", {})),
        "actual_families": _plain(parent.get("actual_families", {})),
    }


def _bounded_projection(result: Mapping[str, Any]) -> dict[str, Any]:
    parents = result.get("parents", ())
    if not isinstance(parents, Sequence) or isinstance(parents, (str, bytes, bytearray)):
        parents = ()
    projected_parents = [
        _parent_projection(parent)
        for parent in parents
        if isinstance(parent, Mapping)
    ]
    return {
        "format": _plain(result.get("format")),
        "status": _plain(result.get("status")),
        "manifest_sha256": _plain(result.get("manifest_sha256")),
        "top_counts": {
            "parent_count": _plain(result.get("parent_count")),
            "human_selection_count": _plain(result.get("human_selection_count")),
            "strict_trusted_selected_count": _plain(result.get("strict_trusted_selected_count")),
        },
        "scientific_approval": _plain(result.get("scientific_approval")),
        "m2_registered": _plain(result.get("m2_registered")),
        "reference_verification": _plain(result.get("reference_verification")),
        "parents": projected_parents,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Verify a parent-funnel bundle without modifying it and print a bounded "
            "computed-diagnostic summary."
        )
    )
    parser.add_argument("--bundle", required=True, help="Existing bundle directory")
    parser.add_argument(
        "--expected-manifest-sha256",
        required=True,
        help="Expected lowercase SHA-256 of manifest.json",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = verify_funnel_bundle(
            args.bundle,
            args.expected_manifest_sha256,
        )
        if not isinstance(result, Mapping):
            raise TypeError("VERIFIER_RESULT_NOT_MAPPING")
        output = _bounded_projection(result)
        print(
            json.dumps(
                output,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )
    except Exception as exc:
        print(
            f"VERIFY_FAILED:{type(exc).__name__}:{exc}",
            file=sys.stderr,
            flush=True,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

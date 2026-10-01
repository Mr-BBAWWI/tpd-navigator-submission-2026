#!/usr/bin/env python3
"""Compile verified campaign evidence and optionally run live research agents."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from packages.agents.campaign_evidence import build_campaign_evidence
from packages.agents.research_campaign import load_diagnostic_evidence, run_campaign


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assessment", type=Path, required=True)
    parser.add_argument("--parent-export", type=Path, action="append", required=True)
    parser.add_argument("--strict-receipt", type=Path)
    parser.add_argument(
        "--diagnostic-json", type=Path, action="append", default=[],
        help="repeatable local diagnostic JSON data (maximum 12, 128 KiB each)",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rounds", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--live-api", action="store_true")
    args = parser.parse_args(argv)

    # Compilation verifies every supplied export closure before a provider exists.
    evidence = build_campaign_evidence(
        args.assessment, args.parent_export, strict_receipt_path=args.strict_receipt
    )
    diagnostics = load_diagnostic_evidence(args.diagnostic_json)
    provider = None
    try:
        if args.live_api:
            import httpx
            from packages.agents.provider import DaconProvider, load_key
            provider = DaconProvider(load_key())
            provider._client.timeout = httpx.Timeout(300, connect=15)
        result = run_campaign(
            evidence, args.output, provider=provider, max_rounds=args.rounds,
            supplemental_evidence=diagnostics,
        )
    finally:
        if provider is not None:
            provider.close()
    return 0 if result["status"] != "failed" else 2


if __name__ == "__main__":
    raise SystemExit(main())

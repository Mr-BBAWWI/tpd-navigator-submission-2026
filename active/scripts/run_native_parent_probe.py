"""Run the deposited 9D12 native parent/Br probe experiment."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(_SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_ROOT))

from packages.science.native_parent_docking import run_native_probe


def _seeds(value: str) -> list[int]:
    try:
        values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("seeds must be comma-separated integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("at least one seed is required")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run an exploratory deposited 9D12 receptor-frame parent/Br docking "
            "comparison. This does not produce an official gate result."
        )
    )
    parser.add_argument("--native-cif", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--parent-id", default="SMARCA2-9D12-A1A1P"
    )
    parser.add_argument("--seeds", type=_seeds, default=[23, 41, 61])
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument(
        "--exclude-remote-incomplete-residue",
        action="append",
        default=[],
        metavar="CHAIN:AUTH_SEQ_ID",
        help=(
            "Explicitly omit a named incomplete canonical residue from the docking "
            "receptor only. The residue must be more than 20 A from every native "
            "ligand heavy atom. May be repeated; no residue is omitted automatically."
        ),
    )
    parser.add_argument(
        "--source-export",
        required=True,
        type=Path,
        metavar="PARENT_EXPORT_DIR",
    )
    arguments = parser.parse_args()
    result = run_native_probe(
        native_cif=arguments.native_cif,
        output=arguments.output,
        parent_id=arguments.parent_id,
        seeds=arguments.seeds,
        exhaustiveness=arguments.exhaustiveness,
        source_export=arguments.source_export,
        exclude_remote_incomplete_residues=arguments.exclude_remote_incomplete_residue,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("status") == "completed_with_limits" else 1


if __name__ == "__main__":
    raise SystemExit(main())

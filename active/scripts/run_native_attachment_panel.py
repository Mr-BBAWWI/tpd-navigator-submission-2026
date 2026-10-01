"""Run the exploratory native-9D12 N3 attachment panel."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from packages.science.native_attachment_experiment import run_native_attachment_panel


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Run the source-supported exploratory native 9D12 N3 attachment panel."
    )
    result.add_argument("--native-probe", required=True, type=Path)
    result.add_argument("--native-cif", required=True, type=Path)
    result.add_argument("--source-export", required=True, type=Path)
    result.add_argument("--output", required=True, type=Path)
    result.add_argument("--seed", type=int, default=23)
    result.add_argument("--exhaustiveness", type=int, default=16)
    return result


def main() -> int:
    arguments = parser().parse_args()
    result = run_native_attachment_panel(
        native_probe=arguments.native_probe,
        native_cif=arguments.native_cif,
        source_export=arguments.source_export,
        output=arguments.output,
        seed=arguments.seed,
        exhaustiveness=arguments.exhaustiveness,
    )
    print(json.dumps({
        "status": result["status"],
        "requested_count": result["requested_count"],
        "retained_count": result["retained_count"],
        "output": str(arguments.output.resolve()),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

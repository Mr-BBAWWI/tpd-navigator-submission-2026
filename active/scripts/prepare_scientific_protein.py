#!/usr/bin/env python3
"""Prepare hash-bound protein hydrogen evidence with local PDB2PQR."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.protein_hydrogen_evidence import prepare_evidence


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdb", required=True, type=Path)
    parser.add_argument("--source-cif", required=True, type=Path)
    parser.add_argument("--pdb2pqr", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--pH", "--ph", dest="pH", type=float, default=7.4)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--fixed-heavy", dest="fixed_heavy", action="store_true", default=True,
        help="Use --noopt --nodebump (default).",
    )
    mode.add_argument(
        "--optimized", dest="fixed_heavy", action="store_false",
        help="Request external hydrogen optimization; scientific review remains pending.",
    )
    args = parser.parse_args(argv)
    try:
        result = prepare_evidence(
            args.source_pdb, args.source_cif, args.outdir, args.pdb2pqr,
            pH=args.pH, fixed_heavy=args.fixed_heavy,
        )
    except (OSError, TypeError, ValueError, RuntimeError):
        sys.stderr.write(json.dumps({
            "status": "failed",
            "error": "protein_preparation_rejected",
        }) + "\n")
        return 2
    sys.stdout.write(json.dumps({
        "status": result["status"],
        "evidence": str(args.outdir.resolve() / "protein_hydrogen_evidence.json"),
        "source_heavy_atoms": result["counts"]["source_pdb_heavy_atoms"],
        "prepared_hydrogens": result["counts"]["prepared_hydrogens"],
        "pending_human_review": result["state_flags"]["pending_human_review"],
    }, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

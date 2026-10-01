#!/usr/bin/env python3
"""Run the hash-bound MolGpKa ligand pKa evidence exporter."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.ligand_pka_evidence import main


if __name__ == "__main__":
    raise SystemExit(main())

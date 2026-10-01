#!/usr/bin/env python3
"""Create machine-only assessment revisions from registered protocols."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.science import scientific_assessment


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reassess immutable registered protocol diagnostics without new computation."
    )
    parser.add_argument("--store", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--job-id", action="append", required=True)
    parser.add_argument("--protocol-id", action="append", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--markdown-output")
    return parser


def _new_output(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists() or path.is_symlink():
        raise SystemExit(label + " already exists")
    if not path.parent.is_dir():
        raise SystemExit(label + " parent directory does not exist")
    return path


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    root = Path(args.store).expanduser().resolve()
    if not root.is_dir():
        raise SystemExit("--store must name an existing Store directory")
    index = root / "index.sqlite3"
    if not index.is_file() or index.is_symlink():
        raise SystemExit("--store has no existing index.sqlite3")
    if not args.project:
        raise SystemExit("--project must not be empty")
    if any(not value for value in args.job_id):
        raise SystemExit("--job-id must not be empty")
    if any(not value for value in args.protocol_id):
        raise SystemExit("--protocol-id must not be empty")
    if len(set(args.job_id)) != len(args.job_id):
        raise SystemExit("duplicate --job-id")
    if len(set(args.protocol_id)) != len(args.protocol_id):
        raise SystemExit("duplicate --protocol-id")

    output = _new_output(args.output, "--output")
    markdown_output = (
        _new_output(args.markdown_output, "--markdown-output")
        if args.markdown_output else None
    )
    if markdown_output is not None and markdown_output == output:
        raise SystemExit("--output and --markdown-output must differ")

    store = Store(root)
    service = ScientificAcceptanceService(store, args.project)
    assessments = []
    for job_id in args.job_id:
        assessment, _created = service.create(
            job_id, protocol_ids=list(args.protocol_id)
        )
        assessments.append(assessment)

    output_value = assessments[0] if len(assessments) == 1 else assessments
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(
            output_value, ensure_ascii=False, sort_keys=True,
            indent=2, allow_nan=False,
        ) + "\n")
    if markdown_output is not None:
        markdown = "\n".join(
            scientific_assessment.report_markdown(assessment)
            for assessment in assessments
        )
        with markdown_output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

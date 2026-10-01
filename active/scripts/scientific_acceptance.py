#!/usr/bin/env python3
"""Local operator CLI for scientific assessment and frozen ternary runs."""
from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.contracts import ContractError, parse_json
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.science import novel_ternary
from scripts import run_novel_ternary

MAX_EXPORT_BYTES = 256 * 1024 * 1024


def _json(value) -> None:
    print(json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False))


def _new_directory(path: Path) -> Path:
    path = path.resolve()
    if path.exists():
        raise ContractError("SCIENTIFIC_OUTPUT_DIRECTORY_EXISTS")
    path.mkdir(parents=True)
    return path


def _service(args) -> ScientificAcceptanceService:
    root = args.store.resolve()
    if not (root / "index.sqlite3").is_file():
        raise ContractError("SCIENTIFIC_EXISTING_STORE_REQUIRED")
    return ScientificAcceptanceService(Store(root), args.project)


def _session(service: ScientificAcceptanceService, reviewer_id: str) -> str:
    key = os.environ.get("TPD_REVIEW_ACCESS_KEY")
    if key is None:
        key = getpass.getpass("TPD review access key: ")
    return service.auth.login(reviewer_id, key)


def _authenticated(args, service):
    if not getattr(args, "reviewer_id", None):
        raise ContractError("SCIENTIFIC_REVIEWER_ID_REQUIRED")
    return _session(service, args.reviewer_id)


def _write_bounded(directory: Path, files: dict[str, bytes]) -> None:
    total = 0
    for name, raw in files.items():
        if not isinstance(raw, bytes):
            raise ContractError("SCIENTIFIC_EXPORT_BYTES_REQUIRED")
        total += len(raw)
        if total > MAX_EXPORT_BYTES:
            raise ContractError("SCIENTIFIC_EXPORT_TOO_LARGE")
        target = directory / name
        if target.parent != directory or target.name != name:
            raise ContractError("SCIENTIFIC_EXPORT_FILENAME")
        with target.open("xb") as handle:
            handle.write(raw)


def _verified_snapshot(service, job_id: str):
    with service.store.db() as db:
        job, result, binding, source = service._verified_result(
            db, job_id, require_current=True
        )
        policy = service._policy(db, job_id, result)
        archive_ref = result["files"]["json"]
        archive_raw = service.port.read(archive_ref)
    return job, result, binding, source, policy, archive_raw


def _check_plan_binding(service, job_id: str, plan: dict):
    novel_ternary.verify_plan(plan)
    _, result, binding, _, policy, _ = _verified_snapshot(service, job_id)
    job_binding = plan.get("bindings", {}).get("job", {})
    if job_binding.get("project") != service.project or job_binding.get("job_id") != job_id:
        raise ContractError("SCIENTIFIC_TERNARY_JOB_BINDING")
    if (job_binding.get("input_sha256") != binding["input_sha256"] or
            job_binding.get("result_sha256") != binding["result_sha256"]):
        raise ContractError("SCIENTIFIC_TERNARY_SOURCE_BINDING")
    candidate_id = plan.get("candidate_graph", {}).get("candidate_id")
    candidate = run_novel_ternary._candidate(result, candidate_id)
    graph = novel_ternary._mapped_graph(candidate)
    if (graph["graph_sha256"] != job_binding.get("candidate_graph_sha256") or
            graph["graph_sha256"] != plan["candidate_graph"].get("graph_sha256")):
        raise ContractError("SCIENTIFIC_TERNARY_GRAPH_BINDING")
    return policy


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--store", type=Path, required=True)
    result.add_argument("--project", default="local-research")
    commands = result.add_subparsers(dest="command", required=True)

    assess = commands.add_parser("assess")
    assess.add_argument("--job-id", required=True)
    assess.add_argument("--reviewer-id")

    compute = commands.add_parser("compute")
    compute.add_argument("--job-id", required=True)
    compute.add_argument("--reviewer-id", required=True)

    policy = commands.add_parser("policy")
    policy.add_argument("--job-id", required=True)

    export = commands.add_parser("export")
    export.add_argument("--assessment-id", required=True)
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--reviewer-id", required=True)

    statement = commands.add_parser("statement")
    statement.add_argument("--text", required=True)
    statement.add_argument("--locator", required=True)
    statement.add_argument("--reviewer-id", required=True)

    prepare = commands.add_parser("prepare-ternary")
    prepare.add_argument("--job-id", required=True)
    prepare.add_argument("--candidate-id", required=True)
    prepare.add_argument("--target-cif", required=True)
    prepare.add_argument("--target-chain", required=True)
    prepare.add_argument("--e3-cif", required=True)
    prepare.add_argument("--e3-chain", required=True)
    prepare.add_argument("--seeds", required=True, type=run_novel_ternary._seeds)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--msa-mode", choices=("single_sequence", "server"),
                         default="single_sequence")

    execute = commands.add_parser("run-novel")
    execute.add_argument("--job-id", required=True)
    execute.add_argument("--plan", type=Path, required=True)
    execute.add_argument("--output", type=Path, required=True)
    execute.add_argument("--boltz-executable", type=Path, required=True)
    execute.add_argument("--checkpoint", type=Path, required=True)
    execute.add_argument("--cache", type=Path, required=True)
    execute.add_argument("--timeout", type=run_novel_ternary._positive_timeout,
                         default=1800.0)

    imported = commands.add_parser("import-ternary")
    imported.add_argument("--job-id", required=True)
    imported.add_argument("--plan", type=Path, required=True)
    imported.add_argument("--run-directory", type=Path, required=True)
    imported.add_argument("--reviewer-id", required=True)

    protein_h = commands.add_parser("import-protein-h")
    protein_h.add_argument("--job-id", required=True)
    protein_h.add_argument("--evidence-json", type=Path, required=True)
    protein_h.add_argument("--reviewer-id", required=True)

    batch = commands.add_parser("batch-plan")
    batch.add_argument("--job-id", required=True)
    batch.add_argument("--output", type=Path, required=True)
    batch.add_argument("--reviewer-id")
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    service = _service(args)

    if args.command == "assess":
        assessment, created = service.create(args.job_id, reviewer_id=args.reviewer_id)
        _json({"created": created, "assessment": assessment})
        return 0

    if args.command == "compute":
        token = _authenticated(args, service)
        try:
            _json(service.compute(args.job_id, token))
        finally:
            service.auth.logout(token)
        return 0

    if args.command == "policy":
        _json(service.policy(args.job_id))
        return 0

    if args.command == "statement":
        token = _authenticated(args, service)
        try:
            _json(service.create_statement(args.text, args.locator, token))
        finally:
            service.auth.logout(token)
        return 0

    if args.command == "export":
        output = _new_directory(args.output)
        token = _authenticated(args, service)
        try:
            packet = service.export(args.assessment_id, token)
            files = {
                "assessment.json": service.port.read(packet["json_ref"]),
                "assessment.md": service.port.read(packet["report_ref"]),
            }
            _write_bounded(output, files)
            _json({"assessment_id": args.assessment_id, "output": str(output),
                   "bytes": sum(map(len, files.values()))})
        except BaseException:
            shutil.rmtree(output, ignore_errors=True)
            raise
        finally:
            service.auth.logout(token)
        return 0

    if args.command == "prepare-ternary":
        output = _new_directory(args.output)
        try:
            _, _, binding, _, _, archive_raw = _verified_snapshot(service, args.job_id)
            design_path = output / "design-result.json"
            design_path.write_bytes(archive_raw)
            plan_path = output / "plan.json"
            cli = [
                "--design-json", str(design_path),
                "--candidate-id", args.candidate_id,
                "--target-cif", args.target_cif,
                "--target-chain", args.target_chain,
                "--e3-cif", args.e3_cif,
                "--e3-chain", args.e3_chain,
                "--job-id", args.job_id,
                "--project", service.project,
                "--output", str(plan_path),
                "--seeds", ",".join(str(seed) for seed in args.seeds),
                "--msa-mode", args.msa_mode,
                "--prepare-only",
            ]
            run_novel_ternary.main(cli)
            plan = json.loads(plan_path.read_text(encoding="utf-8"))
            _check_plan_binding(service, args.job_id, plan)
            _json({"status": "prepared", "plan": str(plan_path),
                   "result_sha256": binding["result_sha256"],
                   "operator_action": (
                       "Run: scripts/scientific_acceptance.py --store <store> "
                       f"--project {service.project} run-novel --job-id {args.job_id} "
                       f"--plan {plan_path} --output <new-run-directory> "
                       "--boltz-executable <path> --checkpoint <path> --cache <path>"
                   )})
        except BaseException:
            shutil.rmtree(output, ignore_errors=True)
            raise
        return 0

    if args.command == "run-novel":
        plan_path = args.plan.resolve()
        if not plan_path.is_file() or plan_path.is_symlink():
            raise ContractError("SCIENTIFIC_TERNARY_PLAN_FILE")
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        _check_plan_binding(service, args.job_id, plan)
        output = args.output.resolve()
        if output.exists():
            raise ContractError("SCIENTIFIC_OUTPUT_DIRECTORY_EXISTS")
        try:
            receipt = novel_ternary.run_plan(
                plan,
                output_root=output,
                executable=args.boltz_executable,
                checkpoint=args.checkpoint,
                cache=args.cache,
                timeout=args.timeout,
            )
            _json({"status": "finished", "run_directory": str(output),
                   "receipt": str(output / "receipt.json"), "summary": receipt})
        except BaseException:
            # Preserve every worker receipt, log, and partial output for failure review.
            raise
        return 0

    if args.command == "import-ternary":
        token = _authenticated(args, service)
        try:
            _json(service.import_ternary(
                args.job_id, args.plan, args.run_directory, token
            ))
        finally:
            service.auth.logout(token)
        return 0

    if args.command == "import-protein-h":
        token = _authenticated(args, service)
        try:
            _json(service.import_protein_hydrogens(
                args.job_id, args.evidence_json, token
            ))
        finally:
            service.auth.logout(token)
        return 0

    if args.command == "batch-plan":
        output = _new_directory(args.output)
        try:
            from packages.science.dual_e3 import DEFAULT_PARENT_ID, catalog
            data = catalog()["reference_parents"]
            available = [DEFAULT_PARENT_ID]
            for row in data.get("available", []):
                if (isinstance(row, dict) and isinstance(row.get("id"), str) and
                        (row.get("ready") is True or row.get("technical_ready") is True)):
                    available.append(row["id"])
            selected = []
            authenticated = False
            if args.reviewer_id:
                token = _authenticated(args, service)
                try:
                    service.auth.session(token)
                    authenticated = True
                    with service.store.db() as db:
                        _, result, _, _ = service._verified_result(db, args.job_id)
                        policy = service._policy(db, args.job_id, result)
                    selected = [row["parent_id"] for row in
                                policy.get("decisions", {}).get("parent_funnel", [])]
                finally:
                    service.auth.logout(token)
            requests = [{
                "result_id": "design:SMARCA2",
                "operation": "design_panel",
                "parameters": {"parent_id": parent_id},
                "manual_submission_required": True,
                "auto_approved": False,
                "auto_executed": False,
            } for parent_id in selected]
            document = {
                "format": "scientific-batch-plan/1",
                "job_id": args.job_id,
                "catalog_available_known_parents": sorted(set(available)),
                "authenticated": authenticated,
                "stored_expert_parent_funnel": selected,
                "new_job_requests": requests,
                "condition_absence": None if selected else
                    "No authenticated active parent_funnel decision; no job requests invented.",
            }
            raw = json.dumps(document, ensure_ascii=False, sort_keys=True,
                             allow_nan=False).encode("utf-8")
            _write_bounded(output, {"batch-plan.json": raw})
            _json({"output": str(output), "request_count": len(requests),
                   "auto_approved": False, "auto_executed": False})
        except BaseException:
            shutil.rmtree(output, ignore_errors=True)
            raise
        return 0

    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())

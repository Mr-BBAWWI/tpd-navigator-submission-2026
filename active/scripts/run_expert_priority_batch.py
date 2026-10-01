#!/usr/bin/env python3
"""Prepare or serially run the fixed 3-analog × dual-E3 × 5-seed exploration."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.store import Store
from packages.science import novel_ternary
from packages.science.novel_ternary import prepare_plan, run_plan, verify_plan, write_plan
from scripts.run_novel_ternary import _candidate

ANALOG_IDS = ("W-c2afc5e73c1a", "W-4c0a639c0a41", "W-80f8f4a11b5d")
E3_ORDER = ("VHL", "CRBN")
DEFAULT_SEEDS = (23, 41, 61, 79, 97)
LINKER_ID = "alkyl_c6"


def _digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                         ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def _seeds(value: str) -> list[int]:
    try:
        result = [int(item.strip()) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("seeds must be comma-separated integers") from error
    if len(result) != 5 or len(set(result)) != 5 or any(seed < 0 or seed >= 2**32 for seed in result):
        raise argparse.ArgumentTypeError("exactly five unique unsigned seeds are required")
    return result


def _linker(candidate: dict) -> str | None:
    value = candidate.get("linker_id", candidate.get("linker"))
    if isinstance(value, dict):
        value = value.get("id", value.get("linker_id"))
    return value


def select_priority_candidates(document: dict) -> list[dict]:
    rows = document.get("protac_candidates")
    if not isinstance(rows, list):
        raise ValueError("DESIGN_PROTAC_CANDIDATES_LIST_REQUIRED")
    scope = document.get("parent_scope")
    if not isinstance(scope, dict):
        raise ValueError("PRIORITY_PARENT_SCOPE_REQUIRED")
    if scope.get("actual_design_parent_id") != "SMARCA2-FX5":
        raise ValueError("PRIORITY_PARENT_SCOPE_MISMATCH")
    if scope.get("receptor_frame") != "6HAZ chain A":
        raise ValueError("PRIORITY_STRUCTURAL_FRAME_MISMATCH")
    selected = []
    seen = set()
    graph_hashes = set()
    for analog_id in ANALOG_IDS:
        for e3 in E3_ORDER:
            matches = [row for row in rows if isinstance(row, dict) and
                       row.get("warhead_analog_id") == analog_id and row.get("e3_type") == e3 and
                       _linker(row) == LINKER_ID]
            if len(matches) != 1:
                raise ValueError("PRIORITY_CANDIDATE_NOT_UNIQUE")
            candidate_id = matches[0].get("candidate_id")
            if not isinstance(candidate_id, str) or not candidate_id or candidate_id in seen:
                raise ValueError("PRIORITY_CANDIDATE_ID_INVALID_OR_DUPLICATE")
            checked = _candidate(document, candidate_id)
            if checked.get("e3_type") != e3 or checked.get("warhead_analog_id") != analog_id:
                raise ValueError("PRIORITY_CANDIDATE_BINDING_MISMATCH")
            graph = novel_ternary._mapped_graph(checked)
            if graph["graph_sha256"] in graph_hashes:
                raise ValueError("PRIORITY_CANDIDATE_GRAPH_DUPLICATE")
            graph_hashes.add(graph["graph_sha256"])
            seen.add(candidate_id)
            selected.append(checked)
    if len(selected) != 6 or {item["e3_type"] for item in selected} != set(E3_ORDER):
        raise ValueError("PRIORITY_E3_COVERAGE_INVALID")
    return selected


def _policy_binding(service: ScientificAcceptanceService, policy: dict) -> dict:
    revision = policy.get("revision")
    if type(revision) is not int or revision < 0:
        raise ValueError("POLICY_REVISION_REQUIRED")
    policy_digest = service._policy_digest(policy)
    if not isinstance(policy_digest, str) or not policy_digest:
        raise ValueError("POLICY_DIGEST_REQUIRED")
    return {"module": novel_ternary.__name__,
            "revision": str(revision), "digest": policy_digest}


def _immutable_plan_payload(plan: dict) -> dict:
    value = dict(plan)
    value.pop("created_utc", None)
    value.pop("plan_digest", None)
    return value


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8", newline="\n")
    temporary.replace(path)


def _inside(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=True).relative_to(root.resolve(strict=True))
        return True
    except (OSError, ValueError):
        return False


def _source(path: Path, chain: str, role: str, terms: list[str],
            expected_sha256: str | None = None) -> dict:
    return {"path": str(path), "label_asym_id": chain,
            "expected_sha256": expected_sha256 or novel_ternary.worker.sha256_file(path),
            "role": role, "description_terms": terms}


def _completed(directory: Path, plan: dict, execution_binding: dict) -> bool:
    receipt_path = directory / "receipt.json"
    if not receipt_path.is_file() or receipt_path.is_symlink():
        return False
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    expected = plan["seeds"]
    if receipt.get("plan_digest") != plan["plan_digest"] or receipt.get("completed_count") != len(expected):
        return False
    paths = receipt.get("all_seed_receipts")
    if not isinstance(paths, list) or len(paths) != len(expected):
        return False
    observed = []
    for value in paths:
        path = Path(value)
        if not path.is_file() or path.is_symlink() or not _inside(path, directory):
            return False
        item = json.loads(path.read_text(encoding="utf-8"))
        if item.get("status") != "completed" or item.get("plan_digest") != plan["plan_digest"]:
            return False
        hashes = item.get("hashes", {})
        if hashes.get("checkpoint_sha256") != execution_binding["checkpoint_sha256"] or \
                hashes.get("executable_sha256") != execution_binding["executable_sha256"]:
            return False
        observed.append(item.get("seed"))
        outputs = hashes.get("output_files_sha256")
        if not isinstance(outputs, dict) or not outputs:
            return False
        for output, digest in outputs.items():
            file_path = path.parent / "boltz_output" / output
            if not file_path.is_file() or file_path.is_symlink() or not _inside(file_path, directory):
                return False
            if novel_ternary.worker.sha256_file(file_path) != digest:
                return False
    return observed == expected


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--store", required=True, type=Path)
    result.add_argument("--job-id", required=True)
    result.add_argument("--output", required=True, type=Path)
    result.add_argument("--boltz-executable", type=Path)
    result.add_argument("--checkpoint", type=Path)
    result.add_argument("--cache", type=Path)
    result.add_argument("--seeds", type=_seeds, default=list(DEFAULT_SEEDS))
    result.add_argument("--prepare-only", action="store_true")
    result.add_argument("--resume", action="store_true")
    result.add_argument("--project", default="local-research")
    result.add_argument("--timeout", type=float, default=1800.0)
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if not args.prepare_only and not all((args.boltz_executable, args.checkpoint, args.cache)):
        parser().error("execution requires --boltz-executable, --checkpoint, and --cache")
    target = ROOT / "cases/design_sources/6HAZ.cif"
    e3_paths = {
        "VHL": ROOT / "cases/acceptance_sources/supporting/catalog/source_snapshots/f9a7ab0aca167f758cad8d77b2173cdfbd7ab5edd2cb91bc6dacd6e9901b5d30/6HAY.cif",
        "CRBN": ROOT / "cases/design_sources/6BOY.cif",
    }
    for source_path in (target, *e3_paths.values()):
        novel_ternary.worker.sha256_file(source_path)
    if args.output.exists() and not args.resume:
        raise ValueError("BATCH_OUTPUT_EXISTS")
    if not args.output.exists():
        args.output.mkdir(parents=True)
    service = ScientificAcceptanceService(Store(args.store.resolve()), args.project)
    with service.store.db() as db:
        job, result, binding, source = service._verified_result(db, args.job_id, require_current=True)
        policy = service._policy(db, args.job_id, result)
        archive_raw = service.port.read(result["files"]["json"])
    design = json.loads(archive_raw)
    candidates = select_priority_candidates(design)
    policy_binding = _policy_binding(service, policy)
    progress = args.output / "batch-progress.jsonl"
    manifest_rows = []
    failures = 0
    execution_binding = None
    if not args.prepare_only:
        execution_binding = {
            "project": args.project, "job_id": args.job_id,
            "policy_digest": policy_binding["digest"],
            "policy_revision": policy_binding["revision"],
            "executable_sha256": novel_ternary.worker.sha256_file(args.boltz_executable),
            "checkpoint_sha256": novel_ternary.worker.sha256_file(args.checkpoint),
        }
        binding_path = args.output / "execution-binding.json"
        if binding_path.exists():
            if json.loads(binding_path.read_text(encoding="utf-8")) != execution_binding:
                raise ValueError("EXECUTION_BINDING_MISMATCH")
        else:
            _atomic_json(binding_path, execution_binding)
    for candidate in candidates:
        e3 = candidate["e3_type"]
        name = f"{candidate['warhead_analog_id']}--{e3}"
        directory, plan_path = args.output / name, args.output / (name + ".plan.json")
        graph = novel_ternary._mapped_graph(candidate)
        job_binding = {"project": args.project, "job_id": args.job_id,
                       "input_sha256": binding["input_sha256"],
                       "result_sha256": binding["result_sha256"],
                       "candidate_graph_sha256": graph["graph_sha256"]}
        plan = prepare_plan(candidate,
            target_source=_source(target, "A", "target", ["SMARCA2", "SNF2L2"]),
            e3_source=_source(e3_paths[e3], "B", "e3",
                              ["von Hippel-Lindau", "VHL"] if e3 == "VHL" else ["cereblon", "CRBN"],
                              "f9a7ab0aca167f758cad8d77b2173cdfbd7ab5edd2cb91bc6dacd6e9901b5d30"
                              if e3 == "VHL" else None),
            job_binding=job_binding, policy_binding=policy_binding, seeds=args.seeds,
            msa_mode="single_sequence")
        if plan_path.exists():
            frozen = json.loads(plan_path.read_text(encoding="utf-8"))
            verify_plan(frozen)
            if _immutable_plan_payload(frozen) != _immutable_plan_payload(plan):
                raise ValueError("IMMUTABLE_PLAN_MISMATCH")
            plan = frozen
        else:
            write_plan(plan, plan_path)
        status = "prepared"
        if not args.prepare_only:
            if directory.exists():
                if args.resume and _completed(directory, plan, execution_binding):
                    status = "reused_verified_completed"
                else:
                    status = "attention_required_incomplete_existing_run"
                    failures += 1
            else:
                try:
                    receipt = run_plan(plan, output_root=directory,
                                       executable=args.boltz_executable, checkpoint=args.checkpoint,
                                       cache=args.cache, timeout=args.timeout, accelerator="gpu",
                                       max_msa_seqs=256)
                    status = "completed" if receipt.get("completed_count") == 5 else "completed_with_failures"
                    failures += int(status != "completed")
                except Exception as error:
                    status = "run_exception:" + (str(error) or type(error).__name__)
                    failures += 1
        row = {"analog_id": candidate["warhead_analog_id"], "e3_type": e3,
               "candidate_id": candidate["candidate_id"], "linker_id": LINKER_ID,
               "plan": str(plan_path), "run_directory": str(directory), "status": status}
        manifest_rows.append(row)
        with progress.open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
    manifest = {"format": "tpd-expert-priority-exploration/1", "runs": manifest_rows,
                "seeds": args.seeds, "expected_gpu_runs": 30, "serial": True,
                "msa_mode": "single_sequence",
                "comparison_limitation": "Novel runs have no MSA and are not a matched-MSA baseline against known server-MSA receipts.",
                "followup_frozen_msa": "pending_configurable_protocol",
                "selection_scope": "expert-selected analog and E3 identities; developer representative alkyl_c6 linker",
                "approval": False, "auto_authorization": False,
                "interpretation": "Exploration only; no cross-E3 winner, efficacy, degradation, or scientific approval claim."}
    immutable_batch = {key: value for key, value in manifest.items() if key != "runs"}
    immutable_batch["runs"] = [{key: row[key] for key in
                                ("analog_id", "e3_type", "candidate_id", "linker_id", "plan", "run_directory")}
                               for row in manifest_rows]
    batch_plan_path = args.output / "batch-plan.json"
    if batch_plan_path.exists():
        if json.loads(batch_plan_path.read_text(encoding="utf-8")) != immutable_batch:
            raise ValueError("IMMUTABLE_MANIFEST_MISMATCH")
    else:
        _atomic_json(batch_plan_path, immutable_batch)
    _atomic_json(args.output / "manifest.json", manifest)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

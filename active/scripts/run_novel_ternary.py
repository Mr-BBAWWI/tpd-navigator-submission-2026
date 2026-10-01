#!/usr/bin/env python3
"""Prepare and optionally execute a frozen, reference-free novel ternary plan."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.boltz_worker import sha256_file
from packages.science import novel_ternary
from packages.science.novel_ternary import (
    _mapped_graph,
    prepare_plan,
    run_plan,
    write_plan,
)


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _actual_document(value):
    """Resolve the known stored-result wrapper without accepting arbitrary nesting."""
    if not isinstance(value, dict):
        raise ValueError("DESIGN_DOCUMENT_OBJECT_REQUIRED")
    if "protac_candidates" in value:
        return value
    for key in ("result", "job_result"):
        wrapped = value.get(key)
        if isinstance(wrapped, str):
            try:
                wrapped = json.loads(wrapped)
            except json.JSONDecodeError as error:
                raise ValueError("DESIGN_RESULT_WRAPPER_INVALID") from error
        if isinstance(wrapped, dict) and "protac_candidates" in wrapped:
            return wrapped
    raise ValueError("DESIGN_PROTAC_CANDIDATES_REQUIRED")


def _digest_string(value, code: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(code)
    return value


def _qualified_analog(document: dict, analog_id: str) -> dict:
    analogs = document.get("analogs")
    if not isinstance(analogs, list):
        raise ValueError("DESIGN_ANALOGS_LIST_REQUIRED")
    matches = [item for item in analogs if isinstance(item, dict) and item.get("id") == analog_id]
    if len(matches) != 1:
        raise ValueError("WARHEAD_ANALOG_ID_NOT_UNIQUE")
    analog = matches[0]
    if analog.get("selected") is not True:
        raise ValueError("WARHEAD_ANALOG_NOT_SELECTED")
    if analog.get("pipeline_status") != "qualified":
        raise ValueError("WARHEAD_ANALOG_PIPELINE_NOT_QUALIFIED")
    if analog.get("assembly_eligible") is not True:
        raise ValueError("WARHEAD_ANALOG_NOT_ASSEMBLY_ELIGIBLE")
    if analog.get("parent_redocking_supported") is not True:
        raise ValueError("WARHEAD_ANALOG_PARENT_REDOCKING_NOT_SUPPORTED")
    docking = analog.get("docking")
    if not isinstance(docking, dict):
        raise ValueError("WARHEAD_ANALOG_DOCKING_REQUIRED")
    if docking.get("status") != "completed_with_limits":
        raise ValueError("WARHEAD_ANALOG_DOCKING_STATUS_INVALID")
    if docking.get("pose_preserved") is not True:
        raise ValueError("WARHEAD_ANALOG_DOCKING_POSE_NOT_PRESERVED")
    passing = docking.get("passing_pose_count")
    if type(passing) is not int or passing <= 0:
        raise ValueError("WARHEAD_ANALOG_DOCKING_NO_PASSING_POSE")
    files = docking.get("files")
    poses_sdf = files.get("poses_sdf") if isinstance(files, dict) else None
    if not isinstance(poses_sdf, dict):
        raise ValueError("WARHEAD_ANALOG_DOCKING_POSES_SDF_REF_REQUIRED")
    return {
        "id": analog_id,
        "selected": analog["selected"],
        "pipeline_status": analog["pipeline_status"],
        "assembly_eligible": analog["assembly_eligible"],
        "parent_redocking_supported": analog["parent_redocking_supported"],
        "docking": {
            "status": docking["status"],
            "pose_preserved": docking["pose_preserved"],
            "passing_pose_count": passing,
            "files": {"poses_sdf": json.loads(json.dumps(poses_sdf, allow_nan=False))},
        },
    }


def _candidate(document: dict, candidate_id: str) -> dict:
    candidates = document.get("protac_candidates")
    if not isinstance(candidates, list):
        raise ValueError("DESIGN_PROTAC_CANDIDATES_LIST_REQUIRED")
    matches = [item for item in candidates
               if isinstance(item, dict) and item.get("candidate_id") == candidate_id]
    if len(matches) != 1:
        raise ValueError("CANDIDATE_ID_NOT_UNIQUE")
    candidate = dict(matches[0])
    if candidate.get("assembly_mode") != "pose_supported_hypothesis":
        raise ValueError("CANDIDATE_NOT_POSE_SUPPORTED_HYPOTHESIS")
    analog_id = candidate.get("warhead_analog_id")
    if not isinstance(analog_id, str) or not analog_id.strip():
        raise ValueError("CANDIDATE_WARHEAD_ANALOG_ID_REQUIRED")
    candidate["analog_qualification"] = _qualified_analog(document, analog_id)
    return candidate


def _seeds(value: str) -> list[int]:
    fields = value.split(",")
    if not fields or any(not field.strip() for field in fields):
        raise argparse.ArgumentTypeError("seeds must be non-empty comma-separated integers")
    try:
        result = [int(field.strip()) for field in fields]
    except ValueError as error:
        raise argparse.ArgumentTypeError("seeds must be integers") from error
    if not 0 < len(result) <= 32:
        raise argparse.ArgumentTypeError("between 1 and 32 seeds are required")
    if len(result) != len(set(result)):
        raise argparse.ArgumentTypeError("seeds must be unique")
    if any(type(seed) is not int or seed < 0 or seed >= 2**32 for seed in result):
        raise argparse.ArgumentTypeError("seeds must be unsigned 32-bit integers")
    return result


def _positive_timeout(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("timeout must be numeric") from error
    if result <= 0 or not math.isfinite(result):
        raise argparse.ArgumentTypeError("timeout must be positive and finite")
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--design-json", required=True)
    result.add_argument("--candidate-id", required=True)
    result.add_argument("--target-cif", required=True)
    result.add_argument("--target-chain", required=True)
    result.add_argument("--e3-cif", required=True)
    result.add_argument("--e3-chain", required=True)
    result.add_argument("--job-id", required=True,
                        help="operator-declared stored-job identity; validated by the service store")
    result.add_argument("--project", required=True,
                        help="operator-declared project identity; not live authentication")
    result.add_argument("--output", required=True)
    result.add_argument("--seeds", required=True, type=_seeds)
    result.add_argument("--msa-mode", choices=("single_sequence", "server"),
                        default="single_sequence")
    result.add_argument("--boltz-executable")
    result.add_argument("--checkpoint")
    result.add_argument("--cache")
    result.add_argument("--timeout", type=_positive_timeout, default=1800.0)
    result.add_argument("--policy-json",
                        help="optional service assessment-policy binding JSON; no expert-policy default")
    result.add_argument("--prepare-only", action="store_true")
    return result


def _policy_binding(path: str | None, max_msa_seqs: int = 256) -> dict:
    module_path = Path(novel_ternary.__file__).resolve()
    material = {
        "module": "packages.science.novel_ternary",
        "module_sha256": sha256_file(module_path),
        "run_settings": {"boltz_version": novel_ternary.SUPPORTED_BOLTZ_VERSION,
                         "msa_mode_source": "CLI", "max_msa_seqs": max_msa_seqs,
                         "recycling_steps": 3, "sampling_steps": 200,
                         "diffusion_samples": 1, "no_kernels": True,
                         "use_potentials": False, "affinity_model": False},
        "developer_protocol_version": "novel-ternary-developer-protocol/1",
        "classification": "exploratory_not_expert",
    }
    revision = material["developer_protocol_version"]
    if path is not None:
        policy = json.loads(Path(path).read_bytes())
        if not isinstance(policy, dict):
            raise ValueError("ASSESSMENT_POLICY_OBJECT_REQUIRED")
        required = {"module", "revision", "digest"}
        if set(policy) != required:
            raise ValueError("ASSESSMENT_POLICY_FIELDS_INVALID")
        if policy["module"] != material["module"]:
            raise ValueError("ASSESSMENT_POLICY_MODULE_INVALID")
        if not isinstance(policy["revision"], str) or not policy["revision"]:
            raise ValueError("ASSESSMENT_POLICY_REVISION_INVALID")
        _digest_string(policy["digest"], "ASSESSMENT_POLICY_DIGEST_INVALID")
        return dict(policy)
    digest = _sha_bytes(json.dumps(material, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=False).encode("utf-8"))
    return {"module": material["module"], "revision": revision, "digest": digest}


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    design_path = Path(args.design_json)
    raw = design_path.read_bytes()
    root_document = json.loads(raw)
    design = _actual_document(root_document)
    candidate = _candidate(design, args.candidate_id)
    e3_type = candidate.get("e3_type")
    if e3_type not in {"CRBN", "VHL"}:
        raise ValueError("E3_TYPE_INVALID")

    input_binding = design.get("input_binding")
    if not isinstance(input_binding, dict):
        raise ValueError("DESIGN_INPUT_BINDING_REQUIRED")
    if not isinstance(input_binding.get("format"), str) or not input_binding["format"]:
        raise ValueError("DESIGN_INPUT_BINDING_FORMAT_INVALID")
    source_files = input_binding.get("source_files")
    if source_files is not None and not isinstance(source_files, (list, dict)):
        raise ValueError("DESIGN_INPUT_BINDING_SOURCE_FILES_INVALID")
    input_digest = _digest_string(input_binding.get("digest"),
                                  "DESIGN_INPUT_BINDING_DIGEST_INVALID")
    graph_digest = _mapped_graph(candidate)["graph_sha256"]
    result_digest = _sha_bytes(raw)
    if result_digest == graph_digest:
        raise ValueError("RESULT_AND_CANDIDATE_GRAPH_DIGEST_NOT_DISTINCT")
    job_binding = {
        "project": args.project,
        "job_id": args.job_id,
        "input_sha256": input_digest,
        "result_sha256": result_digest,
        "candidate_graph_sha256": graph_digest,
    }

    target_path, e3_path = Path(args.target_cif), Path(args.e3_cif)
    target_source = {
        "path": str(target_path), "label_asym_id": args.target_chain,
        "expected_sha256": sha256_file(target_path), "role": "target",
        "description_terms": ["SMARCA2", "SNF2L2"],
    }
    e3_source = {
        "path": str(e3_path), "label_asym_id": args.e3_chain,
        "expected_sha256": sha256_file(e3_path), "role": "e3",
        "description_terms": (["cereblon", "CRBN"] if e3_type == "CRBN"
                              else ["von Hippel-Lindau", "VHL"]),
    }
    policy_binding = _policy_binding(args.policy_json)
    plan = prepare_plan(candidate, target_source=target_source, e3_source=e3_source,
                        job_binding=job_binding, policy_binding=policy_binding,
                        seeds=args.seeds, msa_mode=args.msa_mode)
    output = Path(args.output)
    if args.prepare_only:
        record = write_plan(plan, output)
        print(json.dumps({"status": "prepared", "artifact": record}, sort_keys=True))
        return 0
    if not args.boltz_executable or not args.checkpoint or not args.cache:
        parser().error("execution requires --boltz-executable, --checkpoint, and --cache")
    output.parent.mkdir(parents=True, exist_ok=True)
    plan_path = output.parent / (output.name + ".plan.json")
    plan_artifact = write_plan(plan, plan_path)
    receipt = run_plan(plan, output_root=output, executable=Path(args.boltz_executable),
                       checkpoint=Path(args.checkpoint), cache=Path(args.cache),
                       timeout=args.timeout)
    print(json.dumps({"status": "finished", "plan": plan_artifact,
                      "receipt": str(output / "receipt.json"), "summary": receipt},
                     ensure_ascii=False, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

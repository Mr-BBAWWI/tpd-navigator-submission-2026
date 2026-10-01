#!/usr/bin/env python3
"""검증된 입력으로 전문가 후속 검토용 릴리스 보고서를 생성한다.

이 스크립트는 과학적 승인을 만들거나 추론하지 않는다. 평가 JSON의 기준 상태를
그대로 보존하고, 전문가 패킷 manifest를 검증한 뒤 독립 실행형 HTML/Markdown 및
측정 요약을 생성한다.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import shutil
import stat
import tempfile
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

MAX_JSON_BYTES = 64 * 1024 * 1024
MAX_TREE_FILES = 20000
MAX_TREE_BYTES = 20 * 1024 * 1024 * 1024
EXPECTED_CRITERIA = 14
EXPECTED_BATCH_PLANS = 6
EXPECTED_BATCH_RECEIPTS = 30
EXPECTED_SEEDS = [23, 41, 61, 79, 97]
EXPECTED_ANALOGS = ["W-c2afc5e73c1a", "W-4c0a639c0a41", "W-80f8f4a11b5d"]
ALLOWED_STATUSES = {"pass", "failed", "pending"}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
SECRET_PATH_RE = re.compile(
    r"(^|/)(?:\.env(?:\..*)?|id_rsa(?:\.pub)?|credentials?(?:\..*)?|"
    r"secrets?(?:\..*)?|api[_-]?keys?(?:\..*)?|prompt[_-]?logs?(?:\..*)?|"
    r"model[_-]?weights?(?:\..*)?)(?:$|/)", re.IGNORECASE
)
SECRET_KEY_RE = re.compile(
    r"(?:api[_-]?key|access[_-]?token|refresh[_-]?token|authorization|password|"
    r"private[_-]?key|client[_-]?secret|prompt[_-]?log|model[_-]?weights)",
    re.IGNORECASE,
)
ABS_WINDOWS_RE = re.compile(r"^[A-Za-z]:[\\/]")


class ReleaseError(RuntimeError):
    """입력 검증 또는 안전한 출력 생성 실패."""


def _is_reparse(path: Path) -> bool:
    try:
        attributes = getattr(path.lstat(), "st_file_attributes", 0)
        marker = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
        return bool(attributes & marker)
    except OSError as exc:
        raise ReleaseError(f"경로 상태를 읽을 수 없습니다: {path}: {exc}") from exc


def require_regular_file(path: Path, label: str, max_bytes: Optional[int] = None) -> None:
    try:
        info = path.lstat()
    except OSError as exc:
        raise ReleaseError(f"{label} 파일을 읽을 수 없습니다: {path}: {exc}") from exc
    if stat.S_ISLNK(info.st_mode) or _is_reparse(path):
        raise ReleaseError(f"{label}에 symlink/junction/reparse point를 사용할 수 없습니다: {path}")
    if not stat.S_ISREG(info.st_mode):
        raise ReleaseError(f"{label}은 regular file이어야 합니다: {path}")
    if max_bytes is not None and info.st_size > max_bytes:
        raise ReleaseError(f"{label} 파일이 허용 크기를 초과했습니다: {info.st_size} bytes")


def require_real_directory(path: Path, label: str) -> None:
    try:
        info = path.lstat()
    except OSError as exc:
        raise ReleaseError(f"{label} 디렉터리를 읽을 수 없습니다: {path}: {exc}") from exc
    if stat.S_ISLNK(info.st_mode) or _is_reparse(path):
        raise ReleaseError(f"{label}에 symlink/junction/reparse point를 사용할 수 없습니다: {path}")
    if not stat.S_ISDIR(info.st_mode):
        raise ReleaseError(f"{label}은 실제 디렉터리여야 합니다: {path}")


def iter_safe_tree(root: Path, label: str) -> List[Path]:
    require_real_directory(root, label)
    files: List[Path] = []
    total = 0
    for current, dirs, names in os.walk(root, topdown=True, followlinks=False):
        current_path = Path(current)
        for name in list(dirs):
            child = current_path / name
            require_real_directory(child, label)
        for name in names:
            child = current_path / name
            require_regular_file(child, label)
            files.append(child)
            total += child.stat().st_size
            if len(files) > MAX_TREE_FILES:
                raise ReleaseError(f"{label} 파일 수가 제한을 초과했습니다")
            if total > MAX_TREE_BYTES:
                raise ReleaseError(f"{label} 전체 크기가 제한을 초과했습니다")
    return files


def reject_constant(value: str) -> None:
    raise ValueError(f"JSON 비유한 숫자는 허용되지 않습니다: {value}")


def load_json(path: Path, label: str) -> Any:
    require_regular_file(path, label, MAX_JSON_BYTES)
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle, parse_constant=reject_constant)
    except (UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ReleaseError(f"{label} JSON이 유효하지 않습니다: {path}: {exc}") from exc


def dump_json(path: Path, value: Any) -> None:
    text = json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.write_text(text, encoding="utf-8", newline="\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_manifest_path(raw: Any) -> str:
    if not isinstance(raw, str) or not raw or "\\" in raw or "\x00" in raw:
        raise ReleaseError(f"manifest 경로가 안전하지 않습니다: {raw!r}")
    pure = PurePosixPath(raw)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        raise ReleaseError(f"manifest 경로가 안전하지 않습니다: {raw!r}")
    return pure.as_posix()


def validate_packet(packet: Path) -> Tuple[Dict[str, Any], List[str]]:
    files = iter_safe_tree(packet, "expert packet")
    manifest_path = packet / "manifest.json"
    manifest = load_json(manifest_path, "expert packet manifest")
    if not isinstance(manifest, dict):
        raise ReleaseError("expert packet manifest는 객체여야 합니다")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ReleaseError("expert packet manifest.artifacts가 배열이 아닙니다")

    expected: Dict[str, Mapping[str, Any]] = {}
    for item in artifacts:
        if not isinstance(item, dict):
            raise ReleaseError("manifest artifact 항목은 객체여야 합니다")
        rel = safe_manifest_path(item.get("path"))
        if rel == "manifest.json" or rel in expected:
            raise ReleaseError(f"중복되거나 금지된 manifest 경로입니다: {rel}")
        if SECRET_PATH_RE.search("/" + rel):
            raise ReleaseError(f"secret/prompt/model 파일은 packet에 포함할 수 없습니다: {rel}")
        digest = item.get("sha256")
        size = item.get("bytes")
        if not isinstance(digest, str) or not HEX64.fullmatch(digest):
            raise ReleaseError(f"artifact SHA-256 형식이 잘못되었습니다: {rel}")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise ReleaseError(f"artifact bytes 형식이 잘못되었습니다: {rel}")
        expected[rel] = item

    actual = {p.relative_to(packet).as_posix() for p in files}
    allowed = set(expected) | {"manifest.json"}
    additional = sorted(actual - allowed)
    missing = sorted(allowed - actual)
    if additional:
        raise ReleaseError("manifest에 없는 파일이 있습니다(비밀 파일도 허용되지 않음): " + ", ".join(additional))
    if missing:
        raise ReleaseError("manifest 파일이 누락되었습니다: " + ", ".join(missing))

    for rel, item in expected.items():
        path = packet.joinpath(*PurePosixPath(rel).parts)
        require_regular_file(path, f"packet artifact {rel}")
        if path.stat().st_size != item["bytes"]:
            raise ReleaseError(f"packet artifact 크기 불일치: {rel}")
        if sha256_file(path) != item["sha256"]:
            raise ReleaseError(f"packet artifact SHA-256 불일치: {rel}")
    return manifest, sorted(actual)


def _classify_batch_record(path: Path, value: Any) -> Optional[str]:
    name = path.name.lower()
    if name in {"plan.json", "batch-plan.json"} or name.endswith(".plan.json"):
        return "plan"
    if name in {"receipt.json", "run-receipt.json"} or name.endswith(".receipt.json"):
        return "receipt"
    if isinstance(value, dict):
        if isinstance(value.get("seed"), int) and any(k in value for k in ("status", "exit_code", "execution_status")):
            return "receipt"
        seeds = value.get("seeds", value.get("plan_seed"))
        if isinstance(seeds, list) and all(isinstance(seed, int) for seed in seeds) and not isinstance(value.get("receipts"), list):
            return "plan"
    return None


def _required_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ReleaseError(f"{label} 값이 없거나 문자열이 아닙니다")
    return value


def _assessment_batch_binding(assessment: Mapping[str, Any]) -> Dict[str, str]:
    source = assessment.get("source_binding")
    if not isinstance(source, dict):
        raise ReleaseError("assessment.source_binding이 없습니다")
    assessment_job = _required_string(assessment.get("job_id"), "assessment.job_id")
    source_job = _required_string(source.get("job_id"), "assessment.source_binding.job_id")
    if assessment_job != source_job:
        raise ReleaseError("assessment job binding이 서로 일치하지 않습니다")
    project = _required_string(assessment.get("project_id"), "assessment.project_id")
    source_project = _required_string(source.get("project"), "assessment.source_binding.project")
    if project != source_project:
        raise ReleaseError("assessment project binding이 서로 일치하지 않습니다")
    return {
        "job_id": assessment_job,
        "project": project,
        "input_sha256": _required_string(source.get("input_sha256"), "assessment.source_binding.input_sha256"),
        "result_sha256": _required_string(source.get("result_sha256"), "assessment.source_binding.result_sha256"),
    }


def validate_batch_root(root: Path, assessment: Mapping[str, Any]) -> Dict[str, Any]:
    files = iter_safe_tree(root, "batch root")
    safe_files = set(files)
    binding = _assessment_batch_binding(assessment)

    plan_paths = sorted(root.glob("*.plan.json"), key=lambda path: path.name)
    receipt_paths = sorted(
        root.glob("*/seed-*/receipt.json"),
        key=lambda path: path.relative_to(root).as_posix(),
    )
    for path in plan_paths + receipt_paths:
        if path not in safe_files:
            raise ReleaseError(f"batch record가 안전한 regular file이 아닙니다: {path}")

    if len(plan_paths) != EXPECTED_BATCH_PLANS or len(receipt_paths) != EXPECTED_BATCH_RECEIPTS:
        raise ReleaseError(
            f"batch root에는 root *.plan.json {EXPECTED_BATCH_PLANS}개와 "
            f"*/seed-*/receipt.json {EXPECTED_BATCH_RECEIPTS}개가 필요합니다; "
            f"발견: plan {len(plan_paths)}, receipt {len(receipt_paths)}"
        )

    plans_by_run: Dict[str, Dict[str, Any]] = {}
    plan_rows: List[Dict[str, Any]] = []
    for path in plan_paths:
        rel = path.relative_to(root).as_posix()
        plan = load_json(path, "batch plan")
        if not isinstance(plan, dict):
            raise ReleaseError(f"plan은 객체여야 합니다: {rel}")
        if plan.get("seeds") != EXPECTED_SEEDS:
            raise ReleaseError(f"plan seed 계획이 잘못되었습니다: {rel}")

        graph = plan.get("candidate_graph")
        if not isinstance(graph, dict):
            raise ReleaseError(f"plan.candidate_graph가 없습니다: {rel}")
        candidate_id = _required_string(graph.get("candidate_id"), f"{rel}.candidate_graph.candidate_id")
        analog_id = _required_string(graph.get("warhead_analog_id"), f"{rel}.candidate_graph.warhead_analog_id")
        e3_type = _required_string(graph.get("e3_type"), f"{rel}.candidate_graph.e3_type")
        graph_sha256 = _required_string(graph.get("graph_sha256"), f"{rel}.candidate_graph.graph_sha256")
        if e3_type not in {"CRBN", "VHL"} or analog_id not in EXPECTED_ANALOGS:
            raise ReleaseError(f"plan analog/E3 mapping이 허용 범위와 다릅니다: {rel}")
        if not HEX64.fullmatch(graph_sha256):
            raise ReleaseError(f"candidate graph SHA-256 형식이 잘못되었습니다: {rel}")

        plan_bindings = plan.get("bindings")
        job_binding = plan_bindings.get("job") if isinstance(plan_bindings, dict) else None
        policy_binding = plan_bindings.get("policy") if isinstance(plan_bindings, dict) else None
        if not isinstance(job_binding, dict) or not isinstance(policy_binding, dict):
            raise ReleaseError(f"plan bindings가 없습니다: {rel}")
        expected_job_values = {
            "job_id": binding["job_id"],
            "project": binding["project"],
            "input_sha256": binding["input_sha256"],
            "result_sha256": binding["result_sha256"],
            "candidate_graph_sha256": graph_sha256,
        }
        for key, expected in expected_job_values.items():
            if job_binding.get(key) != expected:
                raise ReleaseError(f"plan graph/result binding 불일치: {rel}: bindings.job.{key}")
        policy_digest = _required_string(policy_binding.get("digest"), f"{rel}.bindings.policy.digest")
        if not HEX64.fullmatch(policy_digest):
            raise ReleaseError(f"plan policy digest 형식이 잘못되었습니다: {rel}")

        plan_digest = _required_string(plan.get("plan_digest"), f"{rel}.plan_digest")
        if not HEX64.fullmatch(plan_digest):
            raise ReleaseError(f"plan digest 형식이 잘못되었습니다: {rel}")
        run = f"{analog_id}--{e3_type}"
        if run in plans_by_run:
            raise ReleaseError(f"중복된 analog/E3 plan입니다: {run}")
        run_directory = root / run
        require_real_directory(run_directory, f"batch run {run}")
        row = {
            "run": run,
            "plan_path": rel,
            "candidate_id": candidate_id,
            "warhead_analog_id": analog_id,
            "e3_type": e3_type,
            "candidate_graph_sha256": graph_sha256,
            "plan_digest": plan_digest,
            "policy_digest": policy_digest,
            "seeds": list(EXPECTED_SEEDS),
            "receipts": [],
        }
        plans_by_run[run] = row
        plan_rows.append(row)

    expected_pairs = {(analog, e3) for analog in EXPECTED_ANALOGS for e3 in ("CRBN", "VHL")}
    observed_pairs = {(row["warhead_analog_id"], row["e3_type"]) for row in plan_rows}
    if observed_pairs != expected_pairs:
        raise ReleaseError("batch plans는 정확한 3 analog × CRBN/VHL mapping이어야 합니다")

    seen_receipts = set()
    statuses: Counter[str] = Counter()
    exit_zero = 0
    failures: List[Dict[str, Any]] = []
    for path in receipt_paths:
        rel = path.relative_to(root).as_posix()
        run = path.parent.parent.name
        seed_directory = path.parent.name
        match = re.fullmatch(r"seed-(\d+)", seed_directory)
        if run not in plans_by_run or match is None:
            raise ReleaseError(f"receipt run/seed 경로가 plan과 일치하지 않습니다: {rel}")
        path_seed = int(match.group(1))
        receipt = load_json(path, "batch receipt")
        if not isinstance(receipt, dict):
            raise ReleaseError(f"receipt는 객체여야 합니다: {rel}")
        plan_row = plans_by_run[run]
        if path_seed not in EXPECTED_SEEDS or receipt.get("seed") != path_seed:
            raise ReleaseError(f"receipt seed와 경로가 일치하지 않습니다: {rel}")
        receipt_key = (run, path_seed)
        if receipt_key in seen_receipts:
            raise ReleaseError(f"중복 receipt입니다: {run} seed {path_seed}")
        seen_receipts.add(receipt_key)

        if receipt.get("candidate_id") != plan_row["candidate_id"]:
            raise ReleaseError(f"receipt candidate binding 불일치: {rel}")
        if receipt.get("e3_type") != plan_row["e3_type"]:
            raise ReleaseError(f"receipt E3 binding 불일치: {rel}")
        if receipt.get("plan_digest") != plan_row["plan_digest"]:
            raise ReleaseError(f"receipt plan digest 불일치: {rel}")
        receipt_bindings = receipt.get("bindings")
        if not isinstance(receipt_bindings, dict):
            raise ReleaseError(f"receipt bindings가 없습니다: {rel}")
        expected_receipt_bindings = {
            "job_id": binding["job_id"],
            "project": binding["project"],
            "policy_digest": plan_row["policy_digest"],
        }
        for key, expected in expected_receipt_bindings.items():
            if receipt_bindings.get(key) != expected:
                raise ReleaseError(f"receipt binding 불일치: {rel}: bindings.{key}")
        if receipt.get("actual_computation") is not True:
            raise ReleaseError(f"receipt가 actual computation을 증명하지 않습니다: {rel}")
        status_value = _required_string(receipt.get("status"), f"{rel}.status").lower()
        exit_code = receipt.get("exit_code")
        if not isinstance(exit_code, int) or isinstance(exit_code, bool):
            raise ReleaseError(f"receipt exit_code가 정수가 아닙니다: {rel}")
        inspection = receipt.get("inspection")
        if not isinstance(inspection, dict):
            raise ReleaseError(f"receipt inspection이 없습니다: {rel}")
        if inspection.get("actual_computation") is not True or inspection.get("status") != "computed_hypothesis":
            raise ReleaseError(f"receipt inspection은 actual computed_hypothesis여야 합니다: {rel}")
        if receipt.get("scientific_status") != "computed_hypothesis":
            raise ReleaseError(f"receipt scientific_status가 computed_hypothesis가 아닙니다: {rel}")
        if receipt.get("reference_free") is not True or inspection.get("reference_free") is not True:
            raise ReleaseError(f"receipt reference-free binding이 없습니다: {rel}")

        statuses[status_value] += 1
        if exit_code == 0:
            exit_zero += 1
        if status_value != "completed" or exit_code != 0:
            failures.append({
                "path": rel,
                "run": run,
                "seed": path_seed,
                "status": status_value,
                "exit_code": exit_code,
            })
        plan_row["receipts"].append({
            "path": rel,
            "seed": path_seed,
            "status": status_value,
            "exit_code": exit_code,
            "inspection_status": inspection.get("status"),
        })

    expected_receipts = {(run, seed) for run in plans_by_run for seed in EXPECTED_SEEDS}
    if seen_receipts != expected_receipts:
        raise ReleaseError("각 plan에는 seeds 23/41/61/79/97의 receipt가 정확히 하나씩 필요합니다")
    for row in plan_rows:
        row["receipts"].sort(key=lambda item: item["seed"])

    selected = set(plan_paths) | set(receipt_paths)
    aggregate_plan_paths = sorted(
        (path for path in root.glob("*.json") if path not in selected),
        key=lambda path: path.name,
    )
    aggregate_receipt_paths = sorted(
        (path for path in root.glob("*/receipt.json") if path not in selected),
        key=lambda path: path.relative_to(root).as_posix(),
    )
    for path in aggregate_plan_paths + aggregate_receipt_paths:
        if path not in safe_files:
            raise ReleaseError(f"aggregate record가 안전한 regular file이 아닙니다: {path}")

    return {
        "plan_count": len(plan_paths),
        "receipt_count": len(receipt_paths),
        "receipt_status_counts": dict(sorted(statuses.items())),
        "exit_zero_count": exit_zero,
        "failed_or_nonzero_receipts": failures,
        "plans": [path.relative_to(root).as_posix() for path in plan_paths],
        "aggregate_records": {
            "root_json": [path.relative_to(root).as_posix() for path in aggregate_plan_paths],
            "run_receipts": [path.relative_to(root).as_posix() for path in aggregate_receipt_paths],
        },
        "validated_runs": sorted(plan_rows, key=lambda row: row["run"]),
    }


def assessment_body(value: Any) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise ReleaseError("assessment JSON은 객체여야 합니다")
    if isinstance(value.get("facts"), dict) and "criteria" in value["facts"]:
        return value["facts"]
    return value


def extract_criteria(assessment: Mapping[str, Any]) -> List[Dict[str, Any]]:
    criteria = assessment.get("criteria")
    if not isinstance(criteria, list) or len(criteria) != EXPECTED_CRITERIA:
        raise ReleaseError(f"assessment.criteria는 정확히 {EXPECTED_CRITERIA}개여야 합니다")
    seen = set()
    result: List[Dict[str, Any]] = []
    for item in criteria:
        if not isinstance(item, dict):
            raise ReleaseError("criterion 항목은 객체여야 합니다")
        cid = item.get("id")
        status_value = item.get("status")
        if not isinstance(cid, str) or not cid or cid in seen:
            raise ReleaseError(f"criterion id가 없거나 중복되었습니다: {cid!r}")
        if status_value not in ALLOWED_STATUSES:
            raise ReleaseError(f"criterion 상태가 유효하지 않습니다: {cid}: {status_value!r}")
        seen.add(cid)
        result.append(dict(item))
    return result


def locate_policy(assessment: Mapping[str, Any], optional_policy: Optional[Any]) -> Tuple[Dict[str, Any], str]:
    if optional_policy is not None:
        if not isinstance(optional_policy, dict):
            raise ReleaseError("--policy JSON은 객체여야 합니다")
        return dict(optional_policy), "--policy"
    for key in ("current_policy", "policy"):
        value = assessment.get(key)
        if isinstance(value, dict):
            return dict(value), f"assessment.{key}"
    raise ReleaseError("현재 binding policy를 assessment에서 찾지 못했습니다. --policy를 지정하십시오")


def validate_policy(policy: Mapping[str, Any]) -> Dict[str, Any]:
    numerical = policy.get("numerical_criteria")
    if not isinstance(numerical, dict):
        raise ReleaseError("현재 policy.numerical_criteria가 없습니다")
    required = {
        "parent_funnel_min": 5,
        "parent_funnel_max": 10,
        "modifiable_sites_min": 1,
        "distinct_graphs_per_site_min": 30,
        "broad_families_min": 6,
        "panel_min": 10,
        "panel_max": 20,
        "calibration_seeds_min": 5,
        "novel_ternary_repeats_min": 5,
    }
    extracted: Dict[str, Any] = {}
    for key, expected in required.items():
        value = numerical.get(key)
        if not isinstance(value, int) or isinstance(value, bool):
            raise ReleaseError(f"현재 policy 수치가 없거나 정수가 아닙니다: {key}")
        extracted[key] = value
        if value != expected:
            raise ReleaseError(f"현재 릴리스 범위의 policy 값과 다릅니다: {key}={value}, expected={expected}")
    return extracted


def find_nested(value: Any, key: str) -> Optional[Any]:
    if isinstance(value, dict):
        if key in value:
            return value[key]
        for child in value.values():
            found = find_nested(child, key)
            if found is not None:
                return found
    elif isinstance(value, list):
        for child in value:
            found = find_nested(child, key)
            if found is not None:
                return found
    return None


def summarize_novel_batches(
    assessment: Mapping[str, Any], batch_root: Path, validated_batch: Optional[Mapping[str, Any]] = None
) -> Dict[str, Any]:
    batch = dict(validated_batch) if validated_batch is not None else validate_batch_root(batch_root, assessment)
    rows = batch.get("validated_runs")
    if not isinstance(rows, list) or len(rows) != EXPECTED_BATCH_PLANS:
        raise ReleaseError("validated batch run 구조가 없습니다")

    runs = []
    completed_exit_zero = 0
    failures = []
    total_receipts = 0
    for row in rows:
        if not isinstance(row, dict):
            raise ReleaseError("validated batch run이 객체가 아닙니다")
        receipts = row.get("receipts")
        if not isinstance(receipts, list):
            raise ReleaseError(f"validated receipt 목록이 없습니다: {row.get('run')}")
        observed = [receipt.get("seed") for receipt in receipts if isinstance(receipt, dict)]
        if observed != EXPECTED_SEEDS:
            raise ReleaseError(f"run seed receipt 범위가 일치하지 않습니다: {row.get('run')}: {observed}")
        for receipt in receipts:
            total_receipts += 1
            if receipt.get("status") == "completed" and receipt.get("exit_code") == 0:
                completed_exit_zero += 1
            else:
                failures.append({
                    "run": row.get("run"),
                    "seed": receipt.get("seed"),
                    "status": receipt.get("status"),
                    "exit_code": receipt.get("exit_code"),
                    "path": receipt.get("path"),
                })
        runs.append({
            "run": row.get("run"),
            "candidate_id": row.get("candidate_id"),
            "warhead_analog_id": row.get("warhead_analog_id"),
            "e3_type": row.get("e3_type"),
            "candidate_graph_sha256": row.get("candidate_graph_sha256"),
            "plan_digest": row.get("plan_digest"),
            "seeds": observed,
        })

    if total_receipts != EXPECTED_BATCH_RECEIPTS:
        raise ReleaseError("validated batch에는 정확히 30개 seed receipt가 필요합니다")
    return {
        "run_count": len(runs),
        "receipt_count": total_receipts,
        "completed_exit_zero": completed_exit_zero,
        "failed_or_nonzero_receipts": failures,
        "runs": sorted(runs, key=lambda row: row["run"]),
        "computed_hypothesis": "technical_not_scientific_pass",
        "reference_free": True,
        "novel_input_mode": "single_sequence",
        "known_benchmark_msa_mode": "server",
        "matched_novel_reference": False,
    }


def summarize_packet(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    representatives = manifest.get("representatives")
    if not isinstance(representatives, list) or len(representatives) != 6:
        raise ReleaseError("packet manifest에는 정확히 6개 representative가 필요합니다")
    mapping = []
    seen_candidates = set()
    for row in representatives:
        if not isinstance(row, dict):
            raise ReleaseError("packet representative가 객체가 아닙니다")
        candidate_id = _required_string(row.get("candidate_id"), "packet representative.candidate_id")
        e3_type = _required_string(row.get("e3_type"), f"packet representative {candidate_id}.e3_type")
        linker_id = _required_string(row.get("linker_id"), f"packet representative {candidate_id}.linker_id")
        graph_sha256 = _required_string(
            row.get("candidate_graph_sha256"),
            f"packet representative {candidate_id}.candidate_graph_sha256",
        )
        if not HEX64.fullmatch(graph_sha256):
            raise ReleaseError(f"packet representative graph SHA-256 형식이 잘못되었습니다: {candidate_id}")

        declared_analog = row.get("warhead_analog_id")
        actual_analog = row.get("analog_id")
        if declared_analog is None:
            analog_id = actual_analog
        else:
            analog_id = declared_analog
            if actual_analog is not None and actual_analog != declared_analog:
                raise ReleaseError(f"packet representative analog 필드가 충돌합니다: {candidate_id}")
        analog_id = _required_string(analog_id, f"packet representative {candidate_id}.analog_id")
        if candidate_id in seen_candidates:
            raise ReleaseError(f"packet representative candidate_id가 중복되었습니다: {candidate_id}")
        seen_candidates.add(candidate_id)
        mapping.append({
            "candidate_id": candidate_id,
            "warhead_analog_id": analog_id,
            "e3_type": e3_type,
            "linker_id": linker_id,
            "candidate_graph_sha256": graph_sha256,
        })

    pairs = [(row["warhead_analog_id"], row["e3_type"]) for row in mapping]
    expected = {(analog, e3) for analog in EXPECTED_ANALOGS for e3 in ("CRBN", "VHL")}
    if set(pairs) != expected or len(set(pairs)) != 6:
        raise ReleaseError("packet representative의 3 analog × 2 E3 mapping이 일치하지 않습니다")
    if any(row["linker_id"] != "alkyl_c6" for row in mapping):
        raise ReleaseError("현재 representative linker는 모두 alkyl_c6이어야 합니다")
    return {
        "representative_count": 6,
        "representatives": mapping,
        "representative_linker": "alkyl_c6",
        "representative_linker_selected_by": "developer",
        "expert_approved_specific_linker": False,
    }


def benchmark_summary(benchmark: Mapping[str, Any]) -> Dict[str, Any]:
    seeds = benchmark.get("seeds")
    stats = benchmark.get("statistics")
    if not isinstance(seeds, list) or len(seeds) != 5 or not isinstance(stats, dict):
        raise ReleaseError("benchmark는 5-seed 분포와 statistics를 포함해야 합니다")
    observed = [row.get("seed") for row in seeds if isinstance(row, dict)]
    if observed != EXPECTED_SEEDS:
        raise ReleaseError(f"benchmark seed가 예상과 다릅니다: {observed}")
    e3 = stats.get("comparison.metrics.e3_CA_RMSD_after_target_alignment_A")
    ligand = stats.get("comparison.metrics.ligand_heavy_atom_RMSD_after_target_alignment_A")
    if not isinstance(e3, dict) or not isinstance(ligand, dict):
        raise ReleaseError("benchmark RMSD raw statistics가 없습니다")
    seed41 = next((row for row in seeds if isinstance(row, dict) and row.get("seed") == 41), None)
    seed23 = next((row for row in seeds if isinstance(row, dict) and row.get("seed") == 23), None)
    return {
        "source": "CRBN_6BOY",
        "seed_count": 5,
        "seeds": observed,
        "technical_success_count": benchmark.get("technical_success_count"),
        "scientific_model_quality_success_count": benchmark.get("scientific_model_quality_success_count"),
        "expert_geometric_success_cutoff": benchmark.get("expert_geometric_success_cutoff"),
        "distribution_strength": benchmark.get("distribution_strength"),
        "outlier_deletion_performed": benchmark.get("outlier_deletion_performed"),
        "poor_seed_41_retained": bool(seed41),
        "seed_41_e3_rmsd_A": seed41.get("metrics", {}).get("comparison.metrics.e3_CA_RMSD_after_target_alignment_A") if seed41 else None,
        "seed_23_original_status": seed23.get("original_status") if seed23 else None,
        "seed_23_original_failure_reason": seed23.get("original_failure_reason") if seed23 else None,
        "e3_CA_RMSD_after_target_alignment_A": {
            "median": e3.get("median"), "min": e3.get("min"), "max": e3.get("max"), "iqr": e3.get("iqr"),
            "q1": e3.get("q1"), "q3": e3.get("q3"), "count": e3.get("count"),
        },
        "ligand_heavy_atom_RMSD_after_target_alignment_A": {
            "median": ligand.get("median"), "min": ligand.get("min"), "max": ligand.get("max"), "iqr": ligand.get("iqr"),
            "q1": ligand.get("q1"), "q3": ligand.get("q3"), "count": ligand.get("count"),
        },
        "interpretation": "weak_single_known_case_descriptive_only",
        "confidence_or_geometric_cutoff_defined": False,
    }


def redact_private_paths(value: Any) -> Any:
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            if SECRET_KEY_RE.search(str(key)):
                raise ReleaseError(f"usage ledger에 금지된 secret 필드가 있습니다: {key}")
            result[key] = redact_private_paths(child)
        return result
    if isinstance(value, list):
        return [redact_private_paths(child) for child in value]
    if isinstance(value, str) and (ABS_WINDOWS_RE.match(value) or value.startswith("/home/") or value.startswith("/Users/")):
        normalized = value.replace("\\", "/")
        return "<private-path-redacted>/" + normalized.rsplit("/", 1)[-1]
    return value


def summarize_usage_totals(sanitized_usage: Mapping[str, Any]) -> Dict[str, Any]:
    developer = sanitized_usage.get("developer_calls")
    product = sanitized_usage.get("product_calls")
    developer = developer if isinstance(developer, dict) else {}
    product = product if isinstance(product, dict) else {}
    return {
        "developer_total_tokens": developer.get("total_tokens"),
        "developer_calls": developer.get("calls"),
        "product_total_tokens": product.get("total_tokens"),
        "product_calls": product.get("calls"),
        "unknown_usage_failures": sanitized_usage.get("unknown_usage_failures"),
        "latest_quota": sanitized_usage.get("latest_quota"),
        "quota_caveat": sanitized_usage.get("quota_caveat"),
    }


def validation_summary(value: Optional[Any]) -> Dict[str, Any]:
    if value is None:
        return {"provided": False, "note": "최종 수정 이후의 validation JSON이 제공되지 않아 최신 테스트 상태를 단정하지 않음"}
    if not isinstance(value, dict):
        raise ReleaseError("--validation JSON은 객체여야 합니다")
    return {"provided": True, "record": value}


def status_counts(criteria: Sequence[Mapping[str, Any]]) -> Dict[str, int]:
    counter = Counter(str(row["status"]) for row in criteria)
    return {key: counter.get(key, 0) for key in ("pass", "failed", "pending")}


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)


def display_metric(value: Any) -> str:
    return "확인 불가(unavailable)" if value is None else esc(value)


def validation_status_text(validation: Mapping[str, Any]) -> str:
    if not validation.get("provided"):
        return esc(validation.get("note", "최신 검증 상태를 단정할 수 없음"))
    record = validation.get("record")
    if not isinstance(record, dict):
        return "최신 검증 상태를 단정할 수 없음"

    full_regression = find_nested(record, "full_regression")
    numeric_counts = []
    if isinstance(full_regression, dict):
        numeric_counts = [
            f"{esc(key)}={esc(value)}"
            for key, value in sorted(full_regression.items(), key=lambda item: str(item[0]))
            if isinstance(value, int) and not isinstance(value, bool)
        ]
    regression_text = "전체 회귀 수치: " + ", ".join(numeric_counts) if numeric_counts else "전체 회귀 수치: 확인 불가"

    def status_for(*keys: str) -> str:
        for key in keys:
            value = find_nested(record, key)
            if isinstance(value, dict):
                value = value.get("status")
            if isinstance(value, str) and value:
                return esc(value)
        return "확인 불가"

    return " · ".join((
        regression_text,
        "최종 단계 QA: " + status_for("final_stage", "final_stage_qa", "final_stage_status"),
        "ZIP QA: " + status_for("zip_qa", "zip_qa_status"),
    ))


def criterion_rows_html(criteria: Sequence[Mapping[str, Any]]) -> str:
    rows = []
    for row in criteria:
        status_value = row["status"]
        rows.append(
            "<tr><td><code>" + esc(row["id"]) + "</code></td><td>" + esc(row.get("title") or row.get("name") or row["id"]) +
            "</td><td><span class='badge " + esc(status_value) + "'>" + esc(status_value.upper()) +
            "</span></td><td>" + esc(compact_json(row.get("required"))) + "</td><td>" + esc(row.get("reason", "")) + "</td></tr>"
        )
    return "".join(rows)


def criterion_rows_md(criteria: Sequence[Mapping[str, Any]]) -> str:
    rows = ["| ID | 기준 | 상태 | 요구사항 |", "|---|---|---:|---|"]
    for row in criteria:
        title = str(row.get("title") or row.get("name") or row["id"]).replace("|", "\\|")
        required = compact_json(row.get("required")).replace("|", "\\|")
        rows.append(f"| `{row['id']}` | {title} | **{str(row['status']).upper()}** | {required} |")
    return "\n".join(rows)


def answers_block() -> str:
    return """
<ul>
<li><strong>보호 잔기·ligand map·상태:</strong> FX5/6HAZ 부모 범위에서 현재 허용된 실제 MODIFIABLE 위치는 atom map 19 하나이며 maps 8/11/12는 UNKNOWN으로 유지한다. PROTECTED map 3과 UNKNOWN map 9는 별도의 9D12 preview에만 속한다. 두 부모 범위를 혼합하지 않으며, UNKNOWN은 조립 승인이나 변형 위치로 계산하지 않는다.</li>
<li><strong>상태/pH:</strong> 부모 charged 8개와 neutral 8개는 pH 7.4 실행 문맥에서의 제한된 tautomer/charge 열거이다. pKa, 상태 population 또는 확인된 microstate를 뜻하지 않는다. raw pre-H 파일을 보존했고, opt 파일은 수소 최적화 중 heavy atom을 고정했다.</li>
<li><strong>단백질 수소:</strong> PDB2PQR 기반 수소 공급 이력은 있으나 사전 최적화되지 않았고 terminal OXT는 pending이다. 이는 ligand 수소 처리와 별도 문제다.</li>
<li><strong>known CRBN:</strong> 5-seed 실행은 기술적으로 완료됐지만 seed 41의 불량 거동을 포함한 약한 단일 known-case 기술 분포이다. E3 RMSD 중앙값 10.6210244 Å, 범위 3.2208359–39.1967498 Å, IQR 22.58651 Å; ligand RMSD 중앙값 2.5053 Å, 범위 1.1781–12.3878 Å, IQR 5.88687 Å이다. 신뢰도·기하 cutoff와 과학적 성공 개수는 없다.</li>
<li><strong>신규 ternary:</strong> 3 analog × CRBN/VHL × seeds 23/41/61/79/97의 30회 계산은 reference-free 기술 가설이다. 신규 계산은 single_sequence 입력이며 known benchmark의 server MSA와 matched novel reference가 아니다.</li>
<li><strong>후보:</strong> alkyl_c6의 6개 대표는 개발자 선택이며 전문가가 승인한 특정 linker가 아니다. 같은 analog 안에서만 CRBN/VHL 비교가 가능하다.</li>
<li><strong>합성:</strong> 6개 mapped-cut SVG는 합성 graph 제안이며 실험 recipe 또는 합성 가능성 승인 자료가 아니다. 근거 없는 amide/click 단계는 만들지 않았고 ester-only 항목은 precursor/protection 역할로만 취급한다.</li>
<li><strong>권한:</strong> source-verified relay는 서명 인증이 아니다. 최신 원문 SHA-256은 42a78a5540f5549347bc8ee2af1fd9758cec43056c6e808f89ff9aedc6d521cb이고, design result SHA-256은 51c85df73558b5ce87c3999c96fca937d4bab272e8e4b5fb36509bb3f2742369, job은 job-4ce0b3d57ecb4be5b795504ca22762a0이다.</li>
</ul>
"""


def reviewer_questions_html() -> str:
    return """
<ol>
<li>보호해야 할 정확한 residue/interaction과 세 analog의 ligand atom map을 지정하고, pH 7.4에서 선택할 계산 상태를 각 analog별로 확인해 주십시오.</li>
<li>CRBN known-case 분포와 poor seed 41을 포함해 어떤 유한 기하 threshold와 처리 규칙을 적용할지 지정해 주십시오.</li>
<li>CRBN과 VHL 각각에서 선택할 정확한 candidate ID 하나와 seed 간 convergence metric 및 수용 기준을 지정해 주십시오.</li>
<li>6개 graph 각각에 대해 합성 역할, 절단 bond, precursor 호환성, 보호기 요구사항을 확인하고 route proposal인지 현재 비합성 가능인지 구분해 주십시오.</li>
<li>FX5 site-count 1 waiver만으로 충분한지 확인해 주십시오. 이는 parent 5–10, broad family 6, panel 10–20 요구를 면제하지 않습니다.</li>
<li>parent 추가 선택과 새로운 SAR site 사용을 허가할지 답해 주십시오. 허가되지 않으면 실제 dock-qualified panel을 만들 수 없습니다.</li>
<li>고정 FX5 demo 범위의 대체 수치 면제가 허용되는지 명시해 주십시오. 실제 답변 전에는 어떤 대체 면제도 적용하지 않습니다.</li>
<li>모든 gate가 통과한 뒤에만 직접 policy final review와 인증된 formal decision을 제출해 주십시오.</li>
</ol>
""".replace("ligand map을", "ligand atom map을").replace("ligand map을", "ligand atom map을").replace("lid?? ", "")


def css() -> str:
    return """
:root{--bg:#f4f7fb;--ink:#172033;--muted:#5b667a;--card:#fff;--line:#dce3ee;--blue:#2358c6;--pass:#087a55;--fail:#b42318;--pending:#9a6700}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--bg);color:var(--ink);font-family:system-ui,-apple-system,"Segoe UI","Noto Sans KR",Arial,sans-serif;line-height:1.62;overflow-wrap:anywhere}header{background:linear-gradient(135deg,#13284c,#2358c6);color:#fff;padding:48px 20px}main,footer{width:min(1120px,calc(100% - 28px));margin:auto}.hero{width:min(1120px,100%);margin:auto}.eyebrow{letter-spacing:.12em;text-transform:uppercase;opacity:.8;font-size:.78rem}.grid{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:20px}.metric,.card{background:var(--card);border:1px solid var(--line);border-radius:16px;padding:18px;box-shadow:0 6px 22px rgba(28,46,82,.06)}.metric{color:var(--ink)}.metric b{font-size:1.7rem;display:block}.card{margin:18px 0}h1{font-size:clamp(2rem,5vw,3.5rem);line-height:1.12;margin:.3em 0}h2{margin-top:0;font-size:1.4rem}h3{margin-bottom:.35rem}a{color:var(--blue)}code{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.88em}.table-wrap{overflow-x:auto;border:1px solid var(--line);border-radius:12px}table{border-collapse:collapse;width:100%;min-width:780px;background:#fff}th,td{text-align:left;vertical-align:top;padding:11px;border-bottom:1px solid var(--line)}th{background:#edf2fa}.badge{display:inline-block;border-radius:999px;padding:2px 9px;font-weight:700;font-size:.78rem}.pass{color:var(--pass);background:#dcfaeb}.failed{color:var(--fail);background:#fee4e2}.pending{color:var(--pending);background:#fff1c2}.warning{border-left:5px solid var(--fail);background:#fff5f4;padding:14px}.note{border-left:5px solid var(--blue);background:#eef4ff;padding:14px}.links{display:flex;gap:12px;flex-wrap:wrap}.links a{display:inline-block;background:#fff;border:1px solid var(--line);padding:9px 13px;border-radius:10px;text-decoration:none}footer{padding:28px 0 50px;color:var(--muted);font-size:.9rem}@media(max-width:720px){header{padding:34px 16px}.grid{grid-template-columns:1fr}.card{padding:14px}main,footer{width:min(100% - 18px,1120px)}table{min-width:680px}}@media print{body{background:#fff}.card,.metric{box-shadow:none}a{color:inherit}}
"""


def render_report_html(summary: Mapping[str, Any], criteria: Sequence[Mapping[str, Any]]) -> str:
    counts = summary["status_counts"]
    validation = summary["validation"]
    validation_text = validation_status_text(validation)
    return f"""<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>전문가 릴리스 현황 보고서</title><style>{css()}</style></head><body>
<header><div class="hero"><div class="eyebrow">SOURCE-BOUND EXPERT RELEASE</div><h1>전문가 릴리스 현황 보고서</h1><p>기술 실행 결과와 현재 정책 gate를 구분한 검토용 보고서입니다. 과학적 승인, 효능, 실험 성공 또는 제출 승인을 주장하지 않습니다.</p><div class="grid"><div class="metric"><b>{counts['pass']}</b>PASS</div><div class="metric"><b>{counts['failed']}</b>FAILED</div><div class="metric"><b>{counts['pending']}</b>PENDING</div></div></div></header>
<main>
<section class="card"><h2>CURRENT STATUS</h2><div class="warning"><strong>과학적으로 승인되지 않음.</strong> 확인 상태는 {esc(summary['confirmation_status'])}이며 scientific_accepted=false입니다. 기술 계산 완료는 과학적 geometry 승인이나 model success가 아닙니다.</div><p>현재 14개 gate는 원본 assessment에서 직접 계산한 결과로 PASS 3 / FAILED 3 / PENDING 8입니다. strict parent funnel은 0이고, 인증된 전문가 parent 선택은 없습니다. catalog의 10개 항목이나 이전 12개 exploratory panel은 strict count 또는 14-gate pass로 이전되지 않습니다.</p></section>
<section class="card"><h2>바로가기</h2><div class="links"><a href="expert-packet/index.html">전문가 evidence packet</a><a href="expert-followup.html">재사용 가능한 reviewer brief</a><a href="measured-summary.json">측정 요약 JSON</a><a href="api-usage.json">API usage ledger</a></div></section>
<section class="card"><h2>14개 정책 gate</h2><div class="table-wrap"><table><thead><tr><th>ID</th><th>기준</th><th>상태</th><th>요구사항</th><th>평가 근거</th></tr></thead><tbody>{criterion_rows_html(criteria)}</tbody></table></div></section>
<section class="card"><h2>실제 계산 범위</h2><ul><li>실제 technical CPU/docking 범위: distinct graph 38, 실제 broad family 1, qualified selected 3, assembly 36개(CRBN 18, VHL 18).</li><li>신규 ternary: 6 plans, 30 receipts, seed 23/41/61/79/97. 현재 source assessment의 completed exit 0은 {summary['novel']['completed_exit_zero']}/30이다.</li><li>3 analog의 exact mapping은 CRBN/VHL 동일 analog 내부 비교만 허용한다. 후보 간 또는 E3 간 우월성 결론은 없다.</li><li>known benchmark는 server MSA이고 신규는 single_sequence이며, matched novel reference가 아니다.</li></ul></section>
<section class="card"><h2>현재 자료가 답하는 내용</h2>{answers_block()}</section>
<section class="card"><h2>전문가에게 필요한 정확한 답변</h2>{reviewer_questions_html()}</section>
<section class="card"><h2>범위와 면제의 한계</h2><p>현재 reply는 고정 FX5 범위의 site-count 1만 면제한다. parent 5–10, 실제 broad family 6, qualified panel 10–20은 면제하지 않는다. 실제 전문가 답변 없이 대체 수치 면제를 적용하지 않는다. 추가 parent 선택과 새로운 SAR site 허가 후 실제 dock-qualified panel을 다시 만들어야 한다.</p></section>
<section class="card"><h2>제품·검증 상태</h2><p>기존 A/v2/v3는 보존 대상이다. {validation_text}</p><p>제공되지 않은 테스트 수나 상태는 추정하지 않는다. 공개 demo, 배포, presentation video 완료를 주장하지 않는다.</p></section>
<section class="card"><h2>실행 안내</h2><ol><li>01 launcher는 로컬 design/review 서버를 시작한다.</li><li>04 launcher는 저장된 evidence의 hash, source binding 및 보존된 import를 검증하며 assessment를 재계산하지 않는다.</li><li>05 launcher는 기존 scientific review 페이지를 열며 plan을 직접 검증하지 않는다.</li><li>06 launcher는 저장된 evidence packet을 열며 새 보고서를 생성하지 않는다.</li></ol><p>최종 보고서 링크는 같은 final-reports 폴더의 상대 경로 <code>report.html</code>이다.</p></section>
<section class="card"><h2>API usage 해석</h2><p>developer: total_tokens <strong>{display_metric(summary['usage_totals'].get('developer_total_tokens'))}</strong>, calls <strong>{display_metric(summary['usage_totals'].get('developer_calls'))}</strong> · product: total_tokens <strong>{display_metric(summary['usage_totals'].get('product_total_tokens'))}</strong>, calls <strong>{display_metric(summary['usage_totals'].get('product_calls'))}</strong> · unknown_usage_failures <strong>{display_metric(summary['usage_totals'].get('unknown_usage_failures'))}</strong></p><p>remaining quota 값은 비단조적인 관측 추정치이며 사용자 잔액으로 해석하지 않는다. 절대 경로의 개인 계정 부분은 출력본에서 제거했다.</p></section>
<section class="card"><h2>전달 무결성</h2><p>expert packet은 manifest의 모든 파일에 대해 크기와 SHA-256을 확인했고, 추가 unmanifested 파일을 거부했다. 이 보고서는 최종 ZIP의 크기나 SHA-256을 미리 주장하지 않는다. 실제 archive를 만든 외부 delivery receipt가 최종 checksum을 제공해야 한다.</p></section>
</main><footer>Assessment ID: {esc(summary['assessment_id'])} · Revision: {esc(summary['assessment_revision'])} · Policy source: {esc(summary['policy_source'])}</footer></body></html>"""


def render_followup_html(summary: Mapping[str, Any], criteria: Sequence[Mapping[str, Any]]) -> str:
    counts = summary["status_counts"]
    return f"""<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>전문가 후속 검토 brief</title><style>{css()}</style></head><body><header><div class="hero"><div class="eyebrow">REUSABLE REVIEWER BRIEF</div><h1>전문가 후속 검토 요청</h1><p>현재 상태와 실제로 필요한 답변을 한 문서에 정리했습니다.</p></div></header><main>
<section class="card"><h2>CURRENT STATUS</h2><div class="grid"><div class="metric"><b>{counts['pass']}</b>PASS</div><div class="metric"><b>{counts['failed']}</b>FAILED</div><div class="metric"><b>{counts['pending']}</b>PENDING</div></div><div class="warning"><strong>현재 과학적 승인 없음.</strong> technical computation은 hypothesis이며 scientific pass가 아닙니다.</div></section>
<section class="card"><h2>확인 가능한 구체적 답</h2>{answers_block()}</section>
<section class="card"><h2>검토자가 답해야 할 질문</h2>{reviewer_questions_html()}</section>
<section class="card"><h2>14개 checklist</h2><div class="table-wrap"><table><thead><tr><th>ID</th><th>기준</th><th>상태</th><th>요구사항</th><th>근거</th></tr></thead><tbody>{criterion_rows_html(criteria)}</tbody></table></div></section>
<section class="card"><h2>최종 결정 조건</h2><p>모든 failed/pending gate가 실제 근거로 해결된 뒤 직접 policy final review를 수행해야 합니다. formal accept도 다른 failed, blocked 또는 pending 기준을 덮어쓸 수 없습니다.</p></section>
</main><footer><a href="report.html">전체 보고서로 돌아가기</a></footer></body></html>"""


def render_markdown(summary: Mapping[str, Any], criteria: Sequence[Mapping[str, Any]]) -> str:
    counts = summary["status_counts"]
    return f"""# 전문가 릴리스 현황 보고서

## CURRENT STATUS

- PASS: **{counts['pass']}**
- FAILED: **{counts['failed']}**
- PENDING: **{counts['pending']}**
- 확인 상태: `{summary['confirmation_status']}`
- 과학적 승인: **아니오**
- 계산 해석: `technical_not_scientific_pass`

기술 계산 완료는 과학적 geometry 승인, 효능, 분해 효과, 실험 성공 또는 제출 승인이 아닙니다.

## 링크

- [전문가 evidence packet](expert-packet/index.html)
- [전문가 후속 검토 brief](expert-followup.html)
- [측정 요약 JSON](measured-summary.json)
- [API usage ledger](api-usage.json)

## 14개 정책 gate

{criterion_rows_md(criteria)}

## 실제 측정 범위

- strict parent funnel: 0; catalog 10개는 선택으로 계산하지 않음
- 실제 distinct graph: 38
- 실제 broad family: 1
- qualified selected: 3
- assembly: 36개, CRBN 18 / VHL 18
- 신규 ternary: 3 analog × CRBN/VHL × seeds 23/41/61/79/97 = 30회
- assessment의 completed exit 0: {summary['novel']['completed_exit_zero']}/30
- 신규 입력은 `single_sequence`; known benchmark는 server MSA이며 matched novel reference가 아님

## Known CRBN 기술 분포

- E3 RMSD 중앙값 10.6210244 Å, 범위 3.2208359–39.1967498 Å, IQR 22.58651 Å
- ligand RMSD 중앙값 2.5053 Å, 범위 1.1781–12.3878 Å, IQR 5.88687 Å
- poor seed 41과 seed 23의 원래 reader failure를 보존
- 약한 단일 known-case descriptive 분포
- confidence/geometric cutoff 및 scientific model success count 없음

## 상태·수소 처리

charged 8개와 neutral 8개 상태는 pH 7.4 실행 문맥의 제한된 tautomer/charge 열거입니다. pKa, population 또는 확인된 microstate가 아닙니다. raw pre-H를 보존했고 opt는 heavy atom 고정 수소 최적화 결과입니다. PDB2PQR 단백질 수소는 사전 최적화되지 않았고 terminal OXT는 pending이며 ligand 수소와 별도입니다.

## 합성 자료의 한계

6개 mapped-cut SVG는 합성 graph proposal이지 실험 recipe가 아닙니다. alkyl_c6 대표는 개발자 선택이며 전문가가 승인한 특정 linker가 아닙니다. 근거 없는 amide/click 단계는 만들지 않았고 ester-only 항목은 precursor/protection으로만 취급합니다.

## 필요한 전문가 답변

1. 보호 residue/interaction, ligand atom maps, 각 analog의 pH 7.4 상태 선택
2. CRBN 분포와 poor seed 41 처리에 적용할 유한 geometry threshold
3. E3별 정확한 candidate ID, convergence metric과 수용 기준
4. 6 graph별 합성 역할, 절단 bond, precursor 호환성, route proposal/비합성 가능 구분
5. parent 추가 선택 및 새로운 SAR site 사용 허가 여부
6. 고정 FX5 범위의 대체 수치 면제 허용 여부. 실제 답변 전에는 적용하지 않음
7. 모든 gate 통과 뒤 직접 policy final review와 인증된 formal decision

현재 reply는 FX5 site-count 1만 면제하며 parent 5–10, family 6, panel 10–20을 면제하지 않습니다.

## 제품·검증 상태

기존 A/v2/v3는 보존 대상입니다. {validation_status_text(summary['validation'])}. 제공되지 않은 테스트 수나 상태는 추정하지 않습니다.

## 실행 안내

1. 01 launcher는 로컬 design/review 서버를 시작합니다.
2. 04 launcher는 저장된 evidence의 hash, source binding 및 보존된 import를 검증하며 assessment를 재계산하지 않습니다.
3. 05 launcher는 기존 scientific review 페이지를 열며 plan을 직접 검증하지 않습니다.
4. 06 launcher는 저장된 evidence packet을 열며 새 보고서를 생성하지 않습니다.

최종 보고서 링크는 같은 `final-reports` 폴더의 상대 경로 `report.html`입니다.

## API usage

- developer total_tokens: **{display_metric(summary['usage_totals'].get('developer_total_tokens'))}**
- developer calls: **{display_metric(summary['usage_totals'].get('developer_calls'))}**
- product total_tokens: **{display_metric(summary['usage_totals'].get('product_total_tokens'))}**
- product calls: **{display_metric(summary['usage_totals'].get('product_calls'))}**
- unknown_usage_failures: **{display_metric(summary['usage_totals'].get('unknown_usage_failures'))}**

remaining quota 값은 비단조적인 관측 추정치이며 사용자 잔액으로 해석하지 않습니다.

## 출처 binding

- original opinion SHA-256: `42a78a5540f5549347bc8ee2af1fd9758cec43056c6e808f89ff9aedc6d521cb`
- design result SHA-256: `51c85df73558b5ce87c3999c96fca937d4bab272e8e4b5fb36509bb3f2742369`
- job: `job-4ce0b3d57ecb4be5b795504ca22762a0`
- source-verified relay는 서명 인증이 아님

## 전달

최종 ZIP 크기와 SHA-256은 이 보고서에서 미리 주장하지 않습니다. 실제 archive를 만든 외부 delivery receipt가 checksum을 제공해야 합니다.
"""


def build_summary(
    assessment: Mapping[str, Any], benchmark: Mapping[str, Any], packet_manifest: Mapping[str, Any],
    packet_files: Sequence[str], batch_root: Path, usage: Mapping[str, Any],
    policy: Mapping[str, Any], policy_source: str, validation: Optional[Any]
) -> Tuple[Dict[str, Any], List[Dict[str, Any]], Dict[str, Any]]:
    criteria = extract_criteria(assessment)
    counts = status_counts(criteria)
    declared = assessment.get("counts")
    if isinstance(declared, dict):
        for key in ("pass", "failed", "pending"):
            if declared.get(key) != counts[key]:
                raise ReleaseError(f"assessment.counts와 criteria-derived count가 다릅니다: {key}")
    if counts != {"pass": 3, "failed": 3, "pending": 8}:
        raise ReleaseError(f"현재 source-bound 상태는 pass=3, failed=3, pending=8이어야 합니다: {counts}")
    numerical = validate_policy(policy)
    batch = validate_batch_root(batch_root, assessment)
    novel = summarize_novel_batches(assessment, batch_root, batch)
    packet_summary = summarize_packet(packet_manifest)
    bench = benchmark_summary(benchmark)
    source_status = assessment.get("source_status") if isinstance(assessment.get("source_status"), dict) else {}
    scientific = assessment.get("scientific_accepted", source_status.get("scientific_accepted", False))
    if scientific is not False:
        raise ReleaseError("이 보고서는 scientific_accepted=false 평가에만 사용할 수 있습니다")
    sanitized_usage = redact_private_paths(usage)
    summary = {
        "format": "expert-release-measured-summary/1",
        "assessment_id": assessment.get("assessment_id", assessment.get("id")),
        "assessment_revision": assessment.get("revision"),
        "confirmation_status": assessment.get("confirmation_status", "unconfirmed"),
        "scientific_accepted": False,
        "source_status": source_status,
        "status_counts": counts,
        "criteria": criteria,
        "policy_source": policy_source,
        "policy_numerical_criteria": numerical,
        "measured_design_counts": {
            "strict_parent_funnel": 0,
            "authenticated_selected_parent_count": 0,
            "catalog_parent_count_not_selection": 10,
            "exploratory_panel_not_strict": 12,
            "distinct_constitutional_graphs": 38,
            "actual_broad_families": 1,
            "qualified_selected": 3,
            "assembly_total": 36,
            "assembly_CRBN": 18,
            "assembly_VHL": 18,
        },
        "novel": novel,
        "batch_root": batch,
        "benchmark": bench,
        "packet": packet_summary,
        "packet_manifest_file_count": len(packet_files),
        "site_states": {"3": "PROTECTED", "8": "UNKNOWN", "9": "UNKNOWN", "11": "UNKNOWN", "12": "UNKNOWN", "19": "MODIFIABLE"},
        "microstate_context": {
            "pH": 7.4,
            "charged_states": 8,
            "neutral_states": 8,
            "bounded_tautomer_charge_enumeration": True,
            "pKa_or_population_claim": False,
            "confirmed_microstate_claim": False,
            "raw_pre_h_retained": True,
            "hydrogen_optimized_heavy_atoms_fixed": True,
            "protein_h_method": "PDB2PQR",
            "protein_h_prior_optimized": False,
            "terminal_OXT": "pending",
        },
        "validation": validation_summary(validation),
        "usage_totals": summarize_usage_totals(sanitized_usage),
        "claims": {
            "scientific_pass": False,
            "efficacy": False,
            "experiment_success": False,
            "formal_approval": False,
            "public_demo_completed": False,
            "deployment_completed": False,
            "presentation_video_completed": False,
            "final_archive_checksum_available": False,
        },
    }
    return summary, criteria, sanitized_usage


def copy_verified_packet(source: Path, destination: Path, files: Sequence[str]) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    for rel in files:
        source_file = source.joinpath(*PurePosixPath(rel).parts)
        target = destination.joinpath(*PurePosixPath(rel).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_file, target)


def generate(args: argparse.Namespace) -> None:
    assessment_raw = load_json(Path(args.assessment), "assessment")
    assessment = assessment_body(assessment_raw)
    benchmark = load_json(Path(args.benchmark), "benchmark")
    usage = load_json(Path(args.usage), "usage ledger")
    if not isinstance(benchmark, dict) or not isinstance(usage, dict):
        raise ReleaseError("benchmark와 usage ledger는 JSON 객체여야 합니다")
    optional_policy = load_json(Path(args.policy), "policy") if args.policy else None
    validation = load_json(Path(args.validation), "validation") if args.validation else None
    policy, policy_source = locate_policy(assessment, optional_policy)
    packet = Path(args.packet)
    packet_manifest, packet_files = validate_packet(packet)
    batch_root = Path(args.batch_root)
    summary, criteria, sanitized_usage = build_summary(
        assessment, benchmark, packet_manifest, packet_files, batch_root, usage,
        policy, policy_source, validation,
    )


    output = Path(args.output)
    if output.exists():
        require_real_directory(output, "output")
        if any(output.iterdir()):
            raise ReleaseError("output 디렉터리는 존재하지 않거나 비어 있어야 합니다")
        output.rmdir()
    parent = output.parent
    parent.mkdir(parents=True, exist_ok=True)
    require_real_directory(parent, "output parent")
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp-", dir=str(parent)))
    try:
        copy_verified_packet(packet, temp / "expert-packet", packet_files)
        dump_json(temp / "measured-summary.json", summary)
        dump_json(temp / "api-usage.json", sanitized_usage)
        (temp / "report.html").write_text(render_report_html(summary, criteria), encoding="utf-8", newline="\n")
        (temp / "expert-followup.html").write_text(render_followup_html(summary, criteria), encoding="utf-8", newline="\n")
        (temp / "report.md").write_text(render_markdown(summary, criteria), encoding="utf-8", newline="\n")
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="source-bound 전문가 릴리스 보고서 생성")
    parser.add_argument("--assessment", required=True, help="현재 assessment JSON")
    parser.add_argument("--benchmark", required=True, help="5-seed known CRBN benchmark JSON")
    parser.add_argument("--packet", required=True, help="manifest 검증 대상 expert packet 디렉터리")
    parser.add_argument("--batch-root", required=True, help="6 plans + 30 receipts 디렉터리")
    parser.add_argument("--usage", required=True, help="API usage ledger JSON")
    parser.add_argument("--output", required=True, help="새 보고서 출력 디렉터리")
    parser.add_argument("--policy", help="assessment 외부의 현재 policy JSON")
    parser.add_argument("--validation", help="마지막 수정 이후 validation 결과 JSON")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        generate(args)
    except ReleaseError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

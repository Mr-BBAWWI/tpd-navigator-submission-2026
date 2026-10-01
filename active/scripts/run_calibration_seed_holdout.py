"""Run preregistered same-6BOY seed-holdout calibration validation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science import boltz_worker as worker

SEEDS = [101, 127, 149, 173, 197]
PRIOR_SEEDS = [23, 41, 61, 79, 97]
MODEL_INDICES = list(range(5))
SOURCE_SETTINGS = {
    "accelerator": "gpu", "recycling_steps": 3, "sampling_steps": 200,
    "diffusion_samples": 1, "max_parallel_samples": 1, "num_workers": 0,
    "preprocessing_threads": 1, "no_kernels": True, "use_potentials": False,
    "msa_mode": "server", "max_msa_seqs": 256,
}
PROTOCOL_SETTINGS = {
    "accelerator": "gpu", "recycling_steps": 10, "sampling_steps": 200,
    "diffusion_samples": 5, "max_parallel_samples": 1, "num_workers": 0,
    "preprocessing_threads": 1, "no_kernels": True, "use_potentials": False,
    "msa_mode": "cached_processed_offline", "max_msa_seqs": 256,
}
CUTOFFS = {
    "target_CA_RMSD_A": 1.5,
    "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
    "e3_CA_RMSD_after_target_alignment_A": 10.0,
    "contact_jaccard": 0.4,
}
LABEL = "holdout_seed_validation_same_6BOY_development_structure"


def sha256_file(path: Path) -> str:
    return worker.sha256_file(Path(path))


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True,
                         separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def module_hash() -> str:
    return sha256_file(Path(worker.__file__))


def runner_hash() -> str:
    return sha256_file(Path(__file__).resolve())


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")


def _finite_number(value: Any, error: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise worker.BoltzWorkerError(error)
    result = float(value)
    if not math.isfinite(result):
        raise worker.BoltzWorkerError(error)
    return result


def complex_iplddt(value: dict) -> float:
    """Return only the raw interface-local confidence; never fall back to iplddt."""
    if not isinstance(value, dict) or "complex_iplddt" not in value:
        raise worker.BoltzWorkerError("RAW_COMPLEX_IPLDDT_REQUIRED")
    return _finite_number(value["complex_iplddt"], "FINITE_RAW_COMPLEX_IPLDDT_REQUIRED")


def confidence_score(value: dict) -> float:
    if not isinstance(value, dict) or "confidence_score" not in value:
        raise worker.BoltzWorkerError("RAW_CONFIDENCE_SCORE_REQUIRED")
    return _finite_number(value["confidence_score"], "FINITE_RAW_CONFIDENCE_SCORE_REQUIRED")


def _indexed_candidates(models: list[dict], score_function) -> list[tuple[float, int]]:
    candidates: list[tuple[float, int]] = []
    for item in models:
        if not isinstance(item, dict) or "model_index" not in item:
            raise worker.BoltzWorkerError("MODEL_INDEX_REQUIRED")
        index = item["model_index"]
        if isinstance(index, bool) or not isinstance(index, int):
            raise worker.BoltzWorkerError("MODEL_INDEX_INVALID")
        raw = item.get("confidence_raw")
        candidates.append((score_function(raw), index))
    if len(candidates) != 5 or {index for _, index in candidates} != set(MODEL_INDICES):
        raise worker.BoltzWorkerError("ALL_FIVE_UNIQUE_MODELS_REQUIRED")
    return candidates


def select_model(models: list[dict]) -> int:
    """Highest finite raw complex_iplddt, with lower model index as the sole tie break."""
    candidates = _indexed_candidates(models, complex_iplddt)
    return min(candidates, key=lambda item: (-item[0], item[1]))[1]


def select_confidence_baseline(models: list[dict]) -> int:
    candidates = _indexed_candidates(models, confidence_score)
    return min(candidates, key=lambda item: (-item[0], item[1]))[1]


def metrics_pass(metrics: dict) -> bool:
    values: dict[str, float] = {}
    for key in CUTOFFS:
        try:
            value = metrics[key]
        except (KeyError, TypeError):
            return False
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return False
        number = float(value)
        if not math.isfinite(number) or number < 0:
            return False
        values[key] = number
    if values["contact_jaccard"] > 1:
        return False
    return (
        values["target_CA_RMSD_A"] <= 1.5 and
        values["ligand_heavy_atom_RMSD_after_target_alignment_A"] <= 3.0 and
        values["e3_CA_RMSD_after_target_alignment_A"] <= 10.0 and
        values["contact_jaccard"] >= 0.4
    )


def protocol_pass(seed_receipts: list[dict]) -> bool:
    observed = [item.get("seed") for item in seed_receipts]
    return (
        len(seed_receipts) == 5 and set(observed) == set(SEEDS) and
        len(set(observed)) == 5 and
        all(item.get("status") == "success" for item in seed_receipts) and
        sum(item.get("selected_model_pass") is True for item in seed_receipts) >= 3
    )


def _replace_option(command: list[str], option: str, value: str) -> None:
    try:
        index = command.index(option)
    except ValueError as error:
        raise worker.BoltzWorkerError("OFFICIAL_COMMAND_OPTION_MISSING:" + option) from error
    if index + 1 >= len(command):
        raise worker.BoltzWorkerError("OFFICIAL_COMMAND_OPTION_VALUE_MISSING:" + option)
    command[index + 1] = value


def effective_command(base: list[str], help_text: str) -> list[str]:
    command = list(base)
    _replace_option(command, "--recycling_steps", "10")
    _replace_option(command, "--sampling_steps", "200")
    _replace_option(command, "--diffusion_samples", "5")
    _replace_option(command, "--max_parallel_samples", "1")
    _replace_option(command, "--max_msa_seqs", "256")
    if "--use_msa_server" in command:
        command.remove("--use_msa_server")
    if re.search(r"(?<![\w-])--model(?:[\s=,]|$)", help_text) and "--model" not in command:
        command.extend(["--model", "boltz2"])
    lowered = " ".join(command).lower()
    if "potential" in lowered or "template" in lowered or "constraint" in lowered:
        raise worker.BoltzWorkerError("POTENTIAL_TEMPLATE_CONSTRAINT_FORBIDDEN")
    return command


def _model_matches(out_dir: Path, model_index: int, suffix: str,
                   confidence: bool = False) -> list[Path]:
    token = re.compile(rf"(?:^|[_-])model[_-]?{model_index}(?:[_\-.]|$)", re.I)
    return sorted(
        path for path in out_dir.rglob("*" + suffix)
        if path.is_file() and not path.is_symlink() and token.search(path.name)
        and (not confidence or "confidence" in path.name.lower())
    )


def _source_receipt_path(directory: Path) -> Path:
    preferred = directory / "receipt.json"
    if preferred.is_file() and not preferred.is_symlink():
        return preferred
    matches = sorted(
        path for path in directory.rglob("*.json")
        if path.is_file() and not path.is_symlink() and "receipt" in path.name.lower()
        and "selection" not in path.name.lower()
    )
    if len(matches) != 1:
        raise worker.BoltzWorkerError("FROZEN_SOURCE_RECEIPT_MISSING_OR_AMBIGUOUS")
    return matches[0]


def _processed_msa_manifest(directory: Path) -> list[dict]:
    entries = []
    for path in sorted(directory.rglob("*")):
        relative = path.relative_to(directory).as_posix()
        lowered = relative.lower()
        if ("msa" in lowered or "processed" in lowered) and path.is_file() and not path.is_symlink():
            entries.append({"relative_path": relative, "sha256": sha256_file(path),
                            "size_bytes": path.stat().st_size})
    if not entries:
        raise worker.BoltzWorkerError("FROZEN_PROCESSED_MSA_ENTRIES_REQUIRED")
    return entries


def _file_record(path: Path) -> dict:
    return {"path": str(path), "sha256": sha256_file(path)}


def _plan(args: argparse.Namespace, paths: dict[str, Path], help_info: dict,
          environment: dict) -> dict:
    source_receipt = _source_receipt_path(paths["frozen_seed_directory"])
    raw_units = [
        {"seed": seed, "model_index": index, "unit_id": f"seed-{seed}/model-{index}"}
        for seed in SEEDS for index in MODEL_INDICES
    ]
    return {
        "format": "tpd-crbn-seed-holdout-preregistration/1.0",
        "label": LABEL,
        "preregistered_before_inference": True,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "scope": {
            "holdout_seed_validation": True,
            "same_6BOY_development_structure": True,
            "independent_target_or_structure": False,
            "training_holdout": False,
            "known_structure_training_overlap_possible": True,
        },
        "seeds": SEEDS,
        "prior_seed_set": PRIOR_SEEDS,
        "heldout_seed_sets_disjoint": set(SEEDS).isdisjoint(PRIOR_SEEDS),
        "raw_unit_count": 25,
        "all_preregistered_raw_units": raw_units,
        "execution": "serial all five seeds and all five models; no seed cherry-picking; no early stop",
        "raw_receipts": "retain all 25 model records and all five seed receipts, including failures",
        "selection_rule": {
            "field": "raw confidence JSON top-level complex_iplddt",
            "rule": "highest finite complex_iplddt; lower model_index on exact tie",
            "fallback": "none; validated confidence_score or iplddt cannot replace missing complex_iplddt",
            "scientific_basis": "interface-local prediction confidence",
            "reference_geometry_available_to_selector": False,
            "oracle_selection": False,
            "selection_receipt_timing": "written exclusively before any reference comparison",
        },
        "prior_rule_selection_history_posthoc_training_study": {
            "original_protocol": {"successes": 2, "denominator": 5},
            "r10_d5_top_confidence_score": {"successes": 2, "denominator": 5},
            "highest_complex_iplddt_exploratory": {"successes": 3, "denominator": 5},
            "status": "exploratory rule observed on prior five seeds and fixed here before new inference",
        },
        "cutoffs_unchanged": CUTOFFS,
        "aggregate_rule": "all four cutoffs for selected model and at least 3/5 seeds; every failed seed remains in denominator and cannot pass",
        "minimum_successes": 3,
        "denominator": 5,
        "source_settings_required": SOURCE_SETTINGS,
        "effective_protocol_settings": PROTOCOL_SETTINGS,
        "potentials": "off",
        "prediction_reuse": False,
        "frozen_cache_reuse": "official worker _install_frozen_seed_cache from source seed 23; processed MSA only",
        "mutates_sql_or_store": False,
        "reference_use": "evaluator only; never prediction template, potential, constraint, or selector input",
        "sources": {
            "reference": _file_record(paths["reference"]),
            "metadata": _file_record(paths["metadata"]),
            "checkpoint": _file_record(paths["checkpoint"]),
            "executable": _file_record(paths["boltz_executable"]),
            "worker": {"path": str(Path(worker.__file__).resolve()), "sha256": module_hash()},
            "runner": {"path": str(Path(__file__).resolve()), "sha256": runner_hash()},
            "frozen_seed_directory": str(paths["frozen_seed_directory"]),
            "frozen_source_seed": 23,
            "frozen_source_receipt": _file_record(source_receipt),
            "immutable_processed_msa_entries": _processed_msa_manifest(paths["frozen_seed_directory"]),
            "boltz_help_sha256": help_info["sha256"],
            "tool_environment": environment,
            "tool_environment_sha256": _canonical_hash(environment),
        },
        "timeout_seconds": args.timeout_seconds,
    }


def _write_plan(output: Path, plan: dict) -> None:
    plan_path = output / "plan.json"
    write_json_exclusive(plan_path, plan)
    digest = sha256_file(plan_path)
    with (output / "plan.sha256").open("x", encoding="ascii", newline="\n") as stream:
        stream.write(digest + "  plan.json\n")


def _verify_preflight(output: Path, plan: dict, environment: dict) -> None:
    plan_path = output / "plan.json"
    hash_path = output / "plan.sha256"
    expected = hash_path.read_text(encoding="ascii").strip().split()[0]
    if sha256_file(plan_path) != expected:
        raise worker.BoltzWorkerError("PREREGISTRATION_PLAN_CHANGED")
    if json.loads(plan_path.read_text(encoding="utf-8")) != plan:
        raise worker.BoltzWorkerError("PREREGISTRATION_PLAN_CONTENT_CHANGED")
    for name in ("reference", "metadata", "checkpoint", "executable", "worker", "runner",
                 "frozen_source_receipt"):
        record = plan["sources"][name]
        if sha256_file(Path(record["path"])) != record["sha256"]:
            raise worker.BoltzWorkerError("PREFLIGHT_SOURCE_CHANGED:" + name)
    frozen = Path(plan["sources"]["frozen_seed_directory"])
    for entry in plan["sources"]["immutable_processed_msa_entries"]:
        if sha256_file(frozen / entry["relative_path"]) != entry["sha256"]:
            raise worker.BoltzWorkerError("PREFLIGHT_PROCESSED_MSA_CHANGED")
    if _canonical_hash(environment) != plan["sources"]["tool_environment_sha256"]:
        raise worker.BoltzWorkerError("PREFLIGHT_ENVIRONMENT_CHANGED")


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--output", required=True)
    value.add_argument("--frozen-seed-directory", required=True)
    value.add_argument("--reference", default="cases/design_sources/6BOY.cif")
    value.add_argument("--metadata", default="cases/design_sources/crbn_benchmark.json")
    value.add_argument("--boltz-executable", required=True)
    value.add_argument("--checkpoint", required=True)
    value.add_argument("--cache", required=True)
    value.add_argument("--timeout-seconds", type=float, default=1800.0)
    value.add_argument("--prepare-only", action="store_true")
    return value


def _resolve(value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def _run_seed(seed: int, run_dir: Path, paths: dict[str, Path], args: argparse.Namespace,
              help_info: dict, environment: dict) -> dict:
    run_dir.mkdir(exist_ok=False)
    receipt: dict[str, Any] = {
        "seed": seed, "status": "failure", "selected_model_pass": False,
        "confidence_baseline_model_pass": False,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "models": [{"model_index": index, "errors": []} for index in MODEL_INDICES],
    }
    input_yaml = run_dir / "6BOY_CRBN.yaml"
    boltz_out = run_dir / "boltz_output"
    log = run_dir / "inference.log"
    reuse = None
    installed = False
    process_ok = False
    selected = None
    baseline_selected = None
    artifact_valid = {index: False for index in MODEL_INDICES}
    comparison_valid = {index: False for index in MODEL_INDICES}
    fatal_error = None
    try:
        packet = worker.prepare_input(paths["reference"], paths["metadata"], input_yaml, "server")
        if packet["input_sha256"] != worker.FROZEN_CRBN_INPUT_SHA256:
            raise worker.BoltzWorkerError("FROZEN_INPUT_HASH_MISMATCH")
        reuse = worker._install_frozen_seed_cache(
            source_directory=paths["frozen_seed_directory"], destination=boltz_out,
            input_yaml=input_yaml, reference=paths["reference"], checkpoint=paths["checkpoint"],
            executable=paths["boltz_executable"], environment=environment,
            settings=SOURCE_SETTINGS, help_sha256=help_info["sha256"],
        )
        installed = True
        base = worker.build_command(
            executable=paths["boltz_executable"], input_yaml=input_yaml, out_dir=boltz_out,
            cache=paths["cache"], checkpoint=paths["checkpoint"], accelerator="gpu", seed=seed,
            msa_mode="server", max_msa_seqs=256, help_has_seed=True,
        )
        command = effective_command(base, help_info["text"])
        receipt.update({"input": packet, "frozen_cache_reuse": reuse,
                        "effective_protocol_settings": PROTOCOL_SETTINGS,
                        "effective_command": command})
        exit_code, state, elapsed = worker._run_process(command, log, args.timeout_seconds)
        process_ok = exit_code == 0 and state == "completed"
        receipt.update({"exit_code": exit_code, "process_state": state,
                        "elapsed_seconds": elapsed})

        valid_candidates = []
        for item in receipt["models"]:
            index = item["model_index"]
            cifs = _model_matches(boltz_out, index, ".cif")
            confidence_files = _model_matches(boltz_out, index, ".json", confidence=True)
            prediction = cifs[0] if len(cifs) == 1 else None
            confidence_path = confidence_files[0] if len(confidence_files) == 1 else None
            if prediction is None:
                item["errors"].append(f"PREDICTION_MISSING_OR_AMBIGUOUS:{len(cifs)}")
            else:
                item["prediction_path"] = str(prediction)
                item["prediction_sha256"] = sha256_file(prediction)
            if confidence_path is None:
                item["errors"].append(f"CONFIDENCE_MISSING_OR_AMBIGUOUS:{len(confidence_files)}")
            else:
                item["confidence_path"] = str(confidence_path)
                item["confidence_sha256"] = sha256_file(confidence_path)
                try:
                    raw = worker._validate_confidence(confidence_path)
                    item["confidence_raw"] = raw
                    item["complex_iplddt"] = complex_iplddt(raw)
                    item["confidence_score"] = confidence_score(raw)
                    valid_candidates.append(item)
                except Exception as error:
                    item["errors"].append("CONFIDENCE_INVALID:" +
                                          (str(error) or type(error).__name__))
            artifact_valid[index] = prediction is not None and item in valid_candidates

        if len(valid_candidates) == 5:
            selected = select_model(valid_candidates)
            baseline_selected = select_confidence_baseline(valid_candidates)
            selection_receipt = {
                "rule": "highest finite raw complex_iplddt; lower model index tie",
                "selected_model_index": selected,
                "confidence_baseline_selected_model_index": baseline_selected,
                "candidates": [
                    {"model_index": item["model_index"],
                     "complex_iplddt": item["complex_iplddt"],
                     "confidence_score": item["confidence_score"]}
                    for item in valid_candidates
                ],
                "reference_metrics_consulted": False,
                "written_before_comparison": True,
            }
            write_json_exclusive(run_dir / "selection-receipt.json", selection_receipt)
            receipt["selection"] = selection_receipt
        else:
            receipt["selection_error"] = "ALL_FIVE_VALID_RAW_CONFIDENCES_REQUIRED"

        # Reference comparison starts only after the immutable selection receipt exists.
        for item in receipt["models"]:
            index = item["model_index"]
            prediction_path = item.get("prediction_path")
            if prediction_path is None:
                item["passes_cutoffs"] = False
                item["errors"].append("COMPARISON_SKIPPED_NO_PREDICTION")
                continue
            try:
                comparison = worker.compare_structures(paths["reference"], Path(prediction_path), packet)
                item["metrics"] = comparison["metrics"]
                item["passes_cutoffs"] = metrics_pass(item["metrics"])
                comparison_valid[index] = True
            except Exception as error:
                item["passes_cutoffs"] = False
                item["errors"].append("COMPARISON_ERROR:" +
                                      (str(error) or type(error).__name__))
    except Exception as error:
        fatal_error = error
        receipt["failure_reason"] = str(error) or type(error).__name__
        receipt["exception_type"] = type(error).__name__
    finally:
        frozen_ok = False
        if installed:
            try:
                worker._verify_frozen_post_run(reuse, boltz_out)
                receipt["frozen_post_run_validation"] = "success"
                frozen_ok = True
            except Exception as error:
                receipt["frozen_post_run_validation"] = "failure"
                receipt["frozen_post_run_error"] = str(error) or type(error).__name__
        else:
            receipt["frozen_post_run_validation"] = "not_installed"
        if boltz_out.exists():
            try:
                receipt["output_files_sha256"] = worker._hash_tree(boltz_out)
            except Exception as error:
                receipt["output_hash_error"] = str(error) or type(error).__name__
        if log.exists():
            try:
                receipt["inference_log_sha256"] = sha256_file(log)
            except Exception as error:
                receipt["log_hash_error"] = str(error) or type(error).__name__
        complete = (
            fatal_error is None and process_ok and frozen_ok and selected is not None and
            baseline_selected is not None and all(artifact_valid.values()) and
            all(comparison_valid.values())
        )
        if complete:
            selected_item = next(x for x in receipt["models"] if x["model_index"] == selected)
            baseline_item = next(x for x in receipt["models"]
                                 if x["model_index"] == baseline_selected)
            receipt["selected_model_pass"] = selected_item["passes_cutoffs"] is True
            receipt["confidence_baseline_model_pass"] = baseline_item["passes_cutoffs"] is True
            receipt["selected_complex_iplddt"] = selected_item["complex_iplddt"]
            receipt["selected_model_confidence_score"] = selected_item["confidence_score"]
            receipt["confidence_baseline_selected_score"] = baseline_item["confidence_score"]
            receipt["status"] = "success"
        else:
            receipt.setdefault("failure_reason", "INCOMPLETE_OR_INVALID_MODEL_SET")
            receipt["selected_model_pass"] = False
            receipt["confidence_baseline_model_pass"] = False
        receipt["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        write_json_exclusive(run_dir / "receipt.json", receipt)
    return receipt


def _distribution(values: list[float]) -> dict:
    if not values:
        return {"count": 0, "median": None, "range": None, "iqr": None,
                "minimum": None, "maximum": None}
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) == 1:
        q1 = q3 = ordered[0]
    elif len(ordered) % 2:
        q1 = statistics.median(ordered[:midpoint])
        q3 = statistics.median(ordered[midpoint + 1:])
    else:
        q1 = statistics.median(ordered[:midpoint])
        q3 = statistics.median(ordered[midpoint:])
    return {
        "count": len(ordered), "median": statistics.median(ordered),
        "minimum": ordered[0], "maximum": ordered[-1],
        "range": ordered[-1] - ordered[0], "iqr": q3 - q1,
    }


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    paths = {name: _resolve(getattr(args, name)) for name in
             ("output", "frozen_seed_directory", "reference", "metadata",
              "boltz_executable", "checkpoint", "cache")}
    output = paths["output"]
    try:
        if output.exists():
            raise worker.BoltzWorkerError("OUTPUT_ALREADY_EXISTS")
        if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
            raise worker.BoltzWorkerError("TIMEOUT_INVALID")
        for name in ("reference", "metadata", "boltz_executable", "checkpoint"):
            if not paths[name].is_file() or paths[name].is_symlink():
                raise worker.BoltzWorkerError(name.upper() + "_REGULAR_FILE_REQUIRED")
        if not paths["frozen_seed_directory"].is_dir() or paths["frozen_seed_directory"].is_symlink():
            raise worker.BoltzWorkerError("FROZEN_SEED_DIRECTORY_REQUIRED")
        if not set(SEEDS).isdisjoint(PRIOR_SEEDS):
            raise worker.BoltzWorkerError("HELDOUT_SEEDS_NOT_DISJOINT")
        help_info = worker.inspect_help(paths["boltz_executable"])
        if not help_info["has_seed"]:
            raise worker.BoltzWorkerError("OFFICIAL_HELP_SEED_OPTION_REQUIRED")
        environment = worker._tool_environment(paths["boltz_executable"])
        output.mkdir(parents=True, exist_ok=False)
        plan = _plan(args, paths, help_info, environment)
        _write_plan(output, plan)
        _verify_preflight(output, plan, environment)
        if args.prepare_only:
            return 0
        paths["cache"].mkdir(parents=True, exist_ok=True)
        receipts = []
        for seed in SEEDS:
            _verify_preflight(output, plan, environment)
            receipts.append(_run_seed(seed, output / f"seed-{seed}", paths, args,
                                      help_info, environment))
        passed = protocol_pass(receipts)
        successful = [item for item in receipts if item.get("status") == "success"]
        summary = {
            "format": "tpd-crbn-seed-holdout-summary/1.0",
            "label": LABEL,
            "scope": plan["scope"],
            "seed_count": 5,
            "denominator_including_failed_units": 5,
            "successful_seed_count": len(successful),
            "failed_unit_count": 5 - len(successful),
            "failed_units": [item.get("seed") for item in receipts
                             if item.get("status") != "success"],
            "selected_complex_iplddt_results_per_seed": [
                {"seed": item["seed"], "status": item["status"],
                 "model_index": item.get("selection", {}).get("selected_model_index"),
                 "complex_iplddt": item.get("selected_complex_iplddt"),
                 "passes": item.get("selected_model_pass") is True}
                for item in receipts
            ],
            "confidence_score_baseline_results_per_seed": [
                {"seed": item["seed"], "status": item["status"],
                 "model_index": item.get("selection", {}).get(
                     "confidence_baseline_selected_model_index"),
                 "confidence_score": item.get("confidence_baseline_selected_score"),
                 "passes": item.get("confidence_baseline_model_pass") is True}
                for item in receipts
            ],
            "selected_complex_iplddt_pass_count": sum(
                item.get("selected_model_pass") is True for item in receipts),
            "confidence_score_baseline_pass_count": sum(
                item.get("confidence_baseline_model_pass") is True for item in receipts),
            "selected_complex_iplddt_distribution": _distribution([
                item["selected_complex_iplddt"] for item in successful]),
            "selected_model_confidence_score_distribution": _distribution([
                item["selected_model_confidence_score"] for item in successful]),
            "confidence_baseline_score_distribution": _distribution([
                item["confidence_baseline_selected_score"] for item in successful]),
            "all_25_raw_records_preserved": all(
                len(item.get("models", [])) == 5 for item in receipts),
            "selection_receipts_written_before_comparison": all(
                item.get("status") != "success" or
                item.get("selection", {}).get("written_before_comparison") is True
                for item in receipts),
            "cutoffs": CUTOFFS,
            "minimum_successes": 3,
            "protocol_pass": passed,
            "interpretation": "Same-6BOY development-structure seed validation only; not an independent target, independent structure, or training holdout.",
            "seeds": receipts,
        }
        write_json_exclusive(output / "protocol-summary.json", summary)
        return 0 if passed else 1
    except (OSError, ValueError, KeyError, worker.BoltzWorkerError) as error:
        print(f"calibration seed-holdout error: {type(error).__name__}: {error}",
              file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

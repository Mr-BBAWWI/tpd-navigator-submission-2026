"""Run the fixed preregistered 6BOY CRBN development calibration protocol."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science import boltz_worker as worker

SEEDS = [23, 41, 61, 79, 97]
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


def sha256_file(path: Path) -> str:
    return worker.sha256_file(Path(path))


def module_hash() -> str:
    return sha256_file(Path(worker.__file__))


def runner_hash() -> str:
    return sha256_file(Path(__file__).resolve())


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")


def confidence_score(value: dict) -> float:
    score = value.get("confidence_score")
    if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(float(score)):
        raise worker.BoltzWorkerError("FINITE_CONFIDENCE_SCORE_REQUIRED")
    return float(score)


def select_model(models: list[dict]) -> int:
    """Select using model index and confidence only; geometry is deliberately inaccessible."""
    candidates = [(confidence_score(item["confidence"]), int(item["model_index"])) for item in models]
    if len(candidates) != 5 or {index for _, index in candidates} != set(MODEL_INDICES):
        raise worker.BoltzWorkerError("ALL_FIVE_MODELS_REQUIRED")
    return min(candidates, key=lambda item: (-item[0], item[1]))[1]


def metrics_pass(metrics: dict) -> bool:
    values = {}
    for key in CUTOFFS:
        try:
            value = metrics[key]
        except (KeyError, TypeError):
            return False
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return False
        value = float(value)
        if not math.isfinite(value) or value < 0:
            return False
        values[key] = value
    if values["contact_jaccard"] > 1:
        return False
    return (
        values["target_CA_RMSD_A"] <= CUTOFFS["target_CA_RMSD_A"] and
        values["ligand_heavy_atom_RMSD_after_target_alignment_A"] <=
        CUTOFFS["ligand_heavy_atom_RMSD_after_target_alignment_A"] and
        values["e3_CA_RMSD_after_target_alignment_A"] <=
        CUTOFFS["e3_CA_RMSD_after_target_alignment_A"] and
        values["contact_jaccard"] >= CUTOFFS["contact_jaccard"]
    )


def protocol_pass(seed_receipts: list[dict]) -> bool:
    seeds = [item.get("seed") for item in seed_receipts]
    return (
        len(seed_receipts) == len(SEEDS) and
        len(set(seeds)) == len(SEEDS) and
        set(seeds) == set(SEEDS) and
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


def _model_matches(out_dir: Path, model_index: int, suffix: str, confidence: bool = False) -> list[Path]:
    token = re.compile(rf"(?:^|[_-])model[_-]?{model_index}(?:[_\-.]|$)", re.I)
    return [
        path for path in out_dir.rglob("*" + suffix)
        if path.is_file() and not path.is_symlink() and token.search(path.name)
        and (not confidence or "confidence" in path.name.lower())
    ]


def _source_receipt_path(directory: Path) -> Path:
    preferred = directory / "receipt.json"
    if preferred.is_file() and not preferred.is_symlink():
        return preferred
    matches = [
        path for path in directory.rglob("*.json")
        if path.is_file() and not path.is_symlink() and "receipt" in path.name.lower()
        and "selection" not in path.name.lower()
    ]
    if len(matches) != 1:
        raise worker.BoltzWorkerError("FROZEN_SOURCE_RECEIPT_MISSING_OR_AMBIGUOUS")
    return matches[0]


def _plan(args: argparse.Namespace, paths: dict[str, Path], help_info: dict, environment: dict) -> dict:
    source_receipt = _source_receipt_path(paths["frozen_seed_directory"])
    return {
        "format": "tpd-crbn-preregistered-development-protocol/1.0",
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "development_set_not_holdout": True,
        "known_structure_training_overlap_possible": True,
        "scientifically_approved": False,
        "mutates_store_or_acceptance": False,
        "reference_use": "evaluator_only; never prediction template or constraint",
        "seeds": SEEDS,
        "expected_models_per_seed": MODEL_INDICES,
        "execution": "serial_all_five_no_early_stop",
        "selection": "highest finite confidence_score; tie smaller model index; no reference metrics",
        "cutoffs": CUTOFFS,
        "aggregate_rule": "exact seed set; all four cutoffs for selected model AND at least 3/5 seeds; failed seeds cannot pass",
        "preserved_original_baseline": {"geometric_successes": 2, "seed_count": 5},
        "no_mixed_protocol": True,
        "source_settings_required": SOURCE_SETTINGS,
        "effective_protocol_settings": PROTOCOL_SETTINGS,
        "sources": {
            "reference": {"path": str(paths["reference"]), "sha256": sha256_file(paths["reference"])},
            "metadata": {"path": str(paths["metadata"]), "sha256": sha256_file(paths["metadata"])},
            "checkpoint": {"path": str(paths["checkpoint"]), "sha256": sha256_file(paths["checkpoint"])},
            "executable": {"path": str(paths["boltz_executable"]), "sha256": sha256_file(paths["boltz_executable"])},
            "worker_module_sha256": module_hash(),
            "runner_module_sha256": runner_hash(),
            "frozen_seed_directory": str(paths["frozen_seed_directory"]),
            "frozen_source_receipt": {"path": str(source_receipt), "sha256": sha256_file(source_receipt)},
            "boltz_help_sha256": help_info["sha256"],
            "tool_environment": environment,
        },
        "timeout_seconds": args.timeout_seconds,
    }


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
        "seed": seed, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "status": "failure", "selected_model_pass": False,
        "scientifically_approved": False, "models": [],
    }
    input_yaml = run_dir / "6BOY_CRBN.yaml"
    boltz_out = run_dir / "boltz_output"
    log = run_dir / "inference.log"
    reuse = None
    reuse_installed = False
    process_ok = False
    selected = None
    artifact_valid: dict[int, bool] = {}
    comparison_valid: dict[int, bool] = {}
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
        reuse_installed = True
        base = worker.build_command(
            executable=paths["boltz_executable"], input_yaml=input_yaml, out_dir=boltz_out,
            cache=paths["cache"], checkpoint=paths["checkpoint"], accelerator="gpu", seed=seed,
            msa_mode="server", max_msa_seqs=256, help_has_seed=True,
        )
        command = effective_command(base, help_info["text"])
        receipt.update({"input": packet, "frozen_cache_reuse": reuse,
                        "effective_protocol_settings": PROTOCOL_SETTINGS, "effective_command": command})
        exit_code, state, elapsed = worker._run_process(command, log, args.timeout_seconds)
        process_ok = exit_code == 0 and state == "completed"
        receipt.update({"exit_code": exit_code, "process_state": state, "elapsed_seconds": elapsed})

        confidence_candidates = []
        for index in MODEL_INDICES:
            item: dict[str, Any] = {"model_index": index, "errors": []}
            cifs = _model_matches(boltz_out, index, ".cif")
            confidence_files = _model_matches(boltz_out, index, ".json", confidence=True)
            prediction = cifs[0] if len(cifs) == 1 else None
            confidence_path = confidence_files[0] if len(confidence_files) == 1 else None
            if prediction is None:
                item["errors"].append(f"PREDICTION_MISSING_OR_AMBIGUOUS:{len(cifs)}")
            else:
                item["prediction_path"] = str(prediction)
                item["prediction_sha256"] = sha256_file(prediction)
            confidence_valid = False
            if confidence_path is None:
                item["errors"].append(f"CONFIDENCE_MISSING_OR_AMBIGUOUS:{len(confidence_files)}")
            else:
                item["confidence_path"] = str(confidence_path)
                item["confidence_sha256"] = sha256_file(confidence_path)
                try:
                    item["confidence_raw"] = json.loads(confidence_path.read_text(encoding="utf-8"))
                except Exception as error:
                    item["errors"].append("CONFIDENCE_READ_ERROR:" + (str(error) or type(error).__name__))
                try:
                    item["confidence"] = worker._validate_confidence(confidence_path)
                    confidence_score(item["confidence"])
                    confidence_valid = True
                    confidence_candidates.append(item)
                except Exception as error:
                    item["errors"].append("CONFIDENCE_INVALID:" + (str(error) or type(error).__name__))
            artifact_valid[index] = prediction is not None and confidence_valid
            receipt["models"].append(item)

        if len(confidence_candidates) == 5:
            selected = select_model(confidence_candidates)
            selection_receipt = {
                "selected_model_index": selected,
                "candidates": [{"model_index": item["model_index"],
                                "confidence_score": confidence_score(item["confidence"])}
                               for item in confidence_candidates],
            }
            write_json_exclusive(run_dir / "selection-receipt.json", selection_receipt)
            receipt["selection"] = selection_receipt
        else:
            receipt["selection_error"] = "ALL_FIVE_VALID_CONFIDENCES_REQUIRED"

        for item in receipt["models"]:
            index = item["model_index"]
            comparison_valid[index] = False
            prediction_path = item.get("prediction_path")
            if prediction_path is None:
                item["errors"].append("COMPARISON_SKIPPED_NO_PREDICTION")
                continue
            try:
                comparison = worker.compare_structures(paths["reference"], Path(prediction_path), packet)
                item["metrics"] = comparison["metrics"]
                item["passes_cutoffs"] = metrics_pass(item["metrics"])
                comparison_valid[index] = True
            except Exception as error:
                item["passes_cutoffs"] = False
                item["errors"].append("COMPARISON_ERROR:" + (str(error) or type(error).__name__))
    except Exception as error:
        fatal_error = error
        receipt["failure_reason"] = str(error) or type(error).__name__
        receipt["exception_type"] = type(error).__name__
    finally:
        frozen_ok = False
        if reuse_installed:
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
            set(artifact_valid) == set(MODEL_INDICES) and all(artifact_valid.values()) and
            set(comparison_valid) == set(MODEL_INDICES) and all(comparison_valid.values())
        )
        if complete:
            selected_item = next(item for item in receipt["models"] if item["model_index"] == selected)
            receipt["selected_model_pass"] = selected_item.get("passes_cutoffs") is True
            receipt["status"] = "success"
        else:
            receipt["selected_model_pass"] = False
            receipt.setdefault("failure_reason", "INCOMPLETE_OR_INVALID_MODEL_SET")
        receipt["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        write_json_exclusive(run_dir / "receipt.json", receipt)
    return receipt


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
        help_info = worker.inspect_help(paths["boltz_executable"])
        if not help_info["has_seed"]:
            raise worker.BoltzWorkerError("OFFICIAL_HELP_SEED_OPTION_REQUIRED")
        environment = worker._tool_environment(paths["boltz_executable"])
        output.mkdir(parents=True, exist_ok=False)
        plan = _plan(args, paths, help_info, environment)
        write_json_exclusive(output / "plan.json", plan)
        if args.prepare_only:
            return 0
        paths["cache"].mkdir(parents=True, exist_ok=True)
        receipts = [
            _run_seed(seed, output / f"seed-{seed}", paths, args, help_info, environment)
            for seed in SEEDS
        ]
        passed = protocol_pass(receipts)
        summary = {
            "format": "tpd-crbn-preregistered-development-summary/1.0",
            "development_set_not_holdout": True,
            "known_structure_training_overlap_possible": True,
            "scientifically_approved": False,
            "preserved_original_baseline": {"geometric_successes": 2, "seed_count": 5},
            "effective_protocol": "10 recycling / 200 sampling / 5 diffusion; no mixed protocol",
            "seeds": receipts, "seed_count": 5,
            "successful_seed_count": sum(x["status"] == "success" for x in receipts),
            "selected_pass_count": sum(x.get("selected_model_pass") is True for x in receipts),
            "missing_or_failed_seeds": [x.get("seed") for x in receipts if x.get("status") != "success"],
            "all_25_model_metrics_preserved": all(
                len(x.get("models", [])) == 5 and
                all("metrics" in model for model in x.get("models", []))
                for x in receipts
            ),
            "cutoffs": CUTOFFS, "protocol_pass": passed,
            "interpretation": "Evidence proposal only; no store, acceptance, approval, or efficacy mutation.",
        }
        write_json_exclusive(output / "protocol-summary.json", summary)
        return 0 if passed else 1
    except (OSError, ValueError, worker.BoltzWorkerError) as error:
        print(f"calibration protocol error: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

"""Run the preregistered two-candidate novel-MSA ternary protocol."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import inspect
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science import boltz_worker as worker
from packages.science import novel_ternary as novel
from scripts import inspect_ternary_ensemble as ensemble
from scripts.run_calibration_protocol import effective_command

SEEDS = [23, 41, 61, 79, 97]
MODELS = list(range(5))
CANDIDATES = [
    ("D-99b12e64986a", "CRBN", "W-c2afc5e73c1a--CRBN.plan.json"),
    ("D-b39273b7a53b", "VHL", "W-80f8f4a11b5d--VHL.plan.json"),
]
RESULT_NAME = "boltz_results_novel_ternary"
EXPECTED_SETTINGS = {
    "accelerator": "gpu",
    "recycling_steps": 10,
    "sampling_steps": 200,
    "diffusion_samples": 5,
    "max_parallel_samples": 1,
    "max_msa_seqs": 256,
    "use_msa_server": False,
    "msa_pairing_strategy": "greedy",
    "preprocessing_threads": 1,
    "potentials": False,
    "templates": [],
    "constraints": [],
}


class ProtocolError(RuntimeError):
    pass


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha(path: Path) -> str:
    return worker.sha256_file(path)


def digest(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def regular(path: Path, code: str) -> Path:
    if not path.is_file() or path.is_symlink():
        raise ProtocolError(code)
    return path


def resolve(value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def write_json(path: Path, value: Any, exclusive: bool = True) -> None:
    mode = "x" if exclusive else "w"
    with path.open(mode, encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")


def read_json(path: Path) -> dict:
    regular(path, "JSON_REGULAR_FILE_REQUIRED")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ProtocolError("JSON_INVALID:" + str(path)) from error
    if not isinstance(value, dict):
        raise ProtocolError("JSON_OBJECT_REQUIRED:" + str(path))
    return value


def append_progress(output: Path, value: dict) -> None:
    with (output / "progress.jsonl").open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def tree_hashes(root: Path, allowed_roots: Iterable[str] | None = None) -> dict[str, str]:
    if not root.is_dir() or root.is_symlink():
        raise ProtocolError("TREE_DIRECTORY_REQUIRED:" + str(root))
    allowed = set(allowed_roots or [])
    result = {}
    for current, directories, files in os.walk(root, followlinks=False):
        base = Path(current)
        for name in directories:
            path = base / name
            if path.is_symlink():
                raise ProtocolError("TREE_SYMLINK_REJECTED:" + str(path))
        for name in files:
            path = base / name
            if path.is_symlink() or not path.is_file():
                raise ProtocolError("TREE_REGULAR_FILE_REQUIRED:" + str(path))
            key = path.relative_to(root).as_posix()
            if allowed and key.split("/", 1)[0] not in allowed:
                continue
            result[key] = sha(path)
    return dict(sorted(result.items()))


def verify_hashes(root: Path, expected: dict[str, str], allowed_roots=None) -> dict[str, str]:
    actual = tree_hashes(root, allowed_roots)
    if actual != expected:
        raise ProtocolError("IMMUTABLE_TREE_HASH_MISMATCH:" + str(root))
    return actual


def source_plan(batch_root: Path, filename: str) -> Path:
    matches = [p for p in batch_root.rglob(filename) if p.is_file() and not p.is_symlink()]
    if len(matches) != 1:
        raise ProtocolError(f"SOURCE_PLAN_MISSING_OR_AMBIGUOUS:{filename}:{len(matches)}")
    return matches[0]


def validate_baseline(batch_root: Path) -> tuple[dict[str, str], int]:
    plans = [p for p in batch_root.rglob("*.plan.json") if p.is_file() and not p.is_symlink()]
    if len(plans) != 6:
        raise ProtocolError("OLD_PLAN_COUNT_NOT_SIX")
    raw = 0
    for path in plans:
        plan = read_json(path)
        novel.verify_plan(plan)
        run_dir = path.with_name(path.name[:-len(".plan.json")])
        for seed in SEEDS:
            prediction = (
                run_dir / f"seed-{seed}" / "boltz_output" / RESULT_NAME /
                "predictions" / "novel_ternary" / "novel_ternary_model_0.cif"
            )
            regular(prediction, "OLD_BASELINE_RAW_PREDICTION_MISSING")
            raw += 1
    if raw != 30:
        raise ProtocolError("OLD_BASELINE_RAW_COUNT_NOT_30")
    return tree_hashes(batch_root), raw


def make_plan(old: dict) -> dict:
    novel.verify_plan(old)
    plan = copy.deepcopy(old)
    target = plan["sources"]["target"]["canonical_sequence"]
    e3 = plan["sources"]["e3"]["canonical_sequence"]
    smiles = plan["candidate_graph"]["actualmapped_smiles"]
    document = novel._document(target, e3, smiles, "server")
    text = yaml.safe_dump(document, sort_keys=False, allow_unicode=False)
    plan["msa_mode"] = "server"
    plan["boltz_input"]["yaml"] = text
    plan["boltz_input"]["sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
    plan["boltz_input"]["templates"] = []
    plan["boltz_input"]["constraints"] = []
    plan["boltz_input"]["potentials"] = False
    plan.pop("plan_digest", None)
    plan["plan_digest"] = novel._digest(plan)
    novel.verify_plan(plan)
    return plan


def wrapper_source() -> str:
    return '''from __future__ import annotations
import argparse
from pathlib import Path
from boltz.main import process_inputs

p = argparse.ArgumentParser()
p.add_argument("--input-yaml", required=True)
p.add_argument("--out-dir", required=True)
p.add_argument("--cache", required=True)
a = p.parse_args()
cache = Path(a.cache)
process_inputs(
    data=[Path(a.input_yaml)],
    out_dir=Path(a.out_dir),
    ccd_path=cache / "ccd.pkl",
    mol_dir=cache / "mols",
    msa_server_url="https://api.colabfold.com",
    msa_pairing_strategy="greedy",
    max_msa_seqs=256,
    use_msa_server=True,
    boltz2=True,
    preprocessing_threads=1,
)
'''


def wrapper_python(executable: Path) -> Path:
    function = worker._wrapper_python
    parameters = inspect.signature(function).parameters
    value = function(executable) if parameters else function()
    return regular(Path(value), "BOLTZ_WRAPPER_PYTHON_REQUIRED")


def run_preprocessing(
    python: Path, wrapper: Path, input_yaml: Path, result_dir: Path,
    cache: Path, log: Path, timeout: float,
) -> dict:
    command = [
        str(python), str(wrapper), "--input-yaml", str(input_yaml),
        "--out-dir", str(result_dir), "--cache", str(cache),
    ]
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = ""
    started = time.monotonic()
    try:
        with log.open("x", encoding="utf-8", newline="\n") as stream:
            completed = subprocess.run(
                command, stdout=stream, stderr=subprocess.STDOUT, text=True,
                env=environment, timeout=timeout, check=False,
            )
    except subprocess.TimeoutExpired as error:
        raise ProtocolError("MSA_PREPROCESSING_TIMEOUT") from error
    if completed.returncode != 0:
        raise ProtocolError(f"MSA_PREPROCESSING_FAILED:{completed.returncode}")
    return {
        "command": command,
        "elapsed_seconds": time.monotonic() - started,
        "exit_code": completed.returncode,
        "log_sha256": sha(log),
    }


def csv_query(path: Path, expected: str) -> dict:
    regular(path, "MSA_CSV_REQUIRED")
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not {"sequence", "key"}.issubset(reader.fieldnames):
            raise ProtocolError("MSA_CSV_SEQUENCE_KEY_COLUMNS_REQUIRED")
        try:
            row = next(reader)
        except StopIteration as error:
            raise ProtocolError("MSA_CSV_QUERY_ROW_REQUIRED") from error
    sequence = "".join(str(row["sequence"]).split()).replace("-", "").upper()
    if sequence != expected.upper():
        raise ProtocolError("MSA_QUERY_SEQUENCE_MISMATCH:" + path.name)
    return {
        "path": str(path), "key": row["key"], "query_sequence": sequence,
        "query_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "file_sha256": sha(path),
    }


def record_sequences(value: Any) -> list[str]:
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key in {"sequence", "seq"} and isinstance(item, str):
                clean = "".join(item.split()).replace("-", "").upper()
                if clean and set(clean) <= set("ABCDEFGHIKLMNPQRSTVWXYZUO"):
                    found.append(clean)
            else:
                found.extend(record_sequences(item))
    elif isinstance(value, list):
        for item in value:
            found.extend(record_sequences(item))
    return found


def validate_preprocessed(result_dir: Path, plan: dict) -> dict:
    record_path = result_dir / "processed" / "records" / "novel_ternary.json"
    npz = [
        result_dir / "processed" / "msa" / "novel_ternary_0.npz",
        result_dir / "processed" / "msa" / "novel_ternary_1.npz",
    ]
    regular(record_path, "PROCESSED_RECORD_REQUIRED")
    for path in npz:
        regular(path, "PROCESSED_MSA_NPZ_REQUIRED")
    target = plan["sources"]["target"]["canonical_sequence"]
    e3 = plan["sources"]["e3"]["canonical_sequence"]
    csvs = {}
    for index, expected in enumerate((target, e3)):
        matches = [
            path for path in result_dir.rglob(f"novel_ternary_{index}.csv")
            if path.is_file() and not path.is_symlink()
        ]
        if len(matches) != 1:
            raise ProtocolError(f"MSA_CSV_MISSING_OR_AMBIGUOUS:{index}:{len(matches)}")
        csvs[str(index)] = csv_query(matches[0], expected)
    record = read_json(record_path)
    sequences = record_sequences(record)
    if sequences and (target.upper() not in sequences or e3.upper() not in sequences):
        raise ProtocolError("PROCESSED_RECORD_CHAIN_SEQUENCE_MISMATCH")
    hashes = tree_hashes(result_dir, {"msa", "processed"})
    if not hashes or any("predictions/" in key for key in hashes):
        raise ProtocolError("FROZEN_PREPROCESSING_TREE_INVALID")
    return {
        "source_sequence_sha256": {
            "target": hashlib.sha256(target.encode()).hexdigest(),
            "e3": hashlib.sha256(e3.encode()).hexdigest(),
        },
        "msa_csv_queries": csvs,
        "processed_msa_sha256": {str(i): sha(path) for i, path in enumerate(npz)},
        "record_sha256": sha(record_path),
        "frozen_files_sha256": hashes,
        "external_ml_key": None,
    }


def copy_frozen(frozen: Path, destination: Path, expected: dict[str, str]) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=False)
    for name in ("msa", "processed"):
        source = frozen / name
        if source.exists():
            if not source.is_dir() or source.is_symlink():
                raise ProtocolError("FROZEN_COMPONENT_INVALID:" + name)
            shutil.copytree(source, destination / name, symlinks=False)
    actual = tree_hashes(destination, {"msa", "processed"})
    if actual != expected:
        raise ProtocolError("COPIED_FROZEN_HASH_MISMATCH")
    return actual


def confidence_score(path: Path) -> tuple[dict, float]:
    value = worker._validate_confidence(path)
    score = value.get("confidence_score")
    if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
        raise ProtocolError("FINITE_CONFIDENCE_SCORE_REQUIRED")
    return value, float(score)


def model_paths(result_dir: Path, index: int) -> tuple[Path, Path]:
    base = result_dir / "predictions" / "novel_ternary"
    return (
        base / f"novel_ternary_model_{index}.cif",
        base / f"confidence_novel_ternary_model_{index}.json",
    )


def inspect_model(run_dir: Path, index: int, cif: Path, confidence: Path, plan: dict) -> dict:
    view = run_dir / "inspection-views" / f"model-{index}"
    prediction_dir = view / "predictions" / "novel_ternary"
    prediction_dir.mkdir(parents=True, exist_ok=False)
    copied_cif = prediction_dir / cif.name
    copied_confidence = prediction_dir / confidence.name
    shutil.copy2(cif, copied_cif)
    shutil.copy2(confidence, copied_confidence)
    if sha(copied_cif) != sha(cif) or sha(copied_confidence) != sha(confidence):
        raise ProtocolError("INSPECTION_VIEW_COPY_HASH_MISMATCH")
    inspection = novel.assess_novel_outputs(view, plan)
    if not isinstance(inspection, dict):
        raise ProtocolError("INSPECTION_OBJECT_REQUIRED")
    if inspection.get("status") != "computed_hypothesis":
        raise ProtocolError("INSPECTION_NOT_COMPUTED_HYPOTHESIS")
    if inspection.get("actual_computation") is not True:
        raise ProtocolError("INSPECTION_ACTUAL_COMPUTATION_REQUIRED")
    return {
        "view_directory": str(view),
        "view_hashes": tree_hashes(view),
        "prediction": str(copied_cif),
        "confidence": str(copied_confidence),
        "inspection": inspection,
    }


def allowed_manifest_change(key: str) -> bool:
    return key == "processed/manifest.json"


def run_seed(candidate: dict, seed: int, args, paths: dict, help_info: dict) -> dict:
    run_dir = paths["output"] / candidate["slug"] / f"seed-{seed}"
    run_dir.mkdir(parents=True, exist_ok=False)
    result_dir = run_dir / "boltz_output" / RESULT_NAME
    receipt: dict[str, Any] = {
        "candidate_id": candidate["id"], "e3_type": candidate["e3"], "seed": seed,
        "started_utc": now(), "status": "failed", "actual_computation": False,
        "scientific_approved": False, "automatic_approval": False,
        "settings": copy.deepcopy(EXPECTED_SETTINGS), "models": [],
    }
    write_json(run_dir / "settings.json", {
        "settings": EXPECTED_SETTINGS, "seed": seed,
        "plan_digest": candidate["plan"]["plan_digest"],
    })
    input_yaml = run_dir / "novel_ternary.yaml"
    input_yaml.write_text(candidate["plan"]["boltz_input"]["yaml"], encoding="utf-8", newline="\n")
    if sha(input_yaml) != candidate["plan"]["boltz_input"]["sha256"]:
        raise ProtocolError("RUN_INPUT_YAML_HASH_MISMATCH")
    before = copy_frozen(candidate["frozen_result"], result_dir, candidate["frozen_hashes"])
    log = run_dir / "inference.log"
    try:
        base = worker.build_command(
            executable=paths["boltz_executable"], input_yaml=input_yaml,
            out_dir=run_dir / "boltz_output", cache=paths["cache"],
            checkpoint=paths["checkpoint"], accelerator="gpu", seed=seed,
            msa_mode="server", max_msa_seqs=256, help_has_seed=help_info["has_seed"],
        )
        command = effective_command(base, help_info["text"])
        receipt["effective_command"] = command
        code, state, elapsed = worker._run_process(command, log, args.timeout_seconds)
        receipt.update({"exit_code": code, "process_state": state, "elapsed_seconds": elapsed})
        receipt["actual_computation"] = code == 0 and state == "completed"
        receipt["raw_output_files_sha256"] = tree_hashes(run_dir / "boltz_output")

        candidates = []
        for index in MODELS:
            item: dict[str, Any] = {"model_index": index, "status": "invalid", "errors": []}
            cif, confidence = model_paths(result_dir, index)
            if not cif.is_file() or cif.is_symlink():
                item["errors"].append("EXACT_PREDICTION_MISSING")
            if not confidence.is_file() or confidence.is_symlink():
                item["errors"].append("EXACT_CONFIDENCE_MISSING")
            if not item["errors"]:
                try:
                    validated, score = confidence_score(confidence)
                    item.update({
                        "raw_prediction_path": str(cif), "raw_prediction_sha256": sha(cif),
                        "raw_confidence_path": str(confidence), "raw_confidence_sha256": sha(confidence),
                        "confidence": validated, "confidence_score": score, "status": "available",
                    })
                    candidates.append((score, index))
                except Exception as error:
                    item["errors"].append(type(error).__name__ + ":" + str(error))
            receipt["models"].append(item)

        if len(candidates) == 5:
            selected = min(candidates, key=lambda pair: (-pair[0], pair[1]))[1]
            selection = {
                "rule": "highest finite confidence_score; tie lowest model index",
                "selected_model_index": selected,
                "candidates": [{"model_index": i, "confidence_score": s} for s, i in candidates],
            }
            write_json(run_dir / "selection-receipt.json", selection)
            receipt["selection"] = selection
        else:
            selected = None
            receipt["selection_error"] = "ALL_FIVE_MODELS_REQUIRED"

        for item in receipt["models"]:
            if item["status"] != "available":
                continue
            try:
                item["inspection_view"] = inspect_model(
                    run_dir, item["model_index"], Path(item["raw_prediction_path"]),
                    Path(item["raw_confidence_path"]), candidate["plan"],
                )
                item["status"] = "inspected"
            except Exception as error:
                item["errors"].append("INSPECTION_ERROR:" + type(error).__name__ + ":" + str(error))
                item["status"] = "invalid"
        after = tree_hashes(result_dir, {"msa", "processed"})
        changed = sorted(key for key in set(before) | set(after) if before.get(key) != after.get(key))
        forbidden = [key for key in changed if not allowed_manifest_change(key)]
        receipt["processed_reuse_post_run"] = {
            "changed_files": changed, "allowed_official_manifest_rewrites": changed,
            "forbidden_changes": forbidden,
        }
        verify_hashes(candidate["frozen_result"], candidate["frozen_hashes"], {"msa", "processed"})
        if forbidden:
            raise ProtocolError("COPIED_REUSE_MUTATED_OUTSIDE_PROCESSED_MANIFEST")
        if not receipt["actual_computation"]:
            raise ProtocolError("INFERENCE_NOT_COMPLETED")
        if selected is None or any(item["status"] != "inspected" for item in receipt["models"]):
            raise ProtocolError("INCOMPLETE_FIVE_MODEL_SET")
        receipt["selected_model_inspection"] = next(
            item["inspection_view"]["inspection"] for item in receipt["models"]
            if item["model_index"] == selected
        )
        receipt["selected_prediction"] = next(
            item["inspection_view"]["prediction"] for item in receipt["models"]
            if item["model_index"] == selected
        )
        receipt["status"] = "completed"
    except KeyboardInterrupt:
        receipt["failure_reason"] = "CANCELLED_BY_OPERATOR"
        write_json(run_dir / "cancellation-receipt.json", {"seed": seed, "cancelled_utc": now()})
        raise
    except Exception as error:
        receipt["failure_reason"] = type(error).__name__ + ":" + str(error)
        write_json(run_dir / "failure-receipt.json", {
            "seed": seed, "failed_utc": now(), "reason": receipt["failure_reason"],
        })
    finally:
        if log.exists() and log.is_file() and not log.is_symlink():
            receipt["inference_log_sha256"] = sha(log)
        receipt["finished_utc"] = now()
        write_json(run_dir / "receipt.json", receipt)
    return receipt


def get_path(value: Any, keys: list[str]) -> Any:
    current = value
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def numeric_by_name(value: Any, token: str) -> float | None:
    if isinstance(value, dict):
        for key, item in value.items():
            normalized = "".join(ch.lower() for ch in str(key) if ch.isalnum())
            if token in normalized and isinstance(item, (int, float)) and not isinstance(item, bool):
                if math.isfinite(float(item)):
                    return float(item)
        for item in value.values():
            found = numeric_by_name(item, token)
            if found is not None:
                return found
    elif isinstance(value, list):
        for item in value:
            found = numeric_by_name(item, token)
            if found is not None:
                return found
    return None


def iqr(values: list[float]) -> float | None:
    if len(values) != 5 or not all(math.isfinite(value) for value in values):
        return None
    ordered = sorted(values)
    return ordered[3] - ordered[1]


def selected_diagnostics(candidate: dict, receipts: list[dict]) -> dict:
    rows, consistency_rows = [], []
    for receipt in receipts:
        row: dict[str, Any] = {"seed": receipt["seed"], "status": receipt["status"]}
        if receipt.get("status") != "completed":
            row["error"] = receipt.get("failure_reason", "SEED_FAILED")
            rows.append(row)
            continue
        prediction = Path(receipt["selected_prediction"])
        try:
            block = novel._source_block(prediction)
            chains = novel._prediction_chain_mapping(block, candidate["plan"])
            _, _, ligand = novel._ligand_atoms(block, candidate["plan"]["candidate_graph"])
            proteins = {}
            proteins.update(worker._protein_atoms(block, chains["target"], "target"))
            proteins.update(worker._protein_atoms(block, chains["e3"], "e3"))
            row.update({
                "protein_protein": ensemble.protein_contact_signature(proteins),
                "target_warhead_site": ensemble.site_signature(
                    proteins, ligand, candidate["plan"]["candidate_graph"], "target"
                ),
                "e3_recruiter_site": ensemble.site_signature(
                    proteins, ligand, candidate["plan"]["candidate_graph"], "e3"
                ),
                "interface_clashes": ensemble.interface_clashes(proteins),
            })
            consistency_rows.append({
                "seed": receipt["seed"], "inspection": {"prediction": str(prediction)}
            })
        except Exception as error:
            row["error"] = type(error).__name__ + ":" + str(error)
        rows.append(row)
    pairs = []
    valid = [row for row in rows if "protein_protein" in row]
    for position, first in enumerate(valid):
        for second in valid[position + 1:]:
            pairs.append({
                "seed_a": first["seed"], "seed_b": second["seed"],
                "protein_protein_jaccard": ensemble.jaccard(
                    first["protein_protein"]["keys"], second["protein_protein"]["keys"]
                ),
                "target_site_jaccard": ensemble.jaccard(
                    first["target_warhead_site"]["keys"], second["target_warhead_site"]["keys"]
                ),
                "e3_site_jaccard": ensemble.jaccard(
                    first["e3_recruiter_site"]["keys"], second["e3_recruiter_site"]["keys"]
                ),
            })
    sets = {row["seed"]: row["protein_protein"]["keys"] for row in valid}
    clusters = [
        {"threshold": threshold, "memberships": ensemble.complete_link_clusters(sets, threshold)}
        for threshold in (0.3, 0.5, 0.7)
    ]
    try:
        consistency = novel._pairwise_consistency(consistency_rows, candidate["plan"])
        consistency_error = None
    except Exception as error:
        consistency, consistency_error = [], type(error).__name__ + ":" + str(error)
    return {
        "all_five_seeds_including_failures": rows,
        "pairwise_contact_jaccard": pairs,
        "cluster_sensitivity": clusters,
        "pairwise_coordinate_consistency": consistency,
        "pairwise_coordinate_consistency_error": consistency_error,
        "interpretation": "Reference-free descriptive geometry only; not automatic approval.",
    }


def successful_seed_receipts(receipts: list[dict]) -> list[dict]:
    successful = [
        receipt for receipt in receipts
        if receipt.get("status") == "completed"
        and receipt.get("actual_computation") is True
        and receipt.get("settings") == EXPECTED_SETTINGS
        and isinstance(receipt.get("seed"), int)
        and not isinstance(receipt.get("seed"), bool)
    ]
    seeds = [receipt["seed"] for receipt in successful]
    if len(receipts) != 5 or len(successful) != 5 or sorted(seeds) != SEEDS or len(set(seeds)) != 5:
        return []
    return successful


def exact_finite_number(value: Any, minimum: float, maximum: float | None = None) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    if not math.isfinite(result) or result < minimum:
        return None
    if maximum is not None and result > maximum:
        return None
    return result


def expert_metrics(receipts: list[dict]) -> dict:
    iptm, endpoint, positive = [], [], []
    successful = successful_seed_receipts(receipts)
    for receipt in successful:
        inspection = receipt.get("selected_model_inspection")
        if not isinstance(inspection, dict):
            continue
        first = exact_finite_number(
            get_path(inspection, ["model_confidence", "iptm"]), 0.0, 1.0
        )
        second = exact_finite_number(
            get_path(inspection, ["descriptive_metrics", "linker_endpoint_distance_A"]),
            0.0,
        )
        if first is not None:
            iptm.append(first)
        if second is not None:
            endpoint.append(second)
        contacts = get_path(
            inspection,
            ["descriptive_metrics", "contacts_by_protein_role_and_ligand_group"],
        )
        target = get_path(contacts, ["target", "warhead"])
        recruiter = get_path(contacts, ["e3", "recruiter"])

        def positive_contact(value):
            if isinstance(value, bool):
                return False
            if isinstance(value, (int, float)):
                return math.isfinite(float(value)) and float(value) > 0
            if isinstance(value, dict):
                count = value.get("count")
                return (
                    isinstance(count, (int, float)) and not isinstance(count, bool)
                    and math.isfinite(float(count)) and float(count) > 0
                )
            if isinstance(value, list):
                return len(value) > 0
            return False

        positive.append(positive_contact(target) and positive_contact(recruiter))
    iptm_iqr, endpoint_iqr = iqr(iptm), iqr(endpoint)
    passed = (
        len(successful) == 5 and len(iptm) == 5 and len(endpoint) == 5 and
        iptm_iqr is not None and iptm_iqr <= 0.20 and
        endpoint_iqr is not None and endpoint_iqr <= 1.5 and
        len(positive) == 5 and all(positive)
    )
    return {
        "finite_selected_ipTM_values": iptm,
        "ipTM_IQR": iptm_iqr,
        "ipTM_IQR_maximum": 0.20,
        "finite_selected_endpoint_values_A": endpoint,
        "endpoint_IQR_A": endpoint_iqr,
        "endpoint_IQR_maximum_A": 1.5,
        "positive_target_warhead_and_e3_recruiter_contacts": positive,
        "quantitative_criteria_met": passed,
        "automatic_approval": False,
        "qualitative_topology_or_clash_acceptance_used": False,
    }


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--batch-root", required=True)
    value.add_argument("--output", required=True)
    value.add_argument("--boltz-executable", required=True)
    value.add_argument("--checkpoint", required=True)
    value.add_argument("--cache", required=True)
    value.add_argument("--timeout-seconds", type=float, default=1800.0)
    value.add_argument("--msa-timeout-seconds", type=float, default=600.0)
    value.add_argument("--prepare-only", action="store_true")
    value.add_argument("--resume-prepared", action="store_true")
    return value


def prepare(args, paths: dict) -> tuple[dict, list[dict], dict]:
    output = paths["output"]
    if output.exists() or output.is_symlink():
        raise ProtocolError("OUTPUT_ALREADY_EXISTS")
    output.mkdir(parents=True, exist_ok=False)
    baseline_hashes, raw_count = validate_baseline(paths["batch_root"])
    help_info = worker.inspect_help(paths["boltz_executable"])
    if not help_info["has_seed"]:
        raise ProtocolError("OFFICIAL_HELP_SEED_OPTION_REQUIRED")
    environment = worker._tool_environment(paths["boltz_executable"])
    wrapper = output / "preprocess_novel_msa.py"
    wrapper.write_text(wrapper_source(), encoding="utf-8", newline="\n")
    python = wrapper_python(paths["boltz_executable"])
    plans_dir = output / "plans"
    frozen_root = output / "frozen"
    plans_dir.mkdir()
    frozen_root.mkdir()
    candidates = []
    for candidate_id, e3, filename in CANDIDATES:
        old_path = source_plan(paths["batch_root"], filename)
        old = read_json(old_path)
        if old["candidate_graph"]["candidate_id"] != candidate_id or old["candidate_graph"]["e3_type"] != e3:
            raise ProtocolError("SOURCE_PLAN_IDENTITY_MISMATCH")
        plan = make_plan(old)
        slug = f"{candidate_id}--{e3}"
        plan_path = plans_dir / f"{slug}.plan.json"
        write_json(plan_path, plan)
        input_yaml = frozen_root / slug / "novel_ternary.yaml"
        input_yaml.parent.mkdir()
        input_yaml.write_text(plan["boltz_input"]["yaml"], encoding="utf-8", newline="\n")
        result_dir = input_yaml.parent / RESULT_NAME
        log = input_yaml.parent / "preprocessing.log"
        preprocessing = run_preprocessing(
            python, wrapper, input_yaml, result_dir, paths["cache"], log,
            args.msa_timeout_seconds,
        )
        validation = validate_preprocessed(result_dir, plan)
        candidates.append({
            "id": candidate_id, "e3": e3, "slug": slug, "plan": plan,
            "plan_path": plan_path, "frozen_result": result_dir,
            "frozen_hashes": validation["frozen_files_sha256"],
            "preprocessing": preprocessing, "validation": validation,
            "old_plan_path": old_path, "old_plan_sha256": sha(old_path),
        })
    protocol = {
        "format": "tpd-novel-msa-preregistered-protocol/1.0",
        "created_utc": now(), "development_set_not_holdout": True,
        "scientific_approved": False, "automatic_approval": False,
        "baseline": {
            "path": str(paths["batch_root"]), "raw_result_count": raw_count,
            "read_only_files_sha256": baseline_hashes,
            "comparison": "new MSA baseline versus old single-sequence baseline",
            "six_boy_msa_reused": False,
            "reason": "Target is SMARCA2, not BRD4.",
        },
        "seeds": SEEDS, "model_indices": MODELS,
        "expected_raw_model_count": 50, "execution": "serial",
        "selection": "highest finite confidence_score; tie lowest model index; before geometry",
        "settings": EXPECTED_SETTINGS,
        "paths": {key: str(value) for key, value in paths.items()},
        "timeouts": {
            "inference_seconds": args.timeout_seconds,
            "msa_preprocessing_seconds": args.msa_timeout_seconds,
        },
        "checkpoint_sha256": sha(paths["checkpoint"]),
        "executable_sha256": sha(paths["boltz_executable"]),
        "executable_help_sha256": help_info["sha256"],
        "environment": environment,
        "harness_hashes": {
            "runner": sha(Path(__file__).resolve()),
            "worker": sha(Path(worker.__file__)),
            "novel": sha(Path(novel.__file__)),
            "calibration_command": sha(ROOT / "scripts" / "run_calibration_protocol.py"),
            "ensemble_helpers": sha(ROOT / "scripts" / "inspect_ternary_ensemble.py"),
            "preprocessing_wrapper": sha(wrapper),
            "preprocessing_python": sha(python),
        },
        "candidates": [
            {
                "candidate_id": item["id"], "e3_type": item["e3"], "slug": item["slug"],
                "new_plan_path": str(item["plan_path"]), "new_plan_sha256": sha(item["plan_path"]),
                "new_plan_digest": item["plan"]["plan_digest"],
                "old_plan_path": str(item["old_plan_path"]),
                "old_plan_sha256": item["old_plan_sha256"],
                "frozen_result": str(item["frozen_result"]),
                "preprocessing": item["preprocessing"], "msa_validation": item["validation"],
            } for item in candidates
        ],
        "limitations": [
            "No reference geometry is used for model selection.",
            "Geometric diagnostics are descriptive only.",
            "No topology, clash, efficacy, degradation, or binding claim is automatically approved.",
        ],
    }
    protocol["protocol_digest"] = digest(protocol)
    write_json(output / "final-protocol.json", protocol)
    return protocol, candidates, help_info


def resume(args, paths: dict) -> tuple[dict, list[dict], dict]:
    protocol = read_json(paths["output"] / "final-protocol.json")
    stored_digest = protocol.get("protocol_digest")
    check = copy.deepcopy(protocol)
    check.pop("protocol_digest", None)
    if stored_digest != digest(check):
        raise ProtocolError("FINAL_PROTOCOL_DIGEST_MISMATCH")
    if protocol.get("format") != "tpd-novel-msa-preregistered-protocol/1.0":
        raise ProtocolError("RESUME_PROTOCOL_FORMAT_MISMATCH")
    if protocol.get("paths") != {key: str(value) for key, value in paths.items()}:
        raise ProtocolError("RESUME_PATHS_MISMATCH")
    if protocol.get("settings") != EXPECTED_SETTINGS:
        raise ProtocolError("RESUME_SETTINGS_MISMATCH")
    if protocol.get("seeds") != SEEDS or protocol.get("model_indices") != MODELS:
        raise ProtocolError("RESUME_SEED_OR_MODEL_SET_MISMATCH")
    if protocol.get("automatic_approval") is not False or protocol.get("scientific_approved") is not False:
        raise ProtocolError("RESUME_APPROVAL_STATE_INVALID")
    if protocol.get("timeouts") != {
        "inference_seconds": args.timeout_seconds,
        "msa_preprocessing_seconds": args.msa_timeout_seconds,
    }:
        raise ProtocolError("RESUME_TIMEOUTS_MISMATCH")
    if sha(paths["checkpoint"]) != protocol.get("checkpoint_sha256"):
        raise ProtocolError("RESUME_CHECKPOINT_HASH_MISMATCH")
    if sha(paths["boltz_executable"]) != protocol.get("executable_sha256"):
        raise ProtocolError("RESUME_EXECUTABLE_HASH_MISMATCH")
    verify_hashes(paths["batch_root"], protocol["baseline"]["read_only_files_sha256"])
    help_info = worker.inspect_help(paths["boltz_executable"])
    if help_info.get("has_seed") is not True:
        raise ProtocolError("RESUME_OFFICIAL_HELP_SEED_OPTION_REQUIRED")
    if help_info.get("sha256") != protocol.get("executable_help_sha256"):
        raise ProtocolError("RESUME_HELP_HASH_MISMATCH")
    if worker._tool_environment(paths["boltz_executable"]) != protocol.get("environment"):
        raise ProtocolError("RESUME_ENVIRONMENT_MISMATCH")

    wrapper = paths["output"] / "preprocess_novel_msa.py"
    python = wrapper_python(paths["boltz_executable"])
    harness_hashes = {
        "runner": sha(Path(__file__).resolve()),
        "worker": sha(Path(worker.__file__)),
        "novel": sha(Path(novel.__file__)),
        "calibration_command": sha(ROOT / "scripts" / "run_calibration_protocol.py"),
        "ensemble_helpers": sha(ROOT / "scripts" / "inspect_ternary_ensemble.py"),
        "preprocessing_wrapper": sha(regular(wrapper, "RESUME_PREPROCESSING_WRAPPER_REQUIRED")),
        "preprocessing_python": sha(python),
    }
    if harness_hashes != protocol.get("harness_hashes"):
        raise ProtocolError("RESUME_HARNESS_HASH_MISMATCH")

    existing_seed_dirs = [
        path for path in paths["output"].rglob("seed-*")
        if path.is_dir() or path.is_symlink()
    ]
    if existing_seed_dirs:
        raise ProtocolError("RESUME_SEED_DIRECTORY_ALREADY_EXISTS")

    stored_candidates = protocol.get("candidates")
    expected_identities = [
        {
            "candidate_id": candidate_id,
            "e3_type": e3,
            "slug": f"{candidate_id}--{e3}",
        }
        for candidate_id, e3, _ in CANDIDATES
    ]
    if not isinstance(stored_candidates, list) or [
        {
            "candidate_id": item.get("candidate_id"),
            "e3_type": item.get("e3_type"),
            "slug": item.get("slug"),
        }
        for item in stored_candidates if isinstance(item, dict)
    ] != expected_identities:
        raise ProtocolError("RESUME_CANDIDATE_IDENTITIES_MISMATCH")

    candidates = []
    for stored in stored_candidates:
        plan_path = Path(stored["new_plan_path"])
        if sha(plan_path) != stored["new_plan_sha256"]:
            raise ProtocolError("RESUME_PLAN_HASH_MISMATCH")
        plan = read_json(plan_path)
        novel.verify_plan(plan)
        if plan.get("plan_digest") != stored.get("new_plan_digest"):
            raise ProtocolError("RESUME_PLAN_DIGEST_MISMATCH")
        graph = plan.get("candidate_graph", {})
        if (
            graph.get("candidate_id") != stored["candidate_id"]
            or graph.get("e3_type") != stored["e3_type"]
        ):
            raise ProtocolError("RESUME_PLAN_IDENTITY_MISMATCH")
        result = Path(stored["frozen_result"])
        hashes = stored["msa_validation"]["frozen_files_sha256"]
        verify_hashes(result, hashes, {"msa", "processed"})
        candidates.append({
            "id": stored["candidate_id"], "e3": stored["e3_type"],
            "slug": stored["slug"], "plan": plan, "plan_path": plan_path,
            "frozen_result": result, "frozen_hashes": hashes,
        })
    return protocol, candidates, help_info


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    paths = {name: resolve(getattr(args, name)) for name in (
        "batch_root", "output", "boltz_executable", "checkpoint", "cache"
    )}
    try:
        if args.prepare_only and args.resume_prepared:
            raise ProtocolError("PREPARE_ONLY_AND_RESUME_MUTUALLY_EXCLUSIVE")
        for value in (args.timeout_seconds, args.msa_timeout_seconds):
            if not math.isfinite(value) or value <= 0:
                raise ProtocolError("TIMEOUT_INVALID")
        regular(paths["boltz_executable"], "BOLTZ_EXECUTABLE_REQUIRED")
        regular(paths["checkpoint"], "CHECKPOINT_REQUIRED")
        paths["cache"].mkdir(parents=True, exist_ok=True)
        if not (paths["cache"] / "mols").is_dir():
            raise ProtocolError("CACHE_MOLS_DIRECTORY_REQUIRED")
        if args.resume_prepared:
            protocol, candidates, help_info = resume(args, paths)
        else:
            protocol, candidates, help_info = prepare(args, paths)
        if args.prepare_only:
            return 0

        receipts_by_candidate = {}
        for candidate in candidates:
            receipts = []
            candidate_dir = paths["output"] / candidate["slug"]
            candidate_dir.mkdir(exist_ok=True)
            for seed in SEEDS:
                receipt = run_seed(candidate, seed, args, paths, help_info)
                receipts.append(receipt)
                append_progress(paths["output"], {
                    "utc": now(), "candidate_id": candidate["id"], "e3_type": candidate["e3"],
                    "seed": seed, "status": receipt.get("status"),
                    "completed_models": sum(
                        model.get("status") == "inspected" for model in receipt.get("models", [])
                    ),
                })
                verify_hashes(
                    paths["batch_root"], protocol["baseline"]["read_only_files_sha256"]
                )
                verify_hashes(candidate["frozen_result"], candidate["frozen_hashes"], {"msa", "processed"})
            receipts_by_candidate[candidate["id"]] = receipts

        summaries = []
        for candidate in candidates:
            receipts = receipts_by_candidate[candidate["id"]]
            valid_seed_set = successful_seed_receipts(receipts)
            summaries.append({
                "candidate_id": candidate["id"], "e3_type": candidate["e3"],
                "seeds": receipts,
                "actual_completion_count": sum(
                    r.get("status") == "completed" and r.get("actual_computation") is True
                    for r in receipts
                ),
                "complete_exact_seed_set": len(valid_seed_set) == 5,
                "selected_model_inspections": [
                    {"seed": r["seed"], "inspection": r.get("selected_model_inspection")}
                    for r in receipts
                ],
                "old_expert_quantitative_metrics": expert_metrics(receipts),
                "geometric_diagnostics": selected_diagnostics(candidate, receipts),
            })
        raw_models = sum(
            sum(model.get("status") == "inspected" for model in receipt.get("models", []))
            for receipts in receipts_by_candidate.values() for receipt in receipts
        )
        raw_artifacts = sum(
            sum(
                "raw_prediction_sha256" in model and "raw_confidence_sha256" in model
                for model in receipt.get("models", [])
            ) * 2
            for receipts in receipts_by_candidate.values() for receipt in receipts
        )
        complete = raw_models == 50 and raw_artifacts == 100 and all(
            item["actual_completion_count"] == 5 and item["complete_exact_seed_set"] is True
            for item in summaries
        )
        summary = {
            "format": "tpd-novel-msa-protocol-summary/1.0",
            "development_set_not_holdout": True,
            "scientific_approved": False,
            "automatic_approval": False,
            "protocol_digest": protocol["protocol_digest"],
            "candidate_count": 2, "seed_count_per_candidate": 5,
            "models_per_seed": 5, "raw_model_count": raw_models,
            "expected_raw_model_count": 50,
            "raw_model_artifact_file_count": raw_artifacts,
            "expected_raw_model_artifact_file_count": 100,
            "actual_completion_count": sum(item["actual_completion_count"] for item in summaries),
            "all_models_share_frozen_input_msa_settings": complete,
            "original_30_raw_baseline_preserved_read_only": True,
            "candidates": summaries,
            "protocol_complete": complete,
            "interpretation": (
                "New-MSA development evidence only. Geometric diagnostics are descriptive; "
                "no qualitative topology or clash result confers scientific approval."
            ),
        }
        write_json(paths["output"] / "protocol-summary.json", summary)
        verify_hashes(paths["batch_root"], protocol["baseline"]["read_only_files_sha256"])
        return 0 if complete else 1
    except KeyboardInterrupt:
        print("novel MSA protocol cancelled", file=sys.stderr)
        return 130
    except Exception as error:
        output = paths.get("output")
        if output is not None and output.is_dir() and not output.is_symlink():
            failure = output / "setup-failure.json"
            if not failure.exists():
                try:
                    write_json(failure, {
                        "failed_utc": now(), "type": type(error).__name__,
                        "reason": str(error) or type(error).__name__,
                        "fallback_attempted": False, "repeated_polling_attempted": False,
                    })
                except OSError:
                    pass
        print(f"novel MSA protocol error: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

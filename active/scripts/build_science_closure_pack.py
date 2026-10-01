#!/usr/bin/env python3
"""Build a verified, compact science-closure evidence bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SEEDS = (101, 127, 149, 173, 197)
FORBIDDEN = {"msa", "processed", "env", "weights", ".git"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_link_or_junction(path: Path) -> bool:
    if path.is_symlink():
        return True
    is_junction = getattr(path, "is_junction", None)
    if is_junction is None:
        return False
    try:
        return bool(is_junction())
    except OSError:
        return False


def safe_join(root: Path, relative: str | Path, *, must_exist: bool = True) -> Path:
    root = root.absolute()
    text = str(relative)
    rel = PurePosixPath(text.replace("\\", "/"))
    windows_rel = PureWindowsPath(text)
    if (rel.is_absolute() or windows_rel.is_absolute() or windows_rel.drive or
            not rel.parts or any(p in ("", ".", "..") or ":" in p for p in rel.parts)):
        raise ValueError(f"unsafe relative path: {relative}")
    candidate = root.joinpath(*rel.parts)
    if os.path.commonpath((str(root), str(candidate.absolute()))) != str(root):
        raise ValueError(f"path escapes root: {relative}")
    current = root
    if is_link_or_junction(current):
        raise ValueError(f"symlink/junction rejected: {current}")
    for part in rel.parts:
        current /= part
        if current.exists() and is_link_or_junction(current):
            raise ValueError(f"symlink/junction rejected: {current}")
    if must_exist and not candidate.is_file():
        raise ValueError(f"required regular file missing: {candidate}")
    return candidate


def verify_file(path: Path, expected_sha256: str, expected_bytes: int | None = None) -> None:
    if not path.is_file() or is_link_or_junction(path):
        raise ValueError(f"not a safe regular file: {path}")
    if expected_bytes is not None and path.stat().st_size != expected_bytes:
        raise ValueError(f"byte-size mismatch: {path}")
    actual = sha256_file(path)
    if actual != expected_sha256:
        raise ValueError(f"SHA-256 mismatch for {path}: {actual}")


def load_json(path: Path) -> Any:
    def object_without_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON object key: {key}")
            result[key] = value
        return result

    def finite_float(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            raise ValueError(f"non-finite JSON number: {token}")
        return value

    def reject_constant(token: str) -> Any:
        raise ValueError(f"non-finite JSON constant: {token}")

    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=object_without_duplicates,
                         parse_float=finite_float, parse_constant=reject_constant)


def dump_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")


def copy_verified(source: Path, destination: Path, expected_sha256: str | None = None) -> None:
    if not source.is_file() or is_link_or_junction(source):
        raise ValueError(f"unsafe copy source: {source}")
    if expected_sha256 is not None:
        verify_file(source, expected_sha256)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    if sha256_file(destination) != sha256_file(source):
        raise ValueError(f"copy verification failed: {destination}")


def selected_index(value: Any) -> int:
    if type(value) is not int:
        raise ValueError("trusted selector did not return an integer model index")
    return value


def receipt_selected_index(value: Any) -> int:
    if not isinstance(value, dict):
        raise ValueError("selection receipt must be an object")
    index = value.get("selected_model_index")
    if type(index) is not int:
        raise ValueError("selection receipt selected_model_index must be an integer")
    return index


def claimed_output(seed_root: Path, claimed: str) -> tuple[Path, str]:
    output = (seed_root / "boltz_output").absolute()
    text = str(claimed)
    posix_claim = PurePosixPath(text.replace("\\", "/"))
    windows_claim = PureWindowsPath(text)
    native_claim = Path(text)

    if native_claim.is_absolute():
        candidate = native_claim.absolute()
        try:
            relative = candidate.relative_to(output).as_posix()
        except ValueError as exc:
            raise ValueError(f"claimed path is outside this seed's boltz_output: {claimed}") from exc
    else:
        if windows_claim.is_absolute() or windows_claim.drive:
            raise ValueError(f"non-native or drive-qualified claimed path rejected: {claimed}")
        if (not posix_claim.parts or any(part in ("", ".", "..") or ":" in part
                                        for part in posix_claim.parts) or
                posix_claim.parts[0] != "boltz_output"):
            raise ValueError(f"invalid relative claimed output: {claimed}")
        relative = PurePosixPath(*posix_claim.parts[1:]).as_posix()
        if not relative or relative == ".":
            raise ValueError(f"invalid claimed output: {claimed}")
        candidate = seed_root.absolute().joinpath(*posix_claim.parts)

    source = safe_join(output, relative)
    if source.absolute() != candidate or source.resolve() != candidate.resolve():
        raise ValueError(f"claimed output does not resolve to the provided source path: {claimed}")
    try:
        source.resolve().relative_to(output.resolve())
    except ValueError as exc:
        raise ValueError(f"resolved claimed output escapes boltz_output: {claimed}") from exc
    return source, relative


def forbidden(relative: str) -> bool:
    return any(part.lower() in FORBIDDEN for part in PurePosixPath(relative).parts)


def iter_prediction_files(seed_root: Path) -> Iterable[tuple[Path, str]]:
    output = seed_root / "boltz_output"
    if not output.is_dir() or is_link_or_junction(output):
        raise ValueError(f"missing boltz_output: {output}")
    for directory, dirs, files in os.walk(output, followlinks=False):
        base = Path(directory)
        for name in list(dirs):
            child = base / name
            if is_link_or_junction(child):
                raise ValueError(f"symlink/junction rejected: {child}")
        for name in files:
            source = base / name
            rel = source.relative_to(output).as_posix()
            if "/predictions/" in f"/{rel}/":
                if forbidden(rel):
                    raise ValueError(f"forbidden prediction path: {rel}")
                yield safe_join(output, rel), rel


def locate_seed_root(holdout: Path, seed: int) -> Path:
    direct = holdout / f"seed-{seed}"
    if direct.is_dir():
        return direct
    matches = [p for p in holdout.rglob(f"seed-{seed}") if p.is_dir()]
    if len(matches) != 1:
        raise ValueError(f"expected one seed-{seed} directory, found {len(matches)}")
    relative = matches[0].relative_to(holdout).as_posix()
    safe_join(holdout, relative, must_exist=False)
    return matches[0]


def metric_fields(cutoffs: Any) -> list[str]:
    if not isinstance(cutoffs, dict) or len(cutoffs) != 4:
        raise ValueError("trusted CUTOFFS must contain exactly four metric fields")
    return list(cutoffs)


def assert_close_metrics(recorded: dict[str, Any], recomputed: dict[str, Any], fields: list[str]) -> None:
    for key in fields:
        if key not in recorded or key not in recomputed:
            raise ValueError(f"missing cutoff metric: {key}")
        recorded_value, recomputed_value = recorded[key], recomputed[key]
        if (isinstance(recorded_value, bool) or isinstance(recomputed_value, bool) or
                not isinstance(recorded_value, (int, float)) or
                not isinstance(recomputed_value, (int, float))):
            raise ValueError(f"metric must be a non-boolean number: {key}")
        if not math.isfinite(float(recorded_value)) or not math.isfinite(float(recomputed_value)):
            raise ValueError(f"metric must be finite: {key}")
        if abs(float(recorded_value) - float(recomputed_value)) > 1e-7:
            raise ValueError(f"recomputed metric mismatch: {key}")


def recursively_find(value: Any, key: str) -> list[Any]:
    found: list[Any] = []
    if isinstance(value, dict):
        for current, child in value.items():
            if current == key:
                found.append(child)
            found.extend(recursively_find(child, key))
    elif isinstance(value, list):
        for child in value:
            found.extend(recursively_find(child, key))
    return found


def validate_pka(pka: Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    manifest_path = safe_join(pka, "manifest.json")
    manifest = load_json(manifest_path)
    entries = manifest.get("files")
    if not isinstance(entries, list):
        raise ValueError("pKa manifest.json has no files list")
    expected = {"manifest.json"}
    for entry in entries:
        rel = str(entry["path"])
        source = safe_join(pka, rel)
        verify_file(source, str(entry["sha256"]), int(entry["bytes"]))
        expected.add(rel)
    actual = set()
    for directory, dirs, files in os.walk(pka, followlinks=False):
        base = Path(directory)
        for name in dirs:
            if is_link_or_junction(base / name):
                raise ValueError(f"pKa symlink/junction rejected: {base / name}")
        for name in files:
            source = base / name
            if is_link_or_junction(source):
                raise ValueError(f"pKa symlink/junction rejected: {source}")
            actual.add(source.relative_to(pka).as_posix())
    if actual != expected:
        raise ValueError(f"pKa folder differs from manifest: extra={actual-expected}, missing={expected-actual}")
    compact = load_json(safe_join(pka, "compact-summary.json"))
    if manifest.get("scientific_approved") is not False or compact.get("scientific_approved") is not False:
        raise ValueError("pKa evidence must remain scientifically unapproved")
    if manifest.get("population_prediction_performed") is not False or compact.get("population_prediction_performed") is not False:
        raise ValueError("pKa population prediction must be false")
    analogs = compact.get("analog_results")
    if not isinstance(analogs, (list, dict)) or len(analogs) != 3:
        raise ValueError("expected exactly three pKa analog results")
    raw = load_json(safe_join(pka, "raw-predictions.json"))
    predictions = raw.get("predictions") if isinstance(raw, dict) else None
    if not isinstance(predictions, dict) or not predictions:
        raise ValueError("raw.predictions must be an analog_id-to-record dictionary")
    target_count = 0
    for analog_id, prediction in predictions.items():
        if not isinstance(analog_id, str) or not isinstance(prediction, dict):
            raise ValueError("invalid raw prediction record")
        sites = prediction.get("sites")
        if not isinstance(sites, list):
            raise ValueError(f"raw prediction sites must be a list: {analog_id}")
        target_count += len(sites)
    if target_count != 26:
        raise ValueError(f"pKa raw predictions contain {target_count} sites, expected 26")
    missing = compact.get("missing_required_maps")
    if not isinstance(missing, dict) or "W-80f8f4a11b5d" not in missing:
        raise ValueError("missing_required_maps must contain W-80f8f4a11b5d")
    target_missing = missing["W-80f8f4a11b5d"]
    if not isinstance(target_missing, list) or 5001 not in target_missing or "N5001" in target_missing:
        raise ValueError("W-80f8f4a11b5d must record missing integer atom map 5001")
    for analog_id, maps in missing.items():
        if not isinstance(maps, list) or any(type(item) is not int for item in maps):
            raise ValueError(f"missing map values must be integer lists: {analog_id}")
        if analog_id != "W-80f8f4a11b5d" and maps:
            raise ValueError(f"unexpected missing required maps for {analog_id}")
    return manifest, compact, manifest_path


def build(args: argparse.Namespace) -> None:
    from scripts import run_calibration_seed_holdout as h

    holdout, pka, reference, output = map(lambda p: Path(p).absolute(),
                                          (args.holdout_dir, args.pka_dir, args.reference, args.output))
    if output.exists():
        raise ValueError("--output must not exist")
    for root in (holdout, pka):
        if not root.is_dir() or is_link_or_junction(root):
            raise ValueError(f"unsafe input directory: {root}")
    if not reference.is_file() or is_link_or_junction(reference):
        raise ValueError("unsafe reference file")

    plan_path = safe_join(holdout, "plan.json")
    plan_sha_path = safe_join(holdout, "plan.sha256")
    summary_path = safe_join(holdout, "protocol-summary.json")
    plan_hash = sha256_file(plan_path)
    token = plan_sha_path.read_text(encoding="utf-8").split()[0]
    if token != plan_hash:
        raise ValueError("plan.sha256 first token does not match plan.json")
    plan, summary = load_json(plan_path), load_json(summary_path)
    if tuple(h.SEEDS) != SEEDS or tuple(plan.get("seeds", ())) != tuple(h.SEEDS):
        raise ValueError("runner, builder, and preregistered plan seeds differ")
    expected_scope = {
        "holdout_seed_validation": True,
        "independent_target_or_structure": False,
        "known_structure_training_overlap_possible": True,
        "same_6BOY_development_structure": True,
        "training_holdout": False,
    }
    scope = plan.get("scope")
    if (
        type(scope) is not dict
        or set(scope) != set(expected_scope)
        or any(type(scope[key]) is not bool or scope[key] is not expected_scope[key] for key in expected_scope)
        or summary.get("scope") != scope
    ):
        raise ValueError("plan and protocol-summary scope differ or are missing")
    if plan.get("preregistered_before_inference") is not True:
        raise ValueError("preregistered protocol flag is missing or false")
    expected_raw_units = [
        {
            "model_index": model_index,
            "seed": seed,
            "unit_id": f"seed-{seed}/model-{model_index}",
        }
        for seed in SEEDS
        for model_index in range(5)
    ]
    raw_units = plan.get("all_preregistered_raw_units")
    if (
        type(raw_units) is not list
        or raw_units != expected_raw_units
        or any(
            type(unit) is not dict
            or type(unit.get("model_index")) is not int
            or type(unit.get("seed")) is not int
            or type(unit.get("unit_id")) is not str
            for unit in raw_units
        )
    ):
        raise ValueError("all_preregistered_raw_units is missing or invalid")
    expected_reference = plan.get("sources", {}).get("reference", {}).get("sha256")
    if expected_reference != sha256_file(reference):
        raise ValueError("reference SHA-256 differs from plan.sources.reference.sha256")
    if summary.get("cutoffs") != h.CUTOFFS:
        raise ValueError("cutoffs changed from trusted preregistered cutoffs")
    if summary.get("selection_receipts_written_before_comparison") is not True:
        raise ValueError("selection receipts were not recorded before comparison")

    fields = metric_fields(h.CUTOFFS)
    rows, selected_passes, baseline_passes, raw_passes = [], 0, 0, 0
    staging = Path(tempfile.mkdtemp(prefix=output.name + ".", dir=output.parent))
    try:
        holdout_out = staging / "known_crbn_holdout"
        for source, rel in ((plan_path, "plan.json"), (plan_sha_path, "plan.sha256"),
                            (summary_path, "protocol-summary.json")):
            copy_verified(source, holdout_out / rel)
        copy_verified(reference, holdout_out / "reference" / reference.name, expected_reference)

        for seed in SEEDS:
            root = locate_seed_root(holdout, seed)
            receipt_path = safe_join(root, "receipt.json")
            selection_path = safe_join(root, "selection-receipt.json")
            log_path = safe_join(root, "inference.log")
            yaml_path = safe_join(root, "6BOY_CRBN.yaml")
            receipt = load_json(receipt_path)
            if receipt.get("seed") != seed or receipt.get("status") != "success":
                raise ValueError(f"seed {seed} is not a successful receipt")
            models = receipt.get("models")
            if (not isinstance(models, list) or
                    any(not isinstance(model, dict) or type(model.get("model_index")) is not int
                        for model in models) or
                    sorted(model["model_index"] for model in models) != list(range(5))):
                raise ValueError(f"seed {seed} must contain strict integer model IDs 0..4")
            selected = selected_index(h.select_model(models))
            baseline = selected_index(h.select_confidence_baseline(models))
            if selected != receipt_selected_index(receipt.get("selection")):
                raise ValueError(f"seed {seed} trusted selection differs from receipt")
            selection_doc = load_json(selection_path)
            if receipt_selected_index(selection_doc) != selected:
                raise ValueError(f"seed {seed} selection-receipt differs from trusted selection")
            hashes = receipt.get("output_files_sha256")
            if not isinstance(hashes, dict):
                raise ValueError(f"seed {seed} lacks output file hashes")
            if sha256_file(log_path) != receipt.get("inference_log_sha256"):
                raise ValueError(f"seed {seed} inference log hash mismatch")

            prediction_sources = list(iter_prediction_files(root))
            existing_prediction_relatives = {rel for _, rel in prediction_sources}
            hashed_prediction_relatives = {
                str(rel) for rel in hashes
                if "/predictions/" in f"/{str(rel).replace(chr(92), '/')}/"
            }
            if existing_prediction_relatives != hashed_prediction_relatives:
                raise ValueError(
                    f"seed {seed} prediction hash inventory mismatch: "
                    f"missing={hashed_prediction_relatives-existing_prediction_relatives}, "
                    f"unhashed={existing_prediction_relatives-hashed_prediction_relatives}"
                )
            suffix_counts = {
                ".cif": sum(rel.lower().endswith(".cif") for rel in existing_prediction_relatives),
                ".json": sum(rel.lower().endswith(".json") for rel in existing_prediction_relatives),
                ".npz": sum(rel.lower().endswith(".npz") for rel in existing_prediction_relatives),
            }
            if len(existing_prediction_relatives) != 25 or suffix_counts != {".cif": 5, ".json": 5, ".npz": 15}:
                raise ValueError(f"seed {seed} must preserve exactly 5 CIF, 5 confidence JSON, and 15 NPZ files")

            selected_metrics = None
            selected_pass = None
            for model in models:
                index = model["model_index"]
                prediction, prediction_rel = claimed_output(root, model["prediction_path"])
                confidence, confidence_rel = claimed_output(root, model["confidence_path"])
                verify_file(prediction, str(model["prediction_sha256"]))
                verify_file(confidence, str(model["confidence_sha256"]))
                if hashes.get(prediction_rel) != model["prediction_sha256"] or hashes.get(confidence_rel) != model["confidence_sha256"]:
                    raise ValueError(f"seed {seed} output_files_sha256 lacks exact relative keys")
                if load_json(confidence) != model.get("confidence_raw"):
                    raise ValueError(f"seed {seed} model {index} confidence_raw mismatch")
                recomputed = h.worker.compare_structures(reference, prediction, receipt["input"])["metrics"]
                assert_close_metrics(model["metrics"], recomputed, fields)
                passed = bool(h.metrics_pass(recomputed))
                if passed != bool(model.get("passes_cutoffs")):
                    raise ValueError(f"seed {seed} model {index} pass flag mismatch")
                raw_passes += int(passed)
                if index == selected:
                    selected_metrics = {key: recomputed[key] for key in fields}
                    selected_pass = passed
                    selected_passes += int(passed)
                if index == baseline:
                    baseline_passes += int(passed)
            if selected_metrics is None:
                raise ValueError(f"seed {seed} selected model absent")
            rows.append({"seed": seed, "selected_model_index": selected,
                         "confidence_baseline_model_index": baseline, "metrics": selected_metrics,
                         "passes": selected_pass})

            seed_out = holdout_out / f"seed-{seed}"
            for source, rel in ((receipt_path, "receipt.json"), (selection_path, "selection-receipt.json"),
                                (log_path, "inference.log"), (yaml_path, "6BOY_CRBN.yaml")):
                copy_verified(source, seed_out / rel)
            copied_relatives = set()
            for source, rel in prediction_sources:
                expected = hashes.get(rel)
                if not isinstance(expected, str):
                    raise ValueError(f"unhashed prediction output: seed {seed} {rel}")
                copy_verified(source, seed_out / "boltz_output" / rel, expected)
                copied_relatives.add(rel)
            if copied_relatives != hashed_prediction_relatives:
                raise ValueError(f"seed {seed} did not copy the complete prediction inventory")

        expected_summary = {
            "successful_seed_count": 5, "denominator_including_failed_units": 5,
            "selected_complex_iplddt_pass_count": selected_passes,
            "confidence_score_baseline_pass_count": baseline_passes,
        }
        for key, value in expected_summary.items():
            if summary.get(key) != value:
                raise ValueError(f"protocol summary derived count mismatch: {key}")
        if summary.get("failed_unit_count", 0) != 0 or summary.get("all_25_raw_records_preserved") is not True:
            raise ValueError("protocol summary does not preserve all 25 successful raw records")

        pka_manifest, pka_compact, pka_manifest_path = validate_pka(pka)
        pka_out = staging / "ligand_pka"
        for entry in pka_manifest["files"]:
            rel = str(entry["path"])
            copy_verified(safe_join(pka, rel), pka_out / rel, str(entry["sha256"]))
        copy_verified(pka_manifest_path, pka_out / "manifest.json")

        compact = {
            "format": "science-closure-pack/1",
            "scientific_approved": False,
            "human_review_performed": False,
            "known_crbn_holdout": {
                "selected_per_seed_metrics": rows,
                "raw_count": 25,
                "raw_pass_count": raw_passes,
                "selected_pass_count": selected_passes,
                "baseline_pass_count": baseline_passes,
                "plan_sha256": plan_hash,
                "original_summary_sha256": sha256_file(summary_path),
                "scope": summary.get("scope"),
            },
            "ligand_pka": pka_compact,
        }
        dump_json(staging / "compact-summary.json", compact)
        if (staging / "compact-summary.json").stat().st_size > 25 * 1024:
            raise ValueError("compact-summary.json exceeds 25 KB")

        readme = f"""# 과학적 클로저 증거 묶음

이 묶음은 동일한 6BOY BRD4-CRBN-dBET6 개발 구조에 대한 제한된 재현 증거입니다. 원래 5개 시드에서 발견된 규칙을 고정한 뒤, 사전등록된 선택 및 컷오프 규칙으로 새로운 시드 {', '.join(map(str, SEEDS))}만 평가했습니다.

- 성공 시드: 5/5
- 원시 모델 기록: 25, 통과: {raw_passes}
- complex_iplddt 선택 모델 통과: {selected_passes}/5
- confidence_score 기준선 선택 모델 통과: {baseline_passes}/5
- 사후 선택이나 컷오프 변경 없음

이 결과는 독립 표적 검증, 독립 구조 검증, 학습 데이터 홀드아웃 또는 다른 방법보다 우월하다는 주장이 아닙니다. pKa 값은 분리된 warhead에 대한 예측값이며 전체 PROTAC 또는 단백질 결합 상태의 pKa가 아닙니다. 결합된 미세상태 population은 계산하지 않았고 W-80f8f4a11b5d의 말단 질소 원자 맵 5001이 누락되었습니다. 이러한 제한과 필요한 다른 검증 요건이 남아 있고 사람의 과학 검토도 수행되지 않았으므로 scientific_approved=false입니다.

정확한 GPU 재실행에는 별도로 동결된 캐시와 모델 가중치가 필요하며 이 묶음에는 배포되지 않습니다. 명령 문자열은 증거로만 보존되며 실행되지 않습니다.
"""
        (staging / "README.md").write_text(readme, encoding="utf-8")

        files = []
        for source in sorted(p for p in staging.rglob("*") if p.is_file()):
            rel = source.relative_to(staging).as_posix()
            files.append({"path": rel, "bytes": source.stat().st_size, "sha256": sha256_file(source)})
        manifest = {
            "format": "science-closure-pack-manifest/1",
            "scientific_approved": False,
            "human_review_performed": False,
            "code_sha256": sha256_file(Path(__file__)),
            "files": files,
        }
        dump_json(staging / "output-manifest.json", manifest)
        os.replace(staging, output)
    except Exception:
        resolved_parent = output.parent.resolve()
        safe_staging = False
        try:
            resolved_staging = staging.resolve()
            safe_staging = (resolved_staging.parent == resolved_parent and
                            staging.name.startswith(output.name + "."))
        except OSError:
            resolved_staging = staging
        if safe_staging:
            shutil.rmtree(resolved_staging, ignore_errors=True)
        else:
            failure_path = output.parent / f"{output.name}.cleanup-failure.txt"
            failure_path.write_text(f"staging preserved; unsafe cleanup path: {staging}\n", encoding="utf-8")
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--holdout-dir", required=True)
    parser.add_argument("--pka-dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    build(parser.parse_args())


if __name__ == "__main__":
    main()

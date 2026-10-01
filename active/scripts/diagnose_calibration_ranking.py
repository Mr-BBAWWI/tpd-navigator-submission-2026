"""Posthoc development-set diagnosis of fixed calibration ranking rules.

This tool never updates acceptance, protocol_pass, or claims independent validation.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_calibration_protocol import CUTOFFS, metrics_pass

SEEDS = [23, 41, 61, 79, 97]
INDICES = list(range(5))
FIXED_SETTINGS = {
    "recycling_steps": 10,
    "sampling_steps": 200,
    "diffusion_samples": 5,
    "msa_mode": "cached_processed_offline",
    "max_msa_seqs": 256,
}
RULES = [
    ("confidence_score_max", "confidence_score", "max", "confidence"),
    ("iptm_max", "iptm", "max", "confidence"),
    ("protein_iptm_max", "protein_iptm", "max", "confidence"),
    ("ligand_iptm_max", "ligand_iptm", "max", "confidence"),
    ("complex_iplddt_max", "complex_iplddt", "max", "confidence"),
    ("complex_ipde_min", "complex_ipde", "min", "confidence"),
    ("clash_count_min", "prediction_descriptive_clashes", "min", "metrics"),
]
UNIT_FIELDS = {
    "confidence_score", "iptm", "protein_iptm", "ligand_iptm", "complex_iplddt"
}
METRIC_KEYS = list(CUTOFFS)


class DiagnosisError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise DiagnosisError("JSON_READ_FAILED:" + str(path)) from error


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")


def finite_number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DiagnosisError("FINITE_NUMERIC_FIELD_REQUIRED:" + field)
    result = float(value)
    if not math.isfinite(result):
        raise DiagnosisError("FINITE_NUMERIC_FIELD_REQUIRED:" + field)
    if field in UNIT_FIELDS and not 0.0 <= result <= 1.0:
        raise DiagnosisError("FIELD_OUT_OF_NATURAL_RANGE:" + field)
    if field in {"complex_ipde", "prediction_descriptive_clashes"} and result < 0.0:
        raise DiagnosisError("FIELD_OUT_OF_NATURAL_RANGE:" + field)
    return result


def _complete_models(models: list[dict]) -> dict[int, dict]:
    if len(models) != 5:
        raise DiagnosisError("ALL_FIVE_MODELS_REQUIRED")
    result: dict[int, dict] = {}
    for model in models:
        index = model.get("model_index")
        if isinstance(index, bool) or not isinstance(index, int) or index not in INDICES:
            raise DiagnosisError("MODEL_INDEX_SET_INVALID")
        if index in result:
            raise DiagnosisError("MODEL_INDEX_SET_INVALID")
        result[index] = model
    if set(result) != set(INDICES):
        raise DiagnosisError("MODEL_INDEX_SET_INVALID")
    return result


def clash_count(value: Any) -> float:
    if not isinstance(value, list):
        raise DiagnosisError("CLASH_RECORD_LIST_REQUIRED:prediction_descriptive_clashes")
    return finite_number(len(value), "prediction_descriptive_clashes")


def select_rule(models: list[dict], rule_name: str) -> int:
    """Select solely from the fixed allowed prediction-derived fields."""
    model_map = _complete_models(models)
    specifications = {item[0]: item for item in RULES}
    if rule_name not in specifications:
        raise DiagnosisError("UNKNOWN_RULE:" + rule_name)
    _, field, direction, source = specifications[rule_name]
    candidates = []
    for index in INDICES:
        model = model_map[index]
        container = model.get(source)
        if not isinstance(container, dict) or field not in container:
            raise DiagnosisError("RULE_FIELD_UNAVAILABLE:" + field)
        primary = (clash_count(container[field]) if field == "prediction_descriptive_clashes"
                   else finite_number(container[field], field))
        confidence = model.get("confidence")
        if not isinstance(confidence, dict) or "confidence_score" not in confidence:
            raise DiagnosisError("RULE_FIELD_UNAVAILABLE:confidence_score")
        tie = finite_number(confidence["confidence_score"], "confidence_score")
        candidates.append((primary, tie, index))
    if direction == "max":
        return min(candidates, key=lambda item: (-item[0], -item[1], item[2]))[2]
    return min(candidates, key=lambda item: (item[0], -item[1], item[2]))[2]


def _require_hash(model: dict, path_key: str, hash_key: str) -> tuple[Path, str]:
    raw_path = model.get(path_key)
    expected = model.get(hash_key)
    if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
        raise DiagnosisError("ABSOLUTE_RAW_PATH_REQUIRED:" + path_key)
    if not isinstance(expected, str) or len(expected) != 64:
        raise DiagnosisError("MISSING_OR_INVALID_HASH:" + hash_key)
    path = Path(raw_path)
    if not path.is_file() or path.is_symlink():
        raise DiagnosisError("RAW_REGULAR_FILE_REQUIRED:" + path_key)
    return path, expected.lower()


def _verify_source(source: Any, label: str) -> None:
    if not isinstance(source, dict):
        raise DiagnosisError("FROZEN_SOURCE_MISSING:" + label)
    path_value, expected = source.get("path"), source.get("sha256")
    if not isinstance(path_value, str) or not Path(path_value).is_absolute():
        raise DiagnosisError("FROZEN_SOURCE_PATH_INVALID:" + label)
    if not isinstance(expected, str) or len(expected) != 64:
        raise DiagnosisError("FROZEN_SOURCE_HASH_MISSING:" + label)
    path = Path(path_value)
    if not path.is_file() or path.is_symlink() or sha256_file(path) != expected.lower():
        raise DiagnosisError("FROZEN_SOURCE_CHANGED:" + label)


def load_verified_protocol(protocol_root: Path) -> tuple[dict, list[dict]]:
    plan = read_json(protocol_root / "plan.json")
    summary = read_json(protocol_root / "protocol-summary.json")
    if plan.get("seeds") != SEEDS or summary.get("seed_count") != 5:
        raise DiagnosisError("EXACT_SEED_SET_REQUIRED")
    if plan.get("expected_models_per_seed") != INDICES:
        raise DiagnosisError("EXACT_MODEL_SET_REQUIRED")
    settings = plan.get("effective_protocol_settings")
    if not isinstance(settings, dict) or any(settings.get(k) != v for k, v in FIXED_SETTINGS.items()):
        raise DiagnosisError("FIXED_PROTOCOL_SETTINGS_CHANGED")
    sources = plan.get("sources")
    if not isinstance(sources, dict):
        raise DiagnosisError("PLAN_SOURCES_MISSING")
    _verify_source(sources.get("checkpoint"), "checkpoint")
    _verify_source(sources.get("frozen_source_receipt"), "fixed_msa_receipt")

    summary_seeds = summary.get("seeds")
    if not isinstance(summary_seeds, list) or len(summary_seeds) != 5:
        raise DiagnosisError("EXACT_SEED_SET_REQUIRED")
    by_seed = {item.get("seed"): item for item in summary_seeds if isinstance(item, dict)}
    if set(by_seed) != set(SEEDS):
        raise DiagnosisError("EXACT_SEED_SET_REQUIRED")

    receipts = []
    pending_hashes: list[tuple[Path, str, str]] = []
    for seed in SEEDS:
        receipt = read_json(protocol_root / f"seed-{seed}" / "receipt.json")
        if receipt.get("seed") != seed or receipt.get("status") != "success":
            raise DiagnosisError("SUCCESSFUL_SEED_RECEIPT_REQUIRED")
        if receipt.get("frozen_post_run_validation") != "success":
            raise DiagnosisError("FIXED_MSA_POSTRUN_VALIDATION_REQUIRED")
        receipt_models = receipt.get("models")
        summary_models = by_seed[seed].get("models")
        _complete_models(receipt_models if isinstance(receipt_models, list) else [])
        if receipt_models != summary_models:
            raise DiagnosisError("SUMMARY_RECEIPT_MODEL_MISMATCH")
        for model in receipt_models:
            prediction, prediction_hash = _require_hash(
                model, "prediction_path", "prediction_sha256"
            )
            confidence, confidence_hash = _require_hash(
                model, "confidence_path", "confidence_sha256"
            )
            pending_hashes.extend([
                (prediction, prediction_hash, "prediction"),
                (confidence, confidence_hash, "confidence"),
            ])
        receipts.append(receipt)

    # Every raw prediction and confidence artifact is verified before any score is used.
    for path, expected, label in pending_hashes:
        if sha256_file(path) != expected:
            raise DiagnosisError("RAW_HASH_MISMATCH:" + label + ":" + str(path))

    for receipt in receipts:
        for model in receipt["models"]:
            raw_confidence = read_json(Path(model["confidence_path"]))
            recorded_raw = model.get("confidence_raw")
            if recorded_raw is not None and raw_confidence != recorded_raw:
                raise DiagnosisError("RAW_CONFIDENCE_RECEIPT_MISMATCH")
            if not isinstance(model.get("confidence"), dict):
                raise DiagnosisError("VALIDATED_CONFIDENCE_REQUIRED")
            if not isinstance(model.get("metrics"), dict):
                raise DiagnosisError("ALL_MODEL_METRICS_REQUIRED")

    baseline = plan.get("preserved_original_baseline")
    if baseline != {"geometric_successes": 2, "seed_count": 5}:
        raise DiagnosisError("PRESERVED_BASELINE_CHANGED")
    if summary.get("selected_pass_count") != 2 or summary.get("all_25_model_metrics_preserved") is not True:
        raise DiagnosisError("PRESERVED_TOP1_OR_ALL25_CHANGED")
    return plan, receipts


def proposals(receipts: list[dict]) -> dict:
    selections: dict[str, Any] = {}
    for name, field, direction, source in RULES:
        indices: dict[str, Any] = {}
        error = None
        for receipt in receipts:
            try:
                indices[str(receipt["seed"])] = select_rule(receipt["models"], name)
            except DiagnosisError as exc:
                indices[str(receipt["seed"])] = None
                error = str(exc)
        selections[name] = {
            "field": field,
            "direction": direction,
            "source": source,
            "tie_break": "higher confidence_score, then lower model index",
            "selection_indices": indices,
            "valid": error is None,
            "invalid_reason": error,
        }
    return {
        "format": "calibration-ranking-posthoc-proposals/1.0",
        "candidate_rules_fixed_before_reference_evaluation": True,
        "development_set_only": True,
        "no_fitted_weights": True,
        "no_seed_selection": True,
        "no_reference_dependent_selection": True,
        "rules": selections,
    }


def metric_breakdown(metrics: dict) -> dict[str, bool]:
    values = {key: finite_number(metrics.get(key), key) for key in METRIC_KEYS}
    return {
        "target_CA_RMSD_A": values["target_CA_RMSD_A"] <= CUTOFFS["target_CA_RMSD_A"],
        "ligand_heavy_atom_RMSD_after_target_alignment_A":
            values["ligand_heavy_atom_RMSD_after_target_alignment_A"] <=
            CUTOFFS["ligand_heavy_atom_RMSD_after_target_alignment_A"],
        "e3_CA_RMSD_after_target_alignment_A":
            values["e3_CA_RMSD_after_target_alignment_A"] <=
            CUTOFFS["e3_CA_RMSD_after_target_alignment_A"],
        "contact_jaccard": values["contact_jaccard"] >= CUTOFFS["contact_jaccard"],
    }


def evaluate(receipts: list[dict], proposal: dict) -> dict:
    evaluations: dict[str, Any] = {}
    for rule_name, rule in proposal["rules"].items():
        rows = []
        if rule["valid"]:
            for receipt in receipts:
                index = rule["selection_indices"][str(receipt["seed"])]
                model = next(x for x in receipt["models"] if x["model_index"] == index)
                metrics = {key: model["metrics"].get(key) for key in METRIC_KEYS}
                rows.append({
                    "seed": receipt["seed"],
                    "selected_model_index": index,
                    "metrics": metrics,
                    "passes_each_cutoff": metric_breakdown(metrics),
                    "all_four_pass": metrics_pass(metrics),
                })
        evaluations[rule_name] = {
            "specificity": "posthoc development-set diagnostic only",
            "selected_pass_count": (
                sum(row["all_four_pass"] for row in rows) if rule["valid"] else None
            ),
            "seed_count": 5,
            "selections": rows,
            "valid": rule["valid"],
            "error": rule["invalid_reason"],
        }

    geometry = []
    for receipt in receipts:
        models = receipt["models"]
        ordered = sorted(
            models,
            key=lambda model: (
                -finite_number(model["confidence"]["confidence_score"], "confidence_score"),
                model["model_index"],
            ),
        )
        rank = {model["model_index"]: position + 1 for position, model in enumerate(ordered)}
        top_index = ordered[0]["model_index"]
        model_rows = []
        for model in sorted(models, key=lambda item: item["model_index"]):
            breakdown = metric_breakdown(model["metrics"])
            model_rows.append({
                "model_index": model["model_index"],
                "confidence_order_label": "rank_" + str(rank[model["model_index"]]),
                "confidence_score": model["confidence"]["confidence_score"],
                "metric_pass_count": sum(breakdown.values()),
                "all_four_pass": metrics_pass(model["metrics"]),
            })
        geometry.append({
            "seed": receipt["seed"],
            "total_passing_models": sum(row["all_four_pass"] for row in model_rows),
            "top1_model_index": top_index,
            "top1_pass": next(row["all_four_pass"] for row in model_rows
                              if row["model_index"] == top_index),
            "oracle_available": any(row["all_four_pass"] for row in model_rows),
            "oracle_warning": "DIAGNOSTIC ONLY; oracle availability is unavailable in real prediction and is not a recommendation.",
            "models": model_rows,
        })

    baseline_count = evaluations["confidence_score_max"]["selected_pass_count"]
    if baseline_count != 2:
        raise DiagnosisError("BASELINE_2_OF_5_NOT_PRESERVED")
    return {
        "format": "calibration-ranking-posthoc-diagnosis/1.0",
        "development_set_not_independent_validation": True,
        "preserved_baseline": "2/5 top-1 seeds; all 25 model metrics preserved",
        "evaluations": evaluations,
        "selection_bottleneck_geometry": geometry,
        "recommendation": "No winner is recommended from highest post-test pass count. Promising hypotheses require preregistration and new independent cases.",
        "limitations": [
            "This is a posthoc development-set analysis.",
            "Reference metrics are used only after proposal indices are persisted.",
            "Oracle availability cannot be used for real prediction or recommendation.",
            "No acceptance state or original protocol result is updated.",
            "Zero results are approved and no independent-validation claim is made.",
        ],
    }


def render_html(receipts: list[dict], proposal: dict, diagnosis: dict) -> str:
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>Calibration ranking diagnosis</title>",
        "<style>body{font-family:sans-serif}table{border-collapse:collapse;margin:1em 0}",
        "th,td{border:1px solid #999;padding:.3em;text-align:right}th:first-child,td:first-child{text-align:left}</style>",
        "</head><body><h1>Posthoc development-set calibration ranking diagnosis</h1>",
        "<p><strong>Zero results approved.</strong> This is not independent validation and does not update acceptance or the frozen protocol result.</p>",
        "<h2>All confidence values</h2><table><tr><th>Seed/model</th>",
    ]
    confidence_fields = [item[1] for item in RULES if item[3] == "confidence"]
    parts.extend("<th>" + html.escape(field) + "</th>" for field in confidence_fields)
    parts.append("<th>prediction_descriptive_clashes</th></tr>")
    for receipt in receipts:
        for model in sorted(receipt["models"], key=lambda item: item["model_index"]):
            parts.append(f"<tr><td>{receipt['seed']}/{model['model_index']}</td>")
            for field in confidence_fields:
                parts.append("<td>" + html.escape(str(model["confidence"].get(field))) + "</td>")
            clashes = model["metrics"].get("prediction_descriptive_clashes")
            clash_display = len(clashes) if isinstance(clashes, list) else "invalid"
            parts.append("<td>" + html.escape(str(clash_display)) + "</td></tr>")
    parts.append("</table><h2>Selected cutoff metrics</h2>")
    for rule_name, result in diagnosis["evaluations"].items():
        parts.append("<h3>" + html.escape(rule_name) + "</h3><table><tr><th>Seed/model</th>")
        parts.extend("<th>" + html.escape(key) + "</th>" for key in METRIC_KEYS)
        parts.append("<th>all four pass</th></tr>")
        for row in result["selections"]:
            parts.append(f"<tr><td>{row['seed']}/{row['selected_model_index']}</td>")
            parts.extend("<td>" + html.escape(str(row["metrics"][key])) + "</td>"
                         for key in METRIC_KEYS)
            parts.append("<td>" + str(row["all_four_pass"]) + "</td></tr>")
        parts.append("</table>")
    parts.append("<h2>Limitations</h2><ul>")
    parts.extend("<li>" + html.escape(item) + "</li>" for item in diagnosis["limitations"])
    parts.append("</ul></body></html>\n")
    return "".join(parts)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--protocol-root", required=True)
    value.add_argument("--output", required=True)
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    protocol_root = Path(args.protocol_root).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    try:
        if output.exists():
            raise DiagnosisError("OUTPUT_ALREADY_EXISTS")
        if not protocol_root.is_dir() or protocol_root.is_symlink():
            raise DiagnosisError("PROTOCOL_ROOT_DIRECTORY_REQUIRED")
        plan, receipts = load_verified_protocol(protocol_root)
        provenance = {
            "plan": {
                "path": str(protocol_root / "plan.json"),
                "sha256": sha256_file(protocol_root / "plan.json"),
            },
            "summary": {
                "path": str(protocol_root / "protocol-summary.json"),
                "sha256": sha256_file(protocol_root / "protocol-summary.json"),
            },
            "receipts": [
                {
                    "seed": receipt["seed"],
                    "path": str(protocol_root / f"seed-{receipt['seed']}" / "receipt.json"),
                    "sha256": sha256_file(
                        protocol_root / f"seed-{receipt['seed']}" / "receipt.json"
                    ),
                    "raw_hashes": [
                        {
                            "model_index": model["model_index"],
                            "prediction_sha256": model["prediction_sha256"],
                            "confidence_sha256": model["confidence_sha256"],
                        }
                        for model in receipt["models"]
                    ],
                }
                for receipt in receipts
            ],
        }
        output.mkdir(parents=True, exist_ok=False)
        write_json_exclusive(output / "source-provenance.json", provenance)
        fixed_proposals = proposals(receipts)
        # This durable write must precede all reference-cutoff evaluation.
        write_json_exclusive(output / "ranking-proposals.json", fixed_proposals)
        diagnosis = evaluate(receipts, fixed_proposals)
        write_json_exclusive(output / "diagnosis.json", diagnosis)
        with (output / "diagnosis.html").open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(render_html(receipts, fixed_proposals, diagnosis))
        return 0
    except (OSError, DiagnosisError) as error:
        print("calibration ranking diagnosis error: " + type(error).__name__ + ": " + str(error),
              file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

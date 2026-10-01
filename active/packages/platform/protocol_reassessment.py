"""Fast, fail-closed adaptation of immutable protocol diagnostics into assessments."""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

from packages.contracts import ContractError, encoded
from packages.science import scientific_assessment

ROOT = Path(__file__).resolve().parents[2]
MAX_PROTOCOL_IDS = 8
KNOWN_KIND = "known_structure_calibration"
NOVEL_KIND = "novel_reference_free_msa"
SUPPORTED_KINDS = {KNOWN_KIND, NOVEL_KIND}
KNOWN_SELECTED_DENOMINATOR = 5
KNOWN_RAW_DENOMINATOR = 25
KNOWN_REQUIRED_PASS_COUNT = 3
NOVEL_CANDIDATE_DENOMINATOR = 2
NOVEL_SELECTED_PER_CANDIDATE = 5
NOVEL_SELECTED_DENOMINATOR = 10
NOVEL_RAW_PER_CANDIDATE = 25
NOVEL_RAW_DENOMINATOR = 50
NOVEL_IPTM_IQR_MAX = 0.20
NOVEL_ENDPOINT_IQR_MAX_A = 1.5


def _require(condition: bool, code: str) -> None:
    if not condition:
        raise ContractError(code)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any) -> str:
    return _sha(encoded(value))


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value, ensure_ascii=False, allow_nan=False))


def _strict_int(value: Any, code: str, *, minimum: int = 0) -> int:
    _require(type(value) is int and value >= minimum, code)
    return value


def _finite_number(value: Any, code: str) -> float:
    _require(type(value) in (int, float) and math.isfinite(float(value)), code)
    return float(value)


def _strict_bool(value: Any, code: str) -> bool:
    _require(type(value) is bool, code)
    return value


def current_verification_source_hashes() -> dict[str, str]:
    paths = {
        "packages/platform/protocol_evidence.py": ROOT / "packages/platform/protocol_evidence.py",
        "packages/science/protocol_evidence.py": ROOT / "packages/science/protocol_evidence.py",
        "contracts/protocol-evidence/v0.1.0/schema.json":
            ROOT / "contracts/protocol-evidence/v0.1.0/schema.json",
    }
    result = {}
    for name, path in paths.items():
        _require(path.is_file() and not path.is_symlink(),
                 "PROTOCOL_REASSESSMENT_VERIFICATION_SOURCE_MISSING")
        result[name] = _sha(path.read_bytes())
    return result


class _CachedStore:
    """Read-through cache so a shared archive is read from disk only once."""

    def __init__(self, store):
        self._store = store
        self.root = store.root
        self._cache: dict[tuple[str, str, str], bytes] = {}

    def read(self, project: str, ref: dict) -> bytes:
        key = (project, ref.get("artifact_id", ""), ref.get("sha256", ""))
        if key not in self._cache:
            self._cache[key] = self._store.read(project, ref)
        return self._cache[key]


class _RegistryReader:
    """Use the registry's verified row reader without running its constructor."""

    def __init__(self, store, project: str):
        from packages.platform import protocol_evidence

        self._protocol_evidence = protocol_evidence
        self.store = _CachedStore(store)
        self.project = project
        schema = protocol_evidence._strict_loads(
            protocol_evidence.SCHEMA_PATH.read_bytes())
        self.validator = Draft202012Validator(schema)

    def read(self, row) -> dict:
        service = self._protocol_evidence.ProtocolEvidenceService.__new__(
            self._protocol_evidence.ProtocolEvidenceService)
        service.store = self.store
        service.project = self.project
        service.validator = self.validator
        return service._read_row(row)


def _table_exists(db) -> bool:
    return db.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protocol_diagnostics'"
    ).fetchone() is not None


def _candidate_index(result: dict) -> dict[str, dict]:
    rows = result.get("protac_candidates")
    _require(isinstance(rows, list), "PROTOCOL_REASSESSMENT_CANDIDATES")
    indexed = {}
    for row in rows:
        _require(isinstance(row, dict), "PROTOCOL_REASSESSMENT_CANDIDATE")
        candidate_id = row.get("candidate_id")
        if not isinstance(candidate_id, str) or not candidate_id:
            continue
        _require(candidate_id not in indexed,
                 "PROTOCOL_REASSESSMENT_CANDIDATE_DUPLICATE")
        indexed[candidate_id] = row
    return indexed


def _validate_common(record: dict, binding: dict,
                     current_hashes: dict[str, str]) -> None:
    _require(record.get("source_binding") == binding,
             "PROTOCOL_REASSESSMENT_SOURCE_BINDING")
    _require(record.get("project_id") == binding["project"] and
             record.get("job_id") == binding["job_id"],
             "PROTOCOL_REASSESSMENT_JOB_BINDING")
    _require(record.get("scientific_approved") is False and
             record.get("formal_acceptance") is False and
             record.get("gates_affected") is False,
             "PROTOCOL_REASSESSMENT_IMMUTABLE_FLAGS")
    engine = record.get("verification_engine")
    _require(isinstance(engine, dict) and engine.get("version") == 1 and
             engine.get("source_hashes") == current_hashes,
             "PROTOCOL_REASSESSMENT_STALE_VERIFICATION_ENGINE")


def _known_evaluation(record: dict) -> dict:
    summary = record.get("summary")
    verification = record.get("measurement_verification")
    _require(isinstance(summary, dict) and isinstance(verification, dict),
             "PROTOCOL_REASSESSMENT_KNOWN_SHAPE")
    required = {
        "benchmark_scope", "ligand_scope", "raw_pass_count",
        "raw_model_denominator", "selected_pass_count",
        "selected_seed_denominator", "claimed_protocol_pass",
        "failure_count", "selected", "preserved_original_baseline",
    }
    _require(required <= set(summary), "PROTOCOL_REASSESSMENT_KNOWN_SHAPE")
    _require(summary.get("benchmark_scope") ==
             "BRD4-CRBN 6BOY related diagnostic; not the SMARCA2 novel target",
             "PROTOCOL_REASSESSMENT_KNOWN_SCOPE")
    _require(isinstance(summary.get("ligand_scope"), str) and
             bool(summary["ligand_scope"].strip()),
             "PROTOCOL_REASSESSMENT_KNOWN_SCOPE")
    _require(isinstance(summary.get("preserved_original_baseline"), str) and
             bool(summary["preserved_original_baseline"].strip()),
             "PROTOCOL_REASSESSMENT_KNOWN_SCOPE")

    raw_denominator = _strict_int(summary.get("raw_model_denominator"),
                                  "PROTOCOL_REASSESSMENT_KNOWN_RAW_COUNT", minimum=1)
    selected_denominator = _strict_int(summary.get("selected_seed_denominator"),
                                       "PROTOCOL_REASSESSMENT_KNOWN_SELECTED_COUNT", minimum=1)
    raw_pass = _strict_int(summary.get("raw_pass_count"),
                           "PROTOCOL_REASSESSMENT_KNOWN_RAW_COUNT")
    selected_pass = _strict_int(summary.get("selected_pass_count"),
                                "PROTOCOL_REASSESSMENT_KNOWN_SELECTED_COUNT")
    failures = _strict_int(summary.get("failure_count"),
                           "PROTOCOL_REASSESSMENT_KNOWN_FAILURE_COUNT")
    _require(raw_denominator == KNOWN_RAW_DENOMINATOR and
             selected_denominator == KNOWN_SELECTED_DENOMINATOR and
             raw_pass <= raw_denominator and selected_pass <= selected_denominator,
             "PROTOCOL_REASSESSMENT_KNOWN_DENOMINATOR")
    selected = summary.get("selected")
    _require(isinstance(selected, list) and len(selected) == selected_denominator,
             "PROTOCOL_REASSESSMENT_KNOWN_SELECTED_COUNT")
    seeds = []
    derived_selected_pass = 0
    for row in selected:
        _require(isinstance(row, dict), "PROTOCOL_REASSESSMENT_KNOWN_SELECTED_ROW")
        seed = _strict_int(row.get("seed"),
                           "PROTOCOL_REASSESSMENT_KNOWN_SEED")
        _strict_int(row.get("selected_model_index"),
                    "PROTOCOL_REASSESSMENT_KNOWN_MODEL_INDEX")
        _finite_number(row.get("selected_e3_CA_RMSD_after_target_alignment_A"),
                       "PROTOCOL_REASSESSMENT_KNOWN_NONFINITE")
        selected_result = _strict_bool(row.get("selected_pass"),
                                       "PROTOCOL_REASSESSMENT_KNOWN_PASS_TYPE")
        seeds.append(seed)
        derived_selected_pass += int(selected_result)
    _require(len(seeds) == len(set(seeds)),
             "PROTOCOL_REASSESSMENT_KNOWN_DUPLICATE_SEED")
    _require(derived_selected_pass == selected_pass,
             "PROTOCOL_REASSESSMENT_KNOWN_SELECTED_MISMATCH")
    claimed = _strict_bool(summary.get("claimed_protocol_pass"),
                           "PROTOCOL_REASSESSMENT_KNOWN_PASS_TYPE")
    original_protocol_claim = selected_pass >= KNOWN_REQUIRED_PASS_COUNT
    _require(claimed == original_protocol_claim,
             "PROTOCOL_REASSESSMENT_KNOWN_CLAIM_MISMATCH")

    measurement_complete = (
        verification.get("geometry_recomputed") is True and
        verification.get("stored_metrics_only") is False and
        verification.get("claimed_pass_reported_separately") is True
    )
    _require(measurement_complete,
             "PROTOCOL_REASSESSMENT_KNOWN_VERIFICATION_SHAPE")
    trusted_numerical_pass = (
        selected_pass >= KNOWN_REQUIRED_PASS_COUNT and
        failures == 0 and measurement_complete
    )

    return {
        "kind": KNOWN_KIND,
        "target_scope": "BRD4_related_6BOY_diagnostic_not_SMARCA2_generalization",
        "machine_subchecks": {
            "raw_model_accounting": {
                "status": "pass" if raw_denominator == KNOWN_RAW_DENOMINATOR else "failed",
                "observed_pass": raw_pass,
                "observed_denominator": raw_denominator,
                "required_denominator": KNOWN_RAW_DENOMINATOR,
                "note": "raw model passes are not selected seed passes",
            },
            "top_confidence_seed_selection": {
                "status": "pass" if len(seeds) == KNOWN_SELECTED_DENOMINATOR else "failed",
                "observed_unique_seeds": seeds,
                "observed_denominator": selected_denominator,
                "required_denominator": KNOWN_SELECTED_DENOMINATOR,
            },
            "measurement_completion": {
                "status": "pass" if measurement_complete else "failed",
                "geometry_recomputed": verification.get("geometry_recomputed"),
                "stored_metrics_only": verification.get("stored_metrics_only"),
                "claimed_pass_reported_separately":
                    verification.get("claimed_pass_reported_separately"),
            },
            "original_record_protocol_claim": {
                "status": "pass" if claimed else "failed",
                "claimed_protocol_pass": claimed,
                "claim_rule_required_pass": KNOWN_REQUIRED_PASS_COUNT,
                "observed_pass": selected_pass,
                "observed_denominator": selected_denominator,
                "note": "원본 protocol_pass는 실제 구현의 3/5 규칙으로 검증되며 신뢰 cutoff를 변경하지 않습니다.",
            },
            "selected_current_quality": {
                "status": "pass" if trusted_numerical_pass else "failed",
                "observed_pass": selected_pass,
                "observed_denominator": selected_denominator,
                "required_pass": KNOWN_REQUIRED_PASS_COUNT,
                "failure_count": failures,
                "trusted_expert_cutoff": "at_least_3_of_5",
            },
            "declared_method_limits": {
                "status": "pass",
                "benchmark_scope": summary["benchmark_scope"],
                "ligand_scope": summary["ligand_scope"],
                "preserved_original_baseline":
                    summary["preserved_original_baseline"],
                "smarca2_generalization_supported": False,
            },
        },
        "numerical_status": "pass" if trusted_numerical_pass else "failed",
        "scope_status": "scope_review_not_implemented",
        "effective_status": "pending" if trusted_numerical_pass else "failed",
        "reason": (
            f"CURRENT selected top-confidence result is {selected_pass}/{selected_denominator}; "
            f"신뢰된 전문가 목표는 최소 {KNOWN_REQUIRED_PASS_COUNT}/{selected_denominator}이며 "
            f"raw model pass count는 별도로 {raw_pass}/{raw_denominator}입니다. "
            "이 기준은 내부 6BOY screening에만 적용되고 SMARCA2 범위 검토는 구현되지 않았습니다."
        ),
        "scientific_authority": False,
        "human_review_synthesized": False,
    }


def _novel_evaluation(record: dict, result: dict) -> dict:
    summary = record.get("summary")
    verification = record.get("measurement_verification")
    _require(isinstance(summary, dict) and isinstance(verification, dict),
             "PROTOCOL_REASSESSMENT_NOVEL_SHAPE")
    required = {
        "candidate_count", "candidate_denominator",
        "selected_geometry_recomputed_count", "selected_geometry_denominator",
        "raw_model_count", "raw_model_denominator", "failure_count",
        "candidates", "baseline30", "core_raw40", "known_calibration_comparison",
    }
    _require(required <= set(summary), "PROTOCOL_REASSESSMENT_NOVEL_SHAPE")
    candidate_count = _strict_int(summary.get("candidate_count"),
                                  "PROTOCOL_REASSESSMENT_NOVEL_COUNT")
    candidate_denominator = _strict_int(summary.get("candidate_denominator"),
                                        "PROTOCOL_REASSESSMENT_NOVEL_COUNT", minimum=1)
    selected_count = _strict_int(summary.get("selected_geometry_recomputed_count"),
                                 "PROTOCOL_REASSESSMENT_NOVEL_SELECTED_COUNT")
    selected_denominator = _strict_int(summary.get("selected_geometry_denominator"),
                                       "PROTOCOL_REASSESSMENT_NOVEL_SELECTED_COUNT", minimum=1)
    raw_count = _strict_int(summary.get("raw_model_count"),
                            "PROTOCOL_REASSESSMENT_NOVEL_RAW_COUNT")
    raw_denominator = _strict_int(summary.get("raw_model_denominator"),
                                  "PROTOCOL_REASSESSMENT_NOVEL_RAW_COUNT", minimum=1)
    failures = _strict_int(summary.get("failure_count"),
                           "PROTOCOL_REASSESSMENT_NOVEL_FAILURE_COUNT")
    _require(candidate_count == candidate_denominator == NOVEL_CANDIDATE_DENOMINATOR and
             selected_count == selected_denominator == NOVEL_SELECTED_DENOMINATOR and
             raw_count == raw_denominator == NOVEL_RAW_DENOMINATOR,
             "PROTOCOL_REASSESSMENT_NOVEL_DENOMINATOR")
    candidates = summary.get("candidates")
    _require(isinstance(candidates, list) and len(candidates) == candidate_count,
             "PROTOCOL_REASSESSMENT_NOVEL_CANDIDATES")
    actual = _candidate_index(result)
    seen_ids = set()
    seen_e3 = set()
    evaluations = []
    all_quantitative = failures == 0
    for candidate in candidates:
        _require(isinstance(candidate, dict),
                 "PROTOCOL_REASSESSMENT_NOVEL_CANDIDATE")
        candidate_id = candidate.get("candidate_id")
        e3_type = candidate.get("e3_type")
        _require(isinstance(candidate_id, str) and candidate_id and
                 e3_type in {"CRBN", "VHL"} and candidate_id not in seen_ids and
                 e3_type not in seen_e3,
                 "PROTOCOL_REASSESSMENT_NOVEL_CANDIDATE_IDENTITY")
        _require(candidate_id in actual and actual[candidate_id].get("e3_type") == e3_type and
                 scientific_assessment._isomeric_record(actual[candidate_id]) is not None,
                 "PROTOCOL_REASSESSMENT_NOVEL_GRAPH_BINDING")
        seen_ids.add(candidate_id)
        seen_e3.add(e3_type)
        selected_models = candidate.get("selected_models")
        selected_model_denominator = _strict_int(
            candidate.get("selected_model_denominator"),
            "PROTOCOL_REASSESSMENT_NOVEL_SELECTED_COUNT", minimum=1)
        candidate_raw_denominator = _strict_int(
            candidate.get("raw_model_denominator"),
            "PROTOCOL_REASSESSMENT_NOVEL_RAW_COUNT", minimum=1)
        candidate_failures = _strict_int(candidate.get("failure_count"),
                                         "PROTOCOL_REASSESSMENT_NOVEL_FAILURE_COUNT")
        _require(isinstance(selected_models, list) and
                 len(selected_models) == selected_model_denominator == NOVEL_SELECTED_PER_CANDIDATE and
                 candidate_raw_denominator == NOVEL_RAW_PER_CANDIDATE,
                 "PROTOCOL_REASSESSMENT_NOVEL_CANDIDATE_DENOMINATOR")
        seeds = []
        for model in selected_models:
            _require(isinstance(model, dict),
                     "PROTOCOL_REASSESSMENT_NOVEL_SELECTED_MODEL")
            seeds.append(_strict_int(model.get("seed"),
                                     "PROTOCOL_REASSESSMENT_NOVEL_SEED"))
            _strict_int(model.get("selected_model_index"),
                        "PROTOCOL_REASSESSMENT_NOVEL_MODEL_INDEX")
        _require(len(seeds) == len(set(seeds)),
                 "PROTOCOL_REASSESSMENT_NOVEL_DUPLICATE_SEED")
        iptm_iqr = _finite_number(candidate.get("ipTM_IQR"),
                                  "PROTOCOL_REASSESSMENT_NOVEL_NONFINITE")
        endpoint_iqr = _finite_number(candidate.get("endpoint_IQR_A"),
                                      "PROTOCOL_REASSESSMENT_NOVEL_NONFINITE")
        contacts = _strict_bool(candidate.get("positive_fragment_contacts_all_five"),
                                "PROTOCOL_REASSESSMENT_NOVEL_CONTACT_TYPE")
        declared = _strict_bool(candidate.get("quantitative_criteria_met"),
                                "PROTOCOL_REASSESSMENT_NOVEL_PASS_TYPE")
        quantitative = (candidate_failures == 0 and
                        iptm_iqr <= NOVEL_IPTM_IQR_MAX and
                        endpoint_iqr <= NOVEL_ENDPOINT_IQR_MAX_A and contacts)
        _require(declared == quantitative,
                 "PROTOCOL_REASSESSMENT_NOVEL_CLAIM_MISMATCH")
        all_quantitative = all_quantitative and quantitative
        evaluations.append({
            "candidate_id": candidate_id,
            "e3_type": e3_type,
            "selected_unique_seeds": seeds,
            "selected_model_denominator": selected_model_denominator,
            "raw_model_denominator": candidate_raw_denominator,
            "ipTM_IQR": iptm_iqr,
            "endpoint_IQR_A": endpoint_iqr,
            "positive_fragment_contacts": {
                "observed": len(seeds) if contacts else None,
                "required": selected_model_denominator,
                "all_five": contacts,
            },
            "subchecks": {
                "ipTM_IQR": {
                    "status": "pass" if iptm_iqr <= NOVEL_IPTM_IQR_MAX else "failed",
                    "observed": iptm_iqr, "maximum": NOVEL_IPTM_IQR_MAX,
                },
                "endpoint_IQR_A": {
                    "status": "pass" if endpoint_iqr <= NOVEL_ENDPOINT_IQR_MAX_A else "failed",
                    "observed": endpoint_iqr, "maximum": NOVEL_ENDPOINT_IQR_MAX_A,
                },
                "positive_fragment_contacts": {
                    "status": "pass" if contacts else "failed",
                    "observed": len(seeds) if contacts else None,
                    "required": selected_model_denominator,
                    "all_five": contacts,
                },
            },
            "quantitative_status": "pass" if quantitative else "failed",
            "geometric_diagnostics": _copy(candidate.get("geometric_diagnostics")),
            "geometry_interpretation":
                "quantitative diagnostics only; Jaccard clustering is not human topology review",
            "residual_clashes_accepted": False,
        })
    _require(seen_e3 == {"CRBN", "VHL"},
             "PROTOCOL_REASSESSMENT_NOVEL_E3_COVERAGE")
    _require(verification.get("raw_models_geometry_recomputed") == selected_count and
             verification.get("raw_models_integrity_verified") == raw_count,
             "PROTOCOL_REASSESSMENT_NOVEL_VERIFICATION_COUNT")
    return {
        "kind": NOVEL_KIND,
        "target_scope": "SMARCA2",
        "machine_subchecks": {
            "candidate_identity": {"status": "pass", "actual_candidate_ids": sorted(seen_ids),
                                   "required_count": NOVEL_CANDIDATE_DENOMINATOR},
            "raw_integrity": {"status": "pass", "observed": raw_count,
                              "required": raw_denominator},
            "selected_geometry_completion": {"status": "pass", "observed": selected_count,
                                             "required": selected_denominator},
            "quantitative_geometry": {"status": "pass" if all_quantitative else "failed",
                                      "candidates": evaluations},
        },
        "numerical_status": "pass" if all_quantitative else "failed",
        "scope_status": "scope_review_not_implemented",
        "effective_status": "pending" if all_quantitative else "failed",
        "human_selection_required": True,
        "human_geometry_inspection_required": True,
        "reason": (
            "Machine completion is a subcheck only; exact candidate-scoped human selection "
            "and geometry inspection remain required."
        ),
        "scientific_authority": False,
        "human_review_synthesized": False,
    }


def _human_disposition(criterion: dict) -> tuple[str | None, str | None]:
    status = criterion.get("status")
    if status == "blocked":
        return "blocked", criterion.get("reason")

    def rejected(value: Any) -> bool:
        if isinstance(value, dict):
            disposition = value.get("disposition")
            if isinstance(disposition, str) and disposition in {"reject", "rejected"}:
                return True
            if value.get("action") == "reject":
                return True
            if value.get("effective_human_review_status") == "rejected":
                return True
            return any(rejected(item) for item in value.values())
        if isinstance(value, list):
            return any(rejected(item) for item in value)
        return False

    reason = criterion.get("reason")
    explicit_reason = (isinstance(reason, str) and
                       "expert" in reason.lower() and
                       ("reject" in reason.lower() or "rejected" in reason.lower()))
    if status == "failed" and (rejected(criterion.get("observed")) or
                               rejected(criterion.get("evidence")) or explicit_reason):
        return "failed", reason
    return None, None


def _category(criterion: dict) -> str:
    status = criterion.get("status")
    observed = criterion.get("observed")
    evaluation = observed.get("current_protocol") if isinstance(observed, dict) else None

    disposition, _ = _human_disposition(criterion)
    if disposition == "failed":
        return "human_review_required"

    if isinstance(evaluation, dict):
        if evaluation.get("numerical_status") == "failed":
            return "quantitative_result_failed"
        if evaluation.get("scope_status") == "scope_review_not_implemented":
            return "additional_development_needed"
        if evaluation.get("human_review_synthesized") is False:
            return "human_review_required"

    identifier = str(criterion.get("id") or "").lower().replace("-", "_")
    additional_development_ids = {
        "parent_funnel",
        "broad_families",
        "broad_chemical_families",
        "broad_chemotype_families",
        "qualified_panel",
        "qualified_panel_size",
        "panel_size",
    }
    missing_evidence_ids = {
        "core_interactions",
        "core_interaction_requirements",
        "core_interaction_preservation",
        "microstates",
        "microstate_requirements",
        "microstate_population",
    }
    if (identifier in additional_development_ids or
            ("parent" in identifier and "funnel" in identifier) or
            ("broad" in identifier and "famil" in identifier) or
            ("panel" in identifier and
             any(token in identifier for token in ("qualified", "size", "minimum")))):
        return "additional_development_needed"
    if (identifier in missing_evidence_ids or
            "microstate" in identifier or
            ("core" in identifier and "interaction" in identifier)):
        return "evidence_missing"

    if status in {"pending", "failed", "blocked"}:
        if criterion.get("needs_expert") is True:
            return "human_review_required"
        return "evidence_missing"
    return "evidence_missing"


def _recompute_summary(assessment: dict) -> None:
    criteria = assessment.get("criteria", [])
    unresolved = [item for item in criteria if item.get("status") != "pass"]
    nonformal = [item for item in criteria
                 if item.get("id") != "formal_expert_decision"]
    assessment["demo_scope"] = (
        "COMPLETE" if all(item.get("status") == "pass" for item in nonformal)
        else "DEMO/INCOMPLETE"
    )
    assessment["readiness"] = (
        "ready_for_expert_review"
        if not any(item.get("status") in {"failed", "blocked"} for item in nonformal)
        else "needs_evidence"
    )
    assessment["pending_actions"] = [{
        "criterion_id": item.get("id"),
        "category": _category(item),
        "action": item.get("reason"),
    } for item in unresolved]
    assessment["scientific_accepted"] = False


def load(db, store, project: str, job_id: str, protocol_ids: Any,
         binding: dict, result: dict, policy: dict, policy_digest: str) -> dict:
    _require(isinstance(protocol_ids, list),
             "PROTOCOL_REASSESSMENT_PROTOCOL_IDS")
    _require(1 <= len(protocol_ids) <= MAX_PROTOCOL_IDS and
             all(isinstance(value, str) and value for value in protocol_ids) and
             len(set(protocol_ids)) == len(protocol_ids),
             "PROTOCOL_REASSESSMENT_PROTOCOL_IDS")
    _require(_table_exists(db), "PROTOCOL_REASSESSMENT_REGISTRY_MISSING")
    placeholders = ",".join("?" for _ in protocol_ids)
    rows = db.execute(
        f"SELECT * FROM protocol_diagnostics WHERE project=? AND job_id=? "
        f"AND protocol_id IN ({placeholders})",
        (project, job_id, *protocol_ids),
    ).fetchall()
    _require(len(rows) == len(protocol_ids),
             "PROTOCOL_REASSESSMENT_PROTOCOL_NOT_FOUND")
    reader = _RegistryReader(store, project)
    current_hashes = current_verification_source_hashes()
    by_id = {}
    kinds = set()
    loaded = []
    for row in rows:
        record = reader.read(row)
        _validate_common(record, binding, current_hashes)
        protocol_id = record["protocol_id"]
        _require(protocol_id in protocol_ids and protocol_id not in by_id,
                 "PROTOCOL_REASSESSMENT_PROTOCOL_IDENTITY")
        kind = record.get("kind")
        _require(kind in SUPPORTED_KINDS,
                 "PROTOCOL_REASSESSMENT_UNSUPPORTED_PROTOCOL_SHAPE")
        _require(kind not in kinds,
                 "PROTOCOL_REASSESSMENT_AMBIGUOUS_PROTOCOL_KIND")
        kinds.add(kind)
        evaluation = (_known_evaluation(record) if kind == KNOWN_KIND else
                      _novel_evaluation(record, result))
        record_ref = json.loads(row["record_ref"])
        archive_ref = json.loads(row["archive_ref"])
        item = {
            "diagnostic_id": record["id"],
            "protocol_id": protocol_id,
            "kind": kind,
            "record_ref": record_ref,
            "record_sha256": record_ref["sha256"],
            "archive_ref": archive_ref,
            "archive_sha256": archive_ref["sha256"],
            "pack_manifest_sha256": record["pack_manifest_sha256"],
            "verification_engine": _copy(record["verification_engine"]),
            "evaluation": evaluation,
        }
        by_id[protocol_id] = item
        loaded.append(item)
    ordered = [by_id[value] for value in protocol_ids]
    closure = {
        "protocol_ids": protocol_ids,
        "records": ordered,
        "source_binding": _copy(binding),
        "policy_revision": policy.get("revision"),
        "policy_digest": policy_digest,
        "policy_cutoff_fingerprint": _digest({
            "numerical_criteria": policy.get("numerical_criteria"),
            "scope": policy.get("scope"),
            "decisions": policy.get("decisions"),
            "opinion": policy.get("opinion"),
            "followup": policy.get("followup"),
            "revision": policy.get("revision"),
            "active_policy_decision_ids": policy.get("active_policy_decision_ids"),
        }),
        "authorization": "scope_review_not_implemented",
        "needs_protocol_scope_review": True,
        "scientific_authority": False,
    }
    closure["fingerprint"] = _digest(closure)
    return closure


def apply(assessment: dict, closure: dict) -> None:
    selected = {row["kind"]: row for row in closure["records"]}
    criterion_ids = set()
    if KNOWN_KIND in selected:
        criterion_ids.add("known_crbn_calibration")
    if NOVEL_KIND in selected:
        criterion_ids.update({"novel_ternary_repeats", "novel_ternary_geometry"})
    for criterion in assessment.get("criteria", []):
        if criterion.get("id") not in criterion_ids:
            continue
        disposition, disposition_reason = _human_disposition(criterion)
        prior = {
            "status": criterion.get("status"),
            "observed": _copy(criterion.get("observed")),
            "reason": criterion.get("reason"),
        }
        if criterion["id"] == "known_crbn_calibration":
            evaluation = _copy(selected[KNOWN_KIND]["evaluation"])
        else:
            evaluation = _copy(selected[NOVEL_KIND]["evaluation"])
            evaluation["criterion_role"] = (
                "machine_completion_subcheck_with_human_selection_required"
                if criterion["id"] == "novel_ternary_repeats"
                else "quantitative_geometry_subchecks_not_human_topology_review"
            )
        criterion["observed"] = {
            "prior_baseline": prior,
            "current_protocol": evaluation,
        }
        if disposition is not None:
            criterion["status"] = disposition
            criterion["reason"] = disposition_reason
        else:
            criterion["status"] = evaluation["effective_status"]
            criterion["reason"] = evaluation["reason"]
        criterion["needs_expert"] = criterion["status"] != "pass"
    assessment["protocol_evaluation"] = _copy(closure)
    assessment["protocol_evaluation"]["affected_criterion_ids"] = sorted(criterion_ids)
    assessment["protocol_evaluation"]["new_scientific_authority"] = False
    _recompute_summary(assessment)


def freshness(db, service, assessment: dict, policy: dict,
              policy_digest: str) -> tuple[bool, str | None]:
    closure = assessment.get("protocol_evaluation")
    if not isinstance(closure, dict):
        return True, None
    if not _table_exists(db):
        return False, "PROTOCOL_REASSESSMENT_REGISTRY_MISSING"
    expected_policy_fingerprint = _digest({
        "numerical_criteria": policy.get("numerical_criteria"),
        "scope": policy.get("scope"),
        "decisions": policy.get("decisions"),
        "opinion": policy.get("opinion"),
        "followup": policy.get("followup"),
        "revision": policy.get("revision"),
        "active_policy_decision_ids": policy.get("active_policy_decision_ids"),
    })
    if (closure.get("policy_digest") != policy_digest or
            closure.get("policy_cutoff_fingerprint") != expected_policy_fingerprint):
        return False, "PROTOCOL_REASSESSMENT_POLICY_CHANGED"
    try:
        from packages.platform import protocol_evidence

        current_hashes = current_verification_source_hashes()
        records = closure.get("records")
        _require(isinstance(records, list) and records,
                 "PROTOCOL_REASSESSMENT_CLOSURE")
        for expected in records:
            row = db.execute(
                "SELECT * FROM protocol_diagnostics WHERE project=? AND job_id=? AND protocol_id=?",
                (service.project, assessment["job_id"], expected["protocol_id"]),
            ).fetchone()
            _require(row is not None and row["id"] == expected["diagnostic_id"],
                     "PROTOCOL_REASSESSMENT_PROTOCOL_NOT_FOUND")
            record_ref = json.loads(row["record_ref"])
            archive_ref = json.loads(row["archive_ref"])
            _require(record_ref == expected["record_ref"] and
                     archive_ref == expected["archive_ref"],
                     "PROTOCOL_REASSESSMENT_REFERENCE_CHANGED")
            record_raw = service._read_registered_ref(db, record_ref)
            archive_raw = service._read_registered_ref(db, archive_ref)
            _require(_sha(record_raw) == expected["record_sha256"] and
                     _sha(archive_raw) == expected["archive_sha256"],
                     "PROTOCOL_REASSESSMENT_ARTIFACT_CHANGED")
            record = protocol_evidence._strict_loads(record_raw)
            _require(record.get("source_binding") == assessment.get("source_binding") and
                     record.get("protocol_id") == expected["protocol_id"] and
                     record.get("kind") == expected["kind"] and
                     record.get("archive_ref") == archive_ref and
                     record.get("verification_engine", {}).get("source_hashes") == current_hashes and
                     expected.get("verification_engine", {}).get("source_hashes") == current_hashes,
                     "PROTOCOL_REASSESSMENT_STALE_INPUT")
        return True, None
    except (ContractError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        return False, str(error)


def mark_stale(value: dict, reason: str | None) -> None:
    closure = value.get("protocol_evaluation")
    affected = set(closure.get("affected_criterion_ids", [])) if isinstance(closure, dict) else set()
    for criterion in value.get("criteria", []):
        if criterion.get("id") not in affected:
            continue
        if criterion.get("status") not in {"failed", "blocked"}:
            criterion["status"] = "pending"
            criterion["needs_expert"] = True
            criterion["reason"] = "Selected protocol evidence is no longer current or intact."
        observed = criterion.get("observed")
        if isinstance(observed, dict):
            observed["protocol_evidence_current"] = False
            observed["protocol_evidence_freshness_reason"] = reason
    _recompute_summary(value)


__all__ = [
    "MAX_PROTOCOL_IDS", "KNOWN_KIND", "NOVEL_KIND", "load", "apply",
    "freshness", "mark_stale", "current_verification_source_hashes",
]

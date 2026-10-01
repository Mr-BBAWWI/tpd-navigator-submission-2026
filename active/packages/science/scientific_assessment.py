"""Deterministic scientific-readiness ledger without scientific approval authority."""
from __future__ import annotations

import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any

from rdkit import Chem

VERSION = "scientific-acceptance/20261001.3"
_ALLOWED_STATUSES = {"pass", "failed", "pending", "blocked", "not_applicable"}
_DEFAULT_NUMERICAL = {"parent_funnel_min": 5, "parent_funnel_max": 10, "modifiable_sites_min": 2, "distinct_graphs_per_site_min": 30, "broad_families_min": 6, "panel_min": 10, "panel_max": 20, "calibration_seeds_min": 3, "novel_ternary_repeats_min": 2}


def _safe(value: Any) -> Any:
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_safe(item) for item in value]
    return str(value)


def _digest(value: Any) -> str:
    encoded = json.dumps(_safe(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _criterion(identifier, title, status, observed, required, evidence, needs_expert, reason):
    if status not in _ALLOWED_STATUSES:
        raise ValueError(f"invalid criterion status: {status}")
    return {"id": identifier, "title": title, "status": status, "observed": _safe(observed), "required": _safe(required), "evidence": _safe(evidence), "needs_expert": bool(needs_expert), "reason": reason}


def _canonical_graph(smiles: Any) -> str | None:
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in mol.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _canonical_isomeric(smiles: Any) -> str | None:
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _graph_record(record: dict[str, Any]) -> str | None:
    has_mapped = record.get("mapped_smiles") is not None
    has_canonical = record.get("canonical_smiles") is not None
    mapped = _canonical_graph(record.get("mapped_smiles")) if has_mapped else None
    canonical = _canonical_graph(record.get("canonical_smiles")) if has_canonical else None
    if has_mapped and mapped is None or has_canonical and canonical is None:
        return None
    if mapped is not None and canonical is not None and mapped != canonical:
        return None
    mapped_iso = _canonical_isomeric(record.get("mapped_smiles")) if has_mapped else None
    canonical_iso = _canonical_isomeric(record.get("canonical_smiles")) if has_canonical else None
    if mapped_iso is not None and canonical_iso is not None and mapped_iso != canonical_iso:
        return None
    return mapped or canonical


def _isomeric_record(record: dict[str, Any]) -> str | None:
    if _graph_record(record) is None:
        return None
    value = record.get("mapped_smiles") if record.get("mapped_smiles") is not None else record.get("canonical_smiles")
    return _canonical_isomeric(value)


def _actual_analogs(result: dict[str, Any]) -> list[dict[str, Any]]:
    values = result.get("analogs", [])
    return [x for x in values if isinstance(x, dict)] if isinstance(values, list) else []


def _actual_candidates(result: dict[str, Any]) -> list[dict[str, Any]]:
    values = result.get("protac_candidates", result.get("candidates", []))
    return [x for x in values if isinstance(x, dict)] if isinstance(values, list) else []


def _broad_family(value):
    return "ring_modification" if value in {"ring_expansion", "ring_contraction"} else value


def _valid_analog(record: dict[str, Any]) -> bool:
    cheap = record.get("cheap_filter")
    return isinstance(cheap, dict) and cheap.get("valid") is True and _graph_record(record) is not None


def _site_graphs(analogs: list[dict[str, Any]], modifiable: set[int]) -> tuple[dict[int, set[str]], int]:
    result: dict[int, set[str]] = {}
    invalid = 0
    for analog in analogs:
        graph = _graph_record(analog)
        if graph is None:
            invalid += 1
            continue
        if not _valid_analog(analog):
            continue
        maps = analog.get("attachment_site_atom_maps", analog.get("modified_atom_maps", []))
        if not isinstance(maps, list):
            continue
        for atom_map in set(maps):
            if type(atom_map) is int and atom_map in modifiable:
                result.setdefault(atom_map, set()).add(graph)
    return result, invalid


def _families(analogs: list[dict[str, Any]]) -> set[str]:
    families = set()
    for analog in analogs:
        if not _valid_analog(analog):
            continue
        family = _broad_family(analog.get("transformation_class"))
        if isinstance(family, str) and family:
            families.add(family)
    return families


def _trusted_policy(policy: Any) -> tuple[dict[str, Any], bool]:
    if policy is None:
        return {"numerical": dict(_DEFAULT_NUMERICAL), "decisions": {}, "scope": None}, False
    if not isinstance(policy, dict):
        raise TypeError("policy must be a dictionary")
    trusted = policy.get("_validated_by_platform") is True
    numerical = policy.get("numerical_criteria", {})
    if not isinstance(numerical, dict):
        raise ValueError("policy.numerical_criteria must be a dictionary")
    unknown = set(numerical) - set(_DEFAULT_NUMERICAL)
    if unknown:
        raise ValueError(f"unknown numerical policy fields: {sorted(unknown)}")
    values = dict(_DEFAULT_NUMERICAL)
    if trusted:
        for key, value in numerical.items():
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"policy numerical value {key} must be a nonnegative integer")
            if key == "novel_ternary_repeats_min" and value < _DEFAULT_NUMERICAL[key] and policy.get("decisions", {}).get("repeat_minimum_reduction_reviewed") is not True:
                continue
            values[key] = value
    decisions = policy.get("decisions", {}) if trusted and isinstance(policy.get("decisions", {}), dict) else {}
    return {"numerical": values, "scope": policy.get("scope") if trusted else None, "decisions": decisions}, trusted


def _supplements(value: Any) -> dict[str, list[Any]]:
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise TypeError("supplements must be a dictionary")
    result = {}
    for key in ("interactions", "synthesis", "ternary", "calibration"):
        items = value.get(key, [])
        if not isinstance(items, list):
            raise ValueError(f"supplements.{key} must be a list")
        result[key] = items
    return result


def _bound(item: Any, policy_digest: str) -> bool:
    return isinstance(item, dict) and item.get("_validated_by_platform") is True and isinstance(item.get("job_id"), str) and bool(item.get("job_id")) and item.get("policy_digest") == policy_digest


def _interaction_status(items: list[Any], policy_digest: str) -> tuple[str, dict[str, Any]]:
    valid = [x for x in items if _bound(x, policy_digest)]
    unbound = len(items) - len(valid)
    failures, passed = [], 0
    for item in valid:
        failures.extend(item.get("failures", []) if isinstance(item.get("failures"), list) else [])
        reports = item.get("requirements", [])
        if isinstance(reports, list) and reports and all(isinstance(report, dict) and (not report.get("required") or report.get("status") == "preserved" and report.get("computed_pass") is True) for report in reports):
            passed += 1
    evidence = {"bound_reports": len(valid), "unbound_reports": unbound, "passing_reports": passed, "failures": failures}
    if failures:
        return "failed", evidence
    if valid and not unbound and passed == len(valid):
        return "pass", evidence
    return "pending", evidence


def _microstate_status(items: list[Any], decisions: dict[str, Any], policy_digest: str) -> tuple[str, dict[str, Any]]:
    choice = decisions.get("microstate_expert_decision")
    if not isinstance(choice, dict) or choice.get("accepted") is not True:
        return "pending", {"reason": "microstate_expert_decision_required", "pKa_prediction_performed": False}
    selected = choice.get("state_choices")
    if not isinstance(selected, dict):
        return "pending", {"reason": "explicit_state_choices_required", "pKa_prediction_performed": False}
    checked = []
    for item in items:
        if not _bound(item, policy_digest):
            return "pending", {"reason": "unbound_interaction_report", "pKa_prediction_performed": False}
        report_id = item.get("report_id") or item.get("job_id")
        index = selected.get(report_id)
        states = item.get("state_alternatives")
        if type(index) is not int or not isinstance(states, list) or index < 0 or index >= len(states):
            return "pending", {"reason": "selected_state_not_computed", "report_id": report_id, "pKa_prediction_performed": False}
        state = states[index]
        receipt = state.get("hydrogen_receipt") if isinstance(state, dict) else None
        requirements = state.get("requirements") if isinstance(state, dict) else None
        if not isinstance(state, dict) or state.get("status") != "computed" or isinstance(receipt, dict) and receipt.get("requires_review") is True:
            return "pending", {"reason": "selected_state_failed_or_minimizer_requires_review", "report_id": report_id, "pKa_prediction_performed": False}
        if not isinstance(requirements, list) or not all(isinstance(x, dict) and x.get("unambiguous") is True and x.get("pending_missing_protein_hydrogen") is not True for x in requirements):
            return "failed", {"reason": "required_interaction_not_unambiguous_in_selected_state", "report_id": report_id, "pKa_prediction_performed": False}
        checked.append({"report_id": report_id, "state_index": index, "population": "unknown"})
    if not checked:
        return "pending", {"reason": "no_bound_state_reports", "pKa_prediction_performed": False}
    return "pass", {"selected_states": checked, "population": "unknown", "population_prediction_performed": False, "pKa_prediction_performed": False}


def _synthesis_status(items: list[Any], policy_digest: str) -> tuple[str, dict[str, Any]]:
    valid = [x for x in items if _bound(x, policy_digest)]
    unbound = len(items) - len(valid)
    blocked = [x.get("candidate_id") for x in valid if x.get("overall_status") == "blocked"]
    documented = [x.get("candidate_id") for x in valid if x.get("documented_route_for_graph") is True]
    evidence = {"bound_reviewed": len(valid), "unbound": unbound, "blocked": blocked, "documented_exact": documented}
    if blocked:
        return "blocked", evidence
    if valid and not unbound and len(documented) == len(valid):
        return "pass", evidence
    return "pending", evidence


def _ternary_metrics(items: list[Any], candidates: list[dict[str, Any]], policy_digest: str) -> dict[str, Any]:
    actual_ids = {x.get("candidate_id") for x in candidates if isinstance(x.get("candidate_id"), str)}
    outcomes, groups = [], {}
    for item in items:
        if not isinstance(item, dict):
            continue
        outcome = {"candidate_id": item.get("candidate_id"), "e3_type": item.get("e3_type"), "seed": item.get("seed"), "execution_success": item.get("execution_success") is True, "failure": item.get("failure"), "validated_binding": _bound(item, policy_digest)}
        outcomes.append(outcome)
        candidate, e3, seed = item.get("candidate_id"), item.get("e3_type"), item.get("seed")
        successful = _bound(item, policy_digest) and candidate in actual_ids and e3 in {"CRBN", "VHL"} and type(seed) is int and item.get("actual_computation") is True and item.get("execution_success") is True and item.get("reference_free") is True and "invented_RMSD" not in item
        if successful:
            groups.setdefault((candidate, e3), set()).add(seed)
    return {"groups": {f"{candidate}|{e3}": len(seeds) for (candidate, e3), seeds in sorted(groups.items())}, "successful_unique_seeds": {f"{candidate}|{e3}": sorted(seeds) for (candidate, e3), seeds in sorted(groups.items())}, "all_outcomes": outcomes}


def _qualified_panel(analogs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen_ids, seen_isomers, result = set(), set(), []
    for item in analogs:
        isomer = _isomeric_record(item)
        identifier = item.get("id")
        docking = item.get("docking")
        qualifies = (
            item.get("selected") is True
            and item.get("pipeline_status") == "qualified"
            and item.get("assembly_eligible") is True
            and item.get("parent_redocking_supported") is True
            and isinstance(docking, dict)
            and docking.get("status") == "completed_with_limits"
            and docking.get("pose_preserved") is True
            and type(docking.get("passing_pose_count")) is int
            and docking.get("passing_pose_count") > 0
            and _valid_analog(item)
            and isinstance(identifier, str) and bool(identifier)
            and isomer is not None
        )
        if qualifies and identifier not in seen_ids and isomer not in seen_isomers:
            seen_ids.add(identifier)
            seen_isomers.add(isomer)
            result.append(item)
    return result


def _assembly_counts(candidates: list[dict[str, Any]]) -> dict[str, int]:
    result = {"CRBN": 0, "VHL": 0}
    seen = set()
    for item in candidates:
        graph = _graph_record(item)
        key = (item.get("candidate_id"), graph, item.get("e3_type"))
        if item.get("preview") is True or item.get("assembly_status") in {"preview", "failed"} or key in seen or not isinstance(key[0], str) or graph is None or key[2] not in result:
            continue
        seen.add(key)
        result[key[2]] += 1
    return result


def _calibration(result: dict[str, Any], decisions: dict[str, Any], minimum: int) -> tuple[str, dict[str, Any]]:
    calibration = result.get("calibration", {}).get("CRBN", {}) if isinstance(result.get("calibration"), dict) else {}
    receipts = calibration.get("seed_receipts", result.get("calibration_receipt", [])) if isinstance(calibration, dict) else []
    if not isinstance(receipts, list):
        receipts = []
    successful, failures, metrics = set(), [], {}
    for receipt in receipts:
        if not isinstance(receipt, dict):
            failures.append(receipt)
            continue
        seed = receipt.get("seed")
        if receipt.get("execution_success") is True and type(seed) is int:
            successful.add(seed)
            raw = {}
            if isinstance(receipt.get("comparison"), dict):
                raw = receipt["comparison"].get("metrics", {})
            if not isinstance(raw, dict) or not raw:
                raw = receipt.get("inspection", {}).get("comparison", {}).get("metrics", {}) if isinstance(receipt.get("inspection"), dict) else receipt.get("metrics", {})
            if isinstance(raw, dict):
                for name, value in raw.items():
                    if type(value) in (int, float) and math.isfinite(float(value)):
                        metrics.setdefault(name, []).append(float(value))
        else:
            failures.append(receipt)
    distributions = {name: {"values": values, "range": [min(values), max(values)], "median": statistics.median(values)} for name, values in metrics.items() if values}
    criterion = decisions.get("calibration_criterion")
    observed = {"successful_unique_seeds": sorted(successful), "receipt_count": len(receipts), "failures": failures, "metric_distributions": distributions, "known_seed_41_retained": 41 in successful}
    if len(successful) < minimum:
        return "pending", observed
    if not isinstance(criterion, dict) or not isinstance(criterion.get("metric"), str) or criterion.get("threshold") is None or not isinstance(criterion.get("expected_seeds"), list) or not criterion.get("scope"):
        return "pending", observed
    expected = {x for x in criterion["expected_seeds"] if type(x) is int}
    retained = all(any(isinstance(r, dict) and r.get("seed") == seed for r in receipts) for seed in expected)
    values = metrics.get(criterion["metric"], [])
    if not retained or len(expected & successful) < minimum or not values:
        return "pending", observed
    mode = criterion.get("comparison", "max_lte")
    threshold = float(criterion["threshold"])
    passes = max(values) <= threshold if mode == "max_lte" else min(values) >= threshold
    if not passes and criterion.get("expert_accept_failed_distribution") is not True:
        return "failed", observed
    return "pass", observed


def _scope_job_id(scope: Any) -> str | None:
    return scope.get("job_id") if isinstance(scope, dict) and isinstance(scope.get("job_id"), str) and scope.get("job_id") else None


def _coverage_bound(item: Any, policy_digest: str, scope: Any) -> bool:
    if not _bound(item, policy_digest):
        return False
    scoped_job = _scope_job_id(scope)
    if scoped_job is None:
        return True
    if item.get("job_id") != scoped_job:
        return False
    current_job = item.get("current_job_id")
    return current_job is None or current_job == scoped_job


def _report_digest(item: dict[str, Any]) -> str:
    value = item.get("report_digest")
    return value if isinstance(value, str) and value else _digest(item)


def _latest_reports(items: list[Any], policy_digest: str, scope: Any, identity_key: str) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    deduplicated: dict[str, tuple[int, dict[str, Any]]] = {}
    unbound = 0
    for index, item in enumerate(items):
        if not isinstance(item, dict) or not _coverage_bound(item, policy_digest, scope):
            unbound += 1
            continue
        digest = _report_digest(item)
        deduplicated.setdefault(digest, (index, item))
    grouped: dict[str, list[tuple[int, dict[str, Any]]]] = {}
    unidentified = 0
    for index, item in deduplicated.values():
        identity = item.get(identity_key)
        if not isinstance(identity, str) or not identity:
            unidentified += 1
            continue
        grouped.setdefault(identity, []).append((index, item))
    latest = {
        identity: max(
            records,
            key=lambda pair: (
                pair[1].get("report_revision") if type(pair[1].get("report_revision")) is int else -1,
                pair[1].get("completed_at") if isinstance(pair[1].get("completed_at"), str) else "",
                pair[0],
            ),
        )[1]
        for identity, records in grouped.items()
    }
    return latest, {
        "submitted_reports": len(items),
        "bound_unique_reports": len(deduplicated),
        "duplicate_digests": max(0, len(items) - unbound - len(deduplicated)),
        "unbound_reports": unbound,
        "unidentified_reports": unidentified,
    }


def _receipt_requires_review(item: dict[str, Any]) -> bool:
    def requires_review(value: Any) -> bool:
        if not isinstance(value, dict):
            return False
        if value.get("requires_review") is True:
            return True
        return any(requires_review(nested) for nested in value.values())

    def authenticated_protein_review(value: Any) -> bool:
        if not isinstance(value, dict) or value.get("requires_review") is not True:
            return False
        binding = value.get("platform_authenticated_human_review_binding")
        review = value.get("human_review")
        evidence_ref = value.get("evidence_ref")
        return bool(
            value.get("effective_human_review_status") == "accepted"
            and isinstance(binding, dict)
            and binding.get("platform_authenticated") is True
            and isinstance(binding.get("project_id"), str) and binding["project_id"]
            and isinstance(binding.get("job_id"), str) and binding["job_id"]
            and isinstance(binding.get("evidence_id"), str) and binding["evidence_id"]
            and binding.get("evidence_id") == value.get("evidence_id")
            and binding.get("evidence_ref") == evidence_ref
            and isinstance(binding.get("decision_id"), str) and binding["decision_id"]
            and isinstance(binding.get("actor_id"), str) and binding["actor_id"]
            and isinstance(review, dict)
            and review.get("accepted") is True
            and review.get("evidence_ref") == evidence_ref
            and review.get("decision_id") == binding.get("decision_id")
            and isinstance(review.get("actor"), dict)
            and review["actor"].get("id") == binding.get("actor_id")
        )

    containers = [item]
    header = item.get("header")
    if isinstance(header, dict):
        containers.append(header)
    for container in containers:
        for key in ("parent_hydrogen_receipt", "pose_hydrogen_receipt",
                    "hydrogen_receipt"):
            if requires_review(container.get(key)):
                return True
        if requires_review(container.get("hydrogen_receipts")):
            return True
        protein = container.get("protein_hydrogen_receipt")
        if (isinstance(protein, dict) and protein.get("requires_review") is True
                and not authenticated_protein_review(protein)):
            return True
    return False


def _state_pH(state: dict[str, Any], item: dict[str, Any]) -> float | None:
    conditions = state.get("pH_conditions", item.get("pH_conditions"))
    value = conditions.get("pH") if isinstance(conditions, dict) else None
    if value is None:
        value = state.get("pH_context", item.get("pH_context"))
    if type(value) not in (int, float) or not math.isfinite(float(value)):
        return None
    return float(value)


def _population_coverage(
        item: dict[str, Any], analog_id: str, state_index: int,
        state: dict[str, Any], analog_graph: str | None, rows: Any,
        scope: Any) -> tuple[dict[str, Any] | None, str]:
    if not isinstance(rows, list):
        return None, "analog_specific_pKa_or_population_evidence_required"
    matching = [row for row in rows
                if isinstance(row, dict) and row.get("analog_id") == analog_id]
    if not matching:
        return None, "analog_specific_pKa_or_population_evidence_required"
    if len(matching) != 1:
        return None, "duplicate_population_evidence"
    row = matching[0]
    report_id = item.get("report_id") or item.get("job_id")
    state_graph = _canonical_isomeric(
        state.get("mapped_smiles", state.get("canonical_smiles")))
    pH = _state_pH(state, item)
    project_id = scope.get("project_id") if isinstance(scope, dict) else None
    quantity_kind, quantity = row.get("quantity_kind"), row.get("value")
    valid = (
        row.get("_validated_by_platform") is True
        and isinstance(row.get("project_id"), str) and row["project_id"]
        and (project_id is None or row["project_id"] == project_id)
        and row.get("job_id") == item.get("job_id")
        and row.get("analog_id") == analog_id
        and row.get("report_id") == report_id
        and type(row.get("state_index")) is int
        and row["state_index"] == state_index
        and state_graph is not None
        and row.get("state_canonical_isomeric_graph") == state_graph
        and analog_graph is not None
        and row.get("analog_canonical_isomeric_graph") == analog_graph
        and pH is not None
        and type(row.get("pH")) in (int, float)
        and math.isfinite(float(row["pH"]))
        and float(row["pH"]) == pH
        and isinstance(item.get("bound_result_sha256"), str)
        and row.get("bound_result_sha256") == item.get("bound_result_sha256")
        and quantity_kind in {"pKa", "state_population"}
        and type(quantity) in (int, float)
        and math.isfinite(float(quantity))
        and row.get("source_kind") in {"measured", "computed"}
        and isinstance(row.get("source"), str) and bool(row["source"].strip())
        and all(isinstance(row.get(key), str) and bool(row[key].strip())
                for key in ("method", "protocol", "conditions",
                            "uncertainty_description", "scientist_interpretation"))
        and isinstance(row.get("evidence_ref"), dict)
        and isinstance(row.get("review_ref"), dict)
        and row["evidence_ref"] != row["review_ref"]
        and isinstance(row.get("source_locator"), str)
        and row["source_locator"].startswith("/")
        and isinstance(row.get("source_interaction_ref"), dict)
        and isinstance(row.get("source_interaction_sha256"), str)
        and bool(row["source_interaction_sha256"])
    )
    if not valid:
        return None, "population_evidence_binding_mismatch"
    if quantity_kind == "state_population" and not 0.0 <= float(quantity) <= 1.0:
        return None, "state_population_out_of_range"
    return {
        "analog_id": analog_id,
        "report_id": report_id,
        "state_index": state_index,
        "state_canonical_isomeric_graph": state_graph,
        "pH": pH,
        "quantity_kind": quantity_kind,
        "value": float(quantity),
        "source": row["source"],
        "source_kind": row["source_kind"],
        "method": row["method"],
        "uncertainty_description": row["uncertainty_description"],
        "scientist_interpretation": row["scientist_interpretation"],
        "evidence_ref": row["evidence_ref"],
        "source_locator": row["source_locator"],
        "review_ref": row["review_ref"],
        "source_interaction_ref": row["source_interaction_ref"],
        "source_interaction_sha256": row["source_interaction_sha256"],
        "bound_result_sha256": row["bound_result_sha256"],
    }, ""


def _covered_interactions(items: list[Any], expected_analogs: list[str], decisions: dict[str, Any], policy_digest: str, scope: Any, *, require_population_evidence: bool = False, analog_graphs: dict[str, str | None] | None = None) -> tuple[str, dict[str, Any], str, dict[str, Any]]:
    latest, counters = _latest_reports(items, policy_digest, scope, "analog_id")
    expected = set(expected_analogs)
    missing = sorted(expected - set(latest))
    unexpected = sorted(set(latest) - expected)
    core_failures, core_pending, passing = [], [], []
    for analog_id, item in latest.items():
        requirements = item.get("requirements")
        required = [x for x in requirements if isinstance(x, dict) and x.get("required") is True] if isinstance(requirements, list) else []
        failures = item.get("failures") if isinstance(item.get("failures"), list) else []
        receipt_review = _receipt_requires_review(item)
        if receipt_review:
            core_pending.append({"analog_id": analog_id, "reason": "parent_or_pose_hydrogen_receipt_requires_review"})
        elif failures:
            core_failures.append({"analog_id": analog_id, "failures": failures})
        if not required:
            core_pending.append({"analog_id": analog_id, "reason": "nonempty_required_requirements_required"})
            continue
        if not receipt_review:
            for requirement in required:
                if requirement.get("pending_missing_protein_hydrogen") is True or requirement.get("status") == "pending":
                    core_pending.append({"analog_id": analog_id, "requirement_id": requirement.get("id"), "reason": "protein_or_requirement_pending"})
                elif requirement.get("status") != "preserved" or requirement.get("computed_pass") is not True:
                    core_failures.append({"analog_id": analog_id, "requirement_id": requirement.get("id"), "reason": "required_interaction_failed"})
        if not failures and not receipt_review and required and all(x.get("status") == "preserved" and x.get("computed_pass") is True and x.get("pending_missing_protein_hydrogen") is not True for x in required):
            passing.append(analog_id)
    core_evidence = {**counters, "expected_analogs": sorted(expected), "covered_analogs": sorted(set(latest) & expected), "missing_analogs": missing, "unexpected_analogs": unexpected, "passing_analogs": sorted(passing), "failures": core_failures, "pending": core_pending}
    if core_failures:
        core_status = "failed"
    elif not expected or missing or core_pending:
        core_status = "pending"
    else:
        core_status = "pass"

    choice = decisions.get("microstate_expert_decision")
    selected = choice.get("state_choices") if isinstance(choice, dict) and choice.get("accepted") is True else None
    micro_failures, micro_pending, checked = [], [], []
    population_rows = choice.get("population_evidence") if isinstance(choice, dict) else None
    population_coverage, population_pending = [], []
    valid_population_ids: set[str] = set()
    analog_graphs = analog_graphs or {}
    if not isinstance(selected, dict):
        micro_pending.append({"reason": "microstate_expert_decision_required"})
    for analog_id in sorted(expected):
        item = latest.get(analog_id)
        if item is None:
            continue
        report_id = item.get("report_id") or item.get("job_id")
        index = selected.get(analog_id, selected.get(report_id)) if isinstance(selected, dict) else None
        states = item.get("state_alternatives")
        if type(index) is not int or not isinstance(states, list) or index < 0 or index >= len(states):
            micro_pending.append({"analog_id": analog_id, "reason": "selected_state_not_computed"})
            continue
        state = states[index]
        if not isinstance(state, dict) or state.get("status") != "computed":
            micro_pending.append({"analog_id": analog_id, "reason": "selected_state_not_computed"})
            continue
        if _receipt_requires_review(item) or _receipt_requires_review(state):
            micro_pending.append({"analog_id": analog_id, "reason": "selected_state_failed_or_minimizer_requires_review"})
            continue
        pH = _state_pH(state, item)
        if pH is None:
            micro_pending.append({"analog_id": analog_id, "reason": "nonempty_pH_conditions_required"})
            continue
        requirements = state.get("requirements")
        required = [x for x in requirements if isinstance(x, dict) and x.get("required", True) is True] if isinstance(requirements, list) else []
        if not required:
            micro_pending.append({"analog_id": analog_id, "reason": "nonempty_required_requirements_required"})
            continue
        for requirement in required:
            if requirement.get("pending_missing_protein_hydrogen") is True or requirement.get("status") == "pending":
                micro_pending.append({"analog_id": analog_id, "reason": "protein_hydrogen_pending", "requirement_id": requirement.get("id")})
            elif requirement.get("unambiguous") is not True or requirement.get("computed_pass", requirement.get("observed")) is not True:
                micro_failures.append({"analog_id": analog_id, "reason": "required_interaction_not_unambiguous_in_selected_state", "requirement_id": requirement.get("id")})
        state_graph = _canonical_isomeric(
            state.get("mapped_smiles", state.get("canonical_smiles")))
        checked.append({
            "analog_id": analog_id,
            "report_id": report_id,
            "state_index": index,
            "state_canonical_isomeric_graph": state_graph,
            "pH": pH,
        })
        if require_population_evidence:
            covered, reason = _population_coverage(
                item, analog_id, index, state, analog_graphs.get(analog_id),
                population_rows, scope)
            if covered is None:
                population_pending.append({
                    "analog_id": analog_id,
                    "reason": reason,
                    "pH": pH,
                    "source": None,
                    "source_kind": None,
                })
            else:
                valid_population_ids.add(analog_id)
                population_coverage.append(covered)
    if require_population_evidence:
        uncovered = sorted(expected - valid_population_ids)
        recorded = {row.get("analog_id") for row in population_pending}
        population_pending.extend({
            "analog_id": analog_id,
            "reason": "analog_specific_pKa_or_population_evidence_required",
            "pH": None,
            "source": None,
            "source_kind": None,
        } for analog_id in uncovered if analog_id not in recorded)
        micro_pending.extend(population_pending)
        if core_status == "pass" and uncovered:
            core_status = "pending"
            core_evidence["pending"] = list(core_evidence["pending"]) + [
                {"analog_id": analog_id,
                 "reason": "analog_specific_pKa_or_population_evidence_required"}
                for analog_id in uncovered
            ]
    core_evidence["population_evidence_required"] = require_population_evidence
    core_evidence["population_evidence_coverage"] = population_coverage
    core_evidence["population_evidence_missing"] = population_pending
    micro_evidence = {**counters, "expected_analogs": sorted(expected), "missing_analogs": missing, "selected_states": checked, "failures": micro_failures, "pending": micro_pending, "population_evidence_required": require_population_evidence, "population_evidence_coverage": population_coverage, "population_evidence_missing": population_pending, "population_prediction_performed": False, "pKa_prediction_performed": False}
    if micro_failures:
        micro_status = "failed"
    elif not expected or missing or micro_pending:
        micro_status = "pending"
    else:
        micro_status = "pass"
    return core_status, core_evidence, micro_status, micro_evidence


def _candidate_index(candidates: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], set[str]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for candidate in candidates:
        identifier = candidate.get("candidate_id")
        if isinstance(identifier, str) and identifier and _isomeric_record(candidate) is not None:
            grouped.setdefault(identifier, []).append(candidate)
    conflicts = {identifier for identifier, rows in grouped.items() if len({_isomeric_record(row) for row in rows}) != 1 or len({row.get("e3_type") for row in rows}) != 1}
    return {identifier: rows[0] for identifier, rows in grouped.items() if identifier not in conflicts}, conflicts


def _chosen_synthesis_ids(decisions: dict[str, Any]) -> list[str] | None:
    criterion = decisions.get("synthesis_criterion")
    values = criterion.get("chosen_candidate_ids") if isinstance(criterion, dict) else None
    if not isinstance(values, list) or not values or any(not isinstance(x, str) or not x for x in values) or len(set(values)) != len(values):
        return None
    return values


def _covered_synthesis(items: list[Any], candidates: list[dict[str, Any]], decisions: dict[str, Any], policy_digest: str, scope: Any) -> tuple[str, dict[str, Any]]:
    chosen = _chosen_synthesis_ids(decisions)
    index, conflicts = _candidate_index(candidates)
    expected = set(chosen or [])
    latest, counters = _latest_reports(items, policy_digest, scope, "candidate_id")
    missing = sorted(expected - set(latest))
    blocked, failed, documented = [], [], []
    for candidate_id in sorted(expected & set(latest)):
        report, candidate = latest[candidate_id], index.get(candidate_id)
        if candidate is None:
            failed.append({"candidate_id": candidate_id, "reason": "candidate_identity_missing_or_conflicting"})
            continue
        report_isomer = _canonical_isomeric(report.get("canonical_smiles", report.get("mapped_smiles")))
        if report_isomer != _isomeric_record(candidate):
            failed.append({"candidate_id": candidate_id, "reason": "review_graph_does_not_match_actual_candidate"})
            continue
        if report.get("overall_status") == "blocked" or report.get("blocked_new_graph") is True:
            blocked.append(candidate_id)
        elif report.get("documented_route_for_graph") is True:
            documented.append(candidate_id)
        else:
            failed.append({"candidate_id": candidate_id, "reason": "exact_graph_route_not_documented"})
    evidence = {**counters, "chosen_candidate_ids": chosen, "covered_candidate_ids": sorted(expected & set(latest)), "missing_candidate_ids": missing, "candidate_identity_conflicts": sorted(conflicts), "blocked": blocked, "failed": failed, "documented_exact": documented}
    if blocked:
        return "blocked", evidence
    if failed:
        return "failed", evidence
    if chosen is None or missing:
        return "pending", evidence
    return "pass", evidence


def _contains_forbidden_geometry_claim(value: Any) -> bool:
    if isinstance(value, dict):
        for key, item in value.items():
            normalized = str(key).lower()
            if normalized == "invented_rmsd" or normalized in {"constraint", "constraints", "reference", "reference_document", "reference_documents", "reference_structure"}:
                return True
            if _contains_forbidden_geometry_claim(item):
                return True
    elif isinstance(value, list):
        return any(_contains_forbidden_geometry_claim(item) for item in value)
    return False


def _comparison_pass(value: float, threshold: float, mode: str) -> bool:
    return value <= threshold if mode in {"lte", "max_lte"} else value >= threshold


def _covered_ternary(items: list[Any], candidates: list[dict[str, Any]], decisions: dict[str, Any], policy_digest: str, scope: Any, minimum: int) -> tuple[str, dict[str, Any], str, dict[str, Any]]:
    index, conflicts = _candidate_index(candidates)
    criterion = decisions.get("ternary_criterion")
    chosen_lists = criterion.get("chosen_candidates") if isinstance(criterion, dict) else None
    chosen = ({e3: values[0] for e3, values in chosen_lists.items()}
              if isinstance(chosen_lists, dict) and set(chosen_lists) == {"CRBN", "VHL"}
              and all(isinstance(values, list) and len(values) == 1 for values in chosen_lists.values())
              else None)
    expected_seeds_raw = criterion.get("expected_seeds") if isinstance(criterion, dict) else None
    expected_seeds = sorted(set(expected_seeds_raw)) if isinstance(expected_seeds_raw, list) and all(type(x) is int for x in expected_seeds_raw) else []
    valid_selection = isinstance(chosen, dict) and set(chosen) == {"CRBN", "VHL"} and all(isinstance(chosen[e3], str) and chosen[e3] in index and index[chosen[e3]].get("e3_type") == e3 for e3 in ("CRBN", "VHL"))
    valid_seed_policy = len(expected_seeds) >= max(2, minimum)
    outcomes, groups, summaries, failures = [], {}, {}, []
    seen_receipts = set()
    for item in items:
        if not isinstance(item, dict):
            outcomes.append({"invalid_report": _safe(item)})
            continue
        candidate_id, e3, seed = item.get("candidate_id"), item.get("e3_type"), item.get("seed")
        bound = _coverage_bound(item, policy_digest, scope)
        mapped = candidate_id in index and index[candidate_id].get("e3_type") == e3
        forbidden = _contains_forbidden_geometry_claim(item)
        summary = item.get("inspection_summary")
        documented = isinstance(summary, dict) and summary.get("compatibility") in {"compatible", "incompatible"} and isinstance(summary.get("metrics"), dict)
        success = bound and mapped and type(seed) is int and item.get("actual_computation") is True and item.get("execution_success") is True and item.get("reference_free") is True and not forbidden and documented
        receipt_key = (candidate_id, e3, seed, _report_digest(item))
        duplicate = receipt_key in seen_receipts
        seen_receipts.add(receipt_key)
        outcomes.append({"candidate_id": candidate_id, "e3_type": e3, "seed": seed, "execution_success": item.get("execution_success") is True, "failure": item.get("failure"), "validated_binding": bound, "candidate_e3_match": mapped, "forbidden_claim": forbidden, "documented_inspection_summary": documented, "duplicate": duplicate})
        if success and not duplicate:
            groups.setdefault((candidate_id, e3), set()).add(seed)
            summaries[(candidate_id, e3, seed)] = summary
        elif bound and item.get("execution_success") is not True:
            failures.append({"candidate_id": candidate_id, "e3_type": e3, "seed": seed, "failure": item.get("failure")})
    group_counts = {f"{candidate}|{e3}": len(seeds) for (candidate, e3), seeds in sorted(groups.items())}
    group_seeds = {f"{candidate}|{e3}": sorted(seeds) for (candidate, e3), seeds in sorted(groups.items())}
    missing_by_group = {}
    if valid_selection and valid_seed_policy:
        for e3 in ("CRBN", "VHL"):
            candidate_id = chosen[e3]
            missing_by_group[f"{candidate_id}|{e3}"] = sorted(set(expected_seeds) - groups.get((candidate_id, e3), set()))
    evidence = {"chosen_per_e3": chosen, "selection_valid": valid_selection, "expected_seeds": expected_seeds, "seed_policy_valid": valid_seed_policy, "groups": group_counts, "successful_unique_seeds": group_seeds, "missing_expected_seeds": missing_by_group, "candidate_identity_conflicts": sorted(conflicts), "failures": failures, "all_outcomes": outcomes}

    geometry = decisions.get("ternary_geometry_criterion")
    geometry_failures, geometry_pending, evaluated = [], [], []
    if not valid_selection or not valid_seed_policy:
        geometry_pending.append("trusted_chosen_per_e3_and_unique_expected_seeds_required")
    if not isinstance(geometry, dict) or geometry.get("disposition") not in {"accept", "reject"}:
        geometry_pending.append("expert_geometry_disposition_required")
    elif geometry.get("disposition") == "reject":
        geometry_failures.append("expert_geometry_rejected")
    else:
        metric, mode, threshold = geometry.get("metric"), geometry.get("comparison"), geometry.get("threshold")
        if not isinstance(metric, str) or mode not in {"lte", "gte", "max_lte", "min_gte"} or type(threshold) not in (int, float) or not math.isfinite(float(threshold)):
            geometry_pending.append("finite_geometry_metric_threshold_required")
        elif valid_selection and valid_seed_policy:
            for e3 in ("CRBN", "VHL"):
                candidate_id = chosen[e3]
                for seed in expected_seeds:
                    summary = summaries.get((candidate_id, e3, seed))
                    metrics = summary.get("metrics") if isinstance(summary, dict) else None
                    value = metrics.get(metric) if isinstance(metrics, dict) else None
                    if type(value) not in (int, float) or not math.isfinite(float(value)):
                        geometry_pending.append(f"missing_metric:{candidate_id}|{e3}|{seed}")
                    else:
                        passed = _comparison_pass(float(value), float(threshold), mode)
                        evaluated.append({"candidate_id": candidate_id, "e3_type": e3, "seed": seed, "metric": metric, "value": float(value), "passed": passed})
                        if not passed:
                            geometry_failures.append(f"threshold_failed:{candidate_id}|{e3}|{seed}")
    geometry_evidence = {**evidence, "geometry_evaluated": evaluated, "geometry_failures": geometry_failures, "geometry_pending": geometry_pending}
    geometry_status = "failed" if geometry_failures else ("pending" if geometry_pending else "pass")
    complete = valid_selection and valid_seed_policy and all(not missing_by_group.get(f"{chosen[e3]}|{e3}") for e3 in ("CRBN", "VHL"))
    repeat_status = "pass" if complete and geometry_status == "pass" else ("failed" if geometry_status == "failed" else "pending")
    return repeat_status, evidence, geometry_status, geometry_evidence


def _linear_quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def _trusted_followup(policy: Any, policy_trusted: bool) -> dict[str, Any] | None:
    if not policy_trusted or not isinstance(policy, dict):
        return None
    followup = policy.get("followup")
    opinion_scope = policy.get("opinion", {}).get("scope")
    if (not isinstance(followup, dict) or not isinstance(opinion_scope, dict)
            or opinion_scope.get("fx5_bound_scope") is not True):
        return None
    from packages.science.expert_followup import load_followup
    current = load_followup()
    if (followup.get("source_policy_digest") != current.get("source_policy_digest")
            or followup.get("source") != current.get("source")
            or followup.get("rules") != current.get("rules")):
        return None
    return current


def _trusted_calibration_report(items: list[Any], policy_digest: str,
                                scope: Any) -> dict[str, Any] | None:
    eligible = []
    for index, item in enumerate(items):
        if not isinstance(item, dict) or not _coverage_bound(item, policy_digest, scope):
            continue
        measurement = item.get("measurement_binding")
        distribution = item.get("distribution")
        if (item.get("id") != "CRBN_6BOY" or
                item.get("kind") != "source_independent_benchmark" or
                not isinstance(measurement, dict) or
                measurement.get("source_id") != "CRBN_6BOY" or
                not isinstance(measurement.get("module_version"), str) or
                not isinstance(measurement.get("module_fingerprint"), str) or
                not isinstance(distribution, dict) or
                not isinstance(distribution.get("seeds"), list)):
            continue
        eligible.append((index, item))
    if not eligible:
        return None
    return max(eligible, key=lambda pair: (
        pair[1].get("report_revision")
        if type(pair[1].get("report_revision")) is int else -1,
        pair[1].get("completed_at")
        if isinstance(pair[1].get("completed_at"), str) else "",
        pair[0],
    ))[1]


def _core_calibration_metric_name(name: Any) -> str | None:
    if not isinstance(name, str) or not name:
        return None
    prefix = "comparison.metrics."
    core = name[len(prefix):] if name.startswith(prefix) else name
    lowered = core.lower()
    if (lowered.startswith("confidence") or "model_confidence" in lowered or
            lowered in {"plddt", "ptm", "iptm", "ligand_iptm", "protein_iptm",
                        "complex_plddt", "complex_iplddt", "complex_pde",
                        "complex_ipde"}):
        return None
    return core


def _calibration_covered(result: dict[str, Any], decisions: dict[str, Any], minimum: int,
                         supplements: list[Any] | None = None,
                         policy_digest: str | None = None,
                         scope: Any = None) -> tuple[str, dict[str, Any]]:
    report = _trusted_calibration_report(
        supplements or [], policy_digest or "", scope)
    supplemental_source = None
    evidence_refs = None
    if report is not None:
        rows = report["distribution"]["seeds"]
        receipts = []
        for row in rows:
            if not isinstance(row, dict):
                receipts.append(row)
                continue
            technical = row.get("technical_execution_success") is True
            receipts.append({
                "seed": row.get("seed"),
                "execution_success": technical,
                "metrics": row.get("metrics") if isinstance(row.get("metrics"), dict) else {},
                "model_confidence": _safe(row.get("model_confidence")),
                "failure": None if technical else {
                    "execution_status": row.get("execution_status"),
                    "inspection_status": row.get("inspection_status"),
                    "original_status": row.get("original_status"),
                    "original_failure_reason": row.get("original_failure_reason"),
                },
                "technical_execution_success": technical,
                "inspection_status": row.get("inspection_status"),
                "original_status": row.get("original_status"),
                "original_failure_reason": row.get("original_failure_reason"),
            })
        supplemental_source = {
            "id": report.get("id"),
            "kind": report.get("kind"),
            "source": report.get("source"),
            "measurement_binding": report.get("measurement_binding"),
        }
        evidence_refs = report.get("registered_artifacts")
    else:
        calibration = result.get("calibration", {}).get("CRBN", {}) if isinstance(result.get("calibration"), dict) else {}
        receipts = calibration.get("seed_receipts", result.get("calibration_receipt", [])) if isinstance(calibration, dict) else []
        receipts = receipts if isinstance(receipts, list) else []
    criterion = decisions.get("calibration_criterion")
    expected_raw = criterion.get("expected_seeds") if isinstance(criterion, dict) else None
    expected = sorted(set(expected_raw)) if isinstance(expected_raw, list) and all(type(x) is int for x in expected_raw) else []
    metric = criterion.get("metric") if isinstance(criterion, dict) else None

    def receipt_metrics(receipt: dict[str, Any]) -> dict[str, Any]:
        raw = receipt.get("comparison", {}).get("metrics", {}) if isinstance(receipt.get("comparison"), dict) else {}
        if not raw and isinstance(receipt.get("inspection"), dict):
            comparison = receipt["inspection"].get("comparison")
            raw = comparison.get("metrics", {}) if isinstance(comparison, dict) else {}
        if not raw and isinstance(receipt.get("metrics"), dict):
            raw = receipt["metrics"]
        if not isinstance(raw, dict):
            return {}
        normalized = {}
        for name, value in raw.items():
            core = _core_calibration_metric_name(name)
            if core is None:
                continue
            if core in normalized and normalized[core] != value:
                # A prefixed and bare spelling may coexist only when identical.
                normalized[core] = None
            else:
                normalized[core] = value
        return normalized

    by_seed: dict[int, list[dict[str, Any]]] = {}
    invalid_receipts = []
    distribution_rows: dict[str, list[tuple[int, float]]] = {}
    for receipt in receipts:
        if not isinstance(receipt, dict) or type(receipt.get("seed")) is not int:
            invalid_receipts.append(_safe(receipt))
            continue
        seed = receipt["seed"]
        by_seed.setdefault(seed, []).append(receipt)
        if receipt.get("execution_success") is True:
            for name, value in receipt_metrics(receipt).items():
                if type(value) in (int, float) and math.isfinite(float(value)):
                    distribution_rows.setdefault(name, []).append((seed, float(value)))

    distributions = {}
    for name, rows in distribution_rows.items():
        values_for_metric = [value for _, value in rows]
        per_seed: dict[str, list[float]] = {}
        for seed, value in rows:
            per_seed.setdefault(str(seed), []).append(value)
        q1 = _linear_quantile(values_for_metric, .25)
        q3 = _linear_quantile(values_for_metric, .75)
        distributions[name] = {
            "values": values_for_metric,
            "seeds": [seed for seed, _ in rows],
            "by_seed": per_seed,
            "range": [min(values_for_metric), max(values_for_metric)],
            "median": statistics.median(values_for_metric),
            "q1": q1,
            "q3": q3,
            "iqr": q3 - q1,
            "quantile_method": "linear interpolation at (n-1)q",
        }

    values, failures, missing_metrics, conflicts = {}, [], [], []
    for seed in expected:
        rows = by_seed.get(seed, [])
        signatures = []
        for row in rows:
            raw = receipt_metrics(row)
            value = raw.get(metric) if isinstance(metric, str) else None
            normalized = float(value) if type(value) in (int, float) and math.isfinite(float(value)) else None
            signatures.append((row.get("execution_success") is True, normalized, _safe(row.get("failure"))))
        if len(set(json.dumps(signature, sort_keys=True) for signature in signatures)) > 1:
            conflicts.append(seed)
            continue
        if not rows:
            continue
        success, value, failure = signatures[0]
        if not success:
            failures.append({"seed": seed, "failure": failure})
        elif value is None:
            missing_metrics.append(seed)
        else:
            values[seed] = value
    observed = {
        "expected_seeds": expected,
        "receipt_seeds": sorted(by_seed),
        "metric": metric,
        "metric_by_seed": values,
        "metric_distributions": distributions,
        "missing_seeds": sorted(set(expected) - set(by_seed)),
        "missing_metric_seeds": missing_metrics,
        "failed_seeds": failures,
        "duplicate_conflict_seeds": conflicts,
        "invalid_receipts": invalid_receipts,
        "all_outcomes": _safe(receipts),
        "known_seed_41_retained": 41 in by_seed,
        "technical_success_count": sum(
            1 for rows in by_seed.values()
            if any(row.get("execution_success") is True for row in rows)),
        "scientific_model_quality_success_count": None,
        "model_confidence_by_seed": {
            str(seed): [_safe(row.get("model_confidence")) for row in rows
                        if row.get("model_confidence") is not None]
            for seed, rows in sorted(by_seed.items())
            if any(row.get("model_confidence") is not None for row in rows)
        },
        "supplemental_source": _safe(supplemental_source),
        "evidence_refs": _safe(evidence_refs),
    }
    if conflicts or failures:
        return "failed", observed
    if not isinstance(criterion, dict) or len(expected) < minimum or not criterion.get("scope") or not isinstance(metric, str):
        return "pending", observed
    threshold, mode = criterion.get("threshold"), criterion.get("comparison", "max_lte")
    if type(threshold) not in (int, float) or not math.isfinite(float(threshold)) or mode not in {"max_lte", "min_gte", "lte", "gte"}:
        return "pending", observed
    if set(values) != set(expected):
        return "pending", observed
    passes = all(_comparison_pass(values[seed], float(threshold), mode) for seed in expected)
    return ("pass" if passes else "failed"), observed


def compute_engine_fingerprint() -> str:
    base = Path(__file__).resolve().parent
    digest = hashlib.sha256()
    for name in ("expert_followup.py", "interaction_review.py", "scientific_assessment.py", "synthesis_review.py"):
        digest.update(name.encode()); digest.update(b"\0"); digest.update((base / name).read_bytes()); digest.update(b"\0")
    return digest.hexdigest()


def evaluate_design(result, *, supplements=None, policy=None):
    """Create a deterministic evidence ledger; never grant scientific acceptance."""
    if not isinstance(result, dict):
        raise TypeError("result must be a dictionary")
    supplied = _supplements(supplements)
    evaluated_policy, policy_trusted = _trusted_policy(policy)
    numerical, decisions = evaluated_policy["numerical"], evaluated_policy["decisions"]
    policy_digest = _digest(policy)
    followup_source = _trusted_followup(policy, policy_trusted)
    followup_assessment = None
    analogs, candidates = _actual_analogs(result), _actual_candidates(result)
    sites = result.get("sites", {}).get("atoms", []) if isinstance(result.get("sites"), dict) else []
    modifiable = sorted({item.get("atom_map") for item in sites if isinstance(item, dict) and item.get("state") == "MODIFIABLE" and type(item.get("atom_map")) is int and item.get("atom_map") > 0})
    graphs, invalid_graphs = _site_graphs(analogs, set(modifiable))
    graph_counts = {str(key): len(value) for key, value in sorted(graphs.items())}
    families, panel = sorted(_families(analogs)), _qualified_panel(analogs)
    parent_records = decisions.get("parent_funnel", []) if policy_trusted else []
    if not isinstance(parent_records, list):
        parent_records = []
    parent_records = [x for x in parent_records if isinstance(x, dict) and isinstance(x.get("parent_id", x.get("id")), str)]
    parent_ids = [x.get("parent_id", x.get("id")) for x in parent_records]
    criteria = []
    low, high = numerical["parent_funnel_min"], numerical["parent_funnel_max"]
    parent_pass = low <= len(parent_records) <= high
    criteria.append(_criterion("parent_funnel", "부모 화합물 funnel 범위", "pass" if parent_pass else ("pending" if not policy_trusted else "failed"), {"count": len(parent_records), "parent_ids": parent_ids}, {"min": low, "max": high}, parent_records, not parent_pass, "Only server-trusted policy.decisions.parent_funnel rows count; caller result booleans and catalogs are ignored."))
    expert_parent = parent_pass and decisions.get("parent_selection_confirmed") is True and all(x.get("expert_selected") is True for x in parent_records)
    criteria.append(_criterion("expert_parent_selection", "전문가 부모 선택", "pass" if expert_parent else "pending", {"actual_known_parent_ids": parent_ids, "confirmed": expert_parent}, True, decisions, not expert_parent, "Five to ten actual known parent IDs require a trusted expert decision; catalogs never count as expert choices."))
    site_min = numerical["modifiable_sites_min"]
    criteria.append(_criterion("modifiable_sites", "실제 MODIFIABLE 위치 수", "pass" if len(modifiable) >= site_min else "failed", {"count": len(modifiable), "atom_maps": modifiable}, {"minimum": site_min}, sites, False, "Only unique positive integer atom maps with state=MODIFIABLE count; UNKNOWN sites can never qualify."))
    per_site = numerical["distinct_graphs_per_site_min"]
    qualifying = sorted(key for key, values in graphs.items() if key in set(modifiable) and len(values) >= per_site)
    criteria.append(_criterion("distinct_constitutional_graphs", "위치별 distinct constitutional graph", "pass" if len(qualifying) >= site_min else "failed", {"counts_by_site": graph_counts, "qualifying_sites": qualifying, "invalid_graph_records": invalid_graphs}, {"minimum_per_site": per_site, "minimum_qualifying_sites": site_min}, {"actual_analog_count": len(analogs)}, False, "Graph counts intersect MODIFIABLE maps and deduplicate constitutional graphs without maps or stereochemistry."))
    criteria.append(_criterion("actual_broad_families", "실제 적용 broad family", "pass" if len(families) >= numerical["broad_families_min"] else "failed", {"count": len(families), "families": families}, {"minimum": numerical["broad_families_min"]}, {"actual_valid_analogs": sum(_valid_analog(x) for x in analogs)}, False, "Families use design-panel broad-family derivation on actual graph-valid, cheap-filter-valid records."))
    panel_count = len(panel)
    criteria.append(_criterion("qualified_panel", "전문가 검토 panel 크기", "pass" if numerical["panel_min"] <= panel_count <= numerical["panel_max"] else "failed", {"qualified_selected_unique": panel_count}, {"min": numerical["panel_min"], "max": numerical["panel_max"]}, [x.get("id") for x in panel], True, "Only selected qualified, assembly-eligible, parent-redocking-supported records with successful passing docking poses count; previews and duplicate IDs/graphs are excluded."))
    e3_counts = _assembly_counts(candidates)
    criteria.append(_criterion("both_e3_assembly", "CRBN/VHL 실제 조립", "pass" if all(e3_counts.values()) else "failed", e3_counts, {"CRBN": ">=1", "VHL": ">=1"}, {"candidate_count": len(candidates)}, False, "Distinct candidate ID plus canonical graph records count; previews and failed assemblies are excluded."))
    expected_analogs = [x["id"] for x in panel]
    interaction_status, interaction_evidence, micro_status, micro_evidence = _covered_interactions(
        supplied["interactions"], expected_analogs, decisions, policy_digest,
        evaluated_policy["scope"],
        require_population_evidence=followup_source is not None,
        analog_graphs={item["id"]: _isomeric_record(item) for item in panel},
    )
    criteria.append(_criterion("core_interaction_preservation", "필수 core interaction 보존", interaction_status, interaction_evidence, {"expected_analogs": expected_analogs, "coverage": "exact latest bound report per analog", "required_requirements": ">=1", "followup_requires_bound_pKa_or_population": followup_source is not None}, supplied["interactions"], interaction_status != "pass", "Every qualified-panel analog requires its latest policy/job-bound report. Duplicate report digests do not count twice; optional-only requirements, review-required receipts, and protein-hydrogen pending states cannot pass. Active bundled follow-up also requires complete analog/state/graph/pH/result-bound quantitative pKa or state-population evidence."))
    criteria.append(_criterion("microstates_h_direction", "microstate 및 수소 방향", micro_status, micro_evidence, {"expected_analogs": expected_analogs, "explicit_state_choice": True, "nonempty_pH_conditions": True, "followup_requires_bound_pKa_or_population": followup_source is not None}, supplied["interactions"], micro_status != "pass", "Every expected analog requires an explicitly selected computed state, nonempty pH conditions, and nonempty unambiguous required requirements. Under active bundled follow-up, pH-only records, generic citations, and reviewer statements without a separate quantitative source cannot pass."))
    calibration_status, calibration_evidence = _calibration_covered(
        result, decisions, numerical["calibration_seeds_min"],
        supplied["calibration"], policy_digest, evaluated_policy["scope"])
    followup_calibration = None
    if followup_source is not None and "calibration_criterion" not in decisions:
        report = _trusted_calibration_report(
            supplied["calibration"], policy_digest, evaluated_policy["scope"])
        if report is not None:
            from packages.science.expert_followup import evaluate_calibration
            candidate_report = evaluate_calibration(report.get("distribution"))
            if candidate_report.get("complete_expected_seed_set") is True:
                followup_calibration = candidate_report
                calibration_status = candidate_report["status"]
                calibration_evidence = {
                    **calibration_evidence,
                    "expert_followup": candidate_report,
                    "source_provenance": followup_source["source"],
                    "source_policy_digest": followup_source["source_policy_digest"],
                    "bound_report": {
                        "id": report.get("id"), "kind": report.get("kind"),
                        "job_id": report.get("job_id"),
                        "policy_digest": report.get("policy_digest"),
                        "measurement_binding": report.get("measurement_binding"),
                    },
                }
    criteria.append(_criterion("known_crbn_calibration", "known CRBN calibration seed 분포", calibration_status, calibration_evidence, {"minimum_unique_seeds": numerical["calibration_seeds_min"], "trusted_threshold_and_scope": True, "retain_all_outcomes": True, "required_source": {"id": "CRBN_6BOY", "kind": "source_independent_benchmark"}}, {"frozen_result_calibration": result.get("calibration", result.get("calibration_receipt", [])), "supplemental_calibration": supplied["calibration"]}, calibration_status != "pass", "Source-verified technical seeds establish run integrity only. The frozen three-seed result remains immutable; a bound supplemental distribution may expand the observation to five seeds. Poor seed 41 behavior and the original seed 23 reader failure remain visible, and no distribution passes without a trusted expert threshold and scope."))
    repeat_min = numerical["novel_ternary_repeats_min"]
    ternary_status, ternary_evidence, geometry_status, geometry_evidence = _covered_ternary(supplied["ternary"], candidates, decisions, policy_digest, evaluated_policy["scope"], repeat_min)
    followup_novel = None
    if followup_source is not None:
        from packages.science.expert_followup import evaluate_novel
        bound_ternary = [item for item in supplied["ternary"]
                         if _coverage_bound(item, policy_digest, evaluated_policy["scope"])]
        followup_novel = evaluate_novel(bound_ternary, {"candidates": candidates})
        ternary_evidence = {
            **ternary_evidence,
            "expert_followup": followup_novel,
            "source_provenance": followup_source["source"],
            "source_policy_digest": followup_source["source_policy_digest"],
            "coverage_binding": {
                "submitted": len(supplied["ternary"]),
                "bound": len(bound_ternary),
                "unbound": len(supplied["ternary"]) - len(bound_ternary),
            },
        }
        geometry_evidence = {**geometry_evidence, "expert_followup": followup_novel}
        if "ternary_criterion" not in decisions:
            followup_geometry = decisions.get("followup_geometry_review")
            if isinstance(followup_geometry, dict):
                from packages.science.expert_followup import evaluate_geometry_review
                geometry_review = evaluate_geometry_review(
                    followup_novel,
                    followup_geometry,
                    validated_by_platform=policy_trusted,
                )
                geometry_status = geometry_review["status"]
                geometry_evidence = {
                    **geometry_evidence,
                    "authenticated_followup_geometry_review": geometry_review,
                }
                if (geometry_status == "pass" and
                        followup_novel["quantitative_status"] == "pass"):
                    ternary_status = "pass"
                elif (geometry_status == "failed" or
                      followup_novel["quantitative_status"] == "failed"):
                    ternary_status = "failed"
                else:
                    ternary_status = "pending"
            elif followup_novel["quantitative_status"] == "failed":
                ternary_status = "failed"
    criteria.append(_criterion("novel_ternary_repeats", "신규 후보 reference-free ternary 반복", ternary_status, ternary_evidence, {"trusted_chosen_candidate_per_e3": True, "minimum_unique_expected_seeds": max(2, repeat_min), "full_expected_seed_coverage": True, "expert_geometry_evaluated": True}, supplied["ternary"], ternary_status != "pass", "Diagnostic reports for all candidates remain visible, but only trusted chosen-per-E3 candidates with exact candidate/E3 identity, complete unique seed coverage, documented inspection summaries, and evaluated expert geometry can pass. Follow-up quantitative convergence alone never resolves geometry."))
    criteria.append(_criterion("novel_ternary_geometry", "신규 ternary geometry 전문가 평가", geometry_status, geometry_evidence, {"expert_disposition": "accept", "finite_metric_threshold": True, "reference_free": True}, supplied["ternary"], geometry_status != "pass", "Execution success or repeat count is not scientific geometry approval. Client-invented RMSD, constraints, and reference documents are rejected recursively."))
    synthesis_status, synthesis_evidence = _covered_synthesis(supplied["synthesis"], candidates, decisions, policy_digest, evaluated_policy["scope"])
    criteria.append(_criterion("exact_synthesis_review", "정확한 graph 합성 검토", synthesis_status, synthesis_evidence, {"trusted_chosen_candidate_ids": True, "exact_canonical_isomeric_identity": True, "complete_coverage": True}, supplied["synthesis"], synthesis_status != "pass", "Every explicitly chosen candidate requires a latest bound exact-isomeric-graph review. Blocked new graphs remain blocked and cannot inherit routes or yields from another candidate."))
    formal_value = decisions.get("formal_expert_decision") if policy_trusted else None
    formal_status = "pass" if formal_value == "accept" else ("failed" if formal_value == "reject" else "pending")
    criteria.append(_criterion("formal_expert_decision", "공식 전문가 결정", formal_status, formal_value, "server-validated formal expert decision", decisions, formal_status != "pass", "A formal reject is FAILED. Formal accept cannot override any other failed, blocked, or pending criterion."))
    blockers = [x["id"] for x in criteria if x["status"] in {"failed", "blocked"}]
    missing = [x["id"] for x in criteria if x["status"] == "pending"]
    if followup_source is not None:
        followup_assessment = {
            "source": followup_source["source"],
            "source_policy_digest": followup_source["source_policy_digest"],
            "manifest": followup_source["manifest"],
            "quantitative_subchecks": {
                "calibration": followup_calibration,
                "novel": followup_novel,
            },
            "authority": False,
            "scientific_accepted": False,
        }
    summary = result.get("summary", {}) if isinstance(result.get("summary"), dict) else {}
    derived_valid = sum(_valid_analog(x) for x in analogs)
    reported_valid = summary.get("valid_analogs")
    consistency = {"reported_valid_analogs": reported_valid, "derived_valid_analogs": derived_valid, "differing_fields": ["valid_analogs"] if reported_valid is not None and reported_valid != derived_valid else [], "reported_values_trusted": False}
    metrics = {"actual_analog_count": len(analogs), "actual_candidate_count": len(candidates), "expected_interaction_analog_count": len(expected_analogs), "expected_interaction_analog_ids": expected_analogs, "modifiable_site_count": len(modifiable), "modifiable_site_maps": modifiable, "distinct_graph_counts_by_site": graph_counts, "actual_broad_families": families, "selected_panel_count": panel_count, "e3_candidate_counts": e3_counts, "summary_consistency": consistency}
    parent_binding = result.get("parent_scope", {}).get("actual_design_parent_id") if isinstance(result.get("parent_scope"), dict) else None
    assessment = {"format": VERSION, "revision": 0, "confirmation_status": "unconfirmed", "input_digest": _digest(result), "policy_digest": policy_digest, "engine_fingerprint": compute_engine_fingerprint(), "policy_server_validated": policy_trusted, "criteria": criteria, "metrics": metrics, "expert_followup": followup_assessment, "demo_scope": "COMPLETE" if all(item["status"] == "pass" for item in criteria if item["id"] != "formal_expert_decision") else "DEMO/INCOMPLETE", "readiness": "ready_for_expert_review" if not blockers and not missing else "needs_evidence", "scientific_accepted": False, "expert_questions": [{"id": "parent-selection", "parent_id": parent_binding, "question": f"실제 알려진 부모 ID {parent_ids} 중 어떤 5–10개가 과학적 근거로 선택됩니까?"}, {"id": "site-counts", "parent_id": parent_binding, "atom_maps": modifiable, "question": f"MODIFIABLE map {modifiable}와 constitutional graph counts {graph_counts}가 기준을 충족합니까?"}, {"id": "interaction-definition", "parent_id": parent_binding, "question": "각 후보의 필수 ligand map/protein atom interaction, baseline, microstate 및 수소 방향 근거는 무엇입니까?"}], "pending_actions": [{"criterion_id": item["id"], "action": item["reason"]} for item in criteria if item["status"] in {"pending", "failed", "blocked"}], "authority_note": "scientific_accepted is always false; only the platform may record authority after all independent criteria are resolved.", "artifact_resolution_note": "The document service must resolve compute artifacts server-side; browser-provided artifact contents are not authoritative."}
    json.dumps(assessment, ensure_ascii=False, allow_nan=False)
    return assessment


def report_markdown(assessment):
    if not isinstance(assessment, dict) or assessment.get("format") != VERSION:
        raise ValueError("assessment must be an evaluate_design result")
    labels = {"pass": "통과", "failed": "실패", "pending": "근거 대기", "blocked": "차단", "not_applicable": "해당 없음"}
    lines = ["# 과학적 평가 준비도 보고서", "", f"- 스키마: `{assessment['format']}`", f"- 준비도: **{assessment.get('readiness')}**", "- 과학적 승인: **아니오**", "", "## 기준별 결과"]
    for item in assessment.get("criteria", []):
        lines.extend((f"### {item.get('title')} — {labels.get(item.get('status'), item.get('status'))}", f"- 관찰값: `{json.dumps(item.get('observed'), ensure_ascii=False, sort_keys=True)}`", f"- 요구값: `{json.dumps(item.get('required'), ensure_ascii=False, sort_keys=True)}`", f"- 설명: {item.get('reason')}", ""))
    lines.append("## 다음 조치")
    lines.extend(f"- **{x.get('criterion_id')}**: {x.get('action')}" for x in assessment.get("pending_actions", []))
    return "\n".join(lines) + "\n"


__all__ = ["VERSION", "evaluate_design", "report_markdown", "compute_engine_fingerprint"]

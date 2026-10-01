"""Source-bound expert follow-up rules and computational-only evaluators.

The archived reply supplies narrowly scoped screening rules.  This module
verifies that source, exposes its provenance, and computes evidence summaries.
It cannot authenticate an expert, approve geometry, or grant scientific or
formal acceptance.
"""
from __future__ import annotations

import hashlib
import json
import math
import statistics
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

from packages.contracts import ContractError

ROOT = Path(__file__).resolve().parents[2]
VERSION = "source-bound-expert-followup/20261001.2"
SOURCE_RELATIVE_PATH = "cases/expert_opinions/20261001/reply.docx"
SOURCE_NAME = "TPD_Navigator_expert_followup_검토회신안_20261001.docx"
SOURCE_SHA256 = "4132ea99e82d2711ac49e3921c9d095e2cef631ef1a0487a46a7cde45380570f"
SOURCE_BYTES = 41532
EXPECTED_SEEDS = [23, 41, 61, 79, 97]


def _canonical(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _extract_hash(value: Any) -> str:
    return _sha(_canonical(value))


# Exact source extracts used by the canonical rules. They are checked against
# the DOCX after the complete-file byte/hash check.
_EXTRACTS: dict[str, Any] = {
    "/body/6": "반드시 보호할 core interaction: (1) ligand map 17 - ASN1464 OD1 방향성 H-bond를 1차 필수 interaction으로 둔다. (2) FX5 core의 소수성 포즈는 VAL1408, PHE1409, ILE1470 주변 접촉군을 보존 대상으로 둔다. (3) parent core RMSD는 analog docking에서 1.0 Å 이하를 보조 기하 기준으로 사용한다. Tyr1421-map20 H-bond와 PHE/TYR aromatic geometry는 현재 방향성/기여가 불확실하므로 필수 hard gate로 두지 않는다.",
    "/body/8": "pH 7.4 상태: 현재 charged/neutral 16-state 열거는 pKa/population 예측이 아니므로 단일 상태 선택의 근거가 부족하다. 세 analog 모두 source-like aromatic tautomer를 기준 구조로 유지하되, N19 protonation은 최소 neutral/protonated 두 상태를 비교한다. W-80f8f4a11b5d는 terminal amine map5001 protonation 조합도 추가 비교한다. analog-specific pKa 또는 합리적인 state population 근거가 확보되기 전에는 H-bond preservation을 최종 hard-pass로 승격하지 않는다.",
    "/body/13": "제안 threshold: seed-level technical reproduction을 다음과 같이 정의한다: target CA RMSD <= 1.5 Å, ligand heavy-atom RMSD after target alignment <= 3.0 Å, E3 CA RMSD after target alignment <= 10.0 Å, contact Jaccard >= 0.40. E3 RMSD 5 Å 이하는 good, 5-10 Å는 screening-acceptable, 10 Å 초과는 fail로 표시한다.",
    "/body/14": "분포 수용 규칙: 5개 unique seed 중 최소 3개가 위 seed-level 기준을 모두 만족해야 known-case calibration을 screening-use 가능으로 본다. median과 range/IQR, 성공/실패 seed 수를 함께 보고한다. seed 41은 실패 사례로 그대로 남긴다.",
    "/body/15": "주의: 이 cutoff는 6BOY CRBN calibration의 내부 screening rule이며 보편적 구조생물학 cutoff나 효능 기준으로 해석하지 않는다.",
    "/body/17": "CRBN 우선 candidate: D-99b12e64986a (warhead analog W-c2afc5e73c1a, alkyl_c6). 5-seed에서 ipTM과 linker endpoint 분산이 비교적 안정적이며, 현재 세 CRBN 후보 중 재현성 관점의 1차 검토 대상으로 선택한다.",
    "/body/18": "VHL 우선 candidate: D-b39273b7a53b (warhead analog W-80f8f4a11b5d, alkyl_c6). 5-seed confidence/ipTM이 높고 linker endpoint 변동이 작아 VHL branch의 1차 검토 대상으로 선택한다.",
    "/body/19": "convergence metric: (1) 예정된 5 seed 100% 완료, (2) ipTM IQR <= 0.20, (3) linker endpoint distance IQR <= 1.5 Å, (4) target-warhead와 E3-recruiter contact가 각 seed에서 소실되지 않는지 확인, (5) target-E3/protein-ligand severe clash가 반복적으로 증가하지 않는지 확인한다.",
    "/body/20": "수용 기준: 위 기준을 만족하고 5 seed 중 최소 4개가 동일한 qualitative ternary topology를 유지할 때 geometry review candidate로 accept한다. 단, reference-free이므로 RMSD를 새로 만들어 판정하지 않는다. CRBN/VHL raw score를 직접 비교해 E3 winner를 선정하지 않는다.",
    "/body/22": "공통 결정: 6개 모두 graph-space route proposal로 분류한다. 현재 자료에서 unsynthesizable로 판정할 근거는 없으며, exact validated route로 승인하지도 않는다.",
    "/body/26": "결정: 승인. FX5에 대해서만 high-confidence MODIFIABLE site map19 한 곳으로 site-count criterion을 충족한 것으로 처리한다.",
    "/body/27": "제한: 이 waiver는 modifiable_sites gate에만 적용한다. parent 5-10, broad family >=6, qualified panel 10-20 요구를 면제하지 않는다. UNKNOWN map8/11/12를 자동 승격하지 않는다.",
    "/body/29": "parent 추가: 허가. reference_parents에서 ready=true인 9개 parent(SMARCA2-5DKC-5BW, 5DKH-5C0, 7Z78-IF8, 8QJT-VLC, 9D12-A1A1P, 9DU0-A1BB5, 9QAC-A1I5P, 9QAD-A1I45, 9S3R-A1JLV)를 exploratory parent pool로 사용해도 된다. 이 중 최소 5개를 실제 funnel에 포함해 비교한다.",
    "/body/30": "strict 승격 조건: co-crystal identity만으로 strict parent로 보지 않는다. 각 parent에서 atom-specific SAR 또는 문헌 attachment precedent(source + locator + 구조 map)가 확보된 site만 MODIFIABLE로 승격한다.",
    "/body/34": "demo-only 예외: 내부 UI/기술 데모에서는 현재 FX5 1-parent/1-family/3-qualified panel을 표시할 수 있으나 반드시 DEMO/INCOMPLETE 라벨을 붙이고 acceptance gate는 FAILED/PENDING 상태를 그대로 유지한다. 이 데모 예외는 공식 기준을 변경하지 않는다.",
    "/body/37": "조건: parent funnel, broad-family, qualified panel, core interaction, microstate, CRBN calibration, novel ternary geometry, exact synthesis review가 모두 evidence-bound 상태로 해결된 후 최종 policy review를 수행한다. formal accept는 다른 failed/pending gate를 덮어쓰지 않는다.",
}


def _rule(identifier: str, locator: str, data: dict[str, Any]) -> dict[str, Any]:
    return {"id": identifier, "locator": locator,
            "extract_sha256": _extract_hash(_EXTRACTS[locator]), **data}


RULES = [
    _rule("core_interactions", "/body/6", {
        "required": {"ligand_map": 17, "protein_residue": "ASN1464", "protein_atom": "OD1",
                     "hydrophobic_residues": ["VAL1408", "PHE1409", "ILE1470"],
                     "parent_core_RMSD_A_max": 1.0},
        "optional": {"ligand_map": 20, "protein_residue": "TYR1421",
                     "aromatic_geometry": True}}),
    _rule("microstates", "/body/8", {
        "status": "pending", "tautomer": "source-like aromatic",
        "state_set_sizes": {"W-c2afc5e73c1a": 2, "W-4c0a639c0a41": 2,
                            "W-80f8f4a11b5d": 4},
        "pKa_or_population_invented": False}),
    _rule("calibration_seed_thresholds", "/body/13", {
        "case": "6BOY", "expected_seeds": EXPECTED_SEEDS,
        "thresholds": {"target_CA_RMSD_A": 1.5,
                       "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
                       "e3_CA_RMSD_after_target_alignment_A": 10.0,
                       "contact_jaccard": 0.40}}),
    _rule("calibration_distribution", "/body/14", {"minimum_passing_seeds": 3,
                                                     "retain_all_outcomes": True}),
    _rule("calibration_scope", "/body/15", {"scope": "internal_6BOY_screening_only",
                                             "bioacceptance": False}),
    _rule("novel_priority_crbn", "/body/17", {"e3_type": "CRBN",
                                                "candidate_id": "D-99b12e64986a",
                                                "warhead_analog_id": "W-c2afc5e73c1a"}),
    _rule("novel_priority_vhl", "/body/18", {"e3_type": "VHL",
                                              "candidate_id": "D-b39273b7a53b",
                                              "warhead_analog_id": "W-80f8f4a11b5d"}),
    _rule("novel_quantitative_convergence", "/body/19", {
        "expected_seeds": EXPECTED_SEEDS, "iptm_IQR_max": 0.20,
        "linker_endpoint_distance_A_IQR_max": 1.5,
        "contacts_each_seed": ["target.warhead>0", "e3.recruiter>0"],
        "clashes": "descriptive_pending_expert_definition"}),
    _rule("novel_geometry_review", "/body/20", {
        "same_qualitative_topology_min": 4, "geometry_status": "pending",
        "reference_free": True, "invented_RMSD_forbidden": True,
        "cross_e3_raw_score_winner_forbidden": True}),
    _rule("synthesis_proposals", "/body/22", {
        "candidate_ids": ["D-bf12802fe7d0", "D-99b12e64986a", "D-eeff6fc3f276",
                          "D-6672edc7e277", "D-b39273b7a53b", "D-9747b9733188"],
        "count": 6, "status": "route_proposal_synthesis_unverified",
        "exact_route_approved": False, "unsynthesizable": False}),
    _rule("fx5_site_waiver", "/body/26", {"parent": "FX5", "modifiable_atom_map": 19}),
    _rule("fx5_waiver_limits", "/body/27", {"other_numerical_waiver": False,
                                             "unknown_maps": [8, 11, 12]}),
    _rule("exploratory_parents", "/body/29", {
        "allowed": ["SMARCA2-5DKC-5BW", "SMARCA2-5DKH-5C0", "SMARCA2-7Z78-IF8",
                    "SMARCA2-8QJT-VLC", "SMARCA2-9D12-A1A1P", "SMARCA2-9DU0-A1BB5",
                    "SMARCA2-9QAC-A1I5P", "SMARCA2-9QAD-A1I45",
                    "SMARCA2-9S3R-A1JLV"], "minimum_compared": 5}),
    _rule("strict_parent_sites", "/body/30", {"source_locator_atom_mapping_required": True}),
    _rule("demo_only", "/body/34", {"formal_numerical_waiver": False,
                                      "acceptance_gates_unchanged": True}),
    _rule("formal_review", "/body/37", {"status": "pending", "override_failed_pending": False}),
]

MANIFEST = {
    "format": VERSION,
    "manifest_type": "source_bound_relayed_expert_followup",
    "review_date": "2026-10-01",
    "origin": "user-relayed-document",
    "source_name": SOURCE_NAME,
    "authority": False,
    "scientific_accepted": False,
}
SOURCE_POLICY_DIGEST = _sha(_canonical({"manifest": MANIFEST, "rules": RULES,
                                        "source_sha256": SOURCE_SHA256,
                                        "source_bytes": SOURCE_BYTES}))


def _source_file(root: Path) -> Path:
    root = Path(root)
    relative = Path(SOURCE_RELATIVE_PATH)
    if relative.is_absolute() or ".." in relative.parts:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_PATH")
    try:
        strict_root = root.resolve(strict=True)
        candidate = root / relative
        current = candidate
        while current != root and current != current.parent:
            if current.exists() and current.is_symlink():
                raise ContractError("EXPERT_FOLLOWUP_SOURCE_FILE")
            current = current.parent
        path = candidate.resolve(strict=True)
    except ContractError:
        raise
    except OSError as error:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_FILE") from error
    if not path.is_file() or path.is_symlink() or strict_root not in path.parents:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_FILE")
    return path


def _extract_docx(path: Path, locator: str) -> Any:
    parts = locator.strip("/").split("/")
    if len(parts) not in {2, 4} or parts[0] != "body" or not parts[1].isdigit():
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
    body_index = int(parts[1])
    if body_index < 0:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
    try:
        with zipfile.ZipFile(path) as archive:
            xml = archive.read("word/document.xml")
        document = ElementTree.fromstring(xml)
        ns = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
        body = document.find(f"{ns}body")
        blocks = list(body) if body is not None else []
        if body_index >= len(blocks):
            raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
        block = blocks[body_index]
        if len(parts) == 2:
            value = "".join(node.text or "" for node in block.iter(f"{ns}t")).strip()
        else:
            if parts[2] != "row" or not parts[3].isdigit():
                raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
            row_index = int(parts[3])
            rows = block.findall(f"{ns}tr")
            if row_index < 0 or row_index >= len(rows):
                raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
            value = ["".join(node.text or "" for node in cell.iter(f"{ns}t")).strip()
                     for cell in rows[row_index].findall(f"{ns}tc")]
    except ContractError:
        raise
    except (OSError, KeyError, zipfile.BadZipFile, ElementTree.ParseError) as error:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR") from error
    if value == "" or value == []:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_LOCATOR")
    return value


def load_followup(root: Path = ROOT) -> dict[str, Any]:
    """Verify and return the canonical source-bound follow-up policy."""
    path = _source_file(Path(root))
    raw = path.read_bytes()
    if len(raw) != SOURCE_BYTES or _sha(raw) != SOURCE_SHA256:
        raise ContractError("EXPERT_FOLLOWUP_SOURCE_HASH")
    verified = {}
    for locator, expected in _EXTRACTS.items():
        observed = _extract_docx(path, locator)
        if observed != expected or _extract_hash(observed) != _extract_hash(expected):
            raise ContractError("EXPERT_FOLLOWUP_SOURCE_EXTRACT")
        verified[locator] = _extract_hash(observed)
    return {"manifest": MANIFEST, "rules": RULES,
            "source": {"path": SOURCE_RELATIVE_PATH, "name": SOURCE_NAME,
                       "sha256": SOURCE_SHA256, "bytes": SOURCE_BYTES,
                       "verified_extracts": verified},
            "authority": False, "source_policy_digest": SOURCE_POLICY_DIGEST,
            "scientific_accepted": False}


def _linear(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def _summary(values: list[float]) -> dict[str, Any]:
    q1, q3 = _linear(values, .25), _linear(values, .75)
    return {"all_finite_values": values, "count": len(values),
            "median": statistics.median(values) if values else None,
            "range": [min(values), max(values)] if values else None,
            "q1": q1, "q3": q3, "iqr": q3 - q1 if values else None,
            "quantile_method": "linear interpolation at (n-1)q"}


def _safe_raw(value: Any) -> Any:
    """Retain invalid raw values without emitting non-standard JSON numbers."""
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        return {"nonfinite_numeric": repr(value)}
    if isinstance(value, dict):
        return {str(key): _safe_raw(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_raw(item) for item in value]
    return {"python_type": type(value).__name__, "representation": repr(value)}


_CALIBRATION_LIMITS = {
    "target_CA_RMSD_A": ("lte", 1.5),
    "ligand_heavy_atom_RMSD_after_target_alignment_A": ("lte", 3.0),
    "e3_CA_RMSD_after_target_alignment_A": ("lte", 10.0),
    "contact_jaccard": ("gte", 0.40),
}


def _normalized_metrics(raw: Any) -> tuple[dict[str, Any], list[str]]:
    if not isinstance(raw, dict):
        return {}, ["metrics_not_object"]
    normalized: dict[str, Any] = {}
    conflicts = []
    prefix = "comparison.metrics."
    for name, value in raw.items():
        core = name[len(prefix):] if isinstance(name, str) and name.startswith(prefix) else name
        if core in normalized and normalized[core] != value:
            conflicts.append(str(core))
        else:
            normalized[core] = value
    return normalized, sorted(set(conflicts))


def evaluate_calibration(distribution: Any) -> dict[str, Any]:
    """Evaluate only the internal 6BOY five-seed screening rule."""
    rows = distribution.get("seeds") if isinstance(distribution, dict) else None
    rows = rows if isinstance(rows, list) else []
    by_seed: dict[int, list[dict[str, Any]]] = {}
    invalid_rows, unexpected = [], []
    finite_values = {name: [] for name in _CALIBRATION_LIMITS}
    evaluated_rows = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or type(row.get("seed")) is not int:
            invalid_rows.append({"index": index, "row": row})
            continue
        seed = row["seed"]
        by_seed.setdefault(seed, []).append(row)
        if seed not in EXPECTED_SEEDS:
            unexpected.append(seed)
        metrics, conflicts = _normalized_metrics(row.get("metrics"))
        metric_results, seed_pass = {}, row.get("technical_execution_success") is True
        for name, (mode, threshold) in _CALIBRATION_LIMITS.items():
            value = metrics.get(name)
            valid = type(value) in (int, float) and math.isfinite(float(value))
            if valid:
                numeric = float(value)
                finite_values[name].append(numeric)
                domain = 0 <= numeric <= 1 if name == "contact_jaccard" else numeric >= 0
                passed = domain and (numeric <= threshold if mode == "lte" else numeric >= threshold)
            else:
                numeric, domain, passed = None, False, False
            metric_results[name] = {"value": numeric, "finite_numeric_not_bool": valid,
                                    "domain_valid": domain, "comparison": mode,
                                    "threshold": threshold, "passed": passed}
            seed_pass = seed_pass and passed
        seed_pass = seed_pass and not conflicts
        e3 = metric_results["e3_CA_RMSD_after_target_alignment_A"]["value"]
        e3_class = ("invalid" if e3 is None or e3 < 0 else "good" if e3 <= 5
                    else "screening_acceptable" if e3 <= 10 else "fail")
        evaluated_rows.append({"index": index, "seed": seed,
                               "technical_execution_success": row.get("technical_execution_success") is True,
                               "original_status": row.get("original_status"),
                               "original_failure_reason": row.get("original_failure_reason"),
                               "metric_conflicts": conflicts, "metrics": metric_results,
                               "e3_geometry_class": e3_class, "passed_all_metrics": seed_pass})
    duplicates = sorted(seed for seed, values in by_seed.items() if len(values) > 1)
    missing = sorted(set(EXPECTED_SEEDS) - set(by_seed))
    seed_verdicts = []
    passing = 0
    for seed in EXPECTED_SEEDS:
        matches = [row for row in evaluated_rows if row["seed"] == seed]
        passed = len(matches) == 1 and matches[0]["passed_all_metrics"]
        passing += int(passed)
        seed_verdicts.append({"seed": seed, "receipt_count": len(matches),
                              "passed": passed, "evaluations": matches})
    complete = not missing
    fail_closed = bool(duplicates or unexpected or invalid_rows)
    if missing:
        status = "pending"
    elif passing < 3 or fail_closed:
        status = "failed"
    else:
        status = "pass"
    failure_count = sum(1 for verdict in seed_verdicts
                        if verdict["receipt_count"] > 0 and not verdict["passed"])
    return {"format": VERSION, "evaluation": "internal_6BOY_screening",
            "status": status, "expected_seeds": EXPECTED_SEEDS,
            "complete_expected_seed_set": complete, "passing_seed_count": passing,
            "failure_count": failure_count, "missing_count": len(missing),
            "required_passing_seed_count": 3, "seed_verdicts": _safe_raw(seed_verdicts),
            "duplicate_seeds": duplicates, "missing_seeds": missing,
            "unexpected_seeds": sorted(set(unexpected)),
            "invalid_rows": _safe_raw(invalid_rows),
            "statistics": {name: _summary(values) for name, values in finite_values.items()},
            "all_outcomes": _safe_raw(rows), "bioacceptance": False, "authority": False,
            "scientific_accepted": False,
            "scope_note": "Only the internal 6BOY CRBN screening calibration rule; no averaging hides seed failures."}


_PRIORITIES = {
    "CRBN": ("D-99b12e64986a", "W-c2afc5e73c1a"),
    "VHL": ("D-b39273b7a53b", "W-80f8f4a11b5d"),
}


def _candidate_rows(candidates: Any) -> list[dict[str, Any]] | None:
    if candidates is None:
        return None
    values = candidates.get("candidates") if isinstance(candidates, dict) else candidates
    if not isinstance(values, list):
        return []
    return [value for value in values if isinstance(value, dict)]


def evaluate_novel(receipts: Any, candidates: Any = None) -> dict[str, Any]:
    """Evaluate priority-candidate convergence without approving geometry."""
    raw = receipts.get("all_seed_receipts") if isinstance(receipts, dict) else receipts
    raw = raw if isinstance(raw, list) else []
    candidate_rows = _candidate_rows(candidates)
    identity = {}
    identity_failures = []
    if candidate_rows is not None:
        for e3, (candidate_id, warhead_id) in _PRIORITIES.items():
            matches = [row for row in candidate_rows if row.get("candidate_id") == candidate_id]
            valid = len(matches) == 1 and matches[0].get("e3_type") == e3 and \
                matches[0].get("warhead_analog_id") == warhead_id
            identity[e3] = {"candidate_id": candidate_id, "warhead_analog_id": warhead_id,
                            "match_count": len(matches), "valid": valid}
            if not valid:
                identity_failures.append(e3)
    groups: dict[str, list[dict[str, Any]]] = {e3: [] for e3 in _PRIORITIES}
    nonpriority = []
    for index, receipt in enumerate(raw):
        if not isinstance(receipt, dict):
            nonpriority.append({"index": index, "invalid": receipt})
            continue
        e3, candidate_id = receipt.get("e3_type"), receipt.get("candidate_id")
        if e3 in _PRIORITIES and candidate_id == _PRIORITIES[e3][0]:
            groups[e3].append({"index": index, "receipt": receipt})
        else:
            nonpriority.append({"index": index, "candidate_id": candidate_id, "e3_type": e3})
    group_reports = {}
    overall_quantitative = "pass"
    for e3, entries in groups.items():
        by_seed: dict[int, list[dict[str, Any]]] = {}
        evaluations, iptms, endpoints, invalid_priority_rows = [], [], [], []
        for entry in entries:
            receipt = entry["receipt"]
            seed = receipt.get("seed")
            if type(seed) is int:
                by_seed.setdefault(seed, []).append(receipt)
            inspection = receipt.get("inspection")
            confidence = inspection.get("model_confidence") if isinstance(inspection, dict) else None
            descriptive = inspection.get("descriptive_metrics") if isinstance(inspection, dict) else None
            contacts = descriptive.get("contacts_by_protein_role_and_ligand_group") if isinstance(descriptive, dict) else None
            iptm = confidence.get("iptm") if isinstance(confidence, dict) else None
            endpoint = descriptive.get("linker_endpoint_distance_A") if isinstance(descriptive, dict) else None
            target_contact = contacts.get("target", {}).get("warhead") if isinstance(contacts, dict) and isinstance(contacts.get("target"), dict) else None
            e3_contact = contacts.get("e3", {}).get("recruiter") if isinstance(contacts, dict) and isinstance(contacts.get("e3"), dict) else None
            valid_iptm = (type(iptm) in (int, float) and math.isfinite(float(iptm))
                          and 0 <= float(iptm) <= 1)
            valid_endpoint = (type(endpoint) in (int, float)
                              and math.isfinite(float(endpoint)) and endpoint >= 0)
            valid_contacts = (type(target_contact) is int and target_contact > 0
                              and type(e3_contact) is int and e3_contact > 0)
            contradictory_failure = (receipt.get("status") in {"failed", "failure"}
                                     or receipt.get("failure") not in (None, False, ""))
            completed = (receipt.get("execution_success") is True and
                         receipt.get("actual_computation") is True and
                         receipt.get("reference_free") is True and
                         not contradictory_failure)
            valid_seed = type(seed) is int and seed in EXPECTED_SEEDS
            valid = completed and valid_seed and valid_iptm and valid_endpoint and valid_contacts
            if not valid_seed:
                invalid_priority_rows.append({"index": entry["index"], "seed": seed})
            if valid_iptm:
                iptms.append(float(iptm))
            if valid_endpoint:
                endpoints.append(float(endpoint))
            evaluations.append({"seed": seed, "completed": completed,
                                "execution_success": receipt.get("execution_success") is True,
                                "contradictory_failure_status": contradictory_failure,
                                "iptm": iptm, "iptm_domain_valid": valid_iptm,
                                "linker_endpoint_distance_A": endpoint,
                                "target_warhead_contacts": target_contact,
                                "e3_recruiter_contacts": e3_contact,
                                "contacts_present": valid_contacts,
                                "protein_ligand_clashes": descriptive.get("protein_ligand_clashes") if isinstance(descriptive, dict) else None,
                                "target_e3_heavy_atom_clashes": descriptive.get("target_e3_heavy_atom_clashes") if isinstance(descriptive, dict) else None,
                                "valid_quantitative_receipt": valid,
                                "qualitative_topology_annotation": receipt.get("qualitative_topology_annotation")})
        duplicates = sorted(seed for seed, values in by_seed.items() if len(values) > 1)
        missing = sorted(set(EXPECTED_SEEDS) - set(by_seed))
        unexpected = sorted(seed for seed in by_seed if seed not in EXPECTED_SEEDS)
        all_valid = all(len(by_seed.get(seed, [])) == 1 and
                        next(item for item in evaluations if item["seed"] == seed)["valid_quantitative_receipt"]
                        for seed in EXPECTED_SEEDS if seed in by_seed)
        iptm_summary, endpoint_summary = _summary(iptms), _summary(endpoints)
        thresholds_pass = (iptm_summary["count"] == 5 and iptm_summary["iqr"] <= .20 and
                           endpoint_summary["count"] == 5 and endpoint_summary["iqr"] <= 1.5)
        if missing:
            status = "pending"
        elif duplicates or unexpected or invalid_priority_rows or not all_valid or not thresholds_pass:
            status = "failed"
        else:
            status = "pass"
        if status == "failed":
            overall_quantitative = "failed"
        elif status == "pending" and overall_quantitative == "pass":
            overall_quantitative = "pending"
        group_reports[e3] = {"candidate_id": _PRIORITIES[e3][0],
                             "warhead_analog_id": _PRIORITIES[e3][1], "status": status,
                             "missing_seeds": missing, "duplicate_seeds": duplicates,
                             "unexpected_seeds": unexpected,
                             "invalid_priority_rows": _safe_raw(invalid_priority_rows),
                             "receipts": _safe_raw(evaluations),
                             "iptm": iptm_summary, "linker_endpoint_distance_A": endpoint_summary,
                             "thresholds": {"iptm_IQR_max": .20,
                                            "linker_endpoint_distance_A_IQR_max": 1.5}}
    if identity_failures:
        overall_quantitative = "failed"
    return {"format": VERSION, "quantitative_status": overall_quantitative,
            "priority_groups": group_reports, "candidate_identity": identity,
            "candidate_identity_failures": identity_failures,
            "geometry_review": {"status": "pending", "scientific_accepted": False,
                                "required_same_qualitative_topology_seeds": 4,
                                "topology_requires_platform_validated_expert_annotations": True,
                                "severe_clash_recurrence_requires_expert_definition": True,
                                "quantitative_pass_is_not_geometry_approval": True},
            "all_raw_receipts": _safe_raw(raw),
            "nonpriority_receipt_index": _safe_raw(nonpriority),
            "priority_only_used_for_evaluation": True,
            "cross_e3_raw_score_winner_selected": False,
            "authority": False, "scientific_accepted": False}


def _geometry_quantitative_status(quantitative: Any) -> str:
    if not isinstance(quantitative, dict) or quantitative.get("format") != VERSION:
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_RECORD")
    groups = quantitative.get("priority_groups")
    identity = quantitative.get("candidate_identity")
    if (not isinstance(groups, dict) or set(groups) != set(_PRIORITIES) or
            not isinstance(identity, dict) or set(identity) != set(_PRIORITIES)):
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_GROUPS")

    statuses = []
    for e3, (candidate_id, warhead_id) in _PRIORITIES.items():
        group = groups[e3]
        identity_row = identity[e3]
        if (not isinstance(group, dict) or
                group.get("candidate_id") != candidate_id or
                group.get("warhead_analog_id") != warhead_id or
                not isinstance(identity_row, dict) or
                identity_row.get("candidate_id") != candidate_id or
                identity_row.get("warhead_analog_id") != warhead_id or
                identity_row.get("match_count") != 1 or
                identity_row.get("valid") is not True):
            raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_IDENTITY")
        status = group.get("status")
        if status not in {"pass", "failed", "pending"}:
            raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_STATUS")
        statuses.append(status)
        if status != "pass":
            continue

        receipts = group.get("receipts")
        seeds = [row.get("seed") for row in receipts
                 if isinstance(row, dict)] if isinstance(receipts, list) else []
        iptm = group.get("iptm")
        endpoint = group.get("linker_endpoint_distance_A")
        if (not isinstance(receipts, list) or len(receipts) != 5 or
                len(seeds) != 5 or len(set(seeds)) != 5 or
                set(seeds) != set(EXPECTED_SEEDS) or
                group.get("missing_seeds") != [] or
                group.get("duplicate_seeds") != [] or
                group.get("unexpected_seeds") != [] or
                group.get("invalid_priority_rows") != [] or
                not all(row.get("valid_quantitative_receipt") is True
                        for row in receipts if isinstance(row, dict)) or
                not isinstance(iptm, dict) or iptm.get("count") != 5 or
                type(iptm.get("iqr")) not in (int, float) or
                not math.isfinite(float(iptm["iqr"])) or iptm["iqr"] > .20 or
                not isinstance(endpoint, dict) or endpoint.get("count") != 5 or
                type(endpoint.get("iqr")) not in (int, float) or
                not math.isfinite(float(endpoint["iqr"])) or
                endpoint["iqr"] > 1.5):
            raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_COVERAGE")
    derived = ("failed" if "failed" in statuses else
               "pending" if "pending" in statuses else "pass")
    if quantitative.get("quantitative_status") != derived:
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_QUANTITATIVE_STATUS")
    return derived


def evaluate_geometry_review(quantitative: Any, trusted_decision: Any,
                             validated_by_platform: bool = False) -> dict[str, Any]:
    """Evaluate an authenticated review as an internal geometry candidate only."""
    if validated_by_platform is not True:
        raise ContractError(
            "EXPERT_FOLLOWUP_GEOMETRY_PLATFORM_VALIDATION_REQUIRED")
    quantitative_status = _geometry_quantitative_status(quantitative)
    if (not isinstance(trusted_decision, dict) or
            trusted_decision.get("_validated_by_platform") is not True):
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_TRUSTED_DECISION")
    scope = trusted_decision.get("scope")
    source_binding = trusted_decision.get("source_binding")
    source_ref = trusted_decision.get("source_ref")
    actor = trusted_decision.get("actor")
    if (not isinstance(scope, dict) or
            not isinstance(scope.get("project_id"), str) or not scope["project_id"] or
            not isinstance(scope.get("job_id"), str) or not scope["job_id"] or
            not isinstance(source_binding, dict) or
            source_binding.get("project_id") != scope["project_id"] or
            source_binding.get("job_id") != scope["job_id"] or
            source_binding.get("followup_source_policy_digest") !=
                SOURCE_POLICY_DIGEST or
            source_binding.get("source_ref") != source_ref or
            not isinstance(source_binding.get("source_sha256"), str) or
            len(source_binding["source_sha256"]) != 64 or
            not isinstance(source_binding.get("result_sha256"), str) or
            len(source_binding["result_sha256"]) != 64 or
            not isinstance(trusted_decision.get("decision_id"), str) or
            not trusted_decision["decision_id"] or
            not isinstance(actor, dict) or
            not isinstance(actor.get("id"), str) or not actor["id"]):
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_SCOPE_BINDING")

    expected = {e3: candidate for e3, (candidate, _) in _PRIORITIES.items()}
    branches = trusted_decision.get("branches")
    if not isinstance(branches, list) or len(branches) != 2:
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_BRANCHES")
    e3_values = [branch.get("e3_type") for branch in branches
                 if isinstance(branch, dict)]
    candidate_values = [branch.get("candidate_id") for branch in branches
                        if isinstance(branch, dict)]
    if (len(e3_values) != 2 or len(set(e3_values)) != 2 or
            set(e3_values) != set(expected) or
            len(candidate_values) != 2 or len(set(candidate_values)) != 2 or
            set(candidate_values) != set(expected.values())):
        raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_BRANCH_COVERAGE")

    branch_reports = {}
    failures = []
    for branch in branches:
        e3_type = branch.get("e3_type")
        candidate_id = branch.get("candidate_id")
        reviews = branch.get("seed_reviews")
        if (candidate_id != expected.get(e3_type) or
                not isinstance(branch.get("clash_definition"), str) or
                not branch["clash_definition"].strip() or
                not isinstance(branch.get("rationale"), str) or
                not branch["rationale"].strip() or
                not isinstance(reviews, list) or len(reviews) != 5):
            raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_BRANCH")
        seeds = [review.get("seed") for review in reviews
                 if isinstance(review, dict)]
        if (len(seeds) != 5 or len(set(seeds)) != 5 or
                set(seeds) != set(EXPECTED_SEEDS)):
            raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_SEEDS")

        labels: dict[str, int] = {}
        reviewed = []
        for review in reviews:
            label = review.get("topology_label")
            acceptable = review.get("severe_clash_acceptable")
            binding = review.get("evidence_binding")
            evidence_ref = review.get("evidence_ref")
            raw_hashes = binding.get("raw_artifact_sha256") \
                if isinstance(binding, dict) else None
            if (not isinstance(label, str) or not label.strip() or
                    type(acceptable) is not bool or
                    not isinstance(review.get("rationale"), str) or
                    not review["rationale"].strip() or
                    not isinstance(binding, dict) or
                    binding.get("evidence_ref") != evidence_ref or
                    binding.get("project_id") != scope["project_id"] or
                    binding.get("job_id") != scope["job_id"] or
                    binding.get("result_sha256") !=
                        source_binding["result_sha256"] or
                    binding.get("candidate_id") != candidate_id or
                    binding.get("e3_type") != e3_type or
                    binding.get("seed") != review["seed"] or
                    binding.get("actual_computation") is not True or
                    binding.get("execution_success") is not True or
                    binding.get("reference_free") is not True or
                    not isinstance(binding.get("candidate_graph_sha256"), str) or
                    len(binding["candidate_graph_sha256"]) != 64 or
                    not isinstance(binding.get("evidence_sha256"), str) or
                    len(binding["evidence_sha256"]) != 64 or
                    not isinstance(binding.get("artifact_sha256"), str) or
                    len(binding["artifact_sha256"]) != 64 or
                    not isinstance(binding.get("plan_digest"), str) or
                    not binding["plan_digest"] or
                    not isinstance(raw_hashes, dict) or not raw_hashes or
                    not all(isinstance(value, str) and len(value) == 64
                            for value in raw_hashes.values())):
                raise ContractError("EXPERT_FOLLOWUP_GEOMETRY_SEED_REVIEW")
            normalized = label.strip()
            labels[normalized] = labels.get(normalized, 0) + 1
            if acceptable is False:
                failures.append(
                    f"severe_clash_rejected:{candidate_id}|{e3_type}|{review['seed']}")
            reviewed.append({
                "seed": review["seed"],
                "topology_label": normalized,
                "severe_clash_acceptable": acceptable,
                "evidence_binding": _safe_raw(binding),
            })
        maximum_same = max(labels.values()) if labels else 0
        if maximum_same < 4:
            failures.append(
                f"insufficient_same_topology:{candidate_id}|{e3_type}")
        branch_reports[e3_type] = {
            "candidate_id": candidate_id,
            "topology_label_counts": labels,
            "maximum_same_topology_seed_count": maximum_same,
            "required_same_topology_seed_count": 4,
            "all_severe_clashes_acceptable": all(
                review["severe_clash_acceptable"] is True
                for review in reviews),
            "clash_definition": branch["clash_definition"],
            "reviews": reviewed,
        }

    pending = []
    if quantitative_status == "failed":
        failures.insert(0, "quantitative_convergence_failed")
    elif quantitative_status == "pending":
        pending.append("quantitative_convergence_pending")
    status = "failed" if failures else "pending" if pending else "pass"
    return {
        "format": VERSION,
        "evaluation": "authenticated_followup_geometry_review",
        "status": status,
        "disposition": (
            "internal_geometry_review_candidate" if status == "pass" else
            "rejected_internal_geometry_review_candidate"
            if status == "failed" else
            "pending_internal_geometry_review_candidate"),
        "quantitative_status": quantitative_status,
        "branches": branch_reports,
        "failures": failures,
        "pending": pending,
        "decision_id": trusted_decision["decision_id"],
        "actor_id": actor["id"],
        "scope": _safe_raw(scope),
        "source_ref": _safe_raw(source_ref),
        "authority": False,
        "formal_acceptance": False,
        "scientific_accepted": False,
    }


__all__ = ["VERSION", "ROOT", "SOURCE_RELATIVE_PATH", "SOURCE_SHA256",
           "SOURCE_BYTES", "RULES", "load_followup", "evaluate_calibration",
           "evaluate_novel", "evaluate_geometry_review"]

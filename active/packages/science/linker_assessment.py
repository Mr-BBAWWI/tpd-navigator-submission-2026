"""확장 linker 설계 가설과 보수적 그래프·경로 근거 평가.

그래프 길이, RDKit descriptor, 제한된 3D conformer 표본 및 문헌 기록의
완전성을 구분한다. 어떤 결과도 결합, ternary complex 형성, 합성 가능성 또는
효능을 입증하지 않는다.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdDistGeom, rdMolDescriptors

from . import linker_design

VERSION = "linker-assessment/20261001.1"
_LENGTH_LABELS = {"short": "상대 짧음", "medium": "상대 중간", "long": "상대 김"}
_FAMILY_ORDER = ("PEG", "alkyl", "PEG-alkyl", "piperazine", "triazole", "rigid-aromatic")

# 결합 길이의 평형값이 아니라, hard rejection에만 쓰는 의도적으로 넉넉한
# 원자·결합 유형별 상한이다. 실제 가능한 최대 신장 또는 에너지 최소값이 아니다.
_BOND_UPPER_A = {
    "TRIPLE": 1.45,
    "DOUBLE": 1.60,
    "AROMATIC": 1.70,
    "SINGLE": 1.90,
}
_TERMINAL_ATTACHMENT_BOND_UPPER_A = 2.20
_CONTOUR_UNCERTAINTY_A = 1.00
_MAX_CONFORMERS = 12
_ETKDG_SEED = 0x5EED


def _new_template(
    template_id: str,
    family: str,
    length: str,
    smiles: str,
    source_ids: list[str],
    state_review: str,
) -> dict[str, Any]:
    family_label = {
        "PEG": "PEG",
        "alkyl": "Alkyl",
        "PEG-alkyl": "PEG–alkyl",
        "piperazine": "Piperazine",
        "triazole": "Triazole",
        "rigid-aromatic": "Rigid aryl",
    }[family]
    return {
        "id": template_id,
        "label": f"{family_label} · {_LENGTH_LABELS[length]}",
        "family": family,
        "smiles": smiles,
        "rationale": (
            f"{family_label} 계열의 { _LENGTH_LABELS[length] } 비교용 설계 가설입니다. "
            "source_ids는 계열 수준의 문헌 맥락만 나타내며 이 정확한 구조, 합성 경로, "
            "결합 또는 효능의 선례를 뜻하지 않습니다."
        ),
        "source_ids": list(source_ids),
        "state_review": state_review,
    }


_VARIANTS = (
    _new_template("peg_short", "PEG", "short", "[*:1001]CCOCC[*:1002]", ["farnaby2019", "piperazine2022"],
                  "말단 결합 환경과 ether 배치 검토 필요"),
    _new_template("peg_long", "PEG", "long", "[*:1001]CCOCCOCCOCCOCCOCC[*:1002]", ["farnaby2019", "piperazine2022"],
                  "높은 유연성과 극성, 말단 결합 환경 검토 필요"),
    _new_template("alkyl_short", "alkyl", "short", "[*:1001]CCCC[*:1002]", ["piperazine2022"],
                  "비극성 사슬의 길이와 말단 결합 환경 검토 필요"),
    _new_template("alkyl_long", "alkyl", "long", "[*:1001]CCCCCCCCC[*:1002]", ["piperazine2022"],
                  "긴 비극성 사슬의 유연성과 전체 물성 검토 필요"),
    _new_template("peg_alkyl_short", "PEG-alkyl", "short", "[*:1001]CCOCCC[*:1002]", ["piperazine2022"],
                  "비대칭 방향과 말단별 반응 맥락 검토 필요"),
    _new_template("peg_alkyl_long", "PEG-alkyl", "long", "[*:1001]CCOCCCCCCOCC[*:1002]", ["piperazine2022"],
                  "비대칭 방향, 유연성 및 ether 배치 검토 필요"),
    _new_template("piperazine_short", "piperazine", "short", "[*:1001]CN1CCN(C[*:1002])CC1", ["piperazine2022"],
                  "질소의 protonation과 치환 상태 검토 필요"),
    _new_template("piperazine_long", "piperazine", "long", "[*:1001]CCCN1CCN(CCC[*:1002])CC1", ["piperazine2022"],
                  "질소의 protonation, 유연성 및 치환 상태 검토 필요"),
    _new_template("triazole_short", "triazole", "short", "[*:1001]Cn1cc(C[*:1002])nn1", ["triazole2024"],
                  "고리 연결 위치, tautomer/protonation 및 반응 맥락 검토 필요"),
    _new_template("triazole_long", "triazole", "long", "[*:1001]CCCn1cc(CCC[*:1002])nn1", ["triazole2024"],
                  "고리 연결 위치, 유연성 및 반응 맥락 검토 필요"),
    _new_template("aryl_short", "rigid-aromatic", "short", "[*:1001]Cc1ccc(C[*:1002])cc1", ["farnaby2019"],
                  "방향족 attachment geometry와 입체 충돌 검토 필요"),
    _new_template("aryl_long", "rigid-aromatic", "long", "[*:1001]CCCc1ccc(CCC[*:1002])cc1", ["farnaby2019"],
                  "방향족 geometry와 양쪽 사슬 유연성 검토 필요"),
)

_ORIGINAL_LENGTH = {
    "peg3": "medium",
    "alkyl_c6": "medium",
    "peg_alkyl": "medium",
    "piperazine": "medium",
    "triazole": "medium",
    "aryl_para": "medium",
}


def _json_safe(value: Any) -> Any:
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _template_mol(template: dict[str, Any]) -> Chem.Mol:
    if type(template) is not dict:
        raise TypeError("template은 dict여야 합니다")
    smiles = template.get("smiles")
    if type(smiles) is not str or not smiles:
        raise ValueError("LINKER_TEMPLATE_SMILES")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("LINKER_TEMPLATE_INVALID")
    dummies = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    maps = [atom.GetAtomMapNum() for atom in dummies]
    if sorted(maps) != [1001, 1002] or len(set(maps)) != 2:
        raise ValueError("LINKER_TEMPLATE_DUMMY_MAPS")
    if any(atom.GetDegree() != 1 for atom in dummies):
        raise ValueError("LINKER_TEMPLATE_DUMMY_TERMINAL")
    if any(atom.GetNeighbors()[0].GetAtomicNum() == 0 for atom in dummies):
        raise ValueError("LINKER_TEMPLATE_EMPTY_ENDPOINT")
    heavy_maps = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 0]
    if any(atom_map != 0 for atom_map in heavy_maps):
        raise ValueError("LINKER_TEMPLATE_HEAVY_MAPS")
    return mol


def validate_library(data: dict[str, Any]) -> dict[str, Any]:
    """확장 library의 ID, 구조, map 및 family/상대 길이 구성을 검증한다."""
    if type(data) is not dict or type(data.get("templates")) is not list:
        raise ValueError("LINKER_LIBRARY_SCHEMA")
    templates = data["templates"]
    ids = [item.get("id") if isinstance(item, dict) else None for item in templates]
    if any(type(item) is not str or not item for item in ids):
        raise ValueError("LINKER_TEMPLATE_ID")
    if len(ids) != len(set(ids)):
        raise ValueError("LINKER_TEMPLATE_DUPLICATE_ID")
    source_ids = {source.get("id") for source in data.get("sources", []) if isinstance(source, dict)}
    family_lengths: dict[str, set[str]] = {family: set() for family in _FAMILY_ORDER}
    for template in templates:
        _template_mol(template)
        if set(template) != {"id", "label", "family", "smiles", "rationale", "source_ids", "state_review"}:
            raise ValueError("LINKER_TEMPLATE_SHAPE")
        if template["family"] not in family_lengths:
            raise ValueError("LINKER_TEMPLATE_FAMILY")
        if type(template["source_ids"]) is not list or not set(template["source_ids"]) <= source_ids:
            raise ValueError("LINKER_TEMPLATE_SOURCE")
        matching = [key for key, text in _LENGTH_LABELS.items() if text in template["label"]]
        if len(matching) != 1:
            raise ValueError("LINKER_TEMPLATE_RELATIVE_LENGTH")
        family_lengths[template["family"]].add(matching[0])
    if any(lengths != set(_LENGTH_LABELS) for lengths in family_lengths.values()):
        raise ValueError("LINKER_LIBRARY_FAMILY_LENGTH_COVERAGE")
    return {
        "valid": True,
        "template_count": len(templates),
        "families": {family: sorted(family_lengths[family]) for family in _FAMILY_ORDER},
        "version": VERSION,
    }


def expanded_library() -> dict[str, Any]:
    """원본 6개 구조와 ID를 유지하면서 6계열×3 상대 길이를 반환한다."""
    original = linker_design.library()
    originals: list[dict[str, Any]] = []
    precedence: dict[str, Any] = {}
    for source_template in original["templates"]:
        template = copy.deepcopy(source_template)
        length = _ORIGINAL_LENGTH.get(template["id"])
        if length is None:
            raise ValueError("ORIGINAL_LINKER_ID_CHANGED")
        family_label = template["label"].split("·", 1)[0].strip()
        template["label"] = f"{family_label} · {_LENGTH_LABELS[length]}"
        originals.append(template)
        precedence[template["id"]] = {
            "design_status": "설계 가설",
            "family_level_source_ids": list(template["source_ids"]),
            "exact_structure_precedence": "검토되지 않음",
        }
    for template in _VARIANTS:
        precedence[template["id"]] = {
            "design_status": "생성된 설계 가설",
            "family_level_source_ids": list(template["source_ids"]),
            "exact_structure_precedence": "확인된 근거 없음",
        }
    result = {
        "version": VERSION,
        "scope": (
            "상대 길이별 linker 설계 가설입니다. 계열 수준 문헌 맥락과 정확한 구조의 "
            "선례를 분리하며 합성, 결합, ternary geometry 또는 효능을 검증하지 않습니다."
        ),
        "review_source_sha256": original.get("review_source_sha256"),
        "sources": copy.deepcopy(original.get("sources", [])),
        "templates": originals + [copy.deepcopy(item) for item in _VARIANTS],
        "precedence": precedence,
        "provenance": {
            "module_version": VERSION,
            "base_library_version": original.get("version"),
            "original_ids_and_smiles_preserved": True,
        },
    }
    validate_library(result)
    return _json_safe(result)


def _endpoints(mol: Chem.Mol) -> tuple[int, int]:
    by_map = {
        atom.GetAtomMapNum(): atom.GetNeighbors()[0].GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 0
    }
    return by_map[1001], by_map[1002]


def _bond_upper(bond: Chem.Bond) -> float:
    key = "AROMATIC" if bond.GetIsAromatic() else str(bond.GetBondType())
    base = _BOND_UPPER_A.get(key, 2.10)
    atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
    if any(atom.GetAtomicNum() in (15, 16) for atom in atoms):
        base = max(base, 2.25)
    elif any(atom.GetAtomicNum() > 16 for atom in atoms):
        base = max(base, 2.35)
    return base


def _maximum_simple_contour(mol: Chem.Mol, start: int, end: int) -> tuple[float, int]:
    heavy_count = sum(atom.GetAtomicNum() > 0 for atom in mol.GetAtoms())
    if heavy_count > 128:
        raise ValueError("LINKER_GRAPH_TOO_LARGE")
    best_length = -1.0
    best_bonds = -1

    def walk(current: int, visited: set[int], length: float, bonds: int) -> None:
        nonlocal best_length, best_bonds
        if current == end:
            if length > best_length:
                best_length, best_bonds = length, bonds
            return
        atom = mol.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            index = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 0 or index in visited:
                continue
            walk(index, visited | {index}, length + _bond_upper(bond), bonds + 1)

    walk(start, {start}, 0.0, 0)
    if best_length < 0:
        raise ValueError("LINKER_ENDPOINTS_DISCONNECTED")
    return best_length, best_bonds


def _fragment_without_dummies(mol: Chem.Mol) -> tuple[Chem.Mol, tuple[int, int]]:
    endpoints = _endpoints(mol)
    dummy_indices = sorted(
        (atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0), reverse=True
    )
    editable = Chem.RWMol(mol)
    for index in dummy_indices:
        editable.RemoveAtom(index)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)

    def shifted(index: int) -> int:
        return index - sum(dummy < index for dummy in dummy_indices)

    return fragment, (shifted(endpoints[0]), shifted(endpoints[1]))


def _descriptors(fragment: Chem.Mol) -> dict[str, Any]:
    return {
        "heavy_atom_count": int(fragment.GetNumHeavyAtoms()),
        "rotatable_bonds_rdkit_strict": int(Lipinski.NumRotatableBonds(fragment)),
        "molecular_weight_g_mol": round(float(Descriptors.MolWt(fragment)), 4),
        "tpsa_A2": round(float(rdMolDescriptors.CalcTPSA(fragment)), 4),
        "clogp_rdkit": round(float(Crippen.MolLogP(fragment)), 4),
        "formal_charge": int(Chem.GetFormalCharge(fragment)),
        "meaning": "분리된 linker 그래프의 계산 descriptor이며 전체 분자 물성이나 측정값이 아닙니다.",
    }


def _sample_3d(fragment: Chem.Mol, endpoints: tuple[int, int]) -> dict[str, Any]:
    sampled = Chem.AddHs(Chem.Mol(fragment))
    parameters = rdDistGeom.ETKDGv3()
    parameters.randomSeed = _ETKDG_SEED
    parameters.pruneRmsThresh = 0.25
    parameters.numThreads = 1
    try:
        conformer_ids = list(rdDistGeom.EmbedMultipleConfs(sampled, numConfs=_MAX_CONFORMERS, params=parameters))
    except Exception as exc:  # RDKit failure is reported as missing computation, not synthetic evidence.
        return {
            "status": "계산되지 않음",
            "reason": f"ETKDG embedding 실패: {type(exc).__name__}",
            "sample_count": 0,
            "method": "RDKit ETKDGv3",
        }
    distances: list[float] = []
    for conformer_id in conformer_ids:
        conformer = sampled.GetConformer(conformer_id)
        distance = float(conformer.GetAtomPosition(endpoints[0]).Distance(conformer.GetAtomPosition(endpoints[1])))
        if math.isfinite(distance):
            distances.append(distance)
    if not distances:
        return {
            "status": "계산되지 않음",
            "reason": "유효한 conformer 거리 없음",
            "sample_count": 0,
            "method": "RDKit ETKDGv3",
        }
    return {
        "status": "계산됨",
        "sample_count": len(distances),
        "seed": _ETKDG_SEED,
        "conformer_cap": _MAX_CONFORMERS,
        "end_to_end_range_A": {
            "min_observed": round(min(distances), 4),
            "max_observed": round(max(distances), 4),
        },
        "method": "RDKit ETKDGv3; terminal linker heavy atoms",
        "interpretation": "표본에서 관찰된 범위이며 max_observed는 이론적 최대 길이가 아닙니다.",
    }


def _endpoint_hypotheses(endpoint_chemistry: Any) -> dict[str, Any]:
    if endpoint_chemistry is None:
        return {
            "status": "미지",
            "declared": None,
            "proposals": [],
            "interpretation": "endpoint chemistry가 없어 화학적 호환성과 반응 유형을 판단하지 않습니다.",
        }
    if type(endpoint_chemistry) is str:
        names = [endpoint_chemistry.lower()]
    elif type(endpoint_chemistry) is dict:
        names = [str(value).lower() for value in endpoint_chemistry.values() if value is not None]
    elif isinstance(endpoint_chemistry, (list, tuple)):
        names = [str(value).lower() for value in endpoint_chemistry]
    else:
        raise TypeError("endpoint_chemistry는 str, dict, list 또는 None이어야 합니다")
    joined = " ".join(names)
    proposals: list[dict[str, str]] = []
    if any(term in joined for term in ("amine", "아민")) and any(
        term in joined for term in ("carbox", "acid", "activated ester", "카복")
    ):
        proposals.append({"type": "amidation", "status": "반응 유형 가설"})
    if any(term in joined for term in ("azide", "아자이드")) and any(
        term in joined for term in ("alkyne", "알카인")
    ):
        proposals.append({"type": "azide–alkyne click (CuAAC 또는 조건별 대안)", "status": "반응 유형 가설"})
    if any(term in joined for term in ("halide", "bromide", "chloride", "tosylate", "mesylate")) and any(
        term in joined for term in ("amine", "alcohol", "thiol", "nucleophile", "아민", "알코올", "티올")
    ):
        proposals.append({"type": "SN2-type substitution", "status": "반응 유형 가설"})
    return {
        "status": "검토 필요",
        "declared": _json_safe(endpoint_chemistry),
        "proposals": proposals,
        "interpretation": (
            "제안은 명칭 기반 endpoint 호환성 가설일 뿐 chemoselectivity, 보호기, 조건, 수율 또는 "
            "실제 합성 경로를 검증하지 않습니다."
        ),
    }


def assess_linker(
    template: dict[str, Any],
    required_attachment_distance_A: float | None = None,
    endpoint_chemistry: Any = None,
) -> dict[str, Any]:
    """Linker를 평가한다. 보수적 contour 상한 초과만 hard rejection으로 사용한다."""
    mol = _template_mol(template)
    if required_attachment_distance_A is not None:
        if type(required_attachment_distance_A) not in (int, float) or isinstance(required_attachment_distance_A, bool):
            raise TypeError("required_attachment_distance_A는 유한한 양수여야 합니다")
        required_attachment_distance_A = float(required_attachment_distance_A)
        if not math.isfinite(required_attachment_distance_A) or required_attachment_distance_A < 0:
            raise ValueError("required_attachment_distance_A는 유한한 0 이상 값이어야 합니다")

    start, end = _endpoints(mol)
    shortest_path = Chem.GetShortestPath(mol, start, end)
    internal_contour, contour_bonds = _maximum_simple_contour(mol, start, end)
    contour_bound = (
        internal_contour
        + 2.0 * _TERMINAL_ATTACHMENT_BOND_UPPER_A
        + _CONTOUR_UNCERTAINTY_A
    )
    fragment, fragment_endpoints = _fragment_without_dummies(mol)
    descriptors = _descriptors(fragment)
    hard_rejected = (
        required_attachment_distance_A is not None
        and required_attachment_distance_A > contour_bound
    )
    flags: list[str] = []
    if descriptors["rotatable_bonds_rdkit_strict"] >= 8 or contour_bound >= 22.0:
        flags.append("긴/유연한 linker 검토 필요")
    if descriptors["formal_charge"] != 0:
        flags.append("명시적 전하 상태 검토 필요")
    if descriptors["tpsa_A2"] >= 75.0:
        flags.append("높은 linker 극성 검토 필요")
    if any(atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 0 for atom in fragment.GetAtoms()):
        flags.append("질소 protonation 상태 검토 필요")

    if required_attachment_distance_A is None:
        geometry = {
            "status": "미지",
            "reason": (
                "target ternary endpoint geometry가 제공되지 않았습니다. warhead-only receptor docking으로 "
                "ternary endpoint 거리를 추론하지 않습니다."
            ),
        }
    else:
        geometry = {
            "status": "사용자 선언 거리와 contour 상한만 비교",
            "declared_endpoint_distance_A": required_attachment_distance_A,
            "source_scope": "외부에서 선언된 값; 본 모듈이 ternary geometry를 추론하지 않음",
        }

    result = {
        "template_id": template.get("id"),
        "status": "hard_rejected_contour_bound" if hard_rejected else "설계 가설 · 검토 필요",
        "hard_rejected": hard_rejected,
        "hard_rejection_reason": (
            "선언된 endpoint 거리가 terminal attachment bond와 불확실성을 포함한 보수적 contour 상한을 초과함"
            if hard_rejected else None
        ),
        "graph_span": {
            "terminal_heavy_atom_shortest_path_bonds": len(shortest_path) - 1,
            "terminal_heavy_atom_shortest_path_atoms": len(shortest_path),
            "definition": "두 dummy에 인접한 linker heavy atom 사이의 최단 그래프 경로",
        },
        "contour_bound": {
            "max_simple_path_bonds": contour_bonds,
            "internal_element_bond_aware_upper_A": round(internal_contour, 4),
            "terminal_attachment_bonds_included": 2,
            "terminal_attachment_bond_upper_each_A": _TERMINAL_ATTACHMENT_BOND_UPPER_A,
            "uncertainty_A": _CONTOUR_UNCERTAINTY_A,
            "conservative_total_upper_A": round(contour_bound, 4),
            "meaning": "hard rejection용 넉넉한 그래프 상한이며 실제 conformer의 이론적 최대 또는 측정 거리가 아닙니다.",
        },
        "descriptors": descriptors,
        "sampled_3d": _sample_3d(fragment, fragment_endpoints),
        "ternary_geometry": geometry,
        "endpoint_chemistry": _endpoint_hypotheses(endpoint_chemistry),
        "soft_flags": flags,
        "scientific_scope": (
            "그래프 조립, descriptor 및 docking은 효능·승인·합성 가능성 증거가 아닙니다."
        ),
        "provenance": {
            "module_version": VERSION,
            "rdkit_version": getattr(Chem, "__version__", None),
            "template_smiles": template.get("smiles"),
        },
    }
    return _json_safe(result)


def ingest_source_json(source: str | bytes | Path | dict[str, Any] | list[Any]) -> list[dict[str, Any]]:
    """로컬 JSON 문자열/파일/객체를 읽는다. 네트워크 수집은 수행하지 않는다."""
    if isinstance(source, Path):
        data = json.loads(source.read_text(encoding="utf-8"))
    elif isinstance(source, bytes):
        data = json.loads(source.decode("utf-8"))
    elif isinstance(source, str):
        stripped = source.lstrip()
        if stripped.startswith("{") or stripped.startswith("["):
            data = json.loads(source)
        else:
            data = json.loads(Path(source).read_text(encoding="utf-8"))
    elif isinstance(source, (dict, list)):
        data = copy.deepcopy(source)
    else:
        raise TypeError("source는 JSON, 로컬 경로, dict 또는 list여야 합니다")
    if isinstance(data, dict) and "evidence_records" in data:
        data = data["evidence_records"]
    elif isinstance(data, dict) and "records" in data:
        data = data["records"]
    elif isinstance(data, dict):
        data = [data]
    if type(data) is not list or any(type(record) is not dict for record in data):
        raise ValueError("ROUTE_EVIDENCE_JSON_SCHEMA")
    return _json_safe(data)


def _substantive(value: Any) -> bool:
    if value is None or value is False:
        return False
    if isinstance(value, str):
        return bool(value.strip()) and value.strip().lower() not in {
            "unknown", "not reported", "n/a", "none", "미상", "없음"
        }
    if isinstance(value, (list, tuple, dict)):
        return bool(value)
    if type(value) in (int, float):
        return not isinstance(value, bool) and math.isfinite(float(value))
    return False


def _nested(record: dict[str, Any], key: str) -> Any:
    if key in record:
        return record[key]
    route = record.get("route")
    return route.get(key) if isinstance(route, dict) else None


def _valid_url(value: Any) -> bool:
    if type(value) is not str:
        return False
    parsed = urlparse(value)
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def validate_evidence_record(record: dict[str, Any]) -> dict[str, Any]:
    """단일 source record의 실제 기재 필드만 점검하고 누락값을 채우지 않는다."""
    if type(record) is not dict:
        raise TypeError("evidence record는 dict여야 합니다")
    source = record.get("source") if isinstance(record.get("source"), dict) else {}
    source_url = record.get("source_url", source.get("url"))
    locator = record.get("locator", source.get("locator"))
    compound_id = record.get("compound_id")
    steps = _nested(record, "steps")
    materials = _nested(record, "materials")
    yield_value = _nested(record, "yield")
    purification = _nested(record, "purification")
    characterization = _nested(record, "characterization")
    checks = {
        "compound_id": type(compound_id) is str and bool(compound_id.strip()),
        "source_url": _valid_url(source_url),
        "locator": _substantive(locator),
        "steps": isinstance(steps, list) and bool(steps) and all(_substantive(step) for step in steps),
        "materials": isinstance(materials, list) and bool(materials) and all(_substantive(item) for item in materials),
        "yield": _substantive(yield_value),
        "purification": _substantive(purification),
        "characterization": _substantive(characterization),
    }
    missing = [key for key, present in checks.items() if not present]
    return {
        "record_id": record.get("record_id", record.get("id")),
        "compound_id": compound_id,
        "source_url": source_url,
        "locator": locator,
        "field_presence": checks,
        "missing_or_invalid_fields": missing,
        "documentation_complete": not missing,
        "source_record": _json_safe(record),
        "interpretation": "필드 완전성 점검이며 독립적 재현, route validation 또는 합성 가능성 판정이 아닙니다.",
    }


def validate_source_json(source: str | bytes | Path | dict[str, Any] | list[Any]) -> dict[str, Any]:
    records = ingest_source_json(source)
    reports = [validate_evidence_record(record) for record in records]
    identifiers = [report["record_id"] for report in reports if report["record_id"] is not None]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("ROUTE_EVIDENCE_DUPLICATE_RECORD_ID")
    return {
        "valid_json_schema": True,
        "record_count": len(reports),
        "records": reports,
        "provenance": {"module_version": VERSION, "collection": "사용자 제공 로컬 JSON; 네트워크 수집 없음"},
    }


def route_evidence(
    evidence_records: Iterable[dict[str, Any]] | dict[str, Any] | str | bytes | Path,
    requested_compound_id: str,
    generated_candidate: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """정확한 compound ID route와 다른 화합물의 재사용 선례를 분리한다."""
    if type(requested_compound_id) is not str or not requested_compound_id.strip():
        raise ValueError("requested_compound_id는 비어 있지 않은 문자열이어야 합니다")
    if generated_candidate is not None and type(generated_candidate) is not dict:
        raise TypeError("generated_candidate는 dict 또는 None이어야 합니다")
    if isinstance(evidence_records, (str, bytes, Path, dict)):
        records = ingest_source_json(evidence_records)
    else:
        records = list(evidence_records)
        if any(type(record) is not dict for record in records):
            raise TypeError("모든 evidence record는 dict여야 합니다")
    reports = [validate_evidence_record(record) for record in records]
    exact = [report for report in reports if report["compound_id"] == requested_compound_id]
    different = [report for report in reports if report["compound_id"] != requested_compound_id]
    complete_exact = [report for report in exact if report["documentation_complete"]]
    reusable = [
        report for report in different
        if report["source_record"].get("reusable_precedent") is True
    ]
    other_compounds = [report for report in different if report not in reusable]

    transformations: Any = []
    if generated_candidate:
        transformations = generated_candidate.get(
            "novel_transformations", generated_candidate.get("transformations", generated_candidate.get("transformation", []))
        )
        if isinstance(transformations, dict):
            transformations = [transformations]
        elif transformations is None:
            transformations = []
        elif not isinstance(transformations, list):
            transformations = [transformations]

    if complete_exact:
        status = "정확한 compound ID의 필수 route 문서 필드가 존재함"
    elif exact:
        status = "정확한 compound ID 기록은 있으나 route 문서 필드가 불완전함"
    else:
        status = "정확한 compound ID route 근거 없음"
    return _json_safe({
        "requested_compound_id": requested_compound_id,
        "status": status,
        "exact_compound_records": exact,
        "complete_exact_documentation_count": len(complete_exact),
        "reusable_precedent_records": reusable,
        "other_compound_records": other_compounds,
        "novel_transformations": {
            "items": transformations,
            "evidence_status": "별도 검토 필요; 생성 후보 자체는 문헌 route 근거가 아님",
        },
        "conclusion": (
            "다른 source compound는 정확한 route를 충족하지 않습니다. 완전한 필드가 있어도 "
            "독립적 재현성이나 실제 합성 성공을 주장하지 않습니다."
        ),
        "provenance": {"module_version": VERSION, "network_collection": False},
    })


def assessment_summary(
    assessment: dict[str, Any],
    route: dict[str, Any] | None = None,
    sa_score: float | dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Linker descriptor 평가, 선택적 SA descriptor 및 route 기록을 분리해 요약한다."""
    if type(assessment) is not dict:
        raise TypeError("assessment는 dict여야 합니다")
    if route is not None and type(route) is not dict:
        raise TypeError("route는 dict 또는 None이어야 합니다")
    if isinstance(sa_score, bool):
        raise TypeError("sa_score는 bool일 수 없습니다")
    if type(sa_score) in (int, float):
        score = float(sa_score)
        if not math.isfinite(score):
            raise ValueError("sa_score는 유한해야 합니다")
        sa = {
            "value": score,
            "status": "외부 계산 descriptor",
            "interpretation": "SA score만으로 synthesizability를 주장할 수 없습니다.",
        }
    elif isinstance(sa_score, dict):
        sa = {
            "value": _json_safe(sa_score),
            "status": "외부 descriptor 기록",
            "interpretation": "SA 관련 descriptor는 실제 route, 수율 또는 synthesizability 증거가 아닙니다.",
        }
    elif sa_score is None:
        sa = {"value": None, "status": "계산/제공되지 않음"}
    else:
        raise TypeError("sa_score는 숫자, dict 또는 None이어야 합니다")
    return _json_safe({
        "linker_assessment": assessment,
        "descriptor_only_sa": sa,
        "route_evidence": route or {
            "status": "제공되지 않음",
            "interpretation": "route 근거가 없어 synthesizability는 미지입니다.",
        },
        "synthesizability_conclusion": "미지/검토 필요",
        "interpretation": (
            "descriptor, 3D 표본, 반응 유형 가설 및 route 문서 완전성은 서로 다른 근거입니다. "
            "이 요약은 합성 가능성, ternary 결합 또는 효능을 확정하지 않습니다."
        ),
        "provenance": {"module_version": VERSION},
    })


__all__ = [
    "VERSION",
    "expanded_library",
    "validate_library",
    "assess_linker",
    "ingest_source_json",
    "validate_source_json",
    "validate_evidence_record",
    "route_evidence",
    "assessment_summary",
]

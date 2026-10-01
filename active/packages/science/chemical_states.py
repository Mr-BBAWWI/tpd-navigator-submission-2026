"""보수적인 화학 상태 준비와 좌표 기반 상호작용 기술.

이 모듈은 가능한 토토머와 단일 산/염기 공액 상태를 열거하지만 pKa나
개체군을 예측하지 않는다. 기본 pH는 실행 문맥일 뿐 상태 선택 근거가 아니다.
수소 결합은 명시적으로 연결된 D-H 좌표와 방향 기하가 모두 있을 때만
기록한다. 나머지 근접 접촉은 검토 대상으로 남긴다.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem.MolStandardize import rdMolStandardize


_HBOND_DA_MAX = 3.5
_HBOND_HA_MAX = 2.6
_HBOND_ANGLE_MIN = 120.0
_HYDROPHOBIC_MAX = 4.5
_AROMATIC_MAX = 6.0
_CLASH_OVERLAP_MIN = 0.4
_ALLOWED_PDB2PQR = {"pdb2pqr", "pdb2pqr.exe"}

_PROTEIN_DONORS = {
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
    "LYS": {"NZ"}, "SER": {"OG"}, "THR": {"OG1"}, "TRP": {"NE1"},
    "TYR": {"OH"},
}
_PROTEIN_ACCEPTORS = {
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"},
    "GLN": {"OE1"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
    "MET": {"SD"},
}
_PROTEIN_AROMATIC = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HID": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIE": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIP": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
}


def _finite(value, name):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} 값은 유한한 수여야 합니다") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} 값은 유한한 수여야 합니다")
    return result


def _coordinates(mol, name="분자"):
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"{name}는 RDKit Mol이어야 합니다")
    if mol.GetNumConformers() != 1:
        raise ValueError(f"{name}에는 정확히 하나의 conformer가 필요합니다")
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if xyz.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
        raise ValueError(f"{name} 좌표가 없거나 유한하지 않습니다")
    return xyz


def _heavy_map_signature(mol):
    result = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        number = atom.GetAtomMapNum()
        if number <= 0:
            raise ValueError("모든 중원자는 고유한 양의 atom-map 번호가 필요합니다")
        if number in result:
            raise ValueError(f"중복 atom-map 번호: {number}")
        result[number] = (
            atom.GetAtomicNum(), atom.GetIsotope(), int(atom.GetChiralTag())
        )
    if not result:
        raise ValueError("분자에 중원자가 없습니다")
    return result


def _stereo_bonds(mol):
    result = {}
    for bond in mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREONONE:
            continue
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a > 0 and b > 0:
            result[tuple(sorted((a, b)))] = int(bond.GetStereo())
    return result


def _state_key(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def _valid_state(candidate, signature, stereo_bonds):
    try:
        Chem.SanitizeMol(candidate)
        Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
    except Exception:
        return False
    if len(Chem.GetMolFrags(candidate)) != 1:
        return False
    try:
        current = _heavy_map_signature(candidate)
    except ValueError:
        return False
    if current != signature:
        return False
    charge = Chem.GetFormalCharge(candidate)
    if not isinstance(charge, int) or abs(charge) > candidate.GetNumAtoms() + 8:
        return False
    maps = {a.GetAtomMapNum(): a.GetIdx() for a in candidate.GetAtoms() if a.GetAtomicNum() > 1}
    for pair, stereo in stereo_bonds.items():
        bond = candidate.GetBondBetweenAtoms(maps[pair[0]], maps[pair[1]])
        if bond is None or int(bond.GetStereo()) != stereo:
            return False
    return True


def _feature_indices(mol):
    donors = {index for match in Lipinski._HDonors(mol) for index in match}
    acceptors = {index for match in Lipinski._HAcceptors(mol) for index in match}
    return donors, acceptors


def _hydrogen_count(mol):
    """명시적 H 원자와 원자 속성의 H를 중복 없이 센다."""
    mol.UpdatePropertyCache(strict=False)
    explicit_atoms = sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    attached = sum(
        int(atom.GetTotalNumHs(includeNeighbors=False))
        for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    )
    return explicit_atoms + attached


def _change_proton(mol, atom_index, delta):
    """한 원자의 H/전하를 정확히 한 단계 바꾼 사본을 만든다."""
    before_h = _hydrogen_count(Chem.Mol(mol))
    before_charge = Chem.GetFormalCharge(mol)
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(atom_index)
    if delta == 1:
        if atom.GetAtomicNum() not in {7, 8, 15, 16}:
            return None
        # NoImplicit을 켜기 전에 기존 implicit H를 explicit-H 속성으로 보존한다.
        existing = int(atom.GetNumExplicitHs()) + int(atom.GetNumImplicitHs())
        atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        atom.SetNumExplicitHs(existing + 1)
        atom.SetNoImplicit(True)
    elif delta == -1:
        hydrogen_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
        atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        if hydrogen_neighbors:
            # 전하를 먼저 바꾸면 H 원자의 인덱스가 중심 원자보다 앞선 경우도 안전하다.
            rw.RemoveAtom(hydrogen_neighbors[0])
        else:
            existing = int(atom.GetNumExplicitHs()) + int(atom.GetNumImplicitHs())
            if existing <= 0:
                return None
            atom.SetNumExplicitHs(existing - 1)
            atom.SetNoImplicit(True)
    else:
        raise ValueError("delta는 +1 또는 -1이어야 합니다")
    result = rw.GetMol()
    result.UpdatePropertyCache(strict=False)
    if (_hydrogen_count(result) != before_h + delta or
            Chem.GetFormalCharge(result) != before_charge + delta):
        return None
    return result


_ACIDIC_H_PATTERNS = tuple(filter(None, (
    Chem.MolFromSmarts("[O;H1]-[C,S,P](=[O,S])"),
    Chem.MolFromSmarts("[O;H1]-[c]"),
    Chem.MolFromSmarts("[S;H1]"),
    Chem.MolFromSmarts("[nH]"),
)))


def _neutral_acidic_indices(mol):
    result = set()
    for pattern in _ACIDIC_H_PATTERNS:
        result.update(match[0] for match in mol.GetSubstructMatches(pattern))
    return result


def _conjugates(mol):
    donors, acceptors = _feature_indices(mol)
    acidic = _neutral_acidic_indices(mol)
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        symbol = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        if symbol in {"N", "O", "P", "S"} and (index in acceptors or charge < 0):
            candidate = _change_proton(mol, index, 1)
            if candidate is not None:
                yield candidate, "conjugate_acid", atom.GetAtomMapNum()
        has_h = atom.GetTotalNumHs(includeNeighbors=True) > 0
        # 중성 상태는 명시적으로 산성인 작용기만 탈양성자화한다. 양이온의
        # 중성화는 donor typing과 무관하게 H가 있는 양전하 헤테로원자에 허용한다.
        can_deprotonate = charge > 0 or (charge == 0 and index in acidic)
        if symbol in {"N", "O", "P", "S"} and has_h and can_deprotonate:
            candidate = _change_proton(mol, index, -1)
            if candidate is not None:
                yield candidate, "conjugate_base", atom.GetAtomMapNum()


def enumerate_microstates(mapped_mol, max_states=16, pH=7.4):
    """토토머와 명시적 단일 산/염기 공액 상태를 제한적으로 열거한다.

    반환 목록의 첫 항목은 입력 화학 상태의 사본이다. pH는 각 분자 속성에
    문맥으로 기록될 뿐, pKa 또는 개체군이 없으므로 상태를 선택하거나 순위를
    매기지 않는다. 각 결과의 ``chemical_state_origin`` 속성이 생성 경로를 담는다.
    """
    if not isinstance(mapped_mol, Chem.Mol):
        raise TypeError("mapped_mol은 RDKit Mol이어야 합니다")
    if isinstance(max_states, bool) or not isinstance(max_states, int) or max_states < 1:
        raise ValueError("max_states는 1 이상의 정수여야 합니다")
    context_ph = _finite(pH, "pH")
    source = Chem.Mol(mapped_mol)
    signature = _heavy_map_signature(source)
    stereo = _stereo_bonds(source)
    if not _valid_state(source, signature, stereo):
        raise ValueError("입력 분자는 연결되고 sanitize 가능한 유효 상태여야 합니다")

    states = []
    seen = set()

    def add(candidate, origin, parent=None, changed_map=None):
        candidate = Chem.Mol(candidate)
        if not _valid_state(candidate, signature, stereo):
            return False
        key = _state_key(candidate)
        if key in seen:
            return False
        seen.add(key)
        candidate.SetProp("chemical_state_origin", origin)
        candidate.SetProp("chemical_state_pH_context", format(context_ph, ".8g"))
        candidate.SetProp("chemical_state_population", "unknown")
        candidate.SetProp("chemical_state_pKa", "unknown")
        if parent is not None:
            candidate.SetProp("chemical_state_parent_key", parent)
        if changed_map:
            candidate.SetIntProp("chemical_state_changed_atom_map", int(changed_map))
        states.append(candidate)
        return True

    add(source, "source_ionization_state_retained")
    enumerator = rdMolStandardize.TautomerEnumerator()
    try:
        tautomers = enumerator.Enumerate(source)
    except Exception as exc:
        raise ValueError("RDKit 토토머 열거에 실패했습니다") from exc
    for tautomer in tautomers:
        if len(states) >= max_states:
            break
        add(tautomer, "rdkit_tautomer")

    cursor = 0
    while cursor < len(states) and len(states) < max_states:
        parent = states[cursor]
        parent_key = _state_key(parent)
        for candidate, origin, atom_map in _conjugates(parent):
            add(candidate, origin, parent_key, atom_map)
            if len(states) >= max_states:
                break
        cursor += 1
    return states


def optimize_ligand_hydrogens(mol, seed=23):
    """중원자를 고정한 force field로 추가 수소 좌표만 완화한다.

    반환값은 ``(structures, receipt)``이며 structures에는 입력 사본,
    수소 추가 직후 구조, 수소 최적화 구조가 각각 들어간다. receipt는 JSON으로
    직렬화할 수 있다. 이 절차는 결합 효능이나 올바른 양성자화 상태를 뜻하지 않는다.
    """
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed는 정수여야 합니다")
    original_xyz = _coordinates(mol, "mol")
    if mol.GetNumAtoms() == 0:
        raise ValueError("빈 분자는 최적화할 수 없습니다")

    initial = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    initial_xyz = _coordinates(initial, "수소 추가 구조")
    heavy = [a.GetIdx() for a in initial.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) != mol.GetNumHeavyAtoms():
        raise ValueError("수소 추가 중 중원자 수가 바뀌었습니다")
    if not np.array_equal(initial_xyz[heavy], original_xyz[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]):
        conf = initial.GetConformer()
        source_heavy = original_xyz[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]
        for index, point in zip(heavy, source_heavy):
            conf.SetAtomPosition(index, tuple(float(x) for x in point))
        initial_xyz = _coordinates(initial, "수소 추가 구조")

    optimized = Chem.Mol(initial)
    forcefield_name = None
    try:
        from rdkit.Chem import AllChem
        properties = AllChem.MMFFGetMoleculeProperties(optimized, mmffVariant="MMFF94s")
        forcefield = None if properties is None else AllChem.MMFFGetMoleculeForceField(optimized, properties)
        if forcefield is not None:
            forcefield_name = "RDKit MMFF94s"
        else:
            forcefield = AllChem.UFFGetMoleculeForceField(optimized)
            if forcefield is not None:
                forcefield_name = "RDKit UFF"
    except Exception as exc:
        raise ValueError("수소 최적화용 force field를 만들 수 없습니다") from exc
    if forcefield is None:
        raise ValueError("이 분자에 적용 가능한 MMFF94s 또는 UFF 매개변수가 없습니다")
    for index in heavy:
        forcefield.AddFixedPoint(index)
    forcefield.Initialize()
    status = int(forcefield.Minimize(maxIts=500))
    raw_xyz = _coordinates(optimized, "최적화 구조")
    raw_displacement = float(np.max(np.linalg.norm(raw_xyz[heavy] - initial_xyz[heavy], axis=1)))

    # 고정점 구현의 수치 잡음도 downstream 좌표를 바꾸지 않도록 원 좌표를 복원한다.
    conf = optimized.GetConformer()
    for index in heavy:
        point = initial_xyz[index]
        conf.SetAtomPosition(index, tuple(float(x) for x in point))
    final_xyz = _coordinates(optimized, "최종 구조")
    final_displacement = float(np.max(np.linalg.norm(final_xyz[heavy] - initial_xyz[heavy], axis=1)))
    if not np.array_equal(final_xyz[heavy], initial_xyz[heavy]):
        raise RuntimeError("중원자 좌표를 정확히 보존하지 못했습니다")

    structures = {
        "input": Chem.Mol(mol),
        "hydrogens_initial": initial,
        "hydrogens_optimized": optimized,
    }
    coordinate_drift_review = raw_displacement > 1.0e-6
    receipt = {
        "status": ("완료" if status == 0 and not coordinate_drift_review
                   else "반복한도_또는_고정점이탈_검토필요"),
        "requires_review": status != 0 or coordinate_drift_review,
        "forcefield": forcefield_name,
        "seed": seed,
        "seed_used": False,
        "seed_note": "기존 중원자 좌표에 수소를 추가했으므로 무작위 3D 임베딩을 수행하지 않았습니다.",
        "heavy_atom_count": len(heavy),
        "added_hydrogen_count": optimized.GetNumAtoms() - len(heavy),
        "fixed_heavy_atom_indices_zero_based": heavy,
        "raw_forcefield_heavy_max_displacement_A": raw_displacement,
        "raw_forcefield_heavy_displacement_tolerance_A": 1.0e-6,
        "raw_forcefield_fixed_points_within_tolerance": not coordinate_drift_review,
        "final_heavy_max_displacement_A": final_displacement,
        "heavy_coordinates_exactly_preserved": True,
        "heavy_coordinate_preservation_method": "최적화 후 입력 중원자 좌표를 정확히 복원",
        "minimizer_return_code": status,
        "limitations": [
            "수소 배치와 국소 force-field 완화만 수행하며 pKa 또는 개체군을 예측하지 않습니다.",
            "생성 구조는 결합, 효능, 선택성 또는 실험 성공을 입증하지 않습니다.",
        ],
    }
    json.dumps(receipt, ensure_ascii=False, allow_nan=False)
    return structures, receipt


def _protein_id(record):
    return ":".join((
        str(record["label_asym_id"]), str(record["auth_seq_id"]),
        str(record["label_comp_id"]).upper(), str(record["label_atom_id"]).upper(),
    ))


def _clean_protein(records, allow_hydrogen=False, name="protein_atoms"):
    if records is None:
        return []
    if not isinstance(records, (list, tuple)):
        raise TypeError(f"{name}는 원자 사전의 목록이어야 합니다")
    required = {"xyz", "type_symbol", "label_comp_id", "auth_seq_id", "label_atom_id", "label_asym_id"}
    result = []
    seen = set()
    for position, source in enumerate(records):
        if not isinstance(source, dict):
            raise TypeError(f"{name}[{position}]는 사전이어야 합니다")
        missing = required - set(source)
        if missing:
            raise ValueError(f"{name}[{position}] 필드 누락: {sorted(missing)}")
        record = dict(source)
        record["type_symbol"] = str(record["type_symbol"]).strip().title()
        atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(record["type_symbol"])
        if atomic_number <= 0 or (atomic_number == 1 and not allow_hydrogen):
            raise ValueError(f"{name}[{position}] 원소가 허용되지 않습니다")
        xyz = np.asarray(record["xyz"], dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError(f"{name}[{position}] 좌표가 유한한 3차원 값이 아닙니다")
        record["xyz"] = xyz
        identifier = _protein_id(record)
        if identifier in seen:
            raise ValueError(f"중복 단백질 원자 식별자: {identifier}")
        seen.add(identifier)
        result.append(record)
    return result


def _protein_role(record):
    residue = str(record["label_comp_id"]).upper()
    atom = str(record["label_atom_id"]).upper()
    donor = atom == "N" and residue != "PRO"
    acceptor = atom in {"O", "OXT"}
    donor |= atom in _PROTEIN_DONORS.get(residue, set())
    acceptor |= atom in _PROTEIN_ACCEPTORS.get(residue, set())
    uncertain = False
    reasons = []
    if residue in {"HIS", "HID", "HIE", "HIP"} and atom in {"ND1", "NE2"}:
        observed = record.get("observed_attached_hydrogens")
        if observed is not None:
            donor = bool(observed)
            other = record.get("histidine_other_nitrogen_has_hydrogen")
            acceptor = not donor and other is not True
            uncertain = other is None
            if uncertain:
                reasons.append("관측 수소는 있으나 히스티딘 두 질소의 전체 상태가 불완전합니다.")
        elif residue == "HID":
            donor, acceptor = atom == "ND1", atom == "NE2"
        elif residue == "HIE":
            donor, acceptor = atom == "NE2", atom == "ND1"
        elif residue == "HIP":
            donor, acceptor = True, False
        else:
            donor = acceptor = True
            uncertain = True
            reasons.append("HIS 호변이성질체/양성자화 상태가 선언되지 않았습니다.")
    if residue in {"CYS", "CYM", "CYX"} and atom == "SG":
        if residue == "CYM":
            donor, acceptor = False, True
        elif residue == "CYX":
            donor = acceptor = False
        else:
            donor = acceptor = True
            uncertain = True
            reasons.append("CYS 티올/티올레이트/이황화 상태가 확정되지 않았습니다.")
    if atom in {"N", "OXT"} and record.get("terminus_state") in (None, "", "unknown"):
        uncertain = True
        reasons.append("말단 상태가 명시되지 않았습니다.")
    return donor, acceptor, uncertain, reasons


def _declared_charge(record):
    if "formal_charge" in record and record["formal_charge"] not in (None, ""):
        value = _finite(record["formal_charge"], "protein formal_charge")
        if not value.is_integer():
            raise ValueError("protein formal_charge는 정수여야 합니다")
        return int(value), "atom_formal_charge"
    state = str(record.get("residue_state", record.get("protonation_state", ""))).upper()
    residue = str(record["label_comp_id"]).upper()
    atom = str(record["label_atom_id"]).upper()
    positive = {"LYS+": {"NZ"}, "ARG+": {"NE", "NH1", "NH2"}, "HIP": {"ND1", "NE2"}}
    negative = {"ASP-": {"OD1", "OD2"}, "GLU-": {"OE1", "OE2"}, "CYM": {"SG"}}
    effective = state or (residue if residue in {"HIP", "CYM"} else "")
    if atom in positive.get(effective, set()):
        return 1, "declared_residue_state"
    if atom in negative.get(effective, set()):
        return -1, "declared_residue_state"
    return 0, None


def _angle(donor, hydrogen, acceptor):
    first = donor - hydrogen
    second = acceptor - hydrogen
    denominator = np.linalg.norm(first) * np.linalg.norm(second)
    if denominator <= 0 or not math.isfinite(float(denominator)):
        return None
    cosine = float(np.clip(np.dot(first, second) / denominator, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def _ligand_id(atom):
    atom_map = atom.GetAtomMapNum()
    return f"L:{atom_map}" if atom_map > 0 else f"L:index{atom.GetIdx()}"


def _interaction(kind, participants, **values):
    identifier = kind + "|" + "|".join(sorted(participants))
    return {"id": identifier, "kind": kind, "participants": participants, **values}


def _protein_hydrogen_parents(hydrogens, protein_by_id):
    result = {}
    uncertainties = []
    for hydrogen in hydrogens:
        parent = hydrogen.get("parent_atom_id", hydrogen.get("donor_atom_id", hydrogen.get("bonded_to")))
        parent_id = None
        if isinstance(parent, str) and parent in protein_by_id:
            parent_id = parent
        elif isinstance(parent, dict):
            try:
                parent_id = _protein_id(parent)
            except KeyError:
                parent_id = None
        if parent_id not in protein_by_id:
            uncertainties.append({
                "kind": "unassigned_protein_hydrogen",
                "hydrogen": _protein_id(hydrogen),
                "description": "명시적인 단백질 D-H 연결 정보를 확인할 수 없습니다.",
            })
            continue
        distance = float(np.linalg.norm(hydrogen["xyz"] - protein_by_id[parent_id]["xyz"]))
        if distance > 1.35:
            uncertainties.append({
                "kind": "invalid_declared_DH_geometry", "hydrogen": _protein_id(hydrogen),
                "donor": parent_id, "distance_A": distance,
                "description": "선언된 D-H 거리가 보수적 공유결합 범위를 벗어납니다.",
            })
            continue
        result.setdefault(parent_id, []).append(hydrogen)
    return result, uncertainties


def _normal(points):
    centered = points - points.mean(axis=0)
    if len(points) < 3 or np.linalg.matrix_rank(centered) < 2:
        return None
    _, _, vectors = np.linalg.svd(centered, full_matrices=False)
    normal = vectors[-1]
    norm = np.linalg.norm(normal)
    return None if norm == 0 else normal / norm


def interaction_profile(ligand_with_coords, protein_atoms, protein_hydrogens=None):
    """주어진 좌표에서 보수적인 비공유 상호작용 기하를 기술한다."""
    ligand_xyz = _coordinates(ligand_with_coords, "ligand_with_coords")
    protein = _clean_protein(protein_atoms)
    if not protein:
        raise ValueError("protein_atoms에는 하나 이상의 중원자가 필요합니다")
    hydrogens = _clean_protein(protein_hydrogens, allow_hydrogen=True, name="protein_hydrogens")
    if any(Chem.GetPeriodicTable().GetAtomicNumber(r["type_symbol"]) != 1 for r in hydrogens):
        raise ValueError("protein_hydrogens에는 H 또는 D만 허용됩니다")
    protein_by_id = {_protein_id(record): record for record in protein}
    protein_h, uncertainties = _protein_hydrogen_parents(hydrogens, protein_by_id)
    ligand_donors, ligand_acceptors = _feature_indices(ligand_with_coords)
    interactions = []

    ligand_h = {}
    for atom in ligand_with_coords.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        ligand_h[atom.GetIdx()] = [
            n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1
        ]

    protein_roles = {}
    for record in protein:
        identifier = _protein_id(record)
        protein_roles[identifier] = _protein_role(record)
        if protein_roles[identifier][2]:
            uncertainties.append({
                "kind": "uncertain_protein_chemical_role", "protein_atom": identifier,
                "reasons": protein_roles[identifier][3],
            })

    def evaluate(donor_id, donor_xyz, hydrogen_xyzs, acceptor_id, acceptor_xyz,
                 role_uncertain=False):
        da = float(np.linalg.norm(donor_xyz - acceptor_xyz))
        if da > _HBOND_DA_MAX:
            return
        if not hydrogen_xyzs:
            interactions.append(_interaction(
                "possible_hbond_contact", [donor_id, acceptor_id], distance_DA_A=round(da, 4),
                requires_review=True,
                description="공여체/수용체 중원자 근접이지만 명시적인 D-H 방향 정보가 없습니다.",
            ))
            return
        accepted = False
        rejected = []
        for hxyz in hydrogen_xyzs:
            ha = float(np.linalg.norm(hxyz - acceptor_xyz))
            angle = _angle(donor_xyz, hxyz, acceptor_xyz)
            if angle is not None and ha <= _HBOND_HA_MAX and angle >= _HBOND_ANGLE_MIN:
                interactions.append(_interaction(
                    "directional_hbond", [donor_id, acceptor_id], distance_DA_A=round(da, 4),
                    distance_HA_A=round(ha, 4), angle_DHA_deg=round(angle, 3),
                    requires_review=bool(role_uncertain),
                    description="명시적 D-H 연결과 거리/각도 기준을 만족한 방향성 기하입니다.",
                ))
                accepted = True
                break
            rejected.append({"distance_HA_A": round(ha, 4),
                             "angle_DHA_deg": None if angle is None else round(angle, 3)})
        if not accepted:
            interactions.append(_interaction(
                "directional_hbond_rejected", [donor_id, acceptor_id],
                distance_DA_A=round(da, 4), tested_hydrogens=rejected,
                requires_review=True,
                description="공여체/수용체 근접은 있으나 명시적 D-H 방향 기준을 만족하지 않습니다.",
            ))

    for atom in ligand_with_coords.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        index = atom.GetIdx()
        lid = _ligand_id(atom)
        if index in ligand_donors:
            hxyz = [ligand_xyz[h] for h in ligand_h[index]]
            for record in protein:
                pid = _protein_id(record)
                if protein_roles[pid][1]:
                    evaluate(lid, ligand_xyz[index], hxyz, pid, record["xyz"], protein_roles[pid][2])
        if index in ligand_acceptors:
            for record in protein:
                pid = _protein_id(record)
                if protein_roles[pid][0]:
                    evaluate(pid, record["xyz"], [h["xyz"] for h in protein_h.get(pid, [])],
                             lid, ligand_xyz[index], protein_roles[pid][2])

    for atom in ligand_with_coords.GetAtoms():
        if atom.GetAtomicNum() <= 1 or atom.GetFormalCharge() == 0:
            continue
        lid = _ligand_id(atom)
        for record in protein:
            charge, source = _declared_charge(record)
            if charge and charge * atom.GetFormalCharge() < 0:
                distance = float(np.linalg.norm(ligand_xyz[atom.GetIdx()] - record["xyz"]))
                if distance <= 5.0:
                    interactions.append(_interaction(
                        "ionic_contact", [lid, _protein_id(record)], distance_A=round(distance, 4),
                        ligand_formal_charge=atom.GetFormalCharge(), protein_charge=charge,
                        protein_charge_source=source, requires_review=source != "atom_formal_charge",
                        description="실제 원자 형식전하 또는 명시된 잔기 상태에 기반한 반대 전하 근접입니다.",
                    ))

    for atom in ligand_with_coords.GetAtoms():
        if atom.GetAtomicNum() not in {6, 16}:
            continue
        for record in protein:
            if record["type_symbol"] not in {"C", "S"}:
                continue
            distance = float(np.linalg.norm(ligand_xyz[atom.GetIdx()] - record["xyz"]))
            if distance <= _HYDROPHOBIC_MAX:
                interactions.append(_interaction(
                    "hydrophobic_contact", [_ligand_id(atom), _protein_id(record)],
                    distance_A=round(distance, 4), requires_review=False,
                    description="탄소/황 원자의 기하학적 소수성 근접입니다.",
                ))

    periodic = Chem.GetPeriodicTable()
    for atom in ligand_with_coords.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        lr = float(periodic.GetRvdw(atom.GetSymbol()))
        for record in protein:
            pr = float(periodic.GetRvdw(record["type_symbol"]))
            distance = float(np.linalg.norm(ligand_xyz[atom.GetIdx()] - record["xyz"]))
            overlap = lr + pr - distance
            if overlap > _CLASH_OVERLAP_MIN:
                interactions.append(_interaction(
                    "vdw_clash", [_ligand_id(atom), _protein_id(record)],
                    distance_A=round(distance, 4), overlap_A=round(overlap, 4),
                    requires_review=True,
                    description="RDKit van der Waals 반지름 합에 대한 좌표 중첩입니다.",
                ))

    ligand_rings = []
    for ring in ligand_with_coords.GetRingInfo().AtomRings():
        if len(ring) >= 5 and all(ligand_with_coords.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            points = ligand_xyz[list(ring)]
            normal = _normal(points)
            if normal is not None:
                ligand_rings.append((ring, points.mean(axis=0), normal))
    grouped = {}
    for record in protein:
        residue = str(record["label_comp_id"]).upper()
        if str(record["label_atom_id"]).upper() in _PROTEIN_AROMATIC.get(residue, set()):
            key = (str(record["label_asym_id"]), str(record["auth_seq_id"]), residue)
            grouped.setdefault(key, []).append(record)
    for ring, center, normal in ligand_rings:
        for key, records in grouped.items():
            expected = _PROTEIN_AROMATIC[key[2]]
            names = {str(r["label_atom_id"]).upper() for r in records}
            if not expected <= names:
                continue
            points = np.asarray([r["xyz"] for r in records if str(r["label_atom_id"]).upper() in expected])
            pnormal = _normal(points)
            if pnormal is None:
                continue
            distance = float(np.linalg.norm(center - points.mean(axis=0)))
            if distance <= _AROMATIC_MAX:
                angle = float(np.degrees(np.arccos(np.clip(abs(np.dot(normal, pnormal)), 0.0, 1.0))))
                interactions.append(_interaction(
                    "aromatic_contact_geometry",
                    ["L:ring:" + ",".join(_ligand_id(ligand_with_coords.GetAtomWithIdx(i)) for i in ring),
                     ":".join(key)], centroid_distance_A=round(distance, 4),
                    normal_angle_deg=round(angle, 3), requires_review=True,
                    annotation="parallel_like" if angle <= 30 else "edge_or_oblique_geometry",
                    description="방향족 고리 중심/법선 기하이며 안정화 에너지나 결합 기여를 뜻하지 않습니다.",
                ))

    unique = {}
    for item in interactions:
        key = item["id"] + "|" + json.dumps(
            {k: v for k, v in item.items() if k not in {"description"}},
            sort_keys=True, ensure_ascii=False,
        )
        unique[key] = item
    interactions = sorted(unique.values(), key=lambda item: (item["id"], item["kind"]))
    return {
        "interactions": interactions,
        "uncertainties": uncertainties,
        "method": {
            "hydrogen_bond": {
                "D_A_max_A": _HBOND_DA_MAX, "H_A_max_A": _HBOND_HA_MAX,
                "D_H_A_min_deg": _HBOND_ANGLE_MIN,
                "policy": "명시적 D-H 연결과 좌표가 있을 때만 방향성 수소 결합으로 기록",
            },
            "protein_typing": "보수적인 표준 잔기 원자표; HIS/CYS/미지정 말단은 검토 대상",
            "ionic_policy": "원자 formal_charge 또는 명시적 residue_state만 사용",
            "clash": "RDKit van der Waals 반지름 합에서 거리를 뺀 중첩",
        },
        "limitations": [
            "좌표와 토폴로지에 대한 기술이며 결합 친화도, 효능 또는 실험적 상호작용을 입증하지 않습니다.",
            "명시적 수소가 없는 공여체/수용체 근접은 수소 결합으로 명명하지 않습니다.",
            "물, 유전환경, 용매화와 동역학은 평가하지 않습니다.",
        ],
    }


def before_after_report(before_profile, after_profile):
    """두 실제 interaction_profile 결과의 식별자 차이를 보고한다."""
    for name, profile in (("before_profile", before_profile), ("after_profile", after_profile)):
        if not isinstance(profile, dict) or not isinstance(profile.get("interactions"), list):
            raise ValueError(f"{name}은 interaction_profile 결과여야 합니다")
    before = {item["id"]: item for item in before_profile["interactions"]}
    after = {item["id"]: item for item in after_profile["interactions"]}
    lost_ids = sorted(set(before) - set(after))
    gained_ids = sorted(set(after) - set(before))
    retained_ids = sorted(set(before) & set(after))
    return {
        "lost_identifiers": lost_ids,
        "gained_identifiers": gained_ids,
        "retained_identifiers": retained_ids,
        "lost": [before[value] for value in lost_ids],
        "gained": [after[value] for value in gained_ids],
        "uncertainties": {
            "before": list(before_profile.get("uncertainties", [])),
            "after": list(after_profile.get("uncertainties", [])),
        },
        "quantitative_efficacy_score": None,
        "interpretation": "좌표 기반 상호작용 표기의 변화이며 효능 또는 결합 자유에너지 점수가 아닙니다.",
    }


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pdb_heavy(path):
    result = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = line[76:78].strip().title() if len(line) >= 78 else ""
            name = line[12:16].strip()
            if not element:
                element = re.sub(r"[^A-Za-z]", "", name)[:1].title()
            if element in {"H", "D"}:
                continue
            try:
                xyz = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except (ValueError, IndexError) as exc:
                raise ValueError(f"좌표 파일 {line_number}행을 해석할 수 없습니다") from exc
            if not np.isfinite(xyz).all():
                raise ValueError("좌표 파일에 비유한 좌표가 있습니다")
            key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
                   line[17:20].strip(), name)
            if key in result:
                raise ValueError(f"중복 중원자 식별자: {key}")
            result[key] = xyz
    if not result:
        raise ValueError("입력 파일에 중원자가 없습니다")
    return result


def _actual_pka_entries(output_dir, additional_paths=()):
    entries = []
    pattern = re.compile(
        r"^\s*(N\+|C-|[A-Za-z]{3})\s*(-?\d+)?\s+([A-Za-z0-9])\s+(-?\d+(?:\.\d+)?)\b",
        re.IGNORECASE,
    )
    paths = {
        path.resolve() for path in Path(output_dir).glob("*")
        if path.is_file() and path.suffix.lower() in {".pka", ".propka", ".log"}
    }
    paths.update(Path(path).resolve() for path in additional_paths if Path(path).is_file())
    for path in sorted(paths):
        in_summary = False
        for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            upper = line.upper()
            if "SUMMARY OF THIS PREDICTION" in upper or "SUMMARY OF PKA" in upper:
                in_summary = True
                continue
            if not in_summary:
                continue
            match = pattern.match(line)
            if match:
                value = float(match.group(4))
                if math.isfinite(value):
                    entries.append({
                        "residue": match.group(1).upper(), "sequence_number": int(match.group(2)) if match.group(2) else None,
                        "chain": match.group(3), "pKa": value,
                        "source_path": str(path.resolve()), "source_line": line_number,
                        "raw_line": line.rstrip(),
                    })
    return entries


def run_pdb2pqr(input_pdb, output_dir, executable, pH=7.4, timeout=300, preserve_heavy=True):
    """허용된 PDB2PQR 실행 파일로 PROPKA와 수소 최적화를 로컬 실행한다.

    기존 입력과 출력은 덮어쓰지 않는다. 실행 실패, 누락된 의존성 또는 출력
    누락 시 예외로 닫힌다. 중원자 이름 교환/좌표 변화는 자동 수용하지 않고
    ``heavy_change_review`` 상태로 반환한다.
    """
    context_ph = _finite(pH, "pH")
    timeout_value = _finite(timeout, "timeout")
    if not isinstance(preserve_heavy, bool):
        raise ValueError("preserve_heavy는 bool이어야 합니다")
    if timeout_value <= 0:
        raise ValueError("timeout은 양수여야 합니다")
    source = Path(input_pdb).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"입력 PDB를 찾을 수 없습니다: {source}")
    destination_dir = Path(output_dir).expanduser().resolve()
    destination_dir.mkdir(parents=True, exist_ok=True)
    if source.parent == destination_dir and source.suffix.lower() == ".pqr":
        raise ValueError("입력 파일과 출력 위치를 분리해야 합니다")

    executable_text = str(executable)
    candidate = Path(executable_text).expanduser()
    resolved_text = str(candidate.resolve()) if candidate.parent != Path(".") or candidate.is_absolute() else shutil.which(executable_text)
    if not resolved_text:
        raise FileNotFoundError("PDB2PQR 실행 파일을 찾을 수 없습니다")
    resolved_executable = Path(resolved_text).resolve()
    if resolved_executable.name.lower() not in _ALLOWED_PDB2PQR or not resolved_executable.is_file():
        raise ValueError("실행 파일은 로컬 pdb2pqr 또는 pdb2pqr.exe만 허용됩니다")

    output_pqr = destination_dir / f"{source.stem}.pqr"
    output_pdb = destination_dir / f"{source.stem}.prepared.pdb"
    stdout_path = destination_dir / "pdb2pqr.stdout.bin"
    stderr_path = destination_dir / "pdb2pqr.stderr.bin"
    existing_outputs = [path for path in (output_pqr, output_pdb, stdout_path, stderr_path) if path.exists()]
    if existing_outputs:
        raise FileExistsError(f"기존 출력을 덮어쓰지 않습니다: {existing_outputs[0]}")
    original_hash = _sha256(source)
    original_atoms = _pdb_heavy(source)
    existing_pka_outputs = [
        path for path in destination_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".pka", ".propka", ".log"}
    ]
    if existing_pka_outputs:
        raise FileExistsError("기존 pKa 출력을 덮어쓰거나 새 실행 결과와 혼합하지 않습니다")
    command = [
        str(resolved_executable), "--ff=PARSE", f"--with-ph={context_ph:.8g}",
        "--titration-state-method=propka", "--keep-chain", "--drop-water",
    ]
    if preserve_heavy:
        command.extend(("--noopt", "--nodebump"))
    command.extend((f"--pdb-output={output_pdb}", str(source), str(output_pqr)))
    try:
        completed = subprocess.run(
            command, cwd=str(destination_dir), shell=False, check=False,
            capture_output=True, text=False, timeout=timeout_value,
            creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError("PDB2PQR 실행에 실패하거나 제한 시간을 초과했습니다") from exc
    stdout_bytes = completed.stdout or b""
    stderr_bytes = completed.stderr or b""
    stdout_path.write_bytes(stdout_bytes)
    stderr_path.write_bytes(stderr_bytes)
    stdout_text = stdout_bytes.decode("utf-8", errors="replace")
    stderr_text = stderr_bytes.decode("utf-8", errors="replace")
    if (completed.returncode != 0 or not output_pqr.is_file()
            or output_pqr.stat().st_size == 0 or not output_pdb.is_file()
            or output_pdb.stat().st_size == 0):
        message = stderr_text.strip()[-1000:]
        raise RuntimeError(f"PDB2PQR이 유효한 출력을 만들지 못했습니다: {message}")
    if _sha256(source) != original_hash:
        raise RuntimeError("실행 중 원본 입력 파일이 변경되었습니다")

    output_atoms = _pdb_heavy(output_pdb)
    original_keys = set(original_atoms)
    output_keys = set(output_atoms)
    common = sorted(original_keys & output_keys)
    displacements = {
        "|".join(key): float(np.linalg.norm(output_atoms[key] - original_atoms[key]))
        for key in common
    }
    maximum = max(displacements.values(), default=None)
    missing = ["|".join(key) for key in sorted(original_keys - output_keys)]
    extra = ["|".join(key) for key in sorted(output_keys - original_keys)]
    exact = not missing and not extra and maximum == 0.0
    status = "complete" if exact else "heavy_change_review"
    receipt = {
        "status": status,
        "requires_review": not exact,
        "input_path": str(source),
        "output_pqr_path": str(output_pqr.resolve()),
        "output_pdb_path": str(output_pdb.resolve()),
        "input_sha256": original_hash,
        "output_sha256": _sha256(output_pqr),
        "output_pdb_sha256": _sha256(output_pdb),
        "stdout_path": str(stdout_path.resolve()),
        "stdout_sha256": _sha256(stdout_path),
        "stderr_path": str(stderr_path.resolve()),
        "stderr_sha256": _sha256(stderr_path),
        "pH_context": context_ph,
        "method": (
            "PDB2PQR with PARSE force field and PROPKA titration-state method; "
            "--noopt and --nodebump preserve input heavy-atom placement and leave added hydrogen orientations unoptimized"
            if preserve_heavy else
            "PDB2PQR with PARSE force field, PROPKA titration-state method, and default hydrogen optimization"
        ),
        "preserve_heavy_requested": preserve_heavy,
        "executable_path": str(resolved_executable),
        "executable_sha256": _sha256(resolved_executable),
        "command": command,
        "chain_policy": "--keep-chain; chain identifiers are included in heavy-atom identity comparison",
        "heavy_atom_mapping": {
            "input_count": len(original_atoms), "output_count": len(output_atoms),
            "mapped_count": len(common), "missing_from_output": missing,
            "extra_in_output": extra, "maximum_displacement_A": maximum,
            "coordinates_exactly_preserved": exact,
        },
        "pKa_entries_from_actual_output": _actual_pka_entries(
            destination_dir, (stdout_path, stderr_path)
        ),
        "stdout_tail": stdout_text[-2000:],
        "stderr_tail": stderr_text[-2000:],
        "limitations": [
            "PROPKA 출력에 실제로 존재하는 pKa 행만 노출하며 누락값을 추정하지 않습니다.",
            "말단 OXT 같은 추가 중원자, 잔기 원자명 교환 또는 중원자 좌표 변화는 heavy_change_review로 남깁니다.",
            "--noopt 사용 시 추가 수소 방향은 최적화되지 않았으며 방향성 접촉 해석에 검토가 필요합니다.",
            "준비된 구조는 효능, 결합 또는 실험 성공을 뜻하지 않습니다.",
        ],
    }
    json.dumps(receipt, ensure_ascii=False, allow_nan=False)
    return receipt


def read_prepared_protein(receipt):
    """Read the hash-verified fixed-width PDB emitted by :func:`run_pdb2pqr`.

    Hydrogen parent topology is reconstructed only when exactly one heavy atom
    in the same residue lies within 1.35 Å; ambiguous hydrogens remain explicit
    uncertainties rather than being assigned heuristically.
    """
    if not isinstance(receipt, dict):
        raise TypeError("receipt는 run_pdb2pqr 결과 사전이어야 합니다")
    path = Path(receipt.get("output_pdb_path", "")).expanduser().resolve()
    expected_hash = receipt.get("output_pdb_sha256")
    if not path.is_file() or not expected_hash or _sha256(path) != expected_hash:
        raise ValueError("준비된 PDB의 provenance hash 검증에 실패했습니다")
    atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                serial = int(line[6:11])
                xyz = np.asarray((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError as exc:
                raise ValueError(f"준비된 PDB {line_number}행을 해석할 수 없습니다") from exc
            atom_name = line[12:16].strip()
            element = line[76:78].strip().title() if len(line) >= 78 else ""
            if not element:
                element = re.sub(r"[^A-Za-z]", "", atom_name)[:1].title()
            atoms.append({
                "xyz": xyz, "type_symbol": element,
                "label_comp_id": line[17:20].strip(),
                "auth_seq_id": line[22:26].strip(),
                "label_atom_id": atom_name,
                "label_asym_id": line[21:22].strip(),
                "insertion_code": line[26:27].strip(),
                "original_atom_serial": serial,
                "topology_source": "PDB2PQR fixed-width PDB; D-H topology reconstructed geometrically",
            })
    heavy = [record for record in atoms if record["type_symbol"] not in {"H", "D"}]
    hydrogens = [record for record in atoms if record["type_symbol"] in {"H", "D"}]
    if not heavy:
        raise ValueError("준비된 PDB에 중원자가 없습니다")
    uncertainties = []
    attached = {_protein_id(record): [] for record in heavy}
    for hydrogen in hydrogens:
        same_residue = [
            record for record in heavy
            if (record["label_asym_id"], record["auth_seq_id"], record["insertion_code"], record["label_comp_id"])
            == (hydrogen["label_asym_id"], hydrogen["auth_seq_id"], hydrogen["insertion_code"], hydrogen["label_comp_id"])
            and float(np.linalg.norm(record["xyz"] - hydrogen["xyz"])) <= 1.35
        ]
        if len(same_residue) == 1:
            parent_id = _protein_id(same_residue[0])
            hydrogen["parent_atom_id"] = parent_id
            hydrogen["parent_topology_reconstructed"] = True
            attached[parent_id].append(hydrogen["label_atom_id"])
        else:
            uncertainties.append({
                "kind": "ambiguous_reconstructed_DH_topology",
                "hydrogen": _protein_id(hydrogen),
                "candidate_parent_count": len(same_residue),
            })
    for record in heavy:
        record["observed_attached_hydrogens"] = attached[_protein_id(record)]
    by_residue = {}
    for record in heavy:
        key = (record["label_asym_id"], record["auth_seq_id"], record["insertion_code"])
        by_residue.setdefault(key, {})[record["label_atom_id"].upper()] = record
    for residue in by_residue.values():
        if "ND1" in residue and "NE2" in residue:
            residue["ND1"]["histidine_other_nitrogen_has_hydrogen"] = bool(
                residue["NE2"]["observed_attached_hydrogens"]
            )
            residue["NE2"]["histidine_other_nitrogen_has_hydrogen"] = bool(
                residue["ND1"]["observed_attached_hydrogens"]
            )
    return {
        "protein_atoms": heavy,
        "protein_hydrogens": hydrogens,
        "uncertainties": uncertainties,
        "source_path": str(path),
        "source_sha256": expected_hash,
        "topology_policy": "unique same-residue heavy atom within 1.35 A only",
    }

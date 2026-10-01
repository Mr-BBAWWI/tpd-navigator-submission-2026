from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from typing import Any, Callable, Iterable

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from packages.science.mapped_stereo import (
    mapped_tetrahedral_parity,
    same_mapped_tetrahedral_stereo,
)


RULE_CATALOG = [
    {
        "rule_id": "LH_OH",
        "transformation_class": "linker_handle_introduction",
        "tier": 1,
        "name": "hydroxyl_handle",
        "design_hypothesis": "검토 대상 C–H 또는 중성 N–H 위치에 수산기 핸들을 도입해 후속 유도체화와 수소결합 탐색을 가능하게 한다.",
        "scope": "C-H or neutral nonaromatic N-H; one hydrogen is consumed.",
    },
    {
        "rule_id": "LH_NH2",
        "transformation_class": "linker_handle_introduction",
        "tier": 1,
        "name": "amino_handle",
        "design_hypothesis": "검토 대상 C–H 또는 중성 N–H 위치에 아미노 핸들을 도입해 극성 상호작용과 결합 확장 방향을 탐색한다.",
        "scope": "C-H or neutral nonaromatic N-H; positively charged nitrogen is excluded.",
    },
    {
        "rule_id": "LH_AMINOMETHYL",
        "transformation_class": "linker_handle_introduction",
        "tier": 1,
        "name": "aminomethyl_handle",
        "design_hypothesis": "아미노메틸기를 도입해 염기성 중심을 한 결합만큼 외부로 이동시키고 합성 검토가 필요한 출구 벡터를 탐색한다.",
        "scope": "C-H or neutral nonaromatic N-H.",
    },
    {
        "rule_id": "LH_HYDROXYETHYL",
        "transformation_class": "linker_handle_introduction",
        "tier": 2,
        "name": "hydroxyethyl_handle",
        "design_hypothesis": "하이드록시에틸기를 도입해 유연한 친수성 출구 벡터와 말단 수소결합 공여체를 탐색한다.",
        "scope": "C-H or neutral nonaromatic N-H.",
    },
    {
        "rule_id": "LH_ALKYNE",
        "transformation_class": "linker_handle_introduction",
        "tier": 2,
        "name": "terminal_alkyne_handle",
        "design_hypothesis": "말단 알카인을 도입해 선형 출구 벡터와 클릭 화학 기반의 후속 결합 가능성을 평가한다.",
        "scope": "C-H or neutral nonaromatic N-H.",
    },
    {
        "rule_id": "LH_AZIDE",
        "transformation_class": "linker_handle_introduction",
        "tier": 2,
        "name": "azide_handle",
        "design_hypothesis": "아지드 핸들을 도입해 클릭 화학용 반응점을 제공하되 새 원자에 필요한 형식전하는 명시적으로 유지한다.",
        "scope": "C-H or neutral nonaromatic N-H.",
    },
    {
        "rule_id": "HS_AROM_CH_TO_N",
        "transformation_class": "heteroatom_swap",
        "tier": 2,
        "name": "aromatic_CH_to_N",
        "design_hypothesis": "방향족 C–H를 피리딘형 질소로 치환해 골격 형상을 유지하면서 수소결합 수용성과 전자분포를 조절한다.",
        "scope": "Neutral aromatic carbon bearing exactly one hydrogen.",
    },
    {
        "rule_id": "HS_RING_CH2_TO_O",
        "transformation_class": "heteroatom_swap",
        "tier": 2,
        "name": "saturated_ring_CH2_to_O",
        "design_hypothesis": "포화 고리의 CH2를 산소로 치환해 고리 형상을 크게 바꾸지 않고 극성과 수용체 특성을 높인다.",
        "scope": "Degree-two neutral aliphatic ring CH2 with two single bonds.",
    },
    {
        "rule_id": "HS_RING_CH2_TO_NH",
        "transformation_class": "heteroatom_swap",
        "tier": 2,
        "name": "saturated_ring_CH2_to_NH",
        "design_hypothesis": "포화 고리의 CH2를 NH로 치환해 고리 내 염기성 및 수소결합 공여 가능성을 탐색한다.",
        "scope": "Degree-two neutral aliphatic ring CH2 with two single bonds.",
    },
    {
        "rule_id": "RE_INSERT_CH2",
        "transformation_class": "ring_expansion",
        "tier": 2,
        "name": "aliphatic_ring_CH2_insertion",
        "design_hypothesis": "비방향족 고리 결합에 CH2를 삽입해 고리 크기와 결합 벡터를 한 단계 확장한다.",
        "scope": "Nonaromatic single bond in an aliphatic 3-7 membered ring.",
    },
    {
        "rule_id": "RC_REMOVE_CH2",
        "transformation_class": "ring_contraction",
        "tier": 2,
        "name": "aliphatic_ring_CH2_deletion",
        "design_hypothesis": "6원 이상의 비방향족 고리에서 degree-2 CH2를 제거해 더 작은 고리의 형태 제약을 탐색한다.",
        "scope": "Degree-two aliphatic CH2 in a ring of size at least six.",
    },
    {
        "rule_id": "CR_CLOSE_CC",
        "transformation_class": "conformational_restriction",
        "tier": 2,
        "name": "intramolecular_CC_closure",
        "design_hypothesis": "두 검토 대상 지방족 C–H 사이에 C–C 결합을 형성해 4–7원 고리를 만들고 유연성을 제한한다.",
        "scope": "Two nonbonded aliphatic carbon atoms separated by 3-6 bonds.",
    },
    {
        "rule_id": "BI_CARBONYL_O_TO_S",
        "transformation_class": "bioisosteric_replacement",
        "tier": 2,
        "name": "carbonyl_O_to_S",
        "design_hypothesis": "카보닐 산소를 황으로 치환해 유사한 결합 방향성을 유지하면서 분극성과 크기를 조절한다.",
        "scope": "Terminal neutral carbonyl oxygen; both oxygen and carbonyl carbon must be allowed.",
    },
    {
        "rule_id": "EV_RELOCATE_OH",
        "transformation_class": "exit_vector_relocation",
        "tier": 2,
        "name": "same_ring_OH_relocation",
        "design_hypothesis": "동일 방향족 고리 내에서 말단 OH의 위치를 이동해 코어를 유지한 채 출구 벡터 방향을 비교한다.",
        "scope": "Terminal OH moved to an aromatic C-H on the same ring.",
    },
    {
        "rule_id": "EV_RELOCATE_NH2",
        "transformation_class": "exit_vector_relocation",
        "tier": 2,
        "name": "same_ring_NH2_relocation",
        "design_hypothesis": "동일 방향족 고리 내에서 말단 NH2의 위치를 이동해 위치 이성질체의 상호작용 방향을 평가한다.",
        "scope": "Terminal neutral NH2 moved to an aromatic C-H on the same ring.",
    },
    {
        "rule_id": "HB_OH_TO_NH2",
        "transformation_class": "hbond_rewiring",
        "tier": 3,
        "name": "terminal_OH_to_NH2",
        "design_hypothesis": "탐색적으로 말단 OH를 NH2로 바꾸어 수소결합 수용체·공여체 패턴의 재배선을 평가한다.",
        "scope": "Exploratory only; terminal neutral OH.",
    },
    {
        "rule_id": "HB_NH2_TO_OH",
        "transformation_class": "hbond_rewiring",
        "tier": 3,
        "name": "terminal_NH2_to_OH",
        "design_hypothesis": "탐색적으로 말단 NH2를 OH로 바꾸어 염기성과 수소결합 패턴의 변화를 평가한다.",
        "scope": "Exploratory only; terminal neutral NH2.",
    },
]

RULE_CATALOG.extend([
    {'rule_id':'LH_CARBOXYL','transformation_class':'linker_handle_introduction','tier':1,
     'name':'carboxyl_handle','design_hypothesis':'검토 대상 C–H 위치에 카복실기를 도입해 후속 amide 결합 가능성을 탐색한다. 반응 경로와 활성은 미검증이다.',
     'scope':'C-H substitution with carboxyl group; synthesis requires review.'},
    {'rule_id':'LH_PROPARGYL','transformation_class':'linker_handle_introduction','tier':2,
     'name':'propargyl_handle','design_hypothesis':'메틸렌을 사이에 둔 말단 알카인으로 후속 결합 방향과 강성을 탐색한다.',
     'scope':'C-H or neutral N-H; a propargyl group, not a direct N-alkyne.'}
])
_N_SUBSTITUENT_FRAGMENTS = [
    ("METHYL", "*C", "Small alkyl substitution probes minimal steric growth.", "N-alkylation feasibility and retained protonation require review."),
    ("ETHYL", "*CC", "Compact alkyl growth probes a short hydrophobic vector.", "N-alkylation feasibility and metabolic liability require review."),
    ("NPROPYL", "*CCC", "Linear propyl growth samples a longer hydrophobic vector.", "Added lipophilicity and metabolic liability require review."),
    ("ISOPROPYL", "*C(C)C", "Branched alkyl growth samples steric occupancy near the attachment atom.", "Steric compatibility and N-alkylation selectivity require review."),
    ("HYDROXYMETHYL", "*CO", "A hydroxymethyl arm adds a compact terminal hydrogen-bond donor and acceptor.", "Chemical stability and synthesis route require review."),
    ("HYDROXYPROPYL", "*CCCO", "A hydroxypropyl arm extends a terminal alcohol farther from the attachment atom.", "Conformational entropy and permeability effects require review."),
    ("METHOXYMETHYL", "*COC", "A methoxymethyl arm probes a compact ether acceptor.", "Acetal-like metabolic stability and synthesis require review."),
    ("METHOXYETHYL", "*CCOC", "A methoxyethyl arm provides a flexible neutral ether vector.", "Flexibility and oxidative metabolism require review."),
    ("ETHOXYETHYL", "*CCOCC", "An ethoxyethyl arm samples a longer neutral polyether-like vector.", "Flexibility, permeability, and oxidative metabolism require review."),
    ("ACETAMIDE", "*CC(=O)N", "A methylene amide introduces donor/acceptor functionality beyond the anchor.", "Amide geometry, hydrolysis, and synthesis require review."),
    ("N_METHYL_ACETAMIDE", "*CC(=O)NC", "An N-methyl amide probes a less donating polar terminus.", "Amide geometry and synthetic accessibility require review."),
    ("HYDROXYETHYL_AMIDE", "*CC(=O)NCCO", "A hydroxyethyl amide combines a directional amide with a terminal alcohol.", "High polarity, flexibility, and synthesis require review."),
    ("METHYL_ESTER", "*CC(=O)OC", "A methyl ester provides a compact carbonyl-containing vector.", "Ester hydrolysis and prodrug-like behavior require review."),
    ("CARBOXYMETHYL", "*CC(=O)O", "A carboxymethyl arm provides a potential amide-coupling handle.", "Ionization, permeability, and synthetic route require review."),
    ("CARBOXYETHYL", "*CCC(=O)O", "A carboxyethyl arm moves a coupling handle farther from the anchor.", "Ionization, flexibility, and synthesis require review."),
    ("CARBAMATE", "*CCOC(=O)N", "A carbamate-bearing arm probes a neutral directional polar motif.", "Hydrolysis and conformational preferences require review."),
    ("DIMETHYLAMINOETHYL", "*CCN(C)C", "A dimethylaminoethyl arm adds a distal basic center without an N-N bond.", "Protonation, permeability, and over-basicity require review."),
    ("AMINOETHYL", "*CCN", "An aminoethyl arm adds a distal primary amine for further derivatization.", "Multiple protonation states and chemoselectivity require review."),
    ("BENZYL", "*Cc1ccccc1", "A benzyl arm samples a rigid aromatic hydrophobic vector.", "Lipophilicity and oxidative metabolism require review."),
    ("PYRIDYLMETHYL", "*Cc1ccncc1", "A pyridylmethyl arm adds an aromatic acceptor and directional vector.", "Regioisomer choice, basicity, and synthesis require review."),
    ("PYRIMIDINYLMETHYL", "*Cc1nccnc1", "A pyrimidinylmethyl arm probes a compact electron-poor heteroaryl vector.", "Heteroaryl substitution chemistry and basicity require review."),
    ("FURANYLMETHYL", "*Cc1ccoc1", "A furanylmethyl arm adds a compact heteroaryl ether-like acceptor.", "Oxidative and reactive-metabolite liability require review."),
    ("CYCLOPROPYL", "*C1CC1", "Direct cyclopropyl substitution probes compact conformational restriction.", "Ring strain and N-alkylation feasibility require review."),
    ("CYCLOPROPYLMETHYL", "*CC1CC1", "A cyclopropylmethyl arm combines compact rigidity with one rotatable bond.", "Ring-opening metabolism and synthesis require review."),
    ("CYCLOBUTYL", "*C1CCC1", "Direct cyclobutyl substitution samples a compact saturated ring vector.", "Ring strain and steric compatibility require review."),
    ("CYCLOBUTYLMETHYL", "*CC1CCC1", "A cyclobutylmethyl arm moves a constrained hydrophobe away from the anchor.", "Lipophilicity and synthetic accessibility require review."),
    ("CYCLOPENTYL", "*C1CCCC1", "Direct cyclopentyl substitution samples a medium saturated ring vector.", "Steric demand and lipophilicity require review."),
    ("CYCLOHEXYLMETHYL", "*CC1CCCCC1", "A cyclohexylmethyl arm probes a larger constrained hydrophobic vector.", "Molecular-weight and lipophilicity growth require review."),
    ("TETRAHYDROFURANYLMETHYL", "*CC1CCOC1", "A tetrahydrofuranylmethyl arm combines ring constraint with an ether acceptor.", "Stereochemistry, ether metabolism, and synthesis require review."),
    ("MORPHOLINYLMETHYL", "*CC1COCCN1", "A morpholinylmethyl arm adds a constrained polar heterocycle through carbon.", "Protonation, regioselectivity, and synthesis require review."),
    ("PIPERIDINYLMETHYL", "*CC1CCNCC1", "A piperidinylmethyl arm adds a constrained distal basic center through carbon.", "Multiple protonation states and lipophilicity require review."),
    ("CYCLOHEXENYLMETHYL", "*CC1=CCCCC1", "A cyclohexenylmethyl arm samples a rigid unsaturated exit vector.", "Alkene metabolism and isomeric composition require review."),
    ("ALLYL", "*CC=C", "An allyl arm provides a compact unsaturated vector for later chemistry review.", "Alkene reactivity and metabolic oxidation require review."),
    ("DIFLUOROETHYL", "*CC(F)F", "A difluoroethyl arm probes compact polarity and conformational effects.", "Fluorinated building-block availability and lipophilicity require review."),
    ("METHYLSULFONYLETHYL", "*CCS(=O)(=O)C", "A methylsulfonylethyl arm adds a strong neutral polar motif.", "High polarity and synthetic compatibility require review."),
]
for _fragment_id, _fragment_smiles, _benefit, _risk in _N_SUBSTITUENT_FRAGMENTS:
    RULE_CATALOG.append({
        "rule_id": "LH_N_" + _fragment_id,
        "transformation_class": "linker_handle_introduction",
        "tier": 2,
        "name": "neutral_N_substitution_" + _fragment_id.lower(),
        "design_hypothesis": _benefit + " No activity or binding improvement is inferred.",
        "scope": "Cited linker-handle scope at a neutral, nonaromatic, nonamide N-H anchor only.",
        "expected_benefit": _benefit,
        "synthetic_risks": [_risk, "The enumerated graph is a design proposal, not a validated reaction product."],
        "fragment_smiles": _fragment_smiles,
    })
RULE_CATALOG.extend([
    {
        "rule_id": "BI_CARBONYL_TO_SULFONYL",
        "transformation_class": "bioisosteric_replacement",
        "tier": 2,
        "name": "carbonyl_to_sulfonyl_or_sulfonamide",
        "design_hypothesis": "An eligible non-ring carbonyl center is expanded to a sulfonyl center to probe a larger, more polar carbonyl bioisostere.",
        "scope": "Non-ring C(=O) with two single-bond substituents; every affected parent atom must be permitted.",
        "expected_benefit": "Samples a tetrahedral sulfone or sulfonamide geometry in place of a planar carbonyl.",
        "synthetic_risks": ["This is a scaffold replacement rather than a direct reaction prediction.", "Oxidation state, polarity, and synthesis require review."],
    },
    {
        "rule_id": "HB_ETHER_TO_NH",
        "transformation_class": "hbond_rewiring",
        "tier": 2,
        "name": "acyclic_ether_O_to_NH",
        "design_hypothesis": "An eligible acyclic ether oxygen is replaced by NH to rewire acceptor/donor and basicity patterns.",
        "scope": "Neutral acyclic degree-two ether with two single bonds; every affected atom must be permitted.",
        "expected_benefit": "Explores a secondary-amine donor/basic center at an existing ether position.",
        "synthetic_risks": ["Basicity and protonation can change substantially.", "A route-specific scaffold synthesis is required."],
    },
])

_RULE_BY_ID = {rule["rule_id"]: rule for rule in RULE_CATALOG}
for _rule in RULE_CATALOG:
    if _rule['rule_id'] in {'LH_OH','LH_NH2','LH_ALKYNE','LH_AZIDE'}:
        _rule['scope'] = 'Neutral C-H only; N anchors are deliberately excluded.'
        _rule['design_hypothesis'] = _rule['design_hypothesis'].replace('또는 중성 N–H ', '')

_EXTENSION_RULE_IDS = {
    "BI_CARBONYL_O_TO_S",
    "BI_CARBONYL_TO_SULFONYL",
    "EV_RELOCATE_OH",
    "EV_RELOCATE_NH2",
    "HB_OH_TO_NH2",
    "HB_NH2_TO_OH",
    "HB_ETHER_TO_NH",
}
for _rule in RULE_CATALOG:
    _rule.setdefault("rationale", _rule["design_hypothesis"])
    _rule.setdefault("applicability", _rule.get("scope"))
    _rule.setdefault("review_required", True)
    _rule.setdefault(
        "api_note",
        "API-authored 2D graph enumeration only; this record is not expert approval, a reaction prediction, or evidence of activity.",
    )
    _rule["requires_expert_allowed_rule"] = _rule["rule_id"] in _EXTENSION_RULE_IDS


def _safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(k): _safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_safe(v) for v in value]
    return str(value)


def _map_index(mol: Chem.Mol) -> dict[int, int]:
    result: dict[int, int] = {}
    for atom in mol.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if atom_map <= 0:
            raise ValueError("Every parent atom must have a positive atom-map number")
        if atom_map in result:
            raise ValueError(f"Duplicate parent atom-map number: {atom_map}")
        result[atom_map] = atom.GetIdx()
    return result


def _without_maps_smiles(mol: Chem.Mol, *, isomeric: bool = True) -> str:
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=isomeric)


def _site_table(sites: Iterable[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    table: dict[int, dict[str, Any]] = {}
    for site in sites:
        atom_map = int(site["atom_map"])
        if atom_map in table:
            raise ValueError(f"Duplicate site declaration for atom map {atom_map}")
        state = str(site.get("state", "UNKNOWN")).upper()
        if state not in {"PROTECTED", "MODIFIABLE", "UNKNOWN"}:
            state = "UNKNOWN"
        table[atom_map] = {
            "atom_map": atom_map,
            "state": state,
            "evidence": _safe(site.get("evidence")),
        }
    return table


def _status(atom_map: int, table: dict[int, dict[str, Any]]) -> dict[str, Any]:
    return dict(table.get(atom_map, {
        "atom_map": atom_map,
        "state": "UNKNOWN",
        "evidence": None,
    }))


def _consume_h(atom: Chem.Atom) -> None:
    if atom.GetTotalNumHs(includeNeighbors=True) < 1:
        raise ValueError("NO_REMOVABLE_H")
    atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
    atom.SetNoImplicit(True)


def _add_h(atom: Chem.Atom) -> None:
    atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
    atom.SetNoImplicit(True)


def _new_atom(atomic_num: int, atom_map: int, hydrogens: int = 0,
              charge: int = 0) -> Chem.Atom:
    atom = Chem.Atom(atomic_num)
    atom.SetAtomMapNum(atom_map)
    atom.SetFormalCharge(charge)
    atom.SetNumExplicitHs(hydrogens)
    atom.SetNoImplicit(True)
    return atom


def _add_fragment(rw: Chem.RWMol, anchor: int, kind: str,
                  first_map: int) -> list[int]:
    if kind == 'CARBOXYL':
        c=rw.AddAtom(_new_atom(6,first_map));o=rw.AddAtom(_new_atom(8,first_map+1));oh=rw.AddAtom(_new_atom(8,first_map+2,1))
        rw.AddBond(anchor,c,Chem.BondType.SINGLE);rw.AddBond(c,o,Chem.BondType.DOUBLE);rw.AddBond(c,oh,Chem.BondType.SINGLE)
        return [first_map,first_map+1,first_map+2]
    specs = {
        "OH": [(8, 1, 0)],
        "NH2": [(7, 2, 0)],
        "AMINOMETHYL": [(6, 2, 0), (7, 2, 0)],
        "HYDROXYETHYL": [(6, 2, 0), (6, 2, 0), (8, 1, 0)],
        "ALKYNE": [(6, 0, 0), (6, 1, 0)],
        "AZIDE": [(7, 0, 0), (7, 0, 1), (7, 0, -1)],
        "PROPARGYL": [(6, 2, 0), (6, 0, 0), (6, 1, 0)],
    }[kind]
    bond_types = {
        "OH": [Chem.BondType.SINGLE],
        "NH2": [Chem.BondType.SINGLE],
        "AMINOMETHYL": [Chem.BondType.SINGLE, Chem.BondType.SINGLE],
        "HYDROXYETHYL": [Chem.BondType.SINGLE] * 3,
        "ALKYNE": [Chem.BondType.SINGLE, Chem.BondType.TRIPLE],
        "AZIDE": [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.DOUBLE],
        "PROPARGYL": [Chem.BondType.SINGLE, Chem.BondType.SINGLE, Chem.BondType.TRIPLE],
    }[kind]
    added: list[int] = []
    previous = anchor
    for offset, ((atomic_num, hs, charge), bond_type) in enumerate(zip(specs, bond_types)):
        atom_map = first_map + offset
        idx = rw.AddAtom(_new_atom(atomic_num, atom_map, hs, charge))
        rw.AddBond(previous, idx, bond_type)
        previous = idx
        added.append(atom_map)
    return added


def _add_smiles_fragment(rw: Chem.RWMol, anchor: int, fragment_smiles: str,
                         first_map: int) -> list[int]:
    fragment = Chem.MolFromSmiles(fragment_smiles)
    if fragment is None:
        raise ValueError("INVALID_FRAGMENT_TEMPLATE")
    dummies = [atom for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(dummies) != 1 or dummies[0].GetDegree() != 1:
        raise ValueError("FRAGMENT_REQUIRES_ONE_TERMINAL_DUMMY")
    dummy = dummies[0]
    attachment = dummy.GetNeighbors()[0]
    attachment_bond = fragment.GetBondBetweenAtoms(dummy.GetIdx(), attachment.GetIdx())
    index_map: dict[int, int] = {}
    added: list[int] = []
    for atom in fragment.GetAtoms():
        if atom.GetIdx() == dummy.GetIdx():
            continue
        copy = Chem.Atom(atom)
        atom_map = first_map + len(added)
        copy.SetAtomMapNum(atom_map)
        index_map[atom.GetIdx()] = rw.AddAtom(copy)
        added.append(atom_map)
    for bond in fragment.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if dummy.GetIdx() in {begin, end}:
            continue
        rw.AddBond(index_map[begin], index_map[end], bond.GetBondType())
    rw.AddBond(anchor, index_map[attachment.GetIdx()], attachment_bond.GetBondType())
    return added


def _is_nonamide_neutral_nh(atom: Chem.Atom) -> bool:
    if not (atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 0
            and not atom.GetIsAromatic() and atom.GetDegree() <= 2
            and atom.GetTotalNumHs(includeNeighbors=True) >= 1):
        return False
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and any(
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(neighbor).GetAtomicNum() in {8, 16}
            for bond in neighbor.GetBonds()
        ):
            return False
        if neighbor.GetAtomicNum() == 16 and sum(
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(neighbor).GetAtomicNum() == 8
            for bond in neighbor.GetBonds()
        ) >= 1:
            return False
    return True


def _ring_sizes_by_atom(mol: Chem.Mol) -> dict[int, set[int]]:
    result: dict[int, set[int]] = defaultdict(set)
    for ring in mol.GetRingInfo().AtomRings():
        for idx in ring:
            result[idx].add(len(ring))
    return result


def _ring_sizes_by_bond(mol: Chem.Mol) -> dict[int, set[int]]:
    result: dict[int, set[int]] = defaultdict(set)
    for ring in mol.GetRingInfo().BondRings():
        for idx in ring:
            result[idx].add(len(ring))
    return result


def _stereo_signature(mol: Chem.Mol) -> tuple[dict[int, tuple[Any, ...]], dict[tuple[int, int], int]]:
    atoms = {
        atom.GetAtomMapNum(): mapped_tetrahedral_parity(atom)
        for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0
    }
    bonds: dict[tuple[int, int], int] = {}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a > 0 and b > 0:
            bonds[tuple(sorted((a, b)))] = int(bond.GetStereo())
    return atoms, bonds


def _graph_preserved(parent: Chem.Mol, product: Chem.Mol,
                     parent_maps: set[int]) -> bool:
    pidx = _map_index(parent)
    qidx = _map_index(product)
    if not parent_maps.issubset(qidx):
        return False
    for atom_map in parent_maps:
        pa, qa = parent.GetAtomWithIdx(pidx[atom_map]), product.GetAtomWithIdx(qidx[atom_map])
        signature = lambda a: (a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(), a.GetIsAromatic(), a.GetTotalNumHs(),
            mapped_tetrahedral_parity(a), sorted((n.GetAtomMapNum(), n.GetAtomicNum()) for n in a.GetNeighbors()))
        if signature(pa) != signature(qa):
            return False
    for bond in parent.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a not in parent_maps and b not in parent_maps:
            continue
        if a not in qidx or b not in qidx:
            return False
        qb = product.GetBondBetweenAtoms(qidx[a], qidx[b])
        if (qb is None or qb.GetBondType() != bond.GetBondType()
                or int(qb.GetStereo()) != int(bond.GetStereo())):
            return False
    return True


def generate(parent: Chem.Mol, sites: list[dict[str, Any]],
             exploratory: bool = False) -> dict[str, Any]:
    if parent is None:
        raise ValueError("parent must be an RDKit Mol")
    base = Chem.Mol(parent)
    Chem.SanitizeMol(base)
    if any(a.GetAtomicNum() <= 1 for a in base.GetAtoms()):
        raise ValueError('HEAVY_ATOM_PARENT_REQUIRED')
    if len(Chem.GetMolFrags(base)) != 1:
        raise ValueError("parent must contain exactly one connected component")
    pidx = _map_index(base)
    parent_maps = set(pidx)
    table = _site_table(sites)
    if set(table) != parent_maps:
        raise ValueError('SITE_MAP_COVERAGE_MISMATCH')
    parent_canonical = _without_maps_smiles(base)
    parent_constitutional = _without_maps_smiles(base, isomeric=False)
    parent_id = "sha256:" + hashlib.sha256(parent_canonical.encode("utf-8")).hexdigest()
    original_site_map = {
        str(atom_map): table[atom_map]["state"] for atom_map in sorted(parent_maps)
    }
    parent_stereo_atoms, parent_stereo_bonds = _stereo_signature(base)
    initial_map = max(5000, max(parent_maps) + 1)
    ring_atoms = _ring_sizes_by_atom(base)
    ring_bonds = _ring_sizes_by_bond(base)
    candidates: dict[str, list[tuple[set[int], Callable[[Chem.RWMol], list[int]], str, set[int]]]] = defaultdict(list)

    for rule_id, kind in [
        ("LH_OH", "OH"), ("LH_NH2", "NH2"),
        ("LH_AMINOMETHYL", "AMINOMETHYL"),
        ("LH_HYDROXYETHYL", "HYDROXYETHYL"),
        ("LH_ALKYNE", "ALKYNE"), ("LH_AZIDE", "AZIDE"),
        ("LH_CARBOXYL", "CARBOXYL"), ("LH_PROPARGYL", "PROPARGYL"),
    ]:
        for atom_map in sorted(parent_maps):
            atom = base.GetAtomWithIdx(pidx[atom_map])
            eligible_c = atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0
            eligible_n = (atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 0
                          and not atom.GetIsAromatic() and atom.GetDegree() <= 2)
            if eligible_n and kind not in {'AMINOMETHYL', 'HYDROXYETHYL', 'PROPARGYL'}:
                continue  # Avoid treating N-O/N-N/ynamine/azidoamine as conservative handles.
            if not (eligible_c or eligible_n) or atom.GetTotalNumHs(includeNeighbors=True) < 1:
                continue

            def edit(rw: Chem.RWMol, m=atom_map, k=kind) -> list[int]:
                idx = next(a.GetIdx() for a in rw.GetAtoms() if a.GetAtomMapNum() == m)
                _consume_h(rw.GetAtomWithIdx(idx))
                return _add_fragment(rw, idx, k, initial_map)

            candidates[rule_id].append(({atom_map}, edit, f"anchor_map={atom_map}", set()))

    for fragment_id, fragment_smiles, _, _ in _N_SUBSTITUENT_FRAGMENTS:
        rule_id = "LH_N_" + fragment_id
        for atom_map in sorted(parent_maps):
            atom = base.GetAtomWithIdx(pidx[atom_map])
            if not _is_nonamide_neutral_nh(atom):
                continue

            def edit(rw: Chem.RWMol, m=atom_map, smi=fragment_smiles) -> list[int]:
                idx = next(a.GetIdx() for a in rw.GetAtoms() if a.GetAtomMapNum() == m)
                _consume_h(rw.GetAtomWithIdx(idx))
                return _add_smiles_fragment(rw, idx, smi, initial_map)

            candidates[rule_id].append((
                {atom_map}, edit, f"N_anchor_map={atom_map};fragment={fragment_id}", set()
            ))

    for atom_map in sorted(parent_maps):
        atom = base.GetAtomWithIdx(pidx[atom_map])
        if (atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0
                and not atom.GetIsAromatic() and not atom.IsInRing()
                and atom.GetDegree() == 2
                and all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds())):
            neighbors = {n.GetAtomMapNum() for n in atom.GetNeighbors()}

            def edit(rw: Chem.RWMol, m=atom_map) -> list[int]:
                a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                a.SetAtomicNum(7)
                a.SetNumExplicitHs(1)
                a.SetNoImplicit(True)
                return []

            candidates["HB_ETHER_TO_NH"].append((
                {atom_map} | neighbors, edit, f"ether_atom_map={atom_map}", set()
            ))

        if (atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0
                and not atom.GetIsAromatic() and not atom.IsInRing()):
            double_oxygens = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() == 8 and n.GetFormalCharge() == 0
                and base.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
            ]
            single_neighbors = [
                n for n in atom.GetNeighbors()
                if base.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE
            ]
            if len(double_oxygens) == 1 and len(single_neighbors) == 2:
                oxygen_map = double_oxygens[0].GetAtomMapNum()
                neighbor_maps = {n.GetAtomMapNum() for n in single_neighbors}

                def edit(rw: Chem.RWMol, m=atom_map) -> list[int]:
                    center = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                    center.SetAtomicNum(16)
                    center.SetNumExplicitHs(0)
                    center.SetNoImplicit(True)
                    oxygen = rw.AddAtom(_new_atom(8, initial_map, 0, 0))
                    rw.AddBond(center.GetIdx(), oxygen, Chem.BondType.DOUBLE)
                    return [initial_map]

                candidates["BI_CARBONYL_TO_SULFONYL"].append((
                    {atom_map, oxygen_map} | neighbor_maps,
                    edit,
                    f"carbonyl_center_map={atom_map}",
                    set(),
                ))

        if (atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0
                and atom.GetIsAromatic()
                and atom.GetTotalNumHs(includeNeighbors=True) == 1):
            def edit(rw: Chem.RWMol, m=atom_map) -> list[int]:
                a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                a.SetAtomicNum(7)
                a.SetNumExplicitHs(0)
                a.SetNoImplicit(True)
                return []
            candidates["HS_AROM_CH_TO_N"].append(({atom_map}, edit, f"atom_map={atom_map}", set()))

        single_neighbors = all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds())
        ring_ch2 = (atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0
                    and not atom.GetIsAromatic() and atom.IsInRing()
                    and atom.GetDegree() == 2 and single_neighbors
                    and atom.GetTotalNumHs(includeNeighbors=True) == 2)
        if ring_ch2:
            for rule_id, atomic_num, hs in [
                ("HS_RING_CH2_TO_O", 8, 0), ("HS_RING_CH2_TO_NH", 7, 1)
            ]:
                def edit(rw: Chem.RWMol, m=atom_map, z=atomic_num, h=hs) -> list[int]:
                    a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                    a.SetAtomicNum(z)
                    a.SetNumExplicitHs(h)
                    a.SetNoImplicit(True)
                    return []
                candidates[rule_id].append(({atom_map}, edit, f"atom_map={atom_map}", set()))

    for bond in sorted(base.GetBonds(), key=lambda b: b.GetIdx()):
        sizes = ring_bonds.get(bond.GetIdx(), set())
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        simple = sum(bond.GetIdx() in r for r in base.GetRingInfo().BondRings()) == 1
        if (bond.GetBondType() == Chem.BondType.SINGLE and sizes and simple
                and min(sizes) >= 3 and max(sizes) <= 7
                and not a.GetIsAromatic() and not b.GetIsAromatic()
                and a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6):
            ma, mb = a.GetAtomMapNum(), b.GetAtomMapNum()

            def edit(rw: Chem.RWMol, x=ma, y=mb) -> list[int]:
                idx = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms()}
                rw.RemoveBond(idx[x], idx[y])
                n = rw.AddAtom(_new_atom(6, initial_map, 2, 0))
                rw.AddBond(idx[x], n, Chem.BondType.SINGLE)
                rw.AddBond(n, idx[y], Chem.BondType.SINGLE)
                return [initial_map]

            candidates["RE_INSERT_CH2"].append(({ma, mb}, edit, f"bond={ma}-{mb}", set()))

    for atom_map in sorted(parent_maps):
        atom = base.GetAtomWithIdx(pidx[atom_map])
        sizes = ring_atoms.get(atom.GetIdx(), set())
        if not sizes or min(sizes) < 6 or sum(atom.GetIdx() in r for r in base.GetRingInfo().AtomRings()) != 1:
            continue
        if not (atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()
                and atom.GetDegree() == 2
                and atom.GetTotalNumHs(includeNeighbors=True) == 2
                and all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds())):
            continue
        n1, n2 = atom.GetNeighbors()
        if n1.GetIsAromatic() or n2.GetIsAromatic():
            continue
        if base.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx()) is not None:
            continue
        m1, m2 = n1.GetAtomMapNum(), n2.GetAtomMapNum()

        def edit(rw: Chem.RWMol, rm=atom_map, x=m1, y=m2) -> list[int]:
            idx = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms()}
            rw.RemoveAtom(idx[rm])
            idx = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms()}
            rw.AddBond(idx[x], idx[y], Chem.BondType.SINGLE)
            return []

        candidates["RC_REMOVE_CH2"].append(
            ({atom_map, m1, m2}, edit, f"removed_map={atom_map}", {atom_map})
        )

    aliphatic_ch = [
        m for m in sorted(parent_maps)
        if base.GetAtomWithIdx(pidx[m]).GetAtomicNum() == 6
        and not base.GetAtomWithIdx(pidx[m]).GetIsAromatic()
        and base.GetAtomWithIdx(pidx[m]).GetTotalNumHs(includeNeighbors=True) >= 1
    ]
    for pos, ma in enumerate(aliphatic_ch):
        for mb in aliphatic_ch[pos + 1:]:
            ia, ib = pidx[ma], pidx[mb]
            if base.GetBondBetweenAtoms(ia, ib) is not None:
                continue
            distance = len(Chem.GetShortestPath(base, ia, ib)) - 1
            if not 3 <= distance <= 6:
                continue

            def edit(rw: Chem.RWMol, x=ma, y=mb) -> list[int]:
                idx = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms()}
                _consume_h(rw.GetAtomWithIdx(idx[x]))
                _consume_h(rw.GetAtomWithIdx(idx[y]))
                rw.AddBond(idx[x], idx[y], Chem.BondType.SINGLE)
                return []

            candidates["CR_CLOSE_CC"].append(({ma, mb}, edit, f"closure={ma}-{mb}", set()))

    for atom in base.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetFormalCharge() != 0 or atom.GetDegree() != 1:
            continue
        bond = atom.GetBonds()[0]
        carbon = atom.GetNeighbors()[0]
        if bond.GetBondType() == Chem.BondType.DOUBLE and carbon.GetAtomicNum() == 6:
            mo, mc = atom.GetAtomMapNum(), carbon.GetAtomMapNum()

            def edit(rw: Chem.RWMol, m=mo) -> list[int]:
                a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                a.SetAtomicNum(16)
                a.SetNumExplicitHs(0)
                a.SetNoImplicit(True)
                return []

            candidates["BI_CARBONYL_O_TO_S"].append(({mo, mc}, edit, f"carbonyl={mc}-{mo}", set()))

    atom_rings = [set(r) for r in base.GetRingInfo().AtomRings()]
    for hetero in base.GetAtoms():
        source_type = None
        if (hetero.GetAtomicNum() == 8 and hetero.GetFormalCharge() == 0
                and hetero.GetDegree() == 1 and hetero.GetTotalNumHs(includeNeighbors=True) == 1):
            source_type = "EV_RELOCATE_OH"
        elif (hetero.GetAtomicNum() == 7 and hetero.GetFormalCharge() == 0
              and hetero.GetDegree() == 1 and hetero.GetTotalNumHs(includeNeighbors=True) == 2):
            source_type = "EV_RELOCATE_NH2"
        if source_type is None:
            continue
        source = hetero.GetNeighbors()[0]
        if not (source.GetAtomicNum() == 6 and source.GetIsAromatic()):
            continue
        source_rings = [ring for ring in atom_rings if source.GetIdx() in ring]
        for dest in base.GetAtoms():
            if not (dest.GetAtomicNum() == 6 and dest.GetIsAromatic()
                    and dest.GetTotalNumHs(includeNeighbors=True) == 1
                    and any(dest.GetIdx() in ring for ring in source_rings)):
                continue
            mh, ms, md = hetero.GetAtomMapNum(), source.GetAtomMapNum(), dest.GetAtomMapNum()

            source_h_count = source.GetTotalNumHs(includeNeighbors=True)
            dest_h_count = dest.GetTotalNumHs(includeNeighbors=True)

            def edit(
                rw: Chem.RWMol,
                h=mh,
                s=ms,
                d=md,
                source_hs=source_h_count,
                dest_hs=dest_h_count,
            ) -> list[int]:
                idx = {a.GetAtomMapNum(): a.GetIdx() for a in rw.GetAtoms()}
                source_atom = rw.GetAtomWithIdx(idx[s])
                dest_atom = rw.GetAtomWithIdx(idx[d])
                if dest_hs < 1:
                    raise ValueError("NO_REMOVABLE_H")

                # Cache hydrogen counts before changing connectivity. RWMol bond
                # edits invalidate implicit-valence data, so GetTotalNumHs() must
                # not be called between the graph mutation and sanitization.
                source_atom.SetNumExplicitHs(source_hs + 1)
                source_atom.SetNoImplicit(True)
                dest_atom.SetNumExplicitHs(dest_hs - 1)
                dest_atom.SetNoImplicit(True)
                rw.RemoveBond(idx[h], idx[s])
                rw.AddBond(idx[h], idx[d], Chem.BondType.SINGLE)
                return []

            candidates[source_type].append(({mh, ms, md}, edit, f"move={ms}->{md}", set()))

    # Tier-3 templates are always enumerated so exact expert policy can authorize
    # production use. Without exact authorization they remain UNKNOWN and are only
    # emitted when the caller explicitly enables exploratory generation.
    for atom_map in sorted(parent_maps):
        atom = base.GetAtomWithIdx(pidx[atom_map])
        if (atom.GetDegree() == 1 and atom.GetFormalCharge() == 0
                and atom.GetAtomicNum() == 8
                and atom.GetTotalNumHs(includeNeighbors=True) == 1):
            def edit(rw: Chem.RWMol, m=atom_map) -> list[int]:
                a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                a.SetAtomicNum(7)
                a.SetNumExplicitHs(2)
                a.SetNoImplicit(True)
                return []
            candidates["HB_OH_TO_NH2"].append(({atom_map}, edit, f"atom_map={atom_map}", set()))
        if (atom.GetDegree() == 1 and atom.GetFormalCharge() == 0
                and atom.GetAtomicNum() == 7
                and atom.GetTotalNumHs(includeNeighbors=True) == 2):
            def edit(rw: Chem.RWMol, m=atom_map) -> list[int]:
                a = next(x for x in rw.GetAtoms() if x.GetAtomMapNum() == m)
                a.SetAtomicNum(8)
                a.SetNumExplicitHs(1)
                a.SetNoImplicit(True)
                return []
            candidates["HB_NH2_TO_OH"].append(({atom_map}, edit, f"atom_map={atom_map}", set()))

    analogs: list[dict[str, Any]] = []
    rejections: list[dict[str, Any]] = []
    # Deduplicate and count by canonical constitutional graph. Isomeric SMILES are
    # retained for traceability but cannot inflate enumeration through stereo tags.
    seen_constitutional: set[str] = set()
    accepted_by_site: dict[int, int] = defaultdict(int)
    per_site_quota = 64
    global_limit = 400
    rules_in_order = [r["rule_id"] for r in RULE_CATALOG]
    max_depth = max((len(candidates[r]) for r in rules_in_order), default=0)

    for rule_id in rules_in_order:
        if not candidates[rule_id]:
            rejections.append({
                "rule_id": rule_id,
                "reason": "NO_APPLICABLE_SITE",
                "candidate": None,
            })

    stop = False
    for depth in range(max_depth):
        for rule_id in rules_in_order:
            if depth >= len(candidates[rule_id]):
                continue
            touched, editor, label, declared_removed = candidates[rule_id][depth]
            quota_sites = sorted(touched)
            if quota_sites and any(accepted_by_site[m] >= per_site_quota for m in quota_sites):
                rejections.append({
                    "rule_id": rule_id,
                    "reason": "PER_SITE_VALID_ANALOG_QUOTA_REACHED",
                    "candidate": label,
                    "touched_atom_maps": quota_sites,
                    "limit": per_site_quota,
                })
                continue
            statuses = [_status(m, table) for m in sorted(touched)]
            # Permission is evaluated independently at every affected original site.
            # Extension rules require the exact rule ID in expert policy; a broad
            # family permission is intentionally insufficient for production use.
            rule = _RULE_BY_ID[rule_id]
            for s in statuses:
                allowed = (s.get('evidence') or {}).get('SAR', {}).get('allowed_transformations', [])
                rule_class = rule['transformation_class']
                exact_required = bool(rule.get('requires_expert_allowed_rule'))
                permitted = (
                    isinstance(allowed, list)
                    and rule_id in allowed
                    if exact_required
                    else isinstance(allowed, list)
                    and (rule_id in allowed or rule_class in allowed)
                )
                if s['state'] == 'MODIFIABLE' and not permitted:
                    s['state'] = 'UNKNOWN'
                    s['scope_note'] = (
                        'Exact expert allowed rule required'
                        if exact_required
                        else 'Outside cited exact-rule or transformation-class scope'
                    )
            protected = [s["atom_map"] for s in statuses if s["state"] == "PROTECTED"]
            unknown = [s["atom_map"] for s in statuses if s["state"] == "UNKNOWN"]
            if protected:
                rejections.append({
                    "rule_id": rule_id,
                    "reason": "PROTECTED_TOUCHED",
                    "candidate": label,
                    "touched_atom_maps": sorted(touched),
                    "protected_atom_maps": protected,
                    "site_status": statuses,
                })
                continue
            if unknown and not exploratory:
                rejections.append({
                    "rule_id": rule_id,
                    "reason": "UNKNOWN_REQUIRES_EXPLORATORY",
                    "candidate": label,
                    "touched_atom_maps": sorted(touched),
                    "unknown_atom_maps": unknown,
                    "site_status": statuses,
                })
                continue
            try:
                rw = Chem.RWMol(base)
                added = editor(rw)
                product = rw.GetMol()
                Chem.SanitizeMol(product)
                Chem.AssignStereochemistry(product, cleanIt=False, force=True)
                if len(Chem.GetMolFrags(product)) != 1:
                    raise ValueError("MULTICOMPONENT")
                if any(a.GetAtomicNum() == 0 for a in product.GetAtoms()):
                    raise ValueError("DUMMY_ATOM")
                qidx = _map_index(product)
                actual_removed = sorted(parent_maps - set(qidx))
                if set(actual_removed) != declared_removed:
                    raise ValueError("UNDECLARED_ATOM_REMOVAL")
                for atom_map in parent_maps & set(qidx):
                    pa = base.GetAtomWithIdx(pidx[atom_map])
                    qa = product.GetAtomWithIdx(qidx[atom_map])
                    if pa.GetFormalCharge() != qa.GetFormalCharge():
                        raise ValueError("PARENT_FORMAL_CHARGE_CHANGED")
                    if atom_map not in touched and pa.GetNumExplicitHs() != qa.GetNumExplicitHs():
                        raise ValueError("UNTOUCHED_EXPLICIT_H_CHANGED")
                    if not same_mapped_tetrahedral_stereo(pa, qa):
                        raise ValueError("LOCAL_STEREOTAG_CHANGED")
                _, product_bond_stereo = _stereo_signature(product)
                for key, stereo in parent_stereo_bonds.items():
                    if stereo != int(Chem.BondStereo.STEREONONE):
                        if key not in product_bond_stereo or product_bond_stereo[key] != stereo:
                            raise ValueError("STEREOBOND_CHANGED")
                protected_set = {m for m in parent_maps if _status(m, table)['state'] == 'PROTECTED'}
                if not _graph_preserved(base, product, protected_set):
                    raise ValueError('PROTECTED_GRAPH_CHANGED')
                canonical = _without_maps_smiles(product)
                constitutional = _without_maps_smiles(product, isomeric=False)
                if constitutional == parent_constitutional:
                    raise ValueError("UNCHANGED_FROM_PARENT")
                if constitutional in seen_constitutional:
                    raise ValueError("DUPLICATE_CONSTITUTIONAL_GRAPH")
                seen_constitutional.add(constitutional)
                rule = _RULE_BY_ID[rule_id]
                modified = sorted(touched - declared_removed)
                confidence = "low" if unknown or rule["tier"] == 3 else "moderate"
                qualified_for_counts = bool(statuses) and all(
                    status["state"] == "MODIFIABLE" for status in statuses
                )
                record = {
                    "mapped_smiles": Chem.MolToSmiles(
                        product, canonical=True, isomericSmiles=True
                    ),
                    "canonical_smiles": canonical,
                    "constitutional_smiles": constitutional,
                    "parent_id": parent_id,
                    "rule_id": rule_id,
                    "transformation_class": rule["transformation_class"],
                    "tier": rule["tier"],
                    "modified_atom_maps": modified,
                    "removed_atom_maps": actual_removed,
                    "added_atom_maps": sorted(added),
                    "design_hypothesis": rule["design_hypothesis"],
                    "rationale": rule["rationale"],
                    "applicability": rule["applicability"],
                    "review_required": bool(rule["review_required"]),
                    "api_note": rule["api_note"],
                    "expected_benefit": rule.get(
                        "expected_benefit",
                        "The enumerated graph tests the stated structural hypothesis; no biological benefit is asserted.",
                    ),
                    "synthetic_risks": rule.get("synthetic_risks", [
                        "Reaction route, chemoselectivity, stability, and isolated yield require experimental review."
                    ]),
                    "enumeration_scope": rule.get("scope"),
                    "attachment_site_atom_maps": sorted(touched),
                    "original_site_map": dict(original_site_map),
                    "protected_graph_preserved": True,
                    "protected_atom_maps": sorted(protected_set),
                    "pharmacophore_preserved": "review",
                    "preservation_scope": "2D protected graph only; binding features and geometry require pose evaluation",
                    "docking_pose_preserved": "review",
                    "site_status": statuses,
                    "confidence": confidence,
                    "scientific_pending": bool(unknown),
                    "qualified_for_counts": qualified_for_counts,
                }
                analogs.append(record)
                for site_map in quota_sites:
                    accepted_by_site[site_map] += 1
                if len(analogs) >= global_limit:
                    stop = True
                    break
            except Exception as exc:
                rejections.append({
                    "rule_id": rule_id,
                    "reason": str(exc) or exc.__class__.__name__,
                    "candidate": label,
                    "touched_atom_maps": sorted(touched),
                    "site_status": statuses,
                })
        if stop:
            break

    if stop:
        rejections.append({
            "rule_id": None,
            "reason": "GLOBAL_VALID_ANALOG_LIMIT_REACHED",
            "candidate": None,
            "limit": global_limit,
        })

    qualified = [
        record for record in analogs if record.get("qualified_for_counts") is True
    ]
    return {
        "analogs": analogs,
        "rejections": rejections,
        "rule_catalog": [_safe(dict(rule)) for rule in RULE_CATALOG],
        "summary": {
            "counting_basis": "unique canonical constitutional graph; stereochemical labels do not create additional counts",
            "all_passed": {
                "constitutional_graph_count": len(analogs),
                "transformation_families": sorted({
                    record["transformation_class"] for record in analogs
                }),
            },
            "qualified": {
                "constitutional_graph_count": len(qualified),
                "transformation_families": sorted({
                    record["transformation_class"] for record in qualified
                }),
                "criterion": "all affected original sites are MODIFIABLE under applicable exact-rule policy",
            },
        },
    }


def select_diverse(records: list[dict[str, Any]], limit: int = 16) -> list[dict[str, Any]]:
    if limit <= 0 or not records:
        return []
    molecules: list[Chem.Mol] = []
    valid_records: list[dict[str, Any]] = []
    for record in records:
        smiles = record.get("mapped_smiles") or record.get("canonical_smiles")
        mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
        if mol is not None:
            molecules.append(mol)
            valid_records.append(record)
    if not molecules:
        return []

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fingerprints = [generator.GetFingerprint(mol) for mol in molecules]
    distances: list[float] = []
    for i in range(1, len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[i], fingerprints[:i]
        )
        distances.extend(1.0 - similarity for similarity in similarities)

    raw_clusters = Butina.ClusterData(
        distances, len(fingerprints), 0.35, isDistData=True, reordering=True
    )
    clusters = sorted((tuple(cluster) for cluster in raw_clusters), key=lambda c: min(c))
    cluster_for: dict[int, int] = {}
    representatives: list[tuple[int, int]] = []

    for cluster_id, cluster in enumerate(clusters):
        for idx in cluster:
            cluster_for[idx] = cluster_id
        best_idx = min(cluster)
        best_score = -1.0
        for idx in cluster:
            if len(cluster) == 1:
                score = 1.0
            else:
                score = sum(
                    DataStructs.TanimotoSimilarity(fingerprints[idx], fingerprints[j])
                    for j in cluster if j != idx
                ) / (len(cluster) - 1)
            if score > best_score or (score == best_score and idx < best_idx):
                best_idx, best_score = idx, score
        representatives.append((cluster_id, best_idx))

    selected: list[tuple[int, int]] = []
    used_indices: set[int] = set()
    cluster_counts: dict[int, int] = defaultdict(int)
    covered_classes: set[str] = set()
    covered_tiers: set[Any] = set()

    # Keep up to three representatives per real cluster, preferring uncovered classes/tiers.
    representatives = [(cluster_for[i], i) for i in range(len(valid_records))]
    while len(selected) < min(limit, len(representatives)):
        best = None
        best_key = None
        for cluster_id, idx in representatives:
            if idx in used_indices or cluster_counts[cluster_id] >= 3:
                continue
            record = valid_records[idx]
            new_class = record.get("transformation_class") not in covered_classes
            new_tier = record.get("tier") not in covered_tiers
            key = (int(cluster_counts[cluster_id] == 0), int(new_class), int(new_tier), -idx)
            if best_key is None or key > best_key:
                best_key = key
                best = (cluster_id, idx)
        if best is None:
            break
        cluster_id, idx = best
        selected.append(best)
        used_indices.add(idx); cluster_counts[cluster_id] += 1
        covered_classes.add(str(valid_records[idx].get("transformation_class")))
        covered_tiers.add(valid_records[idx].get("tier"))

    output: list[dict[str, Any]] = []
    for cluster_id, idx in selected:
        record = dict(valid_records[idx])
        record["cluster_id"] = cluster_id
        record["cluster_size"] = len(clusters[cluster_id])
        output.append(_safe(record))
    return output

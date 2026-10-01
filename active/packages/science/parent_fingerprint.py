"""Conservative, coordinate-bound parent interaction fingerprint packet.

The pH value is execution context only. No pKa, population, affinity, efficacy,
site-role approval, or protected-site approval is inferred.
"""
from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .chemical_states import enumerate_microstates, interaction_profile, optimize_ligand_hydrogens

VERSION = "parent-fingerprint/20261002.2"
_SENSITIVE = {"directional_hbond", "possible_hbond_contact", "directional_hbond_rejected", "ionic_contact"}
_AROMATIC_RING_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}


def _maps(mol: Chem.Mol, name: str) -> dict[int, int]:
    if not isinstance(mol, Chem.Mol) or mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError(f"{name} requires one actual 3D conformer")
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if xyz.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
        raise ValueError(f"{name} has invalid coordinates")
    result = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        number = atom.GetAtomMapNum()
        if number <= 0 or number in result:
            raise ValueError(f"{name} requires unique positive heavy-atom maps")
        result[number] = atom.GetIdx()
    if not result:
        raise ValueError(f"{name} has no mapped heavy atoms")
    return result


def transfer_heavy_coordinates(state: Chem.Mol, source: Chem.Mol) -> Chem.Mol:
    """Transfer exact mapped heavy coordinates and place unoptimized H at its parent."""
    source_maps = _maps(source, "coordinate source")
    state_maps = {}
    for atom in state.GetAtoms():
        if atom.GetAtomicNum() > 1:
            number = atom.GetAtomMapNum()
            if number <= 0 or number in state_maps:
                raise ValueError("state requires unique positive heavy-atom maps")
            state_maps[number] = atom.GetIdx()
    if set(source_maps) != set(state_maps):
        raise ValueError("state and coordinate source heavy maps differ")
    source_conf = source.GetConformer()
    state = Chem.Mol(state)
    state.RemoveAllConformers()
    conf = Chem.Conformer(state.GetNumAtoms())
    conf.Set3D(True)
    for atom in state.GetAtoms():
        if atom.GetAtomicNum() > 1:
            point = source_conf.GetAtomPosition(source_maps[atom.GetAtomMapNum()])
        else:
            heavy = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
            point = source_conf.GetAtomPosition(source_maps[heavy[0].GetAtomMapNum()]) if len(heavy) == 1 else source_conf.GetAtomPosition(next(iter(source_maps.values())))
        conf.SetAtomPosition(atom.GetIdx(), point)
    state.AddConformer(conf, assignId=True)
    return state


def _protein_id(record: dict[str, Any]) -> str:
    required = ("label_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id")
    if any(key not in record for key in required):
        raise ValueError("protein record lacks exact chain/residue/component/atom identity")
    return ":".join((str(record["label_asym_id"]), str(record["auth_seq_id"]), str(record["label_comp_id"]).upper(), str(record["label_atom_id"]).upper()))


def _protein_index(records: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(records, (list, tuple)) or not records:
        raise ValueError("protein_atoms must be a nonempty exact heavy-atom list")
    result = {}
    for raw in records:
        if not isinstance(raw, dict):
            raise TypeError("protein atom records must be dictionaries")
        record = dict(raw)
        identifier = _protein_id(record)
        if identifier in result:
            raise ValueError(f"duplicate protein atom identity: {identifier}")
        xyz = np.asarray(record.get("xyz"), dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise ValueError(f"invalid coordinates for {identifier}")
        record["xyz"] = xyz
        result[identifier] = record
    return result


def _ligand_maps(participants: Any) -> list[int]:
    result = []
    for value in participants if isinstance(participants, list) else []:
        if not isinstance(value, str):
            continue
        if value.startswith("L:ring:"):
            for token in value[7:].split(","):
                if token.startswith("L:") and token[2:].isdigit():
                    result.append(int(token[2:]))
        elif value.startswith("L:") and value[2:].isdigit():
            result.append(int(value[2:]))
    return sorted(set(result))


def _protein_tokens(participants: Any, protein: dict[str, dict[str, Any]]) -> list[str]:
    exact = []
    for value in participants if isinstance(participants, list) else []:
        if not isinstance(value, str):
            continue
        fields = value.split(":")
        if len(fields) == 4:
            token = f"{fields[0]}:{fields[1]}:{fields[2].upper()}:{fields[3].upper()}"
            if token in protein:
                exact.append(token)
        elif len(fields) == 3:
            chain, seq, comp = fields[0], fields[1], fields[2].upper()
            allowed = _AROMATIC_RING_ATOMS.get(comp, set())
            exact.extend(
                key for key, record in protein.items()
                if key.split(":")[:3] == [chain, seq, comp]
                and str(record["label_atom_id"]).upper() in allowed
            )
    return sorted(set(exact))


def _distance(item: dict[str, Any]) -> float | None:
    for key in ("distance_DA_A", "distance_A", "centroid_distance_A"):
        value = item.get(key)
        if isinstance(value, (int, float)) and np.isfinite(value):
            return float(value)
    return None


def _contacts(profile: dict[str, Any], mol: Chem.Mol, protein: dict[str, dict[str, Any]], protein_h_supplied: bool) -> list[dict[str, Any]]:
    mapping = _maps(mol, "profile molecule")
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    rows = []
    for item in profile.get("interactions", []):
        if not isinstance(item, dict):
            continue
        maps = _ligand_maps(item.get("participants"))
        atoms = _protein_tokens(item.get("participants"), protein)
        kind = str(item.get("kind", "unknown"))
        residue_participant = any(isinstance(v, str) and len(v.split(":")) == 3 for v in item.get("participants", []))
        for ligand_map in maps:
            if ligand_map not in mapping:
                continue
            for protein_id in atoms:
                record = protein[protein_id]
                heavy_distance = float(np.linalg.norm(xyz[mapping[ligand_map]] - record["xyz"]))
                if residue_participant and heavy_distance > 4.0:
                    continue
                directional_available = kind in {"directional_hbond", "directional_hbond_rejected"}
                protein_donor = str(item.get("donor", "")).startswith(protein_id) or item.get("protein_is_donor") is True
                confirmed = kind == "directional_hbond" and (not protein_donor or protein_h_supplied)
                ambiguity = bool(item.get("requires_review")) or kind in {"possible_hbond_contact", "aromatic_contact_geometry"}
                if directional_available and protein_donor and not protein_h_supplied:
                    ambiguity = True
                rows.append({
                    "key": f"{kind}|L:{ligand_map}|{protein_id}",
                    "protein_atom_id": protein_id,
                    "chain": str(record["label_asym_id"]),
                    "auth_seq_id": str(record["auth_seq_id"]),
                    "residue": str(record["label_comp_id"]).upper(),
                    "protein_atom": str(record["label_atom_id"]).upper(),
                    "ligand_atom_map": ligand_map,
                    "interaction_kind": kind,
                    "heavy_distance_A": round(heavy_distance, 4),
                    "profile_distance_A": _distance(item),
                    "directional_hydrogen": ({
                        "distance_HA_A": item.get("distance_HA_A"),
                        "angle_DHA_deg": item.get("angle_DHA_deg"),
                        "tested_hydrogens": item.get("tested_hydrogens"),
                        "confirmed_geometry": confirmed,
                        "ambiguity": "protein donor hydrogen absent" if protein_donor and not protein_h_supplied else None,
                    } if directional_available else None),
                    "protein_hydrogens_supplied": protein_h_supplied,
                    "chemical_role_ambiguity": ambiguity,
                    "source_interaction": item,
                })
    unique = {}
    for row in rows:
        unique[(row["key"], row["heavy_distance_A"])] = row
    return sorted(unique.values(), key=lambda row: (row["protein_atom_id"], row["ligand_atom_map"], row["interaction_kind"]))


def _write_sdf(path: Path, mol: Chem.Mol, properties: dict[str, Any]) -> None:
    copy = Chem.Mol(mol)
    for key, value in properties.items():
        copy.SetProp(str(key), json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else str(value))
    stream = io.StringIO()
    writer = Chem.SDWriter(stream)
    if writer is None:
        raise OSError(f"cannot create SDF writer for: {path}")
    writer.write(copy)
    writer.close()
    path.write_text(stream.getvalue(), encoding="utf-8")


def _svg(mol: Chem.Mol) -> str:
    """Draw a 2D copy labelled with the actual source charge and bound H count."""
    copy = Chem.Mol(mol)
    labels = []
    for atom in copy.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        parts = [f"{atom.GetSymbol()}:{atom.GetAtomMapNum()}"]
        if atom.GetAtomicNum() != 6:
            hydrogen_count = atom.GetTotalNumHs(includeNeighbors=True)
            if hydrogen_count > 0:
                parts.append(f"H{hydrogen_count}")
        formal_charge = atom.GetFormalCharge()
        if formal_charge == 1:
            parts.append("+")
        elif formal_charge == -1:
            parts.append("-")
        elif formal_charge != 0:
            parts.append(f"{formal_charge:+d}")
        label = " ".join(parts)
        atom.SetProp("atomLabel", label)
        labels.append(label)

    copy = Chem.RemoveHs(copy)
    rdDepictor.Compute2DCoords(copy, clearConfs=True)
    drawer = rdMolDraw2D.MolDraw2DSVG(1100, 700)
    drawer.DrawMolecule(copy)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    svg_start = svg.find("<svg")
    svg_open_end = svg.find(">", svg_start)
    if svg_start < 0 or svg_open_end < 0:
        raise ValueError("RDKit did not produce an SVG root element")
    accessibility = (
        "<title>출처 구조의 형식전하·결합 수소 표시, 상태 승인 아님</title>"
        f"<desc>Actual source atom labels: {'; '.join(labels)}</desc>"
    )
    return svg[: svg_open_end + 1] + accessibility + svg[svg_open_end + 1 :]


def _html_report(packet: dict[str, Any]) -> str:
    payload = json.dumps(packet, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    return """<!doctype html><html lang='ko'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>원형 상호작용 근거 패킷</title><style>
*{box-sizing:border-box}body{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1280px;margin:0 auto;padding:clamp(1rem,3vw,2rem);color:#17202a;background:#fff;line-height:1.45}h1,h2,h3{line-height:1.2}h1{font-size:clamp(1.55rem,4vw,2.25rem)}h2{font-size:clamp(1.25rem,3vw,1.65rem);overflow-wrap:anywhere}.warn{background:#fff4d6;border-left:.35rem solid #d49b00;padding:1rem;border-radius:.25rem}.figures{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,28rem),1fr));gap:1rem;margin:1.25rem 0}.figure-card,.panel{border:1px solid #ccd1d1;border-radius:.5rem;background:#fff;padding:1rem;min-width:0}.figure-card{margin:0}.figure-card img{display:block;width:100%;max-width:100%;height:auto;max-height:30rem;object-fit:contain;background:#fff}.figure-card figcaption{font-weight:650;margin-bottom:.75rem}.figure-card a,.links a{overflow-wrap:anywhere}.controls{display:flex;flex-wrap:wrap;gap:.75rem;align-items:center;margin:1.25rem 0}.controls select{max-width:100%;padding:.55rem}.links{display:flex;flex-wrap:wrap;gap:.6rem 1.2rem;margin:.6rem 0}.summary-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,22rem),1fr));gap:1rem;margin:1rem 0}.compact-list{margin:.5rem 0 0;padding-left:1.25rem}.compact-list li{overflow-wrap:anywhere;margin:.2rem 0}.table-scroll{width:100%;overflow-x:auto;-webkit-overflow-scrolling:touch;border:1px solid #ccd1d1;border-radius:.4rem}table{border-collapse:collapse;width:100%;min-width:58rem}th,td{border-right:1px solid #ccd1d1;border-bottom:1px solid #ccd1d1;padding:.5rem;text-align:left;vertical-align:top}th:last-child,td:last-child{border-right:0}tbody tr:last-child td{border-bottom:0}th{background:#eef3f5}.direction{white-space:normal;min-width:13rem}.muted{color:#566573}.pending{font-weight:650;color:#7d4e00}code{font-size:.85rem;overflow-wrap:anywhere}@media(max-width:40rem){body{padding:.8rem}th,td{padding:.4rem}.panel,.figure-card{padding:.75rem}}
</style></head><body><h1>원형 상호작용 근거 패킷</h1><p class='warn'>이 문서는 읽기 전용 근거 생성 결과입니다. 통과 기준, 효능 주장, 생물학적 승인, 전문가 승인, 결합 부위 역할 승인 또는 보호 부위 승인을 의미하지 않습니다. 표시된 pH는 실행 조건일 뿐입니다.</p><section aria-labelledby='drawings-heading'><h2 id='drawings-heading'>생성된 2D 원자 맵 구조</h2><div id='figures' class='figures'></div></section><div class='controls'><label for='state'>열거 상태 선택</label><select id='state'></select></div><h2 id='title'></h2><section class='panel' aria-labelledby='files-heading'><h3 id='files-heading'>선택한 상태의 저장 구조</h3><div id='files' class='links'></div><p class='muted'>원본 파일은 최적화 전 입력으로만 식별됩니다. 파일 이름만으로 수소 원자 포함 여부를 추론하지 않습니다.</p></section><div class='summary-grid'><section class='panel' aria-labelledby='common-heading'><h3 id='common-heading'>열거된 모든 상태에 공통인 접촉</h3><p id='common-count'></p><ul id='common-list' class='compact-list'></ul></section><section class='panel' aria-labelledby='protected-heading'><h3 id='protected-heading'>제안된 보호 상호작용 집합</h3><p class='pending'>전문가 검토 대기 중인 제안일 뿐입니다. 이 패킷은 보호 기준을 확정하지 않습니다.</p><p id='protected-count'></p><ul id='protected-list' class='compact-list'></ul></section></div><h2>원자 수준 접촉</h2><div class='table-scroll' role='region' aria-label='원자 수준 접촉' tabindex='0'><table><thead><tr><th>단백질 원자</th><th>리간드 맵</th><th>상호작용</th><th>중원자 거리 Å</th><th>방향성 기하 진단</th><th>화학적 역할 모호성</th></tr></thead><tbody id='rows'></tbody></table></div><script id='data' type='application/json'>""" + payload + """</script><script>
'use strict';
const d=JSON.parse(document.getElementById('data').textContent);
const stateSelect=document.getElementById('state');
const rows=document.getElementById('rows');
const title=document.getElementById('title');
const filesBox=document.getElementById('files');
const commonCount=document.getElementById('common-count');
const commonList=document.getElementById('common-list');
const protectedCount=document.getElementById('protected-count');
const protectedList=document.getElementById('protected-list');
function add(parent,tag,text,className){const node=document.createElement(tag);if(className)node.className=className;if(text!==undefined)node.textContent=String(text);parent.appendChild(node);return node}
function localPath(value){if(typeof value!=='string')return null;const path=value.trim();if(!path||path.startsWith('//')||path.indexOf(String.fromCharCode(92))!==-1)return null;for(let i=0;i<path.length;i++){const code=path.charCodeAt(i);if(code<32||code===127)return null}const colon=path.indexOf(':');if(colon>0){const first=path.charCodeAt(0);const firstIsLetter=(first>=65&&first<=90)||(first>=97&&first<=122);if(firstIsLetter){let scheme=true;for(let i=1;i<colon;i++){const code=path.charCodeAt(i);const allowed=(code>=65&&code<=90)||(code>=97&&code<=122)||(code>=48&&code<=57)||code===43||code===45||code===46;if(!allowed){scheme=false;break}}if(scheme)return null}}return path}
function addFileLink(parent,path,label){const safe=localPath(path);if(!safe)return;const link=add(parent,'a',label);link.href=safe}
function roleLabel(role){return role==='source_charged_parent'?'전하를 띤 원형':role==='neutral_design_parent'?'중성 설계 원형':String(role)}
function kindLabel(kind){if(kind==='aromatic_contact_geometry')return '방향족 원자 근접성(π-스태킹 확인이 아님)';return String(kind)}
function compactValue(value){if(value===null||value===undefined)return null;if(typeof value==='string'||typeof value==='number'||typeof value==='boolean')return String(value);return null}
function testedHydrogenLines(value){if(!Array.isArray(value)){const scalar=compactValue(value);return scalar===null?[]:['검사한 H: '+scalar]}const lines=[];for(let i=0;i<value.length;i++){const record=value[i];if(!record||typeof record!=='object'||Array.isArray(record)){const scalar=compactValue(record);if(scalar!==null)lines.push('검사한 H '+(i+1)+': '+scalar);continue}const parts=[];let atomId=compactValue(record.hydrogen_atom_id);if(atomId===null)atomId=compactValue(record.atom_id);const ha=compactValue(record.distance_HA_A);const angle=compactValue(record.angle_DHA_deg);if(atomId!==null)parts.push('원자 ID: '+atomId);if(ha!==null)parts.push('H···A: '+ha+' Å');if(angle!==null)parts.push('D–H···A: '+angle+'°');lines.push('검사한 H '+(i+1)+': '+(parts.length?parts.join('; '):'기록됨'))}return lines}
function directionLines(row){const geometry=row.directional_hydrogen;if(!geometry)return ['해당 없음'];const result=[];const ha=compactValue(geometry.distance_HA_A);const angle=compactValue(geometry.angle_DHA_deg);if(ha!==null)result.push('H···A: '+ha+' Å');if(angle!==null)result.push('D–H···A: '+angle+'°');for(const line of testedHydrogenLines(geometry.tested_hydrogens))result.push(line);if(geometry.confirmed_geometry){result.push(row.protein_hydrogens_supplied?'방향성 기하가 관찰됨(진단 전용)':'리간드-H 기하 진단이며 단백질 H는 없음');result.push('생물학적 승인이 아님')}else{result.push('방향성 기하가 확인되지 않음')}if(geometry.ambiguity)result.push(String(geometry.ambiguity));return result}
function drawFigures(){const box=document.getElementById('figures');box.replaceChildren();const paths=d.structure_svgs||{};for(const role of ['source_charged_parent','neutral_design_parent']){const path=localPath(paths[role]);if(!path)continue;const figure=add(box,'figure',undefined,'figure-card');add(figure,'figcaption',roleLabel(role)+' — 생성된 2D 원자 맵 구조');const image=add(figure,'img');image.src=path;image.alt=roleLabel(role)+' 2D 원자 맵 구조';image.loading='eager';addFileLink(figure,path,'생성된 SVG 열기') }if(!box.childNodes.length)add(box,'p','기록된 2D 구조 그림이 없습니다.','muted')}
function fillContactList(target,contacts){target.replaceChildren();for(const row of contacts||[])add(target,'li',row.key||[kindLabel(row.interaction_kind),'L:'+row.ligand_atom_map,row.protein_atom_id].join(' | '));if(!(contacts||[]).length)add(target,'li','기록 없음')}
const states=[];
for(const parent of d.parents||[])for(const state of parent.states||[])states.push({parent:parent,state:state});
for(let i=0;i<states.length;i++){const option=document.createElement('option');option.value=String(i);option.textContent=roleLabel(states[i].parent.parent_role)+' / 상태 '+String(states[i].state.state_index).padStart(2,'0');stateSelect.appendChild(option)}
function draw(){if(!states.length){title.textContent='기록된 열거 상태가 없습니다';return}const selected=states[Number(stateSelect.value)||0];const parent=selected.parent;const state=selected.state;title.textContent=roleLabel(parent.parent_role)+' · 형식 전하 '+state.formal_charge+' · '+state.origin;filesBox.replaceChildren();const structures=state.structures||state.files||{};addFileLink(filesBox,structures.raw_sdf,'최적화 전 원본 SDF(수소 포함 여부를 단정하지 않음)');addFileLink(filesBox,structures.hydrogen_optimized_sdf,'중원자 고정 및 추가 수소 최적화 SDF');if(!filesBox.childNodes.length)add(filesBox,'span','기록된 SDF 경로가 없습니다.','muted');const common=parent.common_across_enumerated_states||[];commonCount.textContent=String(common.length)+'개 접촉';fillContactList(commonList,common);const proposal=parent.proposed_protected_interaction_set||{};const proposed=proposal.contacts||[];protectedCount.textContent=String(proposed.length)+'개 제안 접촉';fillContactList(protectedList,proposed);rows.replaceChildren();for(const row of state.contacts||[]){const tr=document.createElement('tr');for(const value of [row.protein_atom_id,row.ligand_atom_map,kindLabel(row.interaction_kind),row.heavy_distance_A])add(tr,'td',value);const directional=add(tr,'td',undefined,'direction');for(const line of directionLines(row))add(directional,'div',line);add(tr,'td',row.chemical_role_ambiguity?'검토 필요':'모호성 플래그 없음');rows.appendChild(tr)}if(!(state.contacts||[]).length){const tr=document.createElement('tr');const td=add(tr,'td','기록된 원자 수준 접촉이 없습니다.','muted');td.colSpan=6;rows.appendChild(tr)}}
drawFigures();stateSelect.addEventListener('change',draw);draw();
</script></body></html>"""


def build_parent_fingerprint(source_charged: Chem.Mol, neutral_design: Chem.Mol, protein_atoms: Any, *, protein_hydrogens=None, output_dir: str | Path | None = None, pH_context: float = 7.4, max_states: int = 8, seed: int = 23) -> dict[str, Any]:
    """Calculate fingerprints for distinct source-charged and neutral-design parents."""
    if isinstance(max_states, bool) or not isinstance(max_states, int) or not 1 <= max_states <= 16:
        raise ValueError("max_states must be an integer from 1 to 16")
    if not np.isfinite(float(pH_context)):
        raise ValueError("pH_context must be finite")
    protein = _protein_index(protein_atoms)
    protein_h_supplied = bool(protein_hydrogens)
    destination = Path(output_dir).resolve() if output_dir is not None else None
    if destination is not None:
        destination.mkdir(parents=True, exist_ok=False)
    parents = []
    for role, source in (("source_charged_parent", source_charged), ("neutral_design_parent", neutral_design)):
        source_maps = _maps(source, role)
        parent_dir = destination / ("charged" if role.startswith("source") else "neutral") if destination else None
        if parent_dir:
            parent_dir.mkdir()
        states = []
        for index, enumerated in enumerate(enumerate_microstates(Chem.Mol(source), max_states=max_states, pH=pH_context), 1):
            transferred = transfer_heavy_coordinates(enumerated, source)
            structures, receipt = optimize_ligand_hydrogens(transferred, seed=seed)
            raw = structures["input"]
            optimized = structures["hydrogens_optimized"]
            profile = interaction_profile(optimized, list(protein.values()), protein_hydrogens=protein_hydrogens)
            contacts = _contacts(profile, optimized, protein, protein_h_supplied)
            files = None
            if parent_dir:
                raw_name, opt_name = f"state{index:02d}.raw.sdf", f"state{index:02d}.opt.sdf"
                props = {"parent_role": role, "state_index": index, "pH_context_only": pH_context, "pKa": "unknown", "population": "unknown"}
                _write_sdf(parent_dir / raw_name, raw, props)
                _write_sdf(parent_dir / opt_name, optimized, {**props, "optimization_receipt": receipt})
                files = {"raw_sdf": f"{parent_dir.name}/{raw_name}", "hydrogen_optimized_sdf": f"{parent_dir.name}/{opt_name}"}
            states.append({
                "state_index": index,
                "origin": enumerated.GetProp("chemical_state_origin") if enumerated.HasProp("chemical_state_origin") else "unknown",
                "mapped_smiles": Chem.MolToSmiles(enumerated, canonical=True, isomericSmiles=True),
                "formal_charge": int(Chem.GetFormalCharge(enumerated)),
                "heavy_atom_maps": sorted(source_maps),
                "pH_context": float(pH_context), "pKa": "unknown", "population": "unknown",
                "enumeration_exhaustive": False, "status": "calculated",
                "optimization_receipt": receipt, "interaction_profile": profile,
                "contacts": contacts, "files": files,
            })
        key_sets = [{row["key"] for row in state["contacts"]} for state in states]
        common = sorted(set.intersection(*key_sets)) if key_sets else []
        lookup = {row["key"]: row for state in states for row in state["contacts"]}
        common_rows = [lookup[key] for key in common]
        direct = [row for row in common_rows if row["interaction_kind"] not in _SENSITIVE]
        sensitive_variation = []
        all_keys = sorted(set().union(*key_sets)) if key_sets else []
        for key in all_keys:
            row = lookup[key]
            presence = [key in values for values in key_sets]
            if row["interaction_kind"] in _SENSITIVE and not all(presence):
                sensitive_variation.append({"interaction": row, "state_presence": presence, "interpretation": "charge/protonation-sensitive exploratory difference; not a hard gate"})
        parents.append({
            "parent_role": role,
            "source_formal_charge": int(Chem.GetFormalCharge(source)),
            "source_atom_map_identity": sorted(source_maps),
            "states": states,
            "common_across_enumerated_states": common_rows,
            "charge_state_sensitive_exploratory_differences": sensitive_variation,
            "proposed_protected_interaction_set": {
                "status": "calculated_proposal_not_expert_approved",
                "auto_site_role_promotion": False,
                "auto_state_promotion": False,
                "contacts": direct,
                "interpretation": "Actual common direct contacts proposed for review only; fingerprint presence is evidence generation, not a pass criterion.",
            },
        })
    packet = {
        "format": VERSION,
        "parents": parents,
        "protein_atom_count": len(protein),
        "protein_hydrogens_supplied": protein_h_supplied,
        "pH_context": float(pH_context), "pH_is_selection_evidence": False,
        "pKa_prediction_performed": False, "population_prediction_performed": False,
        "enumeration": "truncated_nonexhaustive",
        "scientific_approved": False, "efficacy_claim": None,
        "fingerprint_policy": "presence is evidence generation, not a pass criterion",
        "limitations": [
            "No actual assay pH was supplied; 7.4 is execution context only.",
            "Protein hydrogens were not optimized.",
            "When protein hydrogens are absent, protein-donor direction is not confirmed.",
            "Enumerated states are truncated and nonexhaustive; pKa and populations are unknown.",
            "Contacts do not establish affinity, efficacy, selectivity, or expert approval.",
        ],
    }
    if destination:
        for role, mol, name in (("source_charged_parent", source_charged, "charged-parent.svg"), ("neutral_design_parent", neutral_design, "neutral-parent.svg")):
            (destination / name).write_text(_svg(mol), encoding="utf-8")
        packet["structure_svgs"] = {"source_charged_parent": "charged-parent.svg", "neutral_design_parent": "neutral-parent.svg"}
        (destination / "fingerprint.json").write_text(json.dumps(packet, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
        (destination / "index.html").write_text(_html_report(packet), encoding="utf-8")
    json.dumps(packet, ensure_ascii=False, allow_nan=False)
    return packet


__all__ = ["VERSION", "transfer_heavy_coordinates", "build_parent_fingerprint"]

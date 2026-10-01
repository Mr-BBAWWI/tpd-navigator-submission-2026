"""B-only CPU case preparation. No M2 dispatch, approvals, LLM or GPU execution."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import time
import urllib.request
import xml.etree.ElementTree as ET

import yaml
from jsonschema import Draft202012Validator
from rdkit import Chem

from .molecules import read_ccd, identity, reconstruct, map_warhead, write_molecule, write_sdf
from .structures import ligand_coordinates, attachment_geometry, assembly_proteins


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def collect_sources(config, cache, fetch=False):
    cache.mkdir(parents=True, exist_ok=True)
    result = []
    for source in config['sources']:
        name = source['file']
        if Path(name).name != name or '/' in name or '\\' in name:
            raise ValueError('Source filename must be a basename')
        path = cache / name
        if not path.exists():
            if not fetch:
                raise ValueError(f'Missing source {name}; use --fetch to retrieve public sources')
            if not source['url'].startswith('https://'):
                raise ValueError('Source URL must use HTTPS')
            with urllib.request.urlopen(source['url'], timeout=40) as response:
                data = response.read(20_000_001)
            if len(data) > 20_000_000 or digest(data) != source['sha256']:
                raise ValueError(f'Source changed or exceeded size limit: {name}; inspect before updating case')
            path.write_bytes(data)
        if digest(path.read_bytes()) != source['sha256']:
            raise ValueError(f'Cached source hash mismatch: {name}')
        result.append({**source, 'bytes': path.stat().st_size})
    return result


def evidence(config, cache):
    article = ET.parse(cache/'article.xml').getroot()
    dois = [e.text for e in article.findall(".//article-id[@pub-id-type='doi']")]
    if config['doi'] not in dois:
        raise ValueError('Article DOI does not match selected case')
    annotations = config['literature_annotations']
    locators = [annotations['attachment']['locator']] + [x['locator'] for x in annotations['observations'] if x['locator']]
    for locator in locators:
        if article.find(locator) is None:
            raise ValueError(f'Evidence locator not found: {locator}')
    text = ' '.join(article.itertext())
    for c in config['candidates']:
        if c['paper_name'] not in text:
            raise ValueError('Compound alias not present in the article')
    return {**annotations, 'doi': config['doi'], 'source_file': 'article.xml',
            'source_sha256': digest((cache/'article.xml').read_bytes()),
            'coverage': 'Main-text JATS and figure captions; supplements and figure image panels not reviewed',
            'correction_status': 'not_separately_checked',
            'note': 'Locator and hash verified by code; numerical/chemical attribution needs human review'}


def run_case(config_path, cache, output, fetch=False):
    started = time.perf_counter()
    config_path, cache, output = Path(config_path), Path(cache), Path(output)
    config = json.loads(config_path.read_text(encoding='utf-8'))
    # Refuse to mix a new input with existing reviewable results.
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory is not empty; choose a new run directory')
    sources = collect_sources(config, cache, fetch)
    ev = evidence(config, cache)
    output.mkdir(parents=True, exist_ok=True)
    run_id = 'b-cpu-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    (output/'case-input.json').write_bytes(config_path.read_bytes())
    save_json(output/'sources.json', sources)
    save_json(output/'literature-evidence.json', ev)
    start = config['starting_ligand']
    warhead = read_ccd(cache/(start['ccd']+'.cif'))
    positioned, start_map = ligand_coordinates(cache/(start['pdb']+'.cif'), warhead,
                                               start['ligand_label_asym_id'], start['ccd'])
    geometry = attachment_geometry(cache/(start['pdb']+'.cif'), positioned,
                                    start['ligand_label_asym_id'], start['protein_label_asym_id'],
                                    start['attachment_atom'], config['contact_reporting_radius_A'])
    write_molecule(warhead, output/'starting-ligand')
    save_json(output/'starting-atom-map.json', start_map)
    save_json(output/'attachment-geometry.json', geometry)
    candidates = []
    for c in config['candidates']:
        mol = read_ccd(cache/(c['ccd']+'.cif'))
        experimental, atom_map = ligand_coordinates(cache/(c['pdb']+'.cif'), mol,
                                                   c['ligand_label_asym_id'], c['ccd'])
        proteins, assembly_chains = assembly_proteins(cache/(c['pdb']+'.cif'), c['assembly_id'],
                                                     c['protein_label_asym_ids'])
        if c['ligand_label_asym_id'] not in assembly_chains:
            raise ValueError('Ligand is not in the selected biological assembly')
        target_proteins = [p for p in proteins if any(
            r['pdbx_db_accession'] == config['target']['uniprot']
            and r['pdbx_db_isoform'] == config['target']['paper_construct_isoform']
            for r in p['database_references'])]
        if len(target_proteins) != 1 or not any(
            [int(r['db_align_beg']), int(r['db_align_end'])] == config['target']['paper_residue_range']
            for r in target_proteins[0]['database_alignments']):
            raise ValueError('Selected target construct does not match case accession/isoform/range')
        mapping = map_warhead(warhead, mol, start['attachment_atom'], c['attachment_atom'])
        assembled, parts = reconstruct(mol, c['cuts'])
        molecule_id = c['id']+':'+digest(identity(assembled)['canonical_isomeric_smiles'].encode())
        write_molecule(assembled, output/c['id'])
        experimental.SetProp('coordinate_kind', 'experimental crystal coordinates; not a model prediction')
        write_sdf(experimental, output/(c['id']+'-experimental.sdf'))
        part_rows = []
        for value in parts:
            p = Chem.MolFromSmiles(value)
            handles = {a.GetAtomMapNum() for a in p.GetAtoms() if a.GetAtomicNum() == 0}
            role = {(1001,): 'warhead', (1002,): 'recruiter', (1001, 1002): 'linker'}.get(tuple(sorted(handles)))
            if role is None:
                raise ValueError('Unexpected fragment attachment handles')
            part_rows.append({'role': role, 'mapped_smiles': value,
                              'source': 'explicit reference CCD bond cuts; not a proposed synthetic reagent'})
        save_json(output/(c['id']+'-parts.json'), {'reference_ccd': c['ccd'], 'cuts': c['cuts'],
                                                  'parts': part_rows, 'reference_identity_match': True})
        save_json(output/(c['id']+'-atom-map.json'), {'molecule_id': molecule_id, 'reference_atoms': atom_map,
                  'assembled_sdf_atom_order': [{'sdf_index_one_based': a.GetIdx()+1,
                     'atom_map': a.GetAtomMapNum(), 'ccd_atom_id': a.GetProp('ccd_atom_id')}
                     for a in assembled.GetAtoms()], 'warhead_mapping': mapping})
        # No reference coordinates/templates/constraints go into this prospective input.
        # MSA is deliberately unresolved: execution configuration is a later decision.
        boltz = {'version': 1, 'sequences': [
            {'protein': {'id': p['label_asym_id'], 'sequence': p['sequence']}} for p in proteins]
            + [{'ligand': {'id': 'L', 'smiles': identity(assembled)['canonical_isomeric_smiles']}}]}
        (output/(c['id']+'-boltz-input.yaml')).write_text(yaml.safe_dump(boltz, sort_keys=False), encoding='utf-8')
        candidates.append({**c, 'molecule_id': molecule_id, 'identity': identity(assembled),
                           'origin': 'known_reference_reconstruction', 'proteins': proteins,
                           'reconstruction': 'serialized fragments rejoined; stereo SMILES and CCD InChIKey agree',
                           'synthesis_review': {'status': 'not_reviewed', 'reason': 'Supplementary synthesis not inspected'},
                           'prediction': {'status': 'not_run', 'input': c['id']+'-boltz-input.yaml',
                                          'msa': 'not_prepared', 'model_version': None,
                                          'approval': 'not_granted_by_this_tool'},
                           'human_review': {'status': 'pending', 'reviewer': None}})
    artifacts = [{'path': p.name, 'sha256': digest(p.read_bytes()), 'bytes': p.stat().st_size}
                 for p in sorted(output.iterdir()) if p.is_file()]
    result = {'format': 'tpd-b-case-preparation/0.1.0-draft', 'case_id': config['case_id'],
              'case_version': config['version'], 'run_id': run_id,
              'data_mode': 'real_public_sources_and_cpu_computation',
              'status': 'cpu_prepared_pending_human_review',
              'input_sha256': digest(config_path.read_bytes()), 'target': config['target'],
              'starting_ligand': {**start, 'identity': identity(warhead)},
              'candidates': candidates, 'artifacts': artifacts,
              'runtime': {'python': platform.python_version(), 'platform': platform.platform(),
                          'packages': {p: importlib.metadata.version(p) for p in ['rdkit', 'gemmi', 'numpy', 'PyYAML']},
                          'elapsed_seconds': round(time.perf_counter()-started, 3), 'gpu_used': False,
                          'remote_paid_compute_used': False},
              'implementation_sha256': {p.name: digest(p.read_bytes()) for p in Path(__file__).parent.glob('*.py')},
              'not_completed': ['named pharmaceutical review', 'A/M2 integration and G1',
                                'Boltz runtime/MSA/GPU execution', 'prediction-versus-experiment comparison',
                                'supplementary synthesis review', 'C03', 'G2/G3 and public replay']}
    lines = ['# SMARCA2 B 담당 첫 실행 결과', '',
             '실제 공개 구조로 CPU 준비를 완료했습니다. 약학 검수·G1·GPU 예측은 미완료입니다.', '',
             '| 관리 ID | 논문 이름 | 원문 화합물 번호 | PDB / CCD | 무거운 원자 | 재조립 |',
             '|---|---|---|---|---:|---|']
    for c in candidates:
        lines.append(f"| {c['id']} | {c['paper_name']} | {c['paper_compound_number']} | {c['pdb']} / {c['ccd']} | {c['identity']['heavy_atom_count']} | 기준 구조와 일치 |")
    lines += ['', f"출발 부착 원자: {start['ccd']} {start['attachment_atom']}.",
              f"선택한 복합체에서 해당 원자의 SASA: {geometry['SASA_A2']:.3f} Å²; 가장 가까운 단백질 무거운 원자: {geometry['nearest_protein_heavy_atom_A']:.3f} Å.",
              'SASA와 거리는 기술적 관측이며 부착 허용 또는 효능 판정 기준이 아닙니다.', '',
              '재조립은 알려진 기준 분자를 지정 결합에서 세 부품으로 나누고 다시 연결한 검사입니다. 새 후보 설계나 합성 경로 검증이 아닙니다.',
              '원문의 PROTAC 1/2는 화합물 번호 2/3입니다. ACBI1(5)의 실험값을 붙이지 않았습니다.',
              '양전하 출발 리간드와 중성 완성 후보의 차이, 대칭 때문에 가능한 두 원자 대응을 기록했습니다.',
              'Boltz 파일은 검토용 입력이며, MSA·모델 버전·승인·실행 환경은 아직 준비되지 않았습니다.', '',
              '상세 인계: [handoff.json](handoff.json), [문헌 관측](literature-evidence.json), [부착점 계산](attachment-geometry.json).',
              '전체 원문과 보충자료를 검수했다는 의미가 아닙니다. C01 DC50의 노출 시간은 미확인으로 보존했습니다.']
    (output/'RESULT.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    result['artifacts'] = [{'path': p.name, 'sha256': digest(p.read_bytes()), 'bytes': p.stat().st_size}
                           for p in sorted(output.iterdir()) if p.is_file()]
    schema = Path(__file__).resolve().parents[2]/'contracts/drafts/b_case_preparation.schema.json'
    Draft202012Validator(json.loads(schema.read_text(encoding='utf-8'))).validate(result)
    save_json(output/'handoff.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, type=Path)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--fetch', action='store_true', help='Download missing pinned public sources only')
    args = parser.parse_args()
    try:
        result = run_case(args.case, args.cache, args.out, args.fetch)
    except Exception as exc:
        parser.exit(1, f'Case preparation failed; no success/approval was recorded: {exc}\n')
    print(json.dumps({'status': result['status'], 'run_id': result['run_id'],
                      'candidates': [c['id'] for c in result['candidates']]}, ensure_ascii=False))


if __name__ == '__main__':
    main()

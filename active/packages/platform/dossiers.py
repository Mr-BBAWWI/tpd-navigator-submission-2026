"""M2-owned, immutable A/B review dossiers. Draft imports never authorize compute.

Only JSON/hash/link consistency is checked here. Chemistry remains B's calculation;
developer-curated paper aliases remain provisional until a named reviewer confirms them.
"""
import hashlib
import json
import re
import uuid
from pathlib import Path

from jsonschema import Draft202012Validator

from packages.contracts import ArtifactReader, ContractError, encoded, now, parse_json
from packages.platform.analysis import AnalysisService

ROOT = Path(__file__).resolve().parents[2]
VERSION = 'candidate-dossier-builder-20260923.2'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, code):
    if not condition:
        raise ContractError(code)


def validate_draft(value):
    schema = json.loads((ROOT / 'contracts/drafts/candidate_dossier.schema.json').read_text('utf-8'))
    require(not list(Draft202012Validator(schema).iter_errors(value)), 'DOSSIER_SCHEMA_REJECTED')


def read_package(directory, case_path):
    """Stage all bytes before persisting; reject duplicate/escaping paths and inconsistencies."""
    directory = Path(directory).resolve(strict=True)
    handoff_raw = (directory / 'handoff.json').read_bytes()
    handoff = parse_json(handoff_raw)
    schema = json.loads((ROOT / 'contracts/drafts/b_case_preparation.schema.json').read_text('utf-8'))
    require(not list(Draft202012Validator(schema).iter_errors(handoff)), 'B_DRAFT_SCHEMA_REJECTED')
    case_raw = Path(case_path).read_bytes()
    case = parse_json(case_raw)
    require(handoff['input_sha256'] == digest(case_raw), 'B_CASE_INPUT_CHANGED')
    require(handoff['case_id'] == case['case_id'] and handoff['target'] == case['target'], 'B_TARGET_MISMATCH')
    staged = {}
    for item in handoff['artifacts']:
        name = item['path']
        require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,100}', name)), 'B_UNSAFE_ARTIFACT_PATH')
        path = (directory / name).resolve(strict=True)
        require(path.parent == directory and path.is_file(), 'B_UNSAFE_ARTIFACT_PATH')
        require(name not in staged, 'B_DUPLICATE_ARTIFACT')
        require(item['bytes'] <= 12_000_000 and path.stat().st_size == item['bytes'], 'B_ARTIFACT_SIZE_MISMATCH')
        raw = path.read_bytes()
        require(digest(raw) == item['sha256'], 'B_ARTIFACT_HASH_MISMATCH')
        staged[name] = raw
    require(parse_json(staged['case-input.json']) == case, 'B_CASE_COPY_MISMATCH')
    sources = parse_json(staged['sources.json'])
    require(len(sources) == len(case['sources']), 'B_SOURCE_MANIFEST_MISMATCH')
    for expected in case['sources']:
        found = [s for s in sources if s['file'] == expected['file']]
        require(len(found) == 1 and all(found[0][k] == expected[k] for k in ('file','url','sha256')), 'B_SOURCE_MANIFEST_MISMATCH')
    cs = handoff['candidates']
    require(len(cs) == len(case['candidates']) and len({c['id'] for c in cs}) == len(cs), 'B_CANDIDATE_SET_MISMATCH')
    for expected in case['candidates']:
        matches = [c for c in cs if c['id'] == expected['id']]
        require(len(matches) == 1, 'B_CANDIDATE_SET_MISMATCH')
        c = matches[0]
        require(all(c.get(k) == v for k,v in expected.items()), 'B_CANDIDATE_IDENTITY_MISMATCH')
        cid, smiles = c['id'], c['identity']['canonical_isomeric_smiles']
        require(c['molecule_id'] == cid + ':' + digest(smiles.encode()), 'B_MOLECULE_ID_MISMATCH')
        require(staged[cid+'.smi'].decode().strip() == smiles, 'B_SMILES_MISMATCH')
        parts, atom_map = parse_json(staged[cid+'-parts.json']), parse_json(staged[cid+'-atom-map.json'])
        require(parts['reference_ccd'] == c['ccd'] and parts['reference_identity_match'] is True, 'B_RECONSTRUCTION_MISMATCH')
        require(parts['cuts'] == c['cuts'] and atom_map['molecule_id'] == c['molecule_id'], 'B_ATOM_MAP_MISMATCH')
        pairs = atom_map['warhead_mapping']['atoms']
        require(any(a['source_ccd_atom'] == case['starting_ligand']['attachment_atom'] and a['candidate_ccd_atom'] == c['attachment_atom'] for a in pairs), 'B_ATTACHMENT_MISMATCH')
        require(c['prediction']['input'] == cid+'-boltz-input.yaml' and c['prediction']['input'] in staged, 'B_PREDICTION_INPUT_MISMATCH')
    require(all(handoff['starting_ligand'].get(k) == v for k,v in case['starting_ligand'].items()), 'B_STARTING_LIGAND_MISMATCH')
    return handoff_raw, handoff, case, staged


class DossierService:
    def __init__(self, store, project):
        self.store, self.project, self.port = store, project, store.scope(project)
        with store.db() as db:
            db.execute('''CREATE TABLE IF NOT EXISTS dossiers (
                id TEXT PRIMARY KEY, project TEXT NOT NULL, run_id TEXT NOT NULL,
                fingerprint TEXT NOT NULL, ref TEXT NOT NULL, UNIQUE(project,run_id,fingerprint))''')

    @staticmethod
    def snapshot(run):
        return {k:run[k] for k in ('revision','bundle_ref','target_ref','target_query_ref','data_mode')}

    def import_case(self, run_id, directory, case_path, links_path, analysis_ids=()):
        run = self.store.get_run(self.project, run_id)
        require(run['data_mode'] == 'real' and run['state'] in {'completed','partial'}, 'DOSSIER_REAL_COLLECTION_REQUIRED')
        hraw, handoff, case, staged = read_package(directory, case_path)
        links_raw = Path(links_path).read_bytes()
        plan = parse_json(links_raw)
        require(plan['case_id'] == case['case_id'], 'LITERATURE_CASE_MISMATCH')
        require(len(set(analysis_ids)) == len(analysis_ids), 'DUPLICATE_ANALYSIS')
        snapshot = self.snapshot(run)
        fingerprint = digest(encoded({'builder_version':VERSION,'snapshot':snapshot,'handoff':digest(hraw),'links':digest(links_raw),'analyses':sorted(analysis_ids)}))
        existing = self._find(run_id, fingerprint)
        if existing:
            return self.view(existing), False
        # Verify the original I1 graph and the resolved target, not just display labels.
        for name in ('target_ref','bundle_ref','target_query_ref'):
            ArtifactReader(self.port.read).graph([run[name]], run['data_mode'])
        target = self.port.json(run['target_ref'])
        require(target['resolution_status'] == 'resolved' and target['selected_candidate_id'] == case['target']['uniprot'], 'DOSSIER_TARGET_MISMATCH')
        bundle = self.port.json(run['bundle_ref'])
        source = next(s for s in case['sources'] if s['file']=='article.xml')
        docs = [d for d in bundle['documents'] if d['source_artifact_ref']['sha256'] == source['sha256'] and
                {'namespace':'doi','value':case['doi']} in d['identifiers']]
        require(len(docs) == 1, 'LITERATURE_SOURCE_VERSION_MISMATCH')
        document = docs[0]
        anchors = {}
        for link in plan['anchors']:
            segs = [s for s in bundle['segments'] if s['locator']['document_id'] == document['document_id'] and s['locator'].get('xpath') == link['xpath']]
            require(len(segs) == 1, 'LITERATURE_LOCATOR_MISMATCH')
            seg = segs[0]
            text = self.port.read(seg['text_artifact_ref']).decode('utf-8')
            require(text.count(link['quote']) == 1, 'LITERATURE_ANCHOR_NOT_UNIQUE')
            require(link['id'] not in anchors, 'LITERATURE_DUPLICATE_ANCHOR')
            start = text.index(link['quote'])
            anchors[link['id']] = {**link,'segment_id':seg['segment_id'],'text_ref':seg['text_artifact_ref'],
                'document_id':document['document_id'],'source_ref':document['source_artifact_ref'],
                'start':start,'end':start+len(link['quote'])}
        candidate_links = {c['id']:c for c in plan['candidates']}
        require(set(candidate_links) == {c['id'] for c in case['candidates']}, 'LITERATURE_CANDIDATE_SET_MISMATCH')
        analysis_service = AnalysisService(self.store, self.project)
        analyses = [analysis_service.view(i) for i in analysis_ids]
        require(all(a['run_id']==run_id and a['state']=='review_ready' and a['current_input'] for a in analyses), 'CURRENT_REVIEW_READY_ANALYSIS_REQUIRED')
        for c in case['candidates']:
            link = candidate_links[c['id']]
            require(link['paper_name'] == c['paper_name'] and link['paper_compound_number'] == c['paper_compound_number'], 'LITERATURE_ALIAS_MISMATCH')
            require(anchors[link['identity_anchor']]['quote'] == f"{c['paper_name']} ( {c['paper_compound_number']} )", 'LITERATURE_ALIAS_NOT_ANCHORED')
            require(anchors[link['structure_anchor']]['quote'] == 'PDB ID: '+c['pdb'], 'LITERATURE_STRUCTURE_NOT_ANCHORED')
        # All preflight checks precede writes. Originals are copied verbatim to scoped immutable artifacts.
        artifact_refs = {name:self.port.put_raw(raw, 'application/json' if name.endswith('.json') else 'application/octet-stream', 'computed') for name,raw in staged.items()}
        handoff_ref = self.port.put_raw(hraw, 'application/json', 'computed')
        link_plan_ref = self.port.put_raw(links_raw, 'application/json', 'computed')
        cards = []
        for c in handoff['candidates']:
            cid, link = c['id'], candidate_links[c['id']]
            evidence = []
            for a in analyses:
                b = a['annotated_bundle']
                labels = {x['compound_id'] for x in b['compounds'] if x['source_label']==c['paper_name'] and x['document_id']==document['document_id']}
                for record in b['evidence_records']:
                    if labels.intersection(record['compound_ids']):
                        evidence.append({'analysis_id':a['id'],'assessment_ref':a['assessment_ref'],'review_ref':a['review_ref'],
                            'evidence_id':record['evidence_id'],'observations':record['observations'],'assay_context':record['assay_context'],
                            'target':record['target'],'basis_kind':record['basis_kind'],'statement':record['statement'],
                            'source_locators':record['source_locators'],'review':record['review'],
                            'association':'exact_document_label_via_curated_alias_pending_review'})
            mapping = parse_json(staged[cid+'-atom-map.json'])['warhead_mapping']
            cards.append({'id':cid,'paper_name':c['paper_name'],'paper_compound_number':c['paper_compound_number'],
                'origin':c['origin'],'pdb':c['pdb'],'ccd':c['ccd'],'molecule_id':c['molecule_id'],'identity':c['identity'],
                'anchors':[link['identity_anchor'],link['structure_anchor'],plan['attachment_anchor']],
                'attachment':{'starting_atom':handoff['starting_ligand']['attachment_atom'],'candidate_atom':c['attachment_atom'],
                    'equivalent_mapping_count':mapping['equivalent_mapping_count'],'normalization':mapping['normalization_for_matching_only']},
                'literature_evidence':evidence,'files':{name:ref for name,ref in artifact_refs.items() if name.startswith(cid)},
                'stages':{'literature':'linked_pending_human_review','structure':'B_cpu_result_imported',
                    'attachment':'B_mapping_and_paper_rationale_linked','assembly':'known_reference_reconstructed_by_B',
                    'gpu':'not_run','pharmacy_review':'pending','g1':'not_approved','report':'internal_dossier_only'},
                'missing':['약학 검수','분자팀 정식 인계와 변경 전후 검토','실행 비용 범위·G1 승인','실제 GPU 계산 및 실험 구조 대조','최종 보고·공개 재생']})
        value = {'format':'tpd-candidate-dossier/0.1.0-draft','builder_version':VERSION,'id':'dossier-'+uuid.uuid4().hex,
            'project_id':self.project,'run_id':run_id,'created_at':now(),'case_id':case['case_id'],'snapshot':snapshot,
            'status':'imported_pending_joint_review','handoff_ref':handoff_ref,'link_plan_ref':link_plan_ref,
            'artifact_refs':artifact_refs,'anchors':anchors,'analysis_ids':list(analysis_ids),'candidates':cards,
            'planned_candidates':[{'id':'C03','origin':'new_linker_hypothesis','status':'not_designed','structure':None}],
            'authority':{'human_approved':False,'dispatch_authorized':False,'formal_phase1_result':False},
            'limitations':['개발자가 첫 사례의 문헌 별칭·위치를 지정했습니다. 범용 자동 구조 식별이 아닙니다.',
                'A는 파일·해시·식별자·원자 대응의 일관성을 확인했습니다. B의 화학 계산을 독립 재실행한 검증이 아닙니다.',
                '문헌 값은 실험 조건과 함께 검토해야 하며 재구성·예측은 분해 효능 입증이 아닙니다.',
                '원문/그림 공개 재생 권리는 별도 검토 전입니다.']}
        validate_draft(value)
        ref = self.port.put_json(value)
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            current = json.loads(db.execute('SELECT body FROM runs WHERE id=? AND project=?',(run_id,self.project)).fetchone()[0])
            require(self.snapshot(current)==snapshot, 'DOSSIER_INPUT_CHANGED')
            other = db.execute('SELECT id FROM dossiers WHERE project=? AND run_id=? AND fingerprint=?',(self.project,run_id,fingerprint)).fetchone()
            if other:
                identifier = other[0]
            else:
                db.execute('INSERT INTO dossiers VALUES(?,?,?,?,?)',(value['id'],self.project,run_id,fingerprint,json.dumps(ref)))
                identifier = value['id']
        return self.view(identifier), identifier==value['id']

    def _find(self, run_id, fingerprint):
        with self.store.db() as db:
            row = db.execute('SELECT id FROM dossiers WHERE project=? AND run_id=? AND fingerprint=?',(self.project,run_id,fingerprint)).fetchone()
        return row[0] if row else None

    def list(self, run_id):
        self.store.get_run(self.project, run_id)
        with self.store.db() as db:
            ids = [r[0] for r in db.execute('SELECT id FROM dossiers WHERE project=? AND run_id=? ORDER BY rowid DESC',(self.project,run_id))]
        return [self.view(i) for i in ids]

    def view(self, identifier):
        with self.store.db() as db:
            row = db.execute('SELECT ref FROM dossiers WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if not row:
            raise KeyError(identifier)
        ref = json.loads(row[0])
        value = self.port.json(ref)
        validate_draft(value)
        require(value['project_id']==self.project and value['id']==identifier, 'DOSSIER_CONTEXT_MISMATCH')
        for r in [value['handoff_ref'],value['link_plan_ref'],*value['artifact_refs'].values()]:
            self.port.read(r)
        for name in ('target_ref','bundle_ref','target_query_ref'):
            ArtifactReader(self.port.read).graph([value['snapshot'][name]], value['snapshot']['data_mode'])
        for a in value['anchors'].values():
            text = self.port.read(a['text_ref']).decode('utf-8')
            self.port.read(a['source_ref'])
            require(text[a['start']:a['end']]==a['quote'], 'DOSSIER_CITATION_CHANGED')
        svc = AnalysisService(self.store,self.project)
        for i in value['analysis_ids']:
            svc.view(i)
        current = self.store.get_run(self.project,value['run_id'])
        value['current_input'] = self.snapshot(current)==value['snapshot']
        value['dossier_ref'] = ref
        return value

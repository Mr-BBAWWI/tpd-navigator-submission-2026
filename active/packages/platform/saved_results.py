"""M2 storage of verified B archives and separately versioned expert opinions.

Local operator imports only. No inference, authority transition, or public release.
"""
import copy
import json
import tempfile
import uuid
from pathlib import Path

from packages.contracts import encoded, now, parse_json, validate_payload
from packages.platform.dossiers import DossierService, digest, require
from packages.platform.handoffs import literature_exchange
from packages.science.handoff import compare_observation, HandoffError
from packages.science.evidence_common import safe_relative, schema_check
from packages.science.quality_evidence import tree_files
from packages.science.review_evidence import read_review_evidence

VERSION = 'a-saved-results/20260926.1'
BASE = 'source/quality/evidence/'
AUTHORITY = {'human_approved': False, 'dispatch_authorized': False,
             'efficacy_claim': 'not_established', 'public_release_ready': False}


def reconciliation(dossier, material):
    exchange = literature_exchange(dossier, material)
    rows = []
    for record in exchange['records']:
        observation = record['observation']
        normalized = copy.deepcopy(observation)
        proposed = []
        if observation['endpoint'] == 'DC 50':
            normalized['endpoint'] = 'DC50'
            proposed.append({'field': 'endpoint', 'original': 'DC 50', 'proposed': 'DC50',
                             'meaning': 'label only; does not establish identical experiments'})
        comparison = compare_observation(material, record['query'], normalized)
        refs = [r for r in material['h1']['observations']
                if r['observation_id'] in comparison['reference_observation_ids']]
        rows.append({'record': record, 'original_comparison': compare_observation(material, record['query'], observation),
                     'label_normalized_comparison': comparison, 'proposals': proposed, 'B_references': refs,
                     'direct_performance_comparison': 'not_enabled',
                     'reason': 'Cell line, exposure time, assay, replicates, Dmax and source location require condition-specific review.'})
    return {'records': rows, 'original_exchange': exchange,
            'missing': [{'compound_id': c['id'], 'endpoint': 'DC50', 'value': None,
                         'status': 'not_identified_in_reviewed_segments'} for c in dossier['candidates']
                        if not any(r['query']['paper_name'] == c['paper_name'] and r['observation']['endpoint'] in ('DC50','DC 50') for r in exchange['records'])]}


def projection(bundle, files, file_refs, dossier):
    collection = parse_json(files[BASE+'collection.json'])
    runs = {r['run_id']: r for r in collection['runs']}
    def reference(prefix, ref):
        if ref is None:
            return None
        name = prefix + ref['path']
        require(name in file_refs and file_refs[name]['sha256'] == ref['sha256'], 'SAVED_REPORT_REFERENCE')
        return file_refs[name]
    candidates = []
    for c in bundle['candidates']:
        samples = []
        predictions = {(p['run_id'],p['rank']): p for p in c['evidence']['prediction']['samples']}
        for s in c['samples']:
            run = runs[s['run_id']]; prediction = predictions[(s['run_id'],s['rank'])]
            comp_ref = prediction['comparison_file']
            comp = parse_json(files[BASE+comp_ref['path']]) if comp_ref else None
            comparisons = comp['comparison'] if comp else None
            mappings = comparisons['ligand']['all_mappings'] if comparisons else []
            # Ranges come from all saved symmetry mappings, never independently chosen best poses.
            parts = {}
            for part in ('warhead','linker','recruiter'):
                vals = [m['parts_rmsd_A'][part] for m in mappings if m['parts_rmsd_A'].get(part) is not None]
                parts[part] = [min(vals),max(vals)] if vals else None
            contact = s['ligand_contacts']['report']
            counts = {m['hydrogen_convention']: m['bad_overlap_pair_count'] for m in contact['models']} if contact else {}
            roles = []
            if contact:
                for m in contact['models']:
                    for direction, result in m['directions'].items():
                        if result['typing']:
                            roles.extend({'convention':m['hydrogen_convention'],'direction':direction,**r}
                                         for r in result['typing']['role_disagreements'])
            samples.append({'run_id':s['run_id'],'rank':s['rank'],'binding':s['binding'],
                'condition_digest':run['condition_digest'],'reported_execution':run['reported_execution'],
                'summary':prediction['summary'],'parts_rmsd_ranges_A':parts,
                'mapping_count':len(mappings),'alignment_method':comparisons['method'] if comparisons else None,
                'part_range_meaning':'All symmetry mappings after the same target alignment; range endpoints may use different mappings.',
                'status':{'structure':s['prior_quality']['structure_status'],
                          'preparation':s['ligand_preparation']['status'],'contacts':s['ligand_contacts']['status']},
                'contact_interpretation':'preliminary','bad_overlap_pairs':counts,'role_disagreements':roles,
                'review_issues':s['review_issues'],
                'reports':{'comparison':reference(BASE,comp_ref),
                           'prior_quality':reference('source/quality/',s['prior_quality']['quality_report']),
                           'preparation':reference('',s['ligand_preparation']['report_ref']),
                           'contacts':reference('',s['ligand_contacts']['report_ref'])}})
        candidates.append({'compound_id':c['compound_id'],'molecule_id':c['molecule_id'],
                           'facts':c['evidence'],'synthesis_assessment':c['synthesis_assessment'],'samples':samples})
    return {'version':VERSION,'summary':bundle['summary'],'candidates':candidates,
            'source_authority':bundle['authority'],'authority':AUTHORITY.copy(),
            'literature':reconciliation(dossier,parse_json(files[BASE+'material.json'])),
            'interpretation':'Known reference calibration cases; descriptive structural reproduction, not prospective generalization or efficacy prediction.'}


class SavedResultsService:
    def __init__(self, store, project):
        self.store, self.project, self.port = store, project, store.scope(project)
        self.dossiers = DossierService(store, project)
        with store.db() as db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS saved_results(id TEXT PRIMARY KEY, project TEXT NOT NULL,
                    dossier_id TEXT NOT NULL, fingerprint TEXT NOT NULL, ref TEXT NOT NULL,
                    UNIQUE(project,dossier_id,fingerprint));
                CREATE TABLE IF NOT EXISTS saved_expert_reviews(id TEXT PRIMARY KEY, project TEXT NOT NULL,
                    result_id TEXT NOT NULL, fingerprint TEXT NOT NULL, ref TEXT NOT NULL,
                    UNIQUE(project,result_id,fingerprint));
            ''')

    def _write(self, db, paths, raw, media='application/json', provenance='computed'):
        identifier = 'a-'+uuid.uuid4().hex
        ref = {'artifact_id':identifier,'version':1,'sha256':digest(raw),'media_type':media,
               'schema_id':'urn:tpd-navigator:raw:1','provenance':provenance}
        validate_payload('ArtifactRef', ref)
        path = self.store.root/'blobs'/identifier
        with path.open('xb') as f:
            paths.append(path); f.write(raw)
        require(digest(path.read_bytes()) == ref['sha256'], 'SAVED_COPY_HASH')
        db.execute('INSERT INTO artifacts VALUES(?,?,?)',(identifier,self.project,json.dumps(ref)))
        return ref

    def _atomic(self, table, scope_column, scope_id, fingerprint, build, validate_current):
        # Artifact registrations and the envelope become visible in the same transaction.
        paths=[]
        try:
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE')
                validate_current(db)
                row=db.execute(f'SELECT id FROM {table} WHERE project=? AND {scope_column}=? AND fingerprint=?',
                               (self.project,scope_id,fingerprint)).fetchone()
                if row: return row[0],False
                value=build(lambda raw,media='application/json',provenance='computed':self._write(db,paths,raw,media,provenance))
                ref=self._write(db,paths,encoded(value))
                db.execute(f'INSERT INTO {table} VALUES(?,?,?,?,?)',
                           (value['id'],self.project,scope_id,fingerprint,json.dumps(ref)))
                return value['id'],True
        except BaseException:
            # Only files exclusively created by this transaction; no recursive removal.
            for path in paths: path.unlink(missing_ok=True)
            raise

    def import_archive(self, dossier_id, directory):
        dossier=self.dossiers.view(dossier_id)
        require(dossier['current_input'], 'SAVED_STALE_DOSSIER')
        require(not self.store.get_run(self.project,dossier['run_id']).get('replay_only'), 'SAVED_REPLAY_READ_ONLY')
        # Freeze bytes first. The authoritative B reader checks the frozen copy, not a racing source.
        try:
            files=tree_files(directory)
            with tempfile.TemporaryDirectory(prefix='tpd-review-') as tmp:
                root=Path(tmp)
                for name,raw in files.items():
                    target=safe_relative(root,name);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
                bundle=read_review_evidence(root)
        except (HandoffError,ValueError,OSError,KeyError) as e:
            require(False,'SAVED_ARCHIVE_INVALID: '+str(e))
        material=parse_json(files[BASE+'material.json'])
        require(bundle['data_mode']=='real' and bundle['case_id']==dossier['case_id'], 'SAVED_CASE_MISMATCH')
        require(material['source_snapshot']['handoff_sha256']==dossier['handoff_ref']['sha256'], 'SAVED_CPU_SNAPSHOT')
        require(files[BASE+'cpu/handoff.json']==self.port.read(dossier['handoff_ref']), 'SAVED_CPU_HANDOFF')
        for name,ref in dossier['artifact_refs'].items():
            require(files.get(BASE+'cpu/'+name)==self.port.read(ref),'SAVED_CPU_ARTIFACT')
        require({c['id']:c['molecule_id'] for c in dossier['candidates']} ==
                {c['compound_id']:c['molecule_id'] for c in bundle['candidates']}, 'SAVED_MOLECULE_MISMATCH')
        fingerprint=digest(encoded({'version':VERSION,'dossier_ref':dossier['dossier_ref'],'manifest':digest(files['manifest.json'])}))
        def current(db):
            run=json.loads(db.execute('SELECT body FROM runs WHERE project=? AND id=?',(self.project,dossier['run_id'])).fetchone()[0])
            require(DossierService.snapshot(run)==dossier['snapshot'] and not run.get('replay_only'), 'SAVED_INPUT_CHANGED')
        def build(put):
            refs={name:put(raw,'application/json' if name.endswith('.json') else 'application/octet-stream') for name,raw in files.items()}
            projected=projection(bundle,files,refs,dossier)
            value={'id':'saved-'+uuid.uuid4().hex,'version':VERSION,'project_id':self.project,
                   'dossier_id':dossier_id,'dossier_ref':dossier['dossier_ref'],'run_id':dossier['run_id'],
                   'snapshot':dossier['snapshot'],'created_at':now(),'bundle_digest':bundle['digest'],
                   'manifest_sha256':digest(files['manifest.json']),'files':refs,'summary':bundle['summary'],
                   'projection_ref':put(encoded(projected)),'authority':AUTHORITY.copy()}
            schema_check(value,'a_saved_results.schema.json')
            return value
        identifier,fresh=self._atomic('saved_results','dossier_id',dossier_id,fingerprint,build,current)
        return self.view(identifier),fresh

    def _record(self, identifier):
        with self.store.db() as db:
            row=db.execute('SELECT ref FROM saved_results WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if row is None: raise KeyError(identifier)
        ref=json.loads(row[0]);value=self.port.json(ref)
        schema_check(value,'a_saved_results.schema.json')
        require(value['project_id']==self.project and value['id']==identifier,'SAVED_CONTEXT')
        dossier=self.dossiers.view(value['dossier_id'])
        require(dossier['dossier_ref']==value['dossier_ref'],'SAVED_DOSSIER_CHANGED')
        value.update(current_input=dossier['current_input'],record_ref=ref)
        return value

    def view(self, identifier):
        record=self._record(identifier)
        value=self.port.json(record['projection_ref'])
        value.update({k:v for k,v in record.items() if k!='files'})
        value['file_count']=len(record['files'])
        value['expert_reviews']=self.reviews(identifier,record)
        return value

    def list(self, dossier_id):
        self.dossiers.view(dossier_id)
        with self.store.db() as db:
            ids=[r[0] for r in db.execute('SELECT id FROM saved_results WHERE project=? AND dossier_id=? ORDER BY rowid DESC',(self.project,dossier_id))]
        return [{k:v for k,v in self.view(i).items() if k in ('id','created_at','summary','bundle_digest','current_input','file_count')} for i in ids]

    def files(self, identifier):
        return self._record(identifier)['files']

    def verify(self, identifier):
        record=self._record(identifier)
        for ref in record['files'].values(): self.port.read(ref)
        self.port.read(record['projection_ref'])
        self.reviews(identifier,record)
        return {'files_verified':len(record['files']),'bundle_digest':record['bundle_digest']}

    def import_review(self, identifier, source, assessment):
        record=self._record(identifier);raw=Path(source).read_bytes()
        schema_check(assessment,'a_expert_response.schema.json')
        require(digest(raw)==assessment['source_sha256'],'REVIEW_SOURCE_HASH')
        require(assessment['bundle_digest']==record['bundle_digest'],'REVIEW_STALE_BUNDLE')
        # Bound original text, not just a free-form summary supplied by the operator.
        from packages.platform.expert_responses import extract_answers
        answers=extract_answers(raw)
        require(len({q['question_id'] for q in assessment['questions']})==len(assessment['questions']),'REVIEW_DUPLICATE_QUESTION')
        require(all(answers.get(q['question_id'])==q['answer_original'] for q in assessment['questions']), 'REVIEW_ANSWER_SOURCE')
        for supplement in assessment['literature_supplements']:
            source=supplement['source']
            source_ref=self.store.reference(self.project,source['artifact_id'])
            require(source_ref['sha256']==source['sha256'],'REVIEW_LITERATURE_VERSION')
            import xml.etree.ElementTree as ET
            nodes=ET.fromstring(self.port.read(source_ref)).findall(source['locator'])
            require(len(nodes)==1 and ' '.join(' '.join(nodes[0].itertext()).split())==source['quote'],'REVIEW_LITERATURE_QUOTE')
        fingerprint=digest(encoded(assessment))
        def build(put):
            return {'id':'review-'+uuid.uuid4().hex,'result_id':identifier,'dossier_ref':record['dossier_ref'],
                    'received_at':now(),'source_ref':put(raw,'application/vnd.openxmlformats-officedocument.wordprocessingml.document','source'),
                    'assessment':assessment,'authority':AUTHORITY.copy()}
        def current(db):
            row=db.execute('SELECT ref FROM saved_results WHERE id=? AND project=?',(identifier,self.project)).fetchone()
            require(row is not None and json.loads(row[0])==record['record_ref'],'REVIEW_INPUT_CHANGED')
        rid,fresh=self._atomic('saved_expert_reviews','result_id',identifier,fingerprint,build,current)
        return next(r for r in self.reviews(identifier,record) if r['id']==rid),fresh

    def reviews(self, identifier, record):
        with self.store.db() as db:
            rows=db.execute('SELECT ref FROM saved_expert_reviews WHERE project=? AND result_id=? ORDER BY rowid',(self.project,identifier)).fetchall()
        result=[]
        for row in rows:
            ref=json.loads(row[0]);v=self.port.json(ref)
            require(v['result_id']==identifier and v['dossier_ref']==record['dossier_ref'] and v['assessment']['bundle_digest']==record['bundle_digest'],'REVIEW_BINDING')
            self.port.read(v['source_ref']);schema_check(v['assessment'],'a_expert_response.schema.json')
            v['record_ref']=ref;result.append(v)
        return result

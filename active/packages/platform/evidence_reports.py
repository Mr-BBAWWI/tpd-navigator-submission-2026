"""M2 stores immutable M6 reports assembled from saved evidence, without inference.

Freshness is a comparison of exact source versions, never a scientific approval.
"""
import json
import uuid
from pathlib import Path

from packages.contracts import encoded, now, validate_payload
from packages.platform.dossiers import DossierService, digest, require
from packages.platform.saved_results import SavedResultsService, AUTHORITY, BASE
from packages.science.evidence_common import schema_check
from packages.reporting.saved_report import assemble, render_markdown

ROOT = Path(__file__).resolve().parents[2]
VERSION = 'saved-evidence-report/20260926.1'
IDENTITY_FILES = ('packages/platform/evidence_reports.py', 'packages/reporting/saved_report.py',
                  'packages/platform/saved_results.py', 'packages/platform/dossiers.py',
                  'packages/platform/handoffs.py', 'packages/platform/expert_responses.py',
                  'contracts/drafts/a_evidence_report.schema.json')


def implementation():
    return {'version': VERSION, 'sha256': {p: digest((ROOT/p).read_bytes()) for p in IDENTITY_FILES}}


def references(value):
    """Enumerate direct and nested provenance links without interpreting their text."""
    found = {}
    aliases = value.get('source_aliases',{})
    def walk(item):
        if isinstance(item, dict):
            if {'artifact_id','version','sha256','media_type','schema_id','provenance'} <= item.keys():
                validate_payload('ArtifactRef', item)
                if item['artifact_id'].startswith('b-material-'):
                    original = item
                    require(original['artifact_id'] in aliases, 'REPORT_UNBOUND_B_SOURCE')
                    item = aliases[original['artifact_id']]
                    require(item['sha256']==original['sha256'], 'REPORT_B_SOURCE_HASH')
                old = found.setdefault(item['artifact_id'], item)
                require(old == item, 'REPORT_CONFLICTING_SOURCE_REF')
            else:
                for v in item.values(): walk(v)
        elif isinstance(item, list):
            for v in item: walk(v)
    walk(value)
    return [found[k] for k in sorted(found)]


def bind_b_sources(body, files):
    """Keep B's original namespace and bind each explicit CPU file path to M2."""
    aliases = {}
    def walk(item):
        if isinstance(item,dict):
            ref = item.get('ref')
            if isinstance(ref,dict) and str(ref.get('artifact_id','')).startswith('b-material-'):
                path = BASE+'cpu/'+item['path']
                local = files.get(path)
                require(local is not None and local['sha256']==ref['sha256'], 'REPORT_B_SOURCE_BINDING')
                old = aliases.setdefault(ref['artifact_id'],local)
                require(old==local, 'REPORT_CONFLICTING_B_SOURCE')
            for v in item.values(): walk(v)
        elif isinstance(item,list):
            for v in item: walk(v)
    walk(body)
    body['source_aliases'] = aliases


class EvidenceReportService:
    def __init__(self, store, project):
        self.store, self.project, self.port = store, project, store.scope(project)
        self.saved = SavedResultsService(store, project)
        with store.db() as db:
            db.execute('''CREATE TABLE IF NOT EXISTS evidence_reports(
                id TEXT PRIMARY KEY, project TEXT NOT NULL, result_id TEXT NOT NULL,
                fingerprint TEXT NOT NULL, ref TEXT NOT NULL, UNIQUE(project,result_id,fingerprint))''')

    def _state(self, db, result_id, dossier_id, run_id):
        def ref(table, identifier):
            row = db.execute(f'SELECT ref FROM {table} WHERE project=? AND id=?', (self.project,identifier)).fetchone()
            if row is None: raise KeyError(identifier)
            return json.loads(row[0])
        def latest(table, column, identifier):
            row = db.execute(f'SELECT id FROM {table} WHERE project=? AND {column}=? ORDER BY rowid DESC LIMIT 1',
                             (self.project,identifier)).fetchone()
            return row[0] if row else None
        row = db.execute('SELECT body FROM runs WHERE project=? AND id=?',(self.project,run_id)).fetchone()
        if row is None: raise KeyError(run_id)
        run = json.loads(row[0])
        return {'run_snapshot':DossierService.snapshot(run), 'replay_only':bool(run.get('replay_only')),
                'dossier_ref':ref('dossiers',dossier_id), 'result_ref':ref('saved_results',result_id),
                'latest_dossier_id':latest('dossiers','run_id',run_id),
                'latest_result_id':latest('saved_results','dossier_id',dossier_id),
                'expert_reviews':[{'id':r[0], 'ref':json.loads(r[1])} for r in db.execute(
                    'SELECT id,ref FROM saved_expert_reviews WHERE project=? AND result_id=? ORDER BY rowid',
                    (self.project,result_id))]}

    def prepare(self, result_id):
        """Return the exact token the local screen must send back to create a report."""
        data = self.saved.view(result_id)
        with self.store.db() as db:
            db.execute('BEGIN')
            state = self._state(db,result_id,data['dossier_id'],data['run_id'])
        require(state['result_ref']==data['record_ref'] and state['dossier_ref']==data['dossier_ref'] and
                state['expert_reviews']==[{'id':r['id'],'ref':r['record_ref']} for r in data['expert_reviews']], 'REPORT_INPUT_CHANGED')
        blocked = []
        if state['run_snapshot'] != data['snapshot']: blocked.append('input_changed')
        if state['latest_dossier_id'] != data['dossier_id']: blocked.append('dossier_superseded')
        if state['latest_result_id'] != result_id: blocked.append('result_superseded')
        if state['replay_only']: blocked.append('replay_read_only')
        binding = {'state':state, 'implementation':implementation()}
        return {'result_id':result_id, 'dossier_id':data['dossier_id'], 'run_id':data['run_id'],
                'input_digest':digest(encoded(binding)), 'binding':binding, 'blocked_reasons':blocked,
                'can_create':not blocked, 'authority':AUTHORITY.copy()}

    def create(self, result_id, expected_digest):
        preparation = self.prepare(result_id)
        require(preparation['can_create'], 'REPORT_CURRENT_WRITABLE_INPUT_REQUIRED')
        require(preparation['input_digest']==expected_digest, 'REPORT_INPUT_CHANGED')
        data = self.saved.view(result_id)
        dossier = self.saved.dossiers.view(data['dossier_id'])
        checked = self.saved.verify(result_id)
        body = assemble(data,dossier)
        bind_b_sources(body,self.saved.files(result_id))
        sources = references(body)
        for source in sources: self.port.read(source)
        binding = preparation['binding']
        # A response arriving during source reads cannot be silently included under an old token.
        require(data['record_ref']==binding['state']['result_ref'] and dossier['dossier_ref']==binding['state']['dossier_ref'] and
                [{'id':r['id'],'ref':r['record_ref']} for r in data['expert_reviews']]==binding['state']['expert_reviews'], 'REPORT_INPUT_CHANGED')
        def current(db):
            require(self._state(db,result_id,data['dossier_id'],data['run_id'])==binding['state'] and
                    implementation()==binding['implementation'], 'REPORT_INPUT_CHANGED')
        def build(put):
            value = {'id':'report-'+uuid.uuid4().hex, 'version':VERSION, 'project_id':self.project,
                     'result_id':result_id, 'dossier_id':data['dossier_id'], 'run_id':data['run_id'], 'created_at':now(),
                     'input_digest':expected_digest, 'binding':binding, 'bundle_digest':data['bundle_digest'],
                     'content_ref':put(encoded(body)), 'source_catalog':sources,
                     'verification':{'bundle_files_verified_at_creation':checked['files_verified'],
                                     'direct_sources_verified_at_creation':len(sources)},
                     'counts':body['counts'], 'status':'draft_pending_human_review', 'authority':AUTHORITY.copy()}
            value['markdown_ref'] = put(render_markdown(body,value).encode('utf-8'),'text/markdown')
            schema_check(value,'a_evidence_report.schema.json')
            return value
        identifier, fresh = self.saved._atomic('evidence_reports','result_id',result_id,expected_digest,build,current)
        return self.view(identifier), fresh

    def view(self, identifier):
        with self.store.db() as db:
            row = db.execute('SELECT ref FROM evidence_reports WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if row is None: raise KeyError(identifier)
        ref = json.loads(row[0]); value = self.port.json(ref)
        schema_check(value,'a_evidence_report.schema.json')
        require(value['id']==identifier and value['project_id']==self.project, 'REPORT_CONTEXT')
        require(digest(encoded(value['binding']))==value['input_digest'], 'REPORT_BINDING_DIGEST')
        content = self.port.json(value['content_ref']); self.port.read(value['markdown_ref'])
        require(references(content)==value['source_catalog'], 'REPORT_SOURCE_CATALOG')
        for source in value['source_catalog']: self.port.read(source)
        with self.store.db() as db:
            db.execute('BEGIN')
            current = self._state(db,value['result_id'],value['dossier_id'],value['run_id'])
        old = value['binding']['state']; reasons = []
        for field, reason in [('run_snapshot','input_changed'),('replay_only','storage_mode_changed'),
            ('dossier_ref','dossier_changed'),('result_ref','result_changed'),
            ('latest_dossier_id','dossier_superseded'),('latest_result_id','result_superseded'),
            ('expert_reviews','expert_opinions_changed')]:
            if old[field]!=current[field]: reasons.append(reason)
        if value['binding']['implementation']!=implementation(): reasons.append('report_implementation_changed')
        value.update(record_ref=ref, content=content, freshness={'current':not reasons, 'reasons':reasons,
                    'review_state':'requires_re_review' if reasons else 'pending_human_review', 'checked_at':now()})
        return value

    def list(self, result_id):
        self.saved._record(result_id)
        with self.store.db() as db:
            ids = [r[0] for r in db.execute('SELECT id FROM evidence_reports WHERE project=? AND result_id=? ORDER BY rowid DESC',
                                          (self.project,result_id))]
        return [{k:v for k,v in self.view(i).items() if k in
                 ('id','created_at','input_digest','counts','freshness','status','markdown_ref','content_ref')} for i in ids]

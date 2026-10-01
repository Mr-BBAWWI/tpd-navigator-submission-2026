"""A consumes B's versioned H1/H2 materials; unfilled execution inputs stay pending."""
import json
import uuid
from pathlib import Path

from packages.contracts import encoded, now, parse_json
from packages.platform.dossiers import DossierService, digest, require
from packages.science.candidate_evidence import build_evidence
from packages.science.evidence_cli import read_export
from packages.science.handoff import match_compound
from packages.science.evidence_common import schema_check

VERSION = 'ab-review-input-20260924.1'


def literature_exchange(dossier, material):
    """Preserve the complete A record; no unit conversion or guessed missing values."""
    records = []
    for candidate in dossier['candidates']:
        query = dict(doi=material['h1']['compounds'][0]['paper']['doi'],
                     paper_name=candidate['paper_name'], compound_number=candidate['paper_compound_number'])
        match = match_compound(material, **query, expected_digest=material['material_digest'])
        require(match['status'] == 'matched' and match['molecule_id'] == candidate['molecule_id'], 'AB_IDENTITY_MISMATCH')
        for evidence in candidate['literature_evidence']:
            locator = evidence['source_locators'][0]
            for index, observation in enumerate(evidence['observations']):
                measurement, context = observation['measurement'], evidence['assay_context']
                value = measurement.get('value')
                if not isinstance(value, (int, float, str, type(None))):
                    value = measurement.get('raw_text')
                records.append({'record_id':evidence['analysis_id']+':'+evidence['evidence_id']+':'+str(index),
                    'query':query, 'source':{'doi':query['doi'], 'sha256':locator['source_artifact_ref']['sha256'],
                        'locator':locator.get('xpath') or locator['block_id'], 'review_status':'unreviewed'},
                    'observation':{'endpoint':observation['endpoint_label'],
                        'kind':'missing' if measurement['kind']=='missing' else 'literature_observation' if evidence['basis_kind']=='experiment' else evidence['basis_kind'],
                        'value':value, 'unit':measurement.get('unit',{}).get('value'), 'relation':measurement.get('relation'),
                        'target':evidence['target']['reported_label']['value'], 'cell_line':context['cell_line']['value'],
                        'assay':context['assay_type']['value'], 'exposure_time_h':None,
                        'locator':locator.get('xpath') or locator['block_id'],
                        'A_original':evidence, 'A_observation':observation,
                        'adapter_note':'Original A conditions/units retained in A_original. No duration or replicate conversion; differences require review.'}})
    return {'format':'tpd-b-literature-exchange/0.1.0-draft','producer':'A LiteratureAssessment via '+VERSION,
            'material_digest':material['material_digest'],'data_mode':'real','records':records}


def pending_inputs():
    # Role assignment only. These entries are not fabricated B results or reviewer identities.
    return [
        {'id':'phase1','owner':'B + A contract mapping','status':'pending','value':None,'reason':'H1/H2 is available; formal Phase1 and same-universe CPU preview are not yet supplied.'},
        {'id':'execution_scope','owner':'A/B','status':'pending','value':None,'reason':'Exact calculation settings and cost scope await agreement.'},
        {'id':'pharmacy_reviewer','owner':'team','status':'pending','value':None,'reason':'Named pharmaceutical reviewer and actual review are not provided.'},
        {'id':'gpu_result','owner':'B','status':'pending','value':None,'reason':'No actual prediction in the current B export.'}]


class HandoffService:
    def __init__(self, store, project):
        self.store, self.project, self.port = store, project, store.scope(project)
        self.dossiers = DossierService(store, project)
        with store.db() as db:
            db.execute('''CREATE TABLE IF NOT EXISTS b_inputs(
                id TEXT PRIMARY KEY, project TEXT NOT NULL, run_id TEXT NOT NULL,
                dossier_id TEXT NOT NULL, fingerprint TEXT NOT NULL, ref TEXT NOT NULL,
                UNIQUE(project,dossier_id,fingerprint))''')

    def import_export(self, dossier_id, directory):
        dossier = self.dossiers.view(dossier_id)
        run = self.store.get_run(self.project, dossier['run_id'])
        require(dossier['current_input'] and not run.get('replay_only'), 'CURRENT_WRITABLE_DOSSIER_REQUIRED')
        directory = Path(directory)
        evidence = read_export(directory)
        manifest = parse_json((directory/'manifest.json').read_bytes())
        files = {name:(directory/name).read_bytes() for name in [*manifest['files'],'manifest.json']}
        # Freeze once, then check every byte actually stored (the reader above may have raced a writer).
        for name, entry in manifest['files'].items():
            require(digest(files[name]) == entry['sha256'] and len(files[name]) == entry['bytes'], 'B_EXPORT_CHANGED')
        require(encoded(evidence) == files['evidence.json'], 'B_EVIDENCE_BYTES_CHANGED')
        material, collection = (parse_json(files[name+'.json']) for name in ('material','collection'))
        require(material['case_id']==dossier['case_id'] and material['source_snapshot']['handoff_sha256']==dossier['handoff_ref']['sha256'], 'B_DOSSIER_SNAPSHOT_MISMATCH')
        require(files['cpu/handoff.json'] == self.port.read(dossier['handoff_ref']), 'B_CPU_HANDOFF_MISMATCH')
        for name, ref in dossier['artifact_refs'].items():
            require(files.get('cpu/'+name) == self.port.read(ref), 'B_CPU_ARTIFACT_MISMATCH')
        # This initial consumer only accepts the no-inference package; future prediction acceptance is explicit.
        require(not collection['runs'] and all(c['prediction']['status']=='not_run' for c in evidence['candidates']), 'PREDICTION_IMPORT_NOT_ENABLED')
        exchange = literature_exchange(dossier, material)
        reconciled = build_evidence(material, collection, exchange)
        fingerprint = digest(encoded({'version':VERSION,'dossier_ref':dossier['dossier_ref'],'manifest':digest(files['manifest.json'])}))
        with self.store.db() as db:
            old = db.execute('SELECT id FROM b_inputs WHERE project=? AND dossier_id=? AND fingerprint=?', (self.project,dossier_id,fingerprint)).fetchone()
        if old:
            return self.view(old[0]), False
        refs = {name:self.port.put_raw(raw, 'application/json' if name.endswith('.json') else 'application/octet-stream', 'computed') for name,raw in files.items()}
        value = {'format':'tpd-ab-review-input/0.1.0-draft','version':VERSION,'id':'binput-'+uuid.uuid4().hex,
            'project_id':self.project,'run_id':dossier['run_id'],'dossier_id':dossier_id,'dossier_ref':dossier['dossier_ref'],
            'snapshot':dossier['snapshot'],'created_at':now(),'material_digest':material['material_digest'],
            'source_evidence_digest':evidence['digest'],'files':refs,'exchange_ref':self.port.put_json(exchange),
            'reconciled_ref':self.port.put_json(reconciled), 'pending_inputs':pending_inputs(),
            'authority':{'human_approved':False,'dispatch_authorized':False,'formal_phase1_result':False}}
        schema_check(value,'a_review_handoff.schema.json')
        ref = self.port.put_json(value)
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            current = json.loads(db.execute('SELECT body FROM runs WHERE id=? AND project=?',(value['run_id'],self.project)).fetchone()[0])
            require(DossierService.snapshot(current)==value['snapshot'] and not current.get('replay_only'), 'B_INPUT_CHANGED')
            old = db.execute('SELECT id FROM b_inputs WHERE project=? AND dossier_id=? AND fingerprint=?',(self.project,dossier_id,fingerprint)).fetchone()
            if old:
                identifier = old[0]
            else:
                db.execute('INSERT INTO b_inputs VALUES(?,?,?,?,?,?)',(value['id'],self.project,value['run_id'],dossier_id,fingerprint,json.dumps(ref)))
                identifier = value['id']
        return self.view(identifier), identifier==value['id']

    def view(self, identifier):
        with self.store.db() as db:
            row = db.execute('SELECT ref FROM b_inputs WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if row is None:
            raise KeyError(identifier)
        ref = json.loads(row[0]); value = self.port.json(ref)
        schema_check(value,'a_review_handoff.schema.json')
        require(value['project_id']==self.project and value['id']==identifier, 'B_INPUT_CONTEXT')
        dossier = self.dossiers.view(value['dossier_id'])
        require(value['dossier_ref']==dossier['dossier_ref'], 'B_INPUT_DOSSIER_CHANGED')
        for item in value['files'].values():
            self.port.read(item)
        self.port.read(value['exchange_ref']); self.port.read(value['reconciled_ref'])
        value.update(current_input=dossier['current_input'], input_ref=ref)
        return value

    def list(self, dossier_id):
        self.dossiers.view(dossier_id)
        with self.store.db() as db:
            ids = [r[0] for r in db.execute('SELECT id FROM b_inputs WHERE project=? AND dossier_id=? ORDER BY rowid DESC',(self.project,dossier_id))]
        return [self.view(i) for i in ids]

    def packet(self, identifier, workflow=None):
        value = self.view(identifier)
        dossier = self.dossiers.view(value['dossier_id'])
        reconciled = self.port.json(value['reconciled_ref'])
        require(workflow is None or workflow['b_input_id']==identifier, 'REVIEW_WORKFLOW_INPUT_MISMATCH')
        packet = {'format':'tpd-human-review-preparation/0.1.0-draft','b_input_id':identifier,
            'input_ref':value['input_ref'],'current_input':value['current_input'],
            'subject':reconciled['shared_hypothesis']['hypothesis_id'],
            'meaning':'START attachment hypothesis with known C01/C02 references; not a GatePacket or approval.',
            'candidates':dossier['candidates'],'shared_hypothesis':reconciled['shared_hypothesis'],
            'literature_reconciliation':reconciled['literature_reconciliation'],
            'workflow':workflow,'pending_inputs':value['pending_inputs'],
            'ready_for_execution':False,'human_approved':False,'dispatch_authorized':False}
        schema_check(packet,'a_review_handoff.schema.json')
        return packet

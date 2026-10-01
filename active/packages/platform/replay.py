"""Closed artifact graph for internal, offline replay into an empty store only."""
import io
import json
import re
import zipfile
import tempfile
from pathlib import Path

from packages.contracts import encoded, parse_json, validate_payload
from packages.platform.dossiers import DossierService, digest, require
from packages.platform.store import Store
from packages.platform.workflows import WorkflowService
from packages.platform.analysis import AnalysisService

FORMAT = 'tpd-internal-replay/0.1.0-draft'
MAX_BYTES, MAX_FILES = 64_000_000, 2500
TABLES = {
    'runs':('id','project','request_key','request_hash','revision','state','body'),
    'analyses':('id','project','run_id','request_key','fingerprint','state','body'),
    'dossiers':('id','project','run_id','fingerprint','ref'),
    'b_inputs':('id','project','run_id','dossier_id','fingerprint','ref'),
    'workflows':('id','project','run_id','request_key','fingerprint','state','body'),
    'workflow_calls':('workflow_id','ordinal','body')}


def references(value):
    if isinstance(value,dict):
        if {'artifact_id','sha256','schema_id','version','media_type','provenance'} <= value.keys():
            if value['artifact_id'].startswith('a-'):
                validate_payload('ArtifactRef',value)
                require(re.fullmatch(r'a-[a-f0-9]{32}',value['artifact_id']) is not None,'REPLAY_ARTIFACT_ID')
                yield value
        else:
            for child in value.values(): yield from references(child)
    elif isinstance(value,list):
        for child in value: yield from references(child)


def row_payloads(records):
    for rows in records.values():
        for row in rows:
            for field in ('body','ref'):
                if field in row: yield parse_json(row[field].encode('utf-8'))


def export_bytes(store, project, run_id):
    AnalysisService(store,project)
    WorkflowService(store,project)
    run=store.get_run(project,run_id)
    require(run['state'] in {'completed','partial'}, 'REPLAY_FINISHED_COLLECTION_REQUIRED')
    records={}
    with store.db() as db:
        db.execute('BEGIN')
        for table in TABLES:
            if table=='runs': rows=db.execute('SELECT * FROM runs WHERE project=? AND id=?',(project,run_id)).fetchall()
            elif table=='workflow_calls':
                rows=db.execute('SELECT c.* FROM workflow_calls c JOIN workflows w ON w.id=c.workflow_id WHERE w.project=? AND w.run_id=? ORDER BY c.workflow_id,c.ordinal',(project,run_id)).fetchall()
            else: rows=db.execute('SELECT * FROM '+table+' WHERE project=? AND run_id=? ORDER BY rowid',(project,run_id)).fetchall()
            records[table]=[dict(r) for r in rows]
    require(not any(r['state'] in {'queued','running'} for table in ('analyses','workflows') for r in records[table]),'REPLAY_ACTIVE_JOB')
    queue=[r for value in row_payloads(records) for r in references(value)]
    blobs,refs={},{}
    while queue:
        ref=queue.pop(); identifier=ref['artifact_id']
        if identifier in refs:
            require(refs[identifier]==ref,'REPLAY_REF_CONFLICT');continue
        raw=store.read(project,ref);refs[identifier]=ref;blobs['blobs/'+identifier]=raw
        require(len(blobs)<=MAX_FILES and sum(map(len,blobs.values()))<=MAX_BYTES,'REPLAY_LIMIT')
        if ref['media_type']=='application/json': queue.extend(references(parse_json(raw)))
    manifest={'format':FORMAT,'project_id':project,'run_id':run_id,'records':records,
        'artifacts':[{'ref':refs[i],'bytes':len(blobs['blobs/'+i])} for i in sorted(refs)],
        'replay_only':True,'public_release_ready':False,'human_approval_transferred':False,
        'notice':'Internal research sources and saved answers. Hashes check integrity, not authenticity or scientific validity.'}
    manifest['digest']=digest(encoded(manifest))
    result=io.BytesIO()
    with zipfile.ZipFile(result,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('manifest.json',encoded(manifest))
        for name,raw in sorted(blobs.items()): archive.writestr(name,raw)
    return result.getvalue(),manifest


def inspect_archive(raw):
    require(len(raw)<=MAX_BYTES,'REPLAY_COMPRESSED_LIMIT')
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        entries=archive.infolist();names=[e.filename for e in entries]
        require(len(names)==len(set(names)) and len(names)<=MAX_FILES+1,'REPLAY_DUPLICATE_OR_LIMIT')
        require(sum(e.file_size for e in entries)<=MAX_BYTES,'REPLAY_SIZE_LIMIT')
        require(all(e.filename=='manifest.json' or re.fullmatch(r'blobs/a-[a-f0-9]{32}',e.filename) for e in entries),'REPLAY_UNSAFE_PATH')
        require(all((e.external_attr>>16)&0o170000 != 0o120000 for e in entries),'REPLAY_SYMLINK')
        manifest=parse_json(archive.read('manifest.json'))
        require(manifest['format']==FORMAT and manifest['replay_only'] is True and manifest['public_release_ready'] is False
                and manifest['human_approval_transferred'] is False,'REPLAY_AUTHORITY')
        require(digest(encoded({k:v for k,v in manifest.items() if k!='digest'}))==manifest['digest'],'REPLAY_MANIFEST_DIGEST')
        refs={};blobs={}
        for item in manifest['artifacts']:
            ref=item['ref'];validate_payload('ArtifactRef',ref);identifier=ref['artifact_id']
            require(re.fullmatch(r'a-[a-f0-9]{32}',identifier) is not None and identifier not in refs,'REPLAY_ARTIFACT_ID')
            data=archive.read('blobs/'+identifier)
            require(len(data)==item['bytes'] and digest(data)==ref['sha256'],'REPLAY_BLOB_HASH')
            refs[identifier]=ref;blobs[identifier]=data
        require(set(names)=={'manifest.json'}|{'blobs/'+i for i in refs},'REPLAY_UNINDEXED_FILES')
    records=manifest['records'];project=manifest['project_id'];run_id=manifest['run_id']
    require(set(records)==set(TABLES) and len(records['runs'])==1,'REPLAY_TABLES')
    require(records['runs'][0]['id']==run_id,'REPLAY_RUN')
    for table,rows in records.items():
        require(all(set(r)==set(TABLES[table]) for r in rows),'REPLAY_ROW_SHAPE')
        for row in rows:
            if table!='workflow_calls':
                require(row['project']==project and (table=='runs' or row['run_id']==run_id),'REPLAY_PROJECT_SCOPE')
            if 'state' in row:
                require(row['state'] not in {'queued','running','resolving','collecting'},'REPLAY_ACTIVE_JOB')
            if 'body' in row and table!='workflow_calls':
                body=parse_json(row['body'].encode())
                require(body['id']==row['id'] and body['project_id']==project and body['state']==row['state'],'REPLAY_ROW_CONTEXT')
                if table=='workflows':
                    require(body['human_approved'] is False and body['dispatch_authorized'] is False,'REPLAY_WORKFLOW_AUTHORITY')
    workflow_ids={r['id'] for r in records['workflows']}
    require(all(r['workflow_id'] in workflow_ids and type(r['ordinal']) is int and r['ordinal']>=0 for r in records['workflow_calls']),'REPLAY_CALL_SCOPE')
    queue=[r for value in row_payloads(records) for r in references(value)];seen=set()
    while queue:
        ref=queue.pop();identifier=ref['artifact_id']
        require(refs.get(identifier)==ref,'REPLAY_MISSING_OR_CONFLICTING_REF')
        if identifier in seen:continue
        seen.add(identifier)
        if ref['media_type']=='application/json':queue.extend(references(parse_json(blobs[identifier])))
    require(seen==set(refs),'REPLAY_UNREACHABLE_ARTIFACT')
    return manifest,blobs


def import_archive(raw, destination, project):
    manifest,blobs=inspect_archive(raw)
    require(manifest['project_id']==project,'REPLAY_WRONG_PROJECT')
    destination=Path(destination)
    require(not destination.exists(),'REPLAY_DESTINATION_EXISTS')
    destination=destination.resolve()
    destination.parent.mkdir(parents=True,exist_ok=True)
    # Staging stays in the destination's parent; failed validation never exposes a partial store.
    with tempfile.TemporaryDirectory(prefix='.tpd-replay-',dir=destination.parent) as scratch:
        staging=Path(scratch).resolve()
        require(staging.parent==destination.parent,'REPLAY_STAGING_BOUNDARY')
        _import_validated(manifest,blobs,staging,project)
        require(not destination.exists(),'REPLAY_DESTINATION_EXISTS')
        staging.rename(destination)
    return Store(destination),manifest


def _import_validated(manifest,blobs,destination,project):
    store=Store(destination); AnalysisService(store,project); WorkflowService(store,project)
    for identifier,data in blobs.items():
        with (store.root/'blobs'/identifier).open('xb') as stream:stream.write(data)
    with store.db() as db:
        db.execute('BEGIN IMMEDIATE')
        for item in manifest['artifacts']:
            ref=item['ref'];db.execute('INSERT INTO artifacts VALUES(?,?,?)',(ref['artifact_id'],project,json.dumps(ref)))
        for table,columns in TABLES.items():
            for source in manifest['records'][table]:
                row=dict(source)
                if table=='runs':
                    body=json.loads(row['body']);body['replay_only']=True;body['replay_source_digest']=manifest['digest'];row['body']=json.dumps(body)
                db.execute('INSERT INTO '+table+' ('+','.join(columns)+') VALUES ('+','.join('?' for _ in columns)+')',[row[c] for c in columns])
    # Validate the existing consumers in the new scope before reporting success.
    dossiers=DossierService(store,project)
    for row in manifest['records']['dossiers']:dossiers.view(row['id'])
    workflows=WorkflowService(store,project)
    for row in manifest['records']['b_inputs']:workflows.handoffs.view(row['id'])
    for row in manifest['records']['workflows']:workflows.view(row['id'])

"""Internal replay of saved structures and workbench jobs; never transfer credentials."""
import io
import json
import re
import tempfile
import zipfile
from pathlib import Path
from packages.contracts import encoded,parse_json,validate_payload
from packages.platform import replay as base
from packages.platform.dossiers import digest,require
from packages.platform.workbench import WorkbenchService
from packages.platform.report_reviews import ReportReviewService
from packages.platform.structure_views import ensure_table
from packages.platform.store import Store

FORMAT='tpd-workbench-replay/0.1.0'
MAX_BYTES=512_000_000
MAX_FILES=10000
TABLES={
 'saved_results':('id','project','dossier_id','fingerprint','ref'),
 'saved_expert_reviews':('id','project','result_id','fingerprint','ref'),
 'evidence_reports':('id','project','result_id','fingerprint','ref'),
 'report_review_requests':('id','project','report_id','fingerprint','ref'),
 'report_review_replacements':('project','parent_id','child_id'),
 'report_review_decisions':('id','project','request_id','revision','actor_id','request_key','payload_hash','ref'),
 'workbench_jobs':('id','project','run_id','result_id','request_key','fingerprint','state','body'),
 'workbench_calls':('job_id','ordinal','body'),
 'structure_sources':('project','result_id','pdb','ref')}

def initialize(store,project):
    WorkbenchService(store,project);ReportReviewService(store,project);ensure_table(store)

def pack(manifest,base_raw,blobs):
    output=io.BytesIO()
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('manifest.json',encoded(manifest));archive.writestr('base.zip',base_raw)
        for identifier,raw in sorted(blobs.items()):archive.writestr('blobs/'+identifier,raw)
    return output.getvalue()

def export_bytes(store,project,run_id):
    initialize(store,project)
    base_raw,base_manifest=base.export_bytes(store,project,run_id)
    _,base_blobs=base.inspect_archive(base_raw)
    records={}
    # Read a single database snapshot. All referenced content is immutable and hash checked.
    with store.db() as db:
        db.execute('BEGIN')
        dossiers={r['id'] for r in base_manifest['records']['dossiers']}
        def select(table,test):
            rows=db.execute('SELECT * FROM '+table+(" WHERE project=?" if table!='workbench_calls' else ''),
                            (project,) if table!='workbench_calls' else ()).fetchall()
            records[table]=[dict(r) for r in rows if test(r)]
            return {r.get('id') for r in records[table]}
        saved=select('saved_results',lambda r:r['dossier_id'] in dossiers)
        select('saved_expert_reviews',lambda r:r['result_id'] in saved)
        reports=select('evidence_reports',lambda r:r['result_id'] in saved)
        requests=select('report_review_requests',lambda r:r['report_id'] in reports)
        select('report_review_replacements',lambda r:r['parent_id'] in requests and r['child_id'] in requests)
        select('report_review_decisions',lambda r:r['request_id'] in requests)
        jobs=select('workbench_jobs',lambda r:r['result_id'] in saved)
        select('workbench_calls',lambda r:r['job_id'] in jobs)
        select('structure_sources',lambda r:r['result_id'] in saved)
    require(not any(r['state'] in {'queued','running'} for r in records['workbench_jobs']),'REPLAY_ACTIVE_JOB')
    refs={a['ref']['artifact_id']:a['ref'] for a in base_manifest['artifacts']};blobs={};total=sum(map(len,base_blobs.values()))
    queue=[r for value in base.row_payloads(records) for r in base.references(value)]
    while queue:
        ref=queue.pop();identifier=ref['artifact_id']
        if identifier in refs:require(refs[identifier]==ref,'REPLAY_REF_CONFLICT');continue
        raw=store.read(project,ref);refs[identifier]=ref;blobs[identifier]=raw;total+=len(raw)
        require(total<=MAX_BYTES and len(refs)<=MAX_FILES,'REPLAY_LIMIT')
        if ref['media_type']=='application/json':queue.extend(base.references(parse_json(raw)))
    manifest={'format':FORMAT,'project_id':project,'run_id':run_id,'base_sha256':digest(base_raw),'records':records,
              'artifacts':[{'ref':refs[i],'bytes':len(b)} for i,b in sorted(blobs.items())],
              'replay_only':True,'human_approval_transferred':False,'public_release_ready':False,
              'notice':'Internal source materials and received opinions; review credentials and sessions excluded.'}
    manifest['digest']=digest(encoded(manifest))
    return pack(manifest,base_raw,blobs),manifest

def inspect_archive(raw):
    require(len(raw)<=MAX_BYTES,'REPLAY_COMPRESSED_LIMIT')
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        entries=archive.infolist();names=[e.filename for e in entries]
        require(len(names)==len(set(names)) and len(names)<=MAX_FILES+2,'REPLAY_DUPLICATE_OR_LIMIT')
        require(sum(e.file_size for e in entries)<=MAX_BYTES,'REPLAY_SIZE_LIMIT')
        require(all(e.filename in {'manifest.json','base.zip'} or re.fullmatch(r'blobs/a-[a-f0-9]{32}',e.filename) for e in entries),'REPLAY_UNSAFE_PATH')
        require(all((e.external_attr>>16)&0o170000!=0o120000 for e in entries),'REPLAY_SYMLINK')
        m=parse_json(archive.read('manifest.json'));base_raw=archive.read('base.zip')
        require(m['format']==FORMAT and m['replay_only'] is True and m['human_approval_transferred'] is False
                and m['public_release_ready'] is False,'REPLAY_AUTHORITY')
        require(digest(encoded({k:v for k,v in m.items() if k!='digest'}))==m['digest'],'REPLAY_MANIFEST_DIGEST')
        require(digest(base_raw)==m['base_sha256'],'REPLAY_BASE_HASH')
        bm,bb=base.inspect_archive(base_raw)
        require(m['project_id']==bm['project_id'] and m['run_id']==bm['run_id'],'REPLAY_BASE_CONTEXT')
        refs={x['ref']['artifact_id']:x['ref'] for x in bm['artifacts']};blobs=dict(bb);extra={}
        for item in m['artifacts']:
            ref=item['ref'];validate_payload('ArtifactRef',ref);identifier=ref['artifact_id']
            require(re.fullmatch(r'a-[a-f0-9]{32}',identifier) is not None and identifier not in refs,'REPLAY_ARTIFACT_ID')
            data=archive.read('blobs/'+identifier)
            require(len(data)==item['bytes'] and digest(data)==ref['sha256'],'REPLAY_BLOB_HASH')
            refs[identifier]=ref;blobs[identifier]=data;extra[identifier]=data
        require(set(names)=={'manifest.json','base.zip'}|{'blobs/'+i for i in extra},'REPLAY_UNINDEXED_FILES')
        require(sum(map(len,blobs.values()))<=MAX_BYTES,'REPLAY_SIZE_LIMIT')
    records=m['records'];project=m['project_id'];run_id=m['run_id']
    require(set(records)==set(TABLES),'REPLAY_TABLES')
    for table,rows in records.items():
        require(all(set(r)==set(TABLES[table]) for r in rows),'REPLAY_ROW_SHAPE')
        for row in rows:
            if table!='workbench_calls':require(row['project']==project,'REPLAY_PROJECT_SCOPE')
    ids=lambda table:{r['id'] for r in records[table]}
    dossier_ids={r['id'] for r in bm['records']['dossiers']}
    for row in records['saved_results']:require(row['dossier_id'] in dossier_ids,'REPLAY_DOSSIER_SCOPE')
    for table in ('saved_expert_reviews','evidence_reports','workbench_jobs','structure_sources'):
        for row in records[table]:require(row['result_id'] in ids('saved_results'),'REPLAY_RESULT_SCOPE')
    for row in records['report_review_requests']:require(row['report_id'] in ids('evidence_reports'),'REPLAY_REPORT_SCOPE')
    for row in records['report_review_decisions']:require(row['request_id'] in ids('report_review_requests'),'REPLAY_REVIEW_SCOPE')
    for row in records['report_review_replacements']:
        require(row['parent_id'] in ids('report_review_requests') and row['child_id'] in ids('report_review_requests'),'REPLAY_REVIEW_SCOPE')
    for row in records['workbench_jobs']:
        body=parse_json(row['body'].encode())
        require(row['state'] not in {'running','queued'} and row['run_id']==run_id,'REPLAY_ACTIVE_JOB')
        require(body['id']==row['id'] and body['project_id']==project and body['run_id']==run_id
                and body['result_id']==row['result_id'] and body['state']==row['state'],'REPLAY_JOB_CONTEXT')
    for row in records['workbench_calls']:
        require(row['job_id'] in ids('workbench_jobs') and type(row['ordinal']) is int and row['ordinal']>0,'REPLAY_CALL_SCOPE')
    queue=[r for value in base.row_payloads(records) for r in base.references(value)];seen=set(bb)
    while queue:
        ref=queue.pop();identifier=ref['artifact_id'];require(refs.get(identifier)==ref,'REPLAY_MISSING_OR_CONFLICTING_REF')
        if identifier in seen:continue
        seen.add(identifier)
        if ref['media_type']=='application/json':queue.extend(base.references(parse_json(blobs[identifier])))
    require(seen==set(refs),'REPLAY_UNREACHABLE_ARTIFACT')
    return m,bm,bb,extra

def import_archive(raw,destination,project):
    m,bm,bb,extra=inspect_archive(raw);destination=Path(destination).resolve()
    require(m['project_id']==project,'REPLAY_WRONG_PROJECT');require(not destination.exists(),'REPLAY_DESTINATION_EXISTS')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.tpd-workbench-',dir=destination.parent) as scratch:
        staging=Path(scratch).resolve();require(staging.parent==destination.parent,'REPLAY_STAGING_BOUNDARY')
        base._import_validated(bm,bb,staging,project);store=Store(staging);initialize(store,project)
        for identifier,data in extra.items():
            with (store.root/'blobs'/identifier).open('xb') as f:f.write(data)
        with store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            for item in m['artifacts']:
                ref=item['ref'];db.execute('INSERT INTO artifacts VALUES(?,?,?)',(ref['artifact_id'],project,json.dumps(ref)))
            for table,columns in TABLES.items():
                for row in m['records'][table]:db.execute('INSERT INTO '+table+' ('+','.join(columns)+') VALUES ('+','.join('?' for _ in columns)+')',[row[c] for c in columns])
        saved=WorkbenchService(store,project).saved
        for row in m['records']['saved_results']:saved.verify(row['id'])
        require(not destination.exists(),'REPLAY_DESTINATION_EXISTS');staging.rename(destination)
    return Store(destination),m

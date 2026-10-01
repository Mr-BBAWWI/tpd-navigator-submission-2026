"""Version-bound human review of saved reports, distinct from G1 compute approval."""
import json
import uuid
from pathlib import Path

from packages.contracts import encoded, now
from packages.platform.dossiers import digest, require
from packages.platform.evidence_reports import EvidenceReportService, implementation as report_implementation
from packages.platform.review_identity import ReviewIdentityService
from packages.platform.saved_results import AUTHORITY
from packages.science.evidence_common import schema_check

ROOT = Path(__file__).resolve().parents[2]
SCOPES = {
    'chemical_state': '원자 역할·프로톤화·수소와 접촉 결과의 해석 범위',
    'structural_comparison': '기준 후보의 구조·부위별 비교와 반복 샘플 해석',
    'literature_conditions': '문헌 값·실험 조건·결측 표시와 출처 연결',
    'synthesis_evidence': '합성 근거의 충족 여부와 미평가 표시',
}
ACKNOWLEDGEMENTS = {
    'scope_only': '판정은 선택한 검토 범위와 고정된 보고서에만 적용됩니다.',
    'no_efficacy': '범위 내 수용은 분해 효능·신규 후보의 타당성 입증이 아닙니다.',
    'no_execution': '이 판정은 G1 승인·GPU 실행·공개 배포를 허용하지 않습니다.',
}


def implementation():
    names = ('packages/platform/report_reviews.py','packages/platform/review_identity.py',
             'contracts/drafts/a_report_review.schema.json')
    return {'version':'report-human-review/20260927.1', 'sha256':{p:digest((ROOT/p).read_bytes()) for p in names}}


def validate(value, kind):
    # One schema owns the request/decision union; the operation is checked explicitly.
    schema_check(value,'a_report_review.schema.json')
    require(value['kind']==kind,'REVIEW_RECORD_KIND')


class ReportReviewService:
    def __init__(self, store, project):
        self.store,self.project,self.port = store,project,store.scope(project)
        self.reports = EvidenceReportService(store,project)
        self.identities = ReviewIdentityService(store,project)
        with store.db() as db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS report_review_requests(
                    id TEXT PRIMARY KEY,project TEXT NOT NULL,report_id TEXT NOT NULL,
                    fingerprint TEXT NOT NULL,ref TEXT NOT NULL,UNIQUE(project,report_id,fingerprint));
                CREATE TABLE IF NOT EXISTS report_review_replacements(
                    project TEXT NOT NULL,parent_id TEXT NOT NULL,child_id TEXT NOT NULL,
                    PRIMARY KEY(project,parent_id));
                CREATE TABLE IF NOT EXISTS report_review_decisions(
                    id TEXT PRIMARY KEY,project TEXT NOT NULL,request_id TEXT NOT NULL,revision INTEGER NOT NULL,
                    actor_id TEXT NOT NULL,request_key TEXT NOT NULL,payload_hash TEXT NOT NULL,ref TEXT NOT NULL,
                    UNIQUE(project,request_id,revision),UNIQUE(project,request_id,actor_id,request_key));
            ''')

    def _request(self, identifier):
        with self.store.db() as db:
            row=db.execute('SELECT ref FROM report_review_requests WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if row is None: raise KeyError(identifier)
        ref=json.loads(row[0]);value=self.port.json(ref);validate(value,'request')
        require(value['id']==identifier and value['project_id']==self.project,'REVIEW_CONTEXT')
        return value,ref

    def _assert_report_current(self, db, report):
        require(self.reports._state(db,report['result_id'],report['dossier_id'],report['run_id'])==report['binding']['state']
                and report_implementation()==report['binding']['implementation']
                and not report['binding']['state']['replay_only'],'REVIEW_REPORT_CHANGED')
        row=db.execute('SELECT ref FROM evidence_reports WHERE project=? AND id=?',(self.project,report['id'])).fetchone()
        require(row is not None and json.loads(row[0])==report['record_ref'],'REVIEW_REPORT_CHANGED')

    def create(self, report_id, expected_report_sha256, scopes, reviewer_id=None, supersedes=None, session_token=None):
        report=self.reports.view(report_id)
        require(report['freshness']['current'],'REVIEW_REPORT_CHANGED')
        require(report['record_ref']['sha256']==expected_report_sha256,'REVIEW_REPORT_CHANGED')
        require(scopes and len(set(scopes))==len(scopes) and set(scopes)<=set(SCOPES),'REVIEW_UNKNOWN_SCOPE')
        scope_ids=sorted(scopes);prior,prior_ref=self._request(supersedes) if supersedes else (None,None)
        if prior:
            require(prior['run_id']==report['run_id'] and sorted(prior['scopes'])==scope_ids,'REVIEW_REPLACEMENT_SCOPE')
        policy=implementation()
        fingerprint=digest(encoded({'report_ref':report['record_ref'],'scopes':scope_ids,'reviewer_id':reviewer_id,
                                    'supersedes_ref':prior_ref,'policy':policy}))
        paths=[];committed=False
        try:
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE');self._assert_report_current(db,report)
                require(implementation()==policy,'REVIEW_POLICY_CHANGED')
                actor=self.identities.actor(db,reviewer_id) if reviewer_id else None
                require(reviewer_id is None or actor and actor['active'],'REVIEWER_UNAVAILABLE')
                creator=self.identities.authenticate(session_token,db) if session_token is not None else None
                if creator:
                    require(actor and creator['id']==actor['id'],'REVIEW_ASSIGNED_ACTOR_REQUIRED')
                    require(not prior or prior['assigned_reviewer'] and prior['assigned_reviewer']['id']==creator['id'],'REVIEW_REPLACEMENT_ACTOR')
                old=db.execute('SELECT id FROM report_review_requests WHERE project=? AND report_id=? AND fingerprint=?',
                               (self.project,report_id,fingerprint)).fetchone()
                if old: identifier=old[0];fresh=False
                else:
                    if prior:
                        row=db.execute('SELECT child_id FROM report_review_replacements WHERE project=? AND parent_id=?',(self.project,prior['id'])).fetchone()
                        require(row is None,'REVIEW_ALREADY_REPLACED')
                    identifier='hreview-'+uuid.uuid4().hex
                    value={'kind':'request','id':identifier,'project_id':self.project,'report_id':report_id,
                           'run_id':report['run_id'],'result_id':report['result_id'],'created_at':now(),
                           'report_ref':report['record_ref'],'report_input_digest':report['input_digest'],
                           'scopes':{k:SCOPES[k] for k in scope_ids},'acknowledgements':ACKNOWLEDGEMENTS.copy(),
                           'assigned_reviewer':{'id':actor['id'],'name':actor['name']} if actor else None,
                           'supersedes_ref':prior_ref,'policy':policy,'authority':AUTHORITY.copy()}
                    value['created_by']=creator
                    validate(value,'request');ref=self.reports.saved._write(db,paths,encoded(value))
                    db.execute('INSERT INTO report_review_requests VALUES(?,?,?,?,?)',(identifier,self.project,report_id,fingerprint,json.dumps(ref)))
                    if prior: db.execute('INSERT INTO report_review_replacements VALUES(?,?,?)',(self.project,prior['id'],identifier))
                    fresh=True
            committed=True
            return self.view(identifier),fresh
        except BaseException:
            # Only uncommitted blobs are removed. A later view failure cannot erase a committed request.
            if not committed:
                for p in paths:p.unlink(missing_ok=True)
            raise

    def _history(self, db, identifier, request_ref):
        rows=db.execute('SELECT revision,ref FROM report_review_decisions WHERE project=? AND request_id=? ORDER BY revision',(self.project,identifier)).fetchall()
        values=[];previous=None
        for sequence,row in enumerate(rows,1):
            ref=json.loads(row['ref']);value=self.port.json(ref);validate(value,'decision')
            require(row['revision']==value['revision']==sequence and value['project_id']==self.project and
                    value['request_ref']==request_ref and value['request_id']==identifier and value['previous_ref']==previous,'REVIEW_HISTORY_CHAIN')
            value['record_ref']=ref;values.append(value);previous=ref
        return values

    def decide(self, identifier, payload, session_token):
        validate({'kind':'submission',**payload},'submission')
        request,request_ref=self._request(identifier)
        require(payload['expected_request_sha256']==request_ref['sha256'],'REVIEW_REQUEST_CHANGED')
        withdrawing=payload['action']=='withdraw'
        report=None if withdrawing else self.reports.view(request['report_id'])
        if report:
            require(report['freshness']['current'] and report['record_ref']==request['report_ref'],'REVIEW_REPORT_CHANGED')
        paths=[];receipt=None
        try:
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE')
                actor=self.identities.authenticate(session_token,db)
                require(request['assigned_reviewer'] is not None and actor['id']==request['assigned_reviewer']['id'],'REVIEW_ASSIGNED_ACTOR_REQUIRED')
                payload_hash=digest(encoded(payload))
                old=db.execute('SELECT ref,payload_hash FROM report_review_decisions WHERE project=? AND request_id=? AND actor_id=? AND request_key=?',
                               (self.project,identifier,actor['id'],payload['idempotency_key'])).fetchone()
                if old:
                    require(old['payload_hash']==payload_hash,'REVIEW_IDEMPOTENCY_CONFLICT')
                    receipt=self.port.json(json.loads(old['ref']));fresh=False
                else:
                    history=self._history(db,identifier,request_ref)
                    require(len(history)==payload['expected_revision'],'REVIEW_DECISION_CONFLICT')
                    if withdrawing:
                        require(history and history[-1]['action']=='record' and not payload['judgments'] and
                                not payload['acknowledgement_ids'] and payload['withdrawal_reason'].strip(),'REVIEW_WITHDRAWAL_INVALID')
                    else:
                        self._assert_report_current(db,report)
                        require(implementation()==request['policy'],'REVIEW_POLICY_CHANGED')
                        require(not db.execute('SELECT 1 FROM report_review_replacements WHERE project=? AND parent_id=?',(self.project,identifier)).fetchone(),'REVIEW_SUPERSEDED')
                        judgments=payload['judgments'];ids=[j['scope_id'] for j in judgments]
                        require(len(set(ids))==len(ids) and set(ids)==set(request['scopes']),'REVIEW_SCOPE_COVERAGE')
                        require(all(j['reason'].strip() for j in judgments) and payload['withdrawal_reason']=='','REVIEW_REASON_REQUIRED')
                        acks=payload['acknowledgement_ids']
                        require(len(set(acks))==len(acks) and set(acks)<=set(request['acknowledgements']),'REVIEW_ACK_INVALID')
                        if any(j['verdict']=='accept_scope' for j in judgments):
                            require(set(acks)==set(request['acknowledgements']),'REVIEW_ACK_REQUIRED')
                    receipt={'kind':'decision','id':'hdecision-'+uuid.uuid4().hex,'project_id':self.project,
                             'request_id':identifier,'request_ref':request_ref,'revision':len(history)+1,
                             'previous_ref':history[-1]['record_ref'] if history else None,'actor':actor,'recorded_at':now(),
                             'action':payload['action'],'judgments':payload['judgments'],'acknowledgement_ids':payload['acknowledgement_ids'],
                             'withdrawal_reason':payload['withdrawal_reason'],'authority':AUTHORITY.copy()}
                    validate(receipt,'decision');ref=self.reports.saved._write(db,paths,encoded(receipt))
                    db.execute('INSERT INTO report_review_decisions VALUES(?,?,?,?,?,?,?,?)',
                        (receipt['id'],self.project,identifier,receipt['revision'],actor['id'],payload['idempotency_key'],payload_hash,json.dumps(ref)))
                    fresh=True
            return {'decision_id':receipt['id'],'revision':receipt['revision'],'created':fresh}
        except BaseException:
            for p in paths:p.unlink(missing_ok=True)
            raise

    def view(self, identifier):
        request,ref=self._request(identifier)
        report=self.reports.view(request['report_id'])
        require(report['record_ref']==request['report_ref'] and report['input_digest']==request['report_input_digest'],'REVIEW_REPORT_BINDING')
        with self.store.db() as db:
            db.execute('BEGIN')
            history=self._history(db,identifier,ref)
            replacement=db.execute('SELECT child_id FROM report_review_replacements WHERE project=? AND parent_id=?',(self.project,identifier)).fetchone()
            actor=self.identities.actor(db,request['assigned_reviewer']['id']) if request['assigned_reviewer'] else None
            current=self.reports._state(db,report['result_id'],report['dossier_id'],report['run_id'])
        reasons=list(report['freshness']['reasons'])
        if current!=report['binding']['state'] and not reasons: reasons.append('input_changed_during_read')
        if implementation()!=request['policy']:reasons.append('review_policy_changed')
        if replacement:reasons.append('request_superseded')
        if request['assigned_reviewer'] and (not actor or not actor['active']):reasons.append('reviewer_inactive')
        last=history[-1] if history else None
        for entry in history:
            require(request['assigned_reviewer'] is not None and entry['actor']['id']==request['assigned_reviewer']['id'],'REVIEW_HISTORY_ACTOR')
        if reasons:status='requires_re_review'
        elif last and last['action']=='withdraw':status='withdrawn'
        elif not actor:status='unassigned'
        elif not last:status='pending'
        elif any(j['verdict']=='request_changes' for j in last['judgments']):status='changes_requested'
        elif any(j['verdict']=='defer' for j in last['judgments']):status='deferred'
        else:status='accepted_for_scope'
        active=not reasons and last is not None and last['action']=='record'
        request.update(record_ref=ref,history=history,revision=len(history),status=status,recheck_reasons=reasons,
                       superseded_by=replacement[0] if replacement else None,
                       effective_scope_acceptances={k:bool(active and any(j['scope_id']==k and j['verdict']=='accept_scope' for j in last['judgments'])) for k in request['scopes']})
        return request

    def list(self, report_id):
        self.reports.view(report_id)
        with self.store.db() as db:
            ids=[r[0] for r in db.execute('SELECT id FROM report_review_requests WHERE project=? AND report_id=? ORDER BY rowid DESC',(self.project,report_id))]
        return [self.view(i) for i in ids]

    def predecessors(self, report_id):
        report=self.reports.view(report_id)
        with self.store.db() as db:
            rows=db.execute('''SELECT r.id,r.ref FROM report_review_requests r
                LEFT JOIN report_review_replacements x ON r.project=x.project AND r.id=x.parent_id
                WHERE r.project=? AND r.report_id<>? AND x.child_id IS NULL ORDER BY r.rowid DESC''',
                (self.project,report_id)).fetchall()
        return [self.view(r['id']) for r in rows if self.port.json(json.loads(r['ref']))['run_id']==report['run_id']]

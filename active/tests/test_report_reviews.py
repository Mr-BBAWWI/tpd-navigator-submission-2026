"""Synthetic decisions only. No actual reviewer opinion or authorization is created."""
import copy
import json
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import unittest

from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.contracts import ContractError
from packages.platform.report_reviews import ReportReviewService, SCOPES, ACKNOWLEDGEMENTS, implementation
from packages.platform.review_identity import ReviewAuthError, ReviewIdentityService, token_hash
from packages.platform.saved_results import AUTHORITY
import test_evidence_reports as fixtures


class ReportReviewTests(unittest.TestCase):
    def setUp(self):
        fixtures.EvidenceReportTests.setUp(self)
        self.report,_=self.svc.create('s',self.svc.prepare('s')['input_digest'])
        self.reviews=ReportReviewService(self.store,PROJECT);self.reviews.reports=self.svc
        self.auth=self.reviews.identities
        self.key=self.auth.register('test-reviewer','SYNTHETIC REVIEWER')
        self.token=self.auth.login('test-reviewer',self.key)
        self.request,_=self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],list(SCOPES),'test-reviewer')

    def submission(self,revision=0,key='synthetic-key-001',verdict='accept_scope'):
        return {'expected_request_sha256':self.request['record_ref']['sha256'],'expected_revision':revision,
                'idempotency_key':key,'action':'record','judgments':[{'scope_id':k,'verdict':verdict,'reason':'Synthetic test rationale; not a real review.'} for k in SCOPES],
                'acknowledgement_ids':list(ACKNOWLEDGEMENTS),'withdrawal_reason':''}

    def withdraw(self,revision=1):
        p=self.submission(revision,'synthetic-withdraw');p.update(action='withdraw',judgments=[],acknowledgement_ids=[],withdrawal_reason='Synthetic withdrawal')
        return self.reviews.decide(self.request['id'],p,self.token)

    def test_identity_login_expiry_logout_and_scope(self):
        for identifier,key in [('test-reviewer','wrong-key'),('unknown',self.key)]:
            with self.assertRaises(ReviewAuthError):self.auth.login(identifier,key)
        other=ReviewIdentityService(self.store,'other')
        with self.assertRaises(ReviewAuthError):other.session(self.token)
        with self.store.db() as db:
            row=db.execute('SELECT key_hash FROM review_actors WHERE id=?',('test-reviewer',)).fetchone()
            self.assertNotEqual(row[0],self.key)
            db.execute('UPDATE review_sessions SET expires_at=0 WHERE token_hash=?',(token_hash(self.token),))
        with self.assertRaises(ReviewAuthError):self.auth.session(self.token)
        token=self.auth.login('test-reviewer',self.key);self.auth.logout(token)
        with self.assertRaises(ReviewAuthError):self.auth.session(token)

    def test_acceptance_is_scoped_and_server_attributed(self):
        result=self.reviews.decide(self.request['id'],self.submission(),self.token)
        view=self.reviews.view(self.request['id']);decision=view['history'][0]
        self.assertEqual(result['revision'],1);self.assertEqual(view['status'],'accepted_for_scope')
        self.assertTrue(all(view['effective_scope_acceptances'].values()))
        self.assertEqual(decision['actor']['name'],'SYNTHETIC REVIEWER');self.assertTrue(decision['recorded_at'])
        self.assertEqual(view['authority'],AUTHORITY);self.assertEqual(decision['authority'],AUTHORITY)
        self.assertNotIn(self.key,json.dumps(view));self.assertNotIn(self.token,json.dumps(view))

    def test_scope_coverage_acknowledgements_and_actor_forgery(self):
        mutations=[('judgments',self.submission()['judgments'][:-1]),('acknowledgement_ids',[]),
                   ('actor',{'id':'forged'}),('expected_request_sha256','a'*64)]
        for field,value in mutations:
            p=self.submission();p[field]=value
            with self.assertRaises((ValueError,ContractError)):self.reviews.decide(self.request['id'],p,self.token)
        self.assertEqual(self.reviews.view(self.request['id'])['revision'],0)

    def test_other_reviewer_or_unassigned_request_cannot_decide(self):
        key=self.auth.register('test-other','SYNTHETIC OTHER');token=self.auth.login('test-other',key)
        with self.assertRaisesRegex(ContractError,'ASSIGNED_ACTOR'):self.reviews.decide(self.request['id'],self.submission(),token)
        row,_=self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],list(SCOPES))
        p=self.submission();p['expected_request_sha256']=row['record_ref']['sha256']
        with self.assertRaisesRegex(ContractError,'ASSIGNED_ACTOR'):self.reviews.decide(row['id'],p,self.token)

    def test_duplicate_and_conflicting_retries(self):
        p=self.submission();first=self.reviews.decide(self.request['id'],p,self.token)
        count=len(list((self.store.root/'blobs').iterdir()))
        retry=self.reviews.decide(self.request['id'],p,self.token)
        self.assertEqual(first['decision_id'],retry['decision_id']);self.assertFalse(retry['created'])
        self.assertEqual(len(list((self.store.root/'blobs').iterdir())),count)
        p['judgments'][0]['reason']='different'
        with self.assertRaisesRegex(ContractError,'IDEMPOTENCY'):self.reviews.decide(self.request['id'],p,self.token)

    def test_concurrent_decisions_have_one_winner(self):
        def submit(index):
            try:return self.reviews.decide(self.request['id'],self.submission(key='concurrent-'+str(index)),self.token)
            except ContractError as error:return str(error)
        with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(submit,range(2)))
        self.assertEqual(sum(isinstance(r,dict) for r in results),1)
        self.assertIn('REVIEW_DECISION_CONFLICT',results)
        self.assertEqual(self.reviews.view(self.request['id'])['revision'],1)

    def test_revision_withdrawal_and_resubmission_preserve_chain(self):
        self.reviews.decide(self.request['id'],self.submission(),self.token)
        self.withdraw()
        view=self.reviews.view(self.request['id']);self.assertEqual(view['status'],'withdrawn')
        self.assertFalse(any(view['effective_scope_acceptances'].values()))
        self.assertEqual(view['history'][1]['previous_ref'],view['history'][0]['record_ref'])
        self.reviews.decide(self.request['id'],self.submission(2,'revision-three','request_changes'),self.token)
        view=self.reviews.view(self.request['id']);self.assertEqual(view['status'],'changes_requested')
        self.assertEqual(view['revision'],3)
        retry=self.withdraw(1);self.assertFalse(retry['created']);self.assertEqual(retry['revision'],2)
        p=self.submission(1,'stale-new-request')
        with self.assertRaisesRegex(ContractError,'DECISION_CONFLICT'):self.reviews.decide(self.request['id'],p,self.token)

    def test_changed_input_clears_effective_acceptance_but_allows_withdrawal(self):
        self.reviews.decide(self.request['id'],self.submission(),self.token)
        run={**self.run,'revision':2}
        with self.store.db() as db:db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps(run),'r'))
        view=self.reviews.view(self.request['id'])
        self.assertEqual(view['status'],'requires_re_review');self.assertFalse(any(view['effective_scope_acceptances'].values()))
        with self.assertRaisesRegex(ContractError,'REPORT_CHANGED'):self.reviews.decide(self.request['id'],self.submission(1,'old-input'),self.token)
        self.assertEqual(self.withdraw()['revision'],2)

    def test_policy_change_and_disable_invalidate_without_erasing(self):
        self.reviews.decide(self.request['id'],self.submission(),self.token)
        identity=implementation();identity['version']='synthetic changed policy'
        with patch('packages.platform.report_reviews.implementation',return_value=identity):
            self.assertIn('review_policy_changed',self.reviews.view(self.request['id'])['recheck_reasons'])
            with self.assertRaisesRegex(ContractError,'POLICY_CHANGED'):self.reviews.decide(self.request['id'],self.submission(1,'changed-policy'),self.token)
        self.auth.disable('test-reviewer')
        self.assertIn('reviewer_inactive',self.reviews.view(self.request['id'])['recheck_reasons'])
        with self.assertRaises(ReviewAuthError):self.auth.session(self.token)

    def test_replacement_scope_and_foreign_actor_protected(self):
        with self.assertRaisesRegex(ContractError,'REPLACEMENT_SCOPE'):
            self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],['chemical_state'],'test-reviewer',self.request['id'])
        key=self.auth.register('test-other','SYNTHETIC OTHER');token=self.auth.login('test-other',key)
        with self.assertRaisesRegex(ContractError,'REPLACEMENT_ACTOR'):
            self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],list(SCOPES),'test-other',self.request['id'],token)
        new,_=self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],list(SCOPES),'test-reviewer',self.request['id'],self.token)
        old=self.reviews.view(self.request['id']);self.assertEqual(old['superseded_by'],new['id'])
        self.assertEqual(new['supersedes_ref'],self.request['record_ref'])
        with self.assertRaisesRegex(ContractError,'SUPERSEDED'):self.reviews.decide(self.request['id'],self.submission(),self.token)
        retry,fresh=self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],list(SCOPES),'test-reviewer',self.request['id'],self.token)
        self.assertFalse(fresh);self.assertEqual(retry['id'],new['id'])

    def test_new_report_can_continue_previous_request(self):
        self.reviews.decide(self.request['id'],self.submission(),self.token)
        self.run['revision']=2;self.data['snapshot']=self.run.copy()
        with self.store.db() as db:db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps(self.run),'r'))
        report,_=self.svc.create('s',self.svc.prepare('s')['input_digest'])
        self.assertEqual([r['id'] for r in self.reviews.predecessors(report['id'])],[self.request['id']])
        successor,_=self.reviews.create(report['id'],report['record_ref']['sha256'],list(SCOPES),'test-reviewer',self.request['id'],self.token)
        self.assertEqual(successor['status'],'pending');self.assertEqual(successor['revision'],0)
        self.assertEqual(self.reviews.predecessors(report['id']),[])

    def test_race_during_decision_checks_state_inside_transaction(self):
        original=self.svc.view
        def racing_view(identifier):
            value=original(identifier)
            with self.store.db() as db:db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps({**self.run,'revision':2}),'r'))
            return value
        with patch.object(self.svc,'view',side_effect=racing_view):
            with self.assertRaisesRegex(ContractError,'REPORT_CHANGED'):self.reviews.decide(self.request['id'],self.submission(),self.token)
        with self.store.db() as db:self.assertEqual(db.execute('SELECT COUNT(*) FROM report_review_decisions').fetchone()[0],0)

    def test_failed_storage_rolls_back_new_blobs(self):
        before=set((self.store.root/'blobs').iterdir());writer=self.svc.saved._write
        def failure(*args,**kwargs):writer(*args,**kwargs);raise OSError('synthetic failure')
        with patch.object(self.svc.saved,'_write',side_effect=failure):
            with self.assertRaises(OSError):self.reviews.decide(self.request['id'],self.submission(),self.token)
        self.assertEqual(set((self.store.root/'blobs').iterdir()),before)
        self.assertEqual(self.reviews.view(self.request['id'])['revision'],0)

    def test_request_failure_cleans_only_uncommitted_files(self):
        before=set((self.store.root/'blobs').iterdir());writer=self.svc.saved._write
        def failure(*args,**kwargs):writer(*args,**kwargs);raise OSError('synthetic failure')
        with patch.object(self.svc.saved,'_write',side_effect=failure):
            with self.assertRaises(OSError):self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],['chemical_state'],'test-reviewer')
        self.assertEqual(set((self.store.root/'blobs').iterdir()),before)
        with patch.object(self.reviews,'view',side_effect=OSError('post-commit view failure')):
            with self.assertRaises(OSError):self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],['chemical_state'],'test-reviewer')
        row,fresh=self.reviews.create(self.report['id'],self.report['record_ref']['sha256'],['chemical_state'],'test-reviewer')
        self.assertFalse(fresh);self.assertEqual(row['status'],'pending')

    def test_history_corruption_is_not_accepted(self):
        self.reviews.decide(self.request['id'],self.submission(),self.token)
        decision=self.reviews.view(self.request['id'])['history'][0]
        path=self.store.root/'blobs'/decision['record_ref']['artifact_id'];path.write_bytes(b'corrupt synthetic record')
        with self.assertRaises(ContractError):self.reviews.view(self.request['id'])

    def test_api_cookie_csrf_and_review_only_boundaries(self):
        app=create_app(self.store.root,enable_worker=False,review_only=True)
        url='/api/report-reviews/'+self.request['id']+'/decisions'
        with TestClient(app) as client:
            self.assertEqual(client.post(url,json=self.submission(),headers={'X-TPD-Local':'1'}).status_code,401)
            self.assertEqual(client.post('/api/review-auth/login',json={'reviewer_id':'test-reviewer','access_key':self.key}).status_code,403)
            response=client.post('/api/review-auth/login',json={'reviewer_id':'test-reviewer','access_key':self.key},headers={'X-TPD-Local':'1'})
            self.assertEqual(response.status_code,200)
            self.assertIn('HttpOnly',response.headers['set-cookie']);self.assertIn('SameSite=strict',response.headers['set-cookie'])
            self.assertNotIn(self.key,response.text)
            self.assertEqual(client.post(url,json=self.submission(),headers={'X-TPD-Local':'1','Origin':'http://other.example'}).status_code,403)
            self.assertEqual(client.post(url,json={**self.submission(),'actor_id':'forged'},headers={'X-TPD-Local':'1'}).status_code,422)
            self.assertEqual(client.post('/api/workflows',json={},headers={'X-TPD-Local':'1'}).status_code,403)
            self.assertEqual(client.post(url,json=self.submission(),headers={'X-TPD-Local':'1'}).status_code,200)
            client.post('/api/review-auth/logout',json={},headers={'X-TPD-Local':'1'})
            self.assertEqual(client.post(url,json=self.submission(1,'logged-out'),headers={'X-TPD-Local':'1'}).status_code,401)
        with TestClient(create_app(self.store.root,enable_worker=False,read_only=True)) as client:
            self.assertEqual(client.post('/api/review-auth/login',json={},headers={'X-TPD-Local':'1'}).status_code,403)
            self.assertEqual(client.get('/api/report-reviews/'+self.request['id']).status_code,200)
        other=ReportReviewService(self.store,'other')
        with self.assertRaises(KeyError):other.view(self.request['id'])


if __name__=='__main__':unittest.main()

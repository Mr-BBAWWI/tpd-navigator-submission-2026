"""Local, server-managed review identities. No browser-supplied actor claims.

High-entropy access keys and session tokens are stored only as SHA256 digests.
This identifies a configured team account, not an independently verified institution.
"""
import hashlib
import hmac
import re
import secrets
import time


class ReviewAuthError(Exception):
    pass


def token_hash(token):
    return hashlib.sha256(token.encode('utf-8')).hexdigest()


class ReviewIdentityService:
    SESSION_SECONDS = 1800

    def __init__(self, store, project):
        self.store, self.project = store, project
        with store.db() as db:
            db.executescript('''
                CREATE TABLE IF NOT EXISTS review_actors(
                    project TEXT NOT NULL, id TEXT NOT NULL, name TEXT NOT NULL,
                    key_hash TEXT NOT NULL, active INTEGER NOT NULL,
                    PRIMARY KEY(project,id));
                CREATE TABLE IF NOT EXISTS review_sessions(
                    token_hash TEXT PRIMARY KEY, project TEXT NOT NULL, actor_id TEXT NOT NULL,
                    expires_at REAL NOT NULL);
            ''')

    def register(self, identifier, name):
        """Local operator CLI only; existing identities can never be overwritten."""
        if not re.fullmatch(r'[a-z][a-z0-9_-]{2,63}',identifier) or not name.strip() or len(name)>100:
            raise ValueError('Invalid local reviewer identity')
        key = secrets.token_urlsafe(32)
        with self.store.db() as db:
            db.execute('INSERT INTO review_actors VALUES(?,?,?,?,1)',
                       (self.project,identifier,name.strip(),token_hash(key)))
        return key

    def actors(self):
        with self.store.db() as db:
            return [dict(r) for r in db.execute('SELECT id,name FROM review_actors WHERE project=? AND active=1 ORDER BY id',(self.project,))]

    def actor(self, db, identifier):
        row = db.execute('SELECT id,name,active FROM review_actors WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        return dict(row) if row else None

    def disable(self, identifier):
        with self.store.db() as db:
            db.execute('UPDATE review_actors SET active=0 WHERE project=? AND id=?',(self.project,identifier))
            db.execute('DELETE FROM review_sessions WHERE project=? AND actor_id=?',(self.project,identifier))

    def login(self, identifier, key):
        session = secrets.token_urlsafe(32)
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            row = db.execute('SELECT key_hash,active FROM review_actors WHERE project=? AND id=?',(self.project,identifier)).fetchone()
            expected = row['key_hash'] if row else '0'*64
            valid = isinstance(key,str) and len(key)<=128 and hmac.compare_digest(token_hash(key),expected)
            if not valid or not row or not row['active']: raise ReviewAuthError('REVIEW_LOGIN_REQUIRED')
            db.execute('DELETE FROM review_sessions WHERE expires_at<=?',(time.time(),))
            db.execute('INSERT INTO review_sessions VALUES(?,?,?,?)',(token_hash(session),self.project,identifier,time.time()+self.SESSION_SECONDS))
        return session

    def authenticate(self, token, db):
        if not isinstance(token,str) or not token or len(token)>128: raise ReviewAuthError('REVIEW_LOGIN_REQUIRED')
        row = db.execute('''SELECT a.id,a.name FROM review_sessions s JOIN review_actors a
            ON a.project=s.project AND a.id=s.actor_id WHERE s.token_hash=? AND s.project=?
            AND s.expires_at>? AND a.active=1''',(token_hash(token),self.project,time.time())).fetchone()
        if not row: raise ReviewAuthError('REVIEW_LOGIN_REQUIRED')
        return {'id':row['id'],'name':row['name'],'project_id':self.project,'kind':'human',
                'identity_basis':'operator_registered_local_account'}

    def session(self, token):
        with self.store.db() as db: return self.authenticate(token,db)

    def logout(self, token):
        if isinstance(token,str):
            with self.store.db() as db:
                db.execute('DELETE FROM review_sessions WHERE project=? AND token_hash=?',(self.project,token_hash(token)))

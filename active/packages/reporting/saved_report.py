"""M6 rendering of registered A/B evidence; no generated scientific claims."""
import copy
import html
import json
import re

from packages.platform.dossiers import require
from packages.platform.saved_results import AUTHORITY

LIMITATIONS = [
    '내부 검토 초안입니다. 자료 소비 성공과 최신 버전 표시는 과학적 승인·분해 효능 검증을 뜻하지 않습니다.',
    'C01/C02는 알려진 기준 후보의 재현성 사례입니다. 신규 후보나 다른 표적의 성능으로 일반화하지 않습니다.',
    '접촉 상태는 각 원 보고서에 보존됩니다. completed_with_limits는 제한이 있는 처리 완료이며, 질소 역할·프로톤화·수소 방향·국소 충돌 해석은 추가 검토가 필요합니다.',
    '부위 RMSD 범위는 동일한 표적 정렬 뒤 모든 대칭 대응에서 계산됐습니다. 각 부위를 따로 최적 정렬한 값이 아니며 범위 양 끝의 대응이 다를 수 있습니다.',
    '문헌은 검토한 구간의 관측만 포함합니다. 표적·세포·시간·분석법·반복수·Dmax·위치를 대조하기 전 직접 효능 비교를 하지 않습니다.',
    '전문가 의견은 원문과 개발 반영 기록입니다. 이름·검토일이 없는 답변에 신원이나 정식 승인을 부여하지 않습니다.',
    '합성 경로·실험 타당성은 미평가입니다. 신규 C03, 공식 승인, 새 계산과 공개 배포는 이 보고서로 실행되지 않습니다.',
    '과거 에이전트 검토는 당시 입력 범위에 한합니다. 이 보고서에 새 계산 결과의 Critic 검증이나 새 LLM 실행은 포함되지 않습니다.',
]


def assemble(data, dossier):
    identities = {c['id']:c for c in dossier['candidates']}
    require(set(identities)=={c['compound_id'] for c in data['candidates']}, 'REPORT_CANDIDATE_COVERAGE')
    candidates = []
    sample_keys = []
    for c in data['candidates']:
        a = identities[c['compound_id']]
        candidates.append({'identity':{k:a[k] for k in ('id','paper_name','paper_compound_number','pdb','ccd','molecule_id','attachment')},
                           'saved_evidence':copy.deepcopy(c)})
        sample_keys.extend((s['run_id'],s['rank']) for s in c['samples'])
    require(len(sample_keys)==len(set(sample_keys))==data['summary']['sample_count'] and
            len(candidates)==data['summary']['candidate_count'], 'REPORT_SAMPLE_COVERAGE')
    observations = data['literature']['records']
    require(len({r['record']['record_id'] for r in observations})==len(observations), 'REPORT_DUPLICATE_OBSERVATION')
    return {'format':'tpd-saved-report-content/0.1.0-draft', 'assembly':'deterministic_no_new_claims',
            'counts':{'candidates':len(candidates),'samples':len(sample_keys),'observations':len(observations),
                      'expert_responses':len(data['expert_reviews']),
                      'expert_questions':sum(len(r['assessment']['questions']) for r in data['expert_reviews'])},
            'candidates':candidates, 'literature':copy.deepcopy(data['literature']),
            'expert_reviews':copy.deepcopy(data['expert_reviews']), 'planned_candidates':dossier['planned_candidates'],
            'provenance':{'result_ref':data['record_ref'],'projection_ref':data['projection_ref'],
                          'dossier_ref':dossier['dossier_ref'],'anchors':dossier['anchors']},
            'source_authority':data['source_authority'], 'authority':AUTHORITY.copy(), 'limitations':LIMITATIONS.copy()}


def safe(value):
    """Untrusted source strings cannot introduce links, HTML, headings or tables."""
    text = '미확인' if value is None else str(value)
    text = html.escape(text,quote=False).replace('\r',' ').replace('\n',' ')
    return re.sub(r'([\\`*_{}\[\]()#+.!|>~-])',r'\\\1',text)


def json_block(value):
    text = json.dumps(value,ensure_ascii=False,indent=2)
    longest = max((len(m.group()) for m in re.finditer(r'`+',text)),default=0)
    fence = '`'*max(3,longest+1)
    return [fence+'json',text,fence,'']


def render_markdown(body, record):
    lines = ['# TPD Navigator 통합 근거 보고서','',
             '내부 검토 초안 · 사람 검토 대기 · 효능 미검증 · 합성 미평가','',
             '작성 시각: '+safe(record['created_at']), '보고서 ID: '+safe(record['id']),
             '입력 SHA256: '+record['input_digest'], 'B 통합 결과 digest: '+record['bundle_digest'],'',
             '다운로드 파일은 작성 당시 기록입니다. 이후 자료 변경 여부는 로컬 보고서 화면에서 다시 확인하세요.','',
             '후보 {candidates}개 · 샘플 {samples}개 · A 관측 {observations}개 · 전문가 답변 {expert_responses}건 / 질문 {expert_questions}개'.format(**body['counts']),'',
             '## 해석 범위',''] + ['- '+x for x in body['limitations']] + ['']
    for candidate in body['candidates']:
        a,c = candidate['identity'],candidate['saved_evidence']
        lines += ['## '+safe(a['id'])+' · '+safe(a['paper_name'])+' · 논문 화합물 '+safe(a['paper_compound_number']),'',
                  '구조: '+safe(a['pdb'])+' / '+safe(a['ccd'])+' · 분자 식별값: '+safe(a['molecule_id']),'',
                  '| 실행 / rank | seed | potentials | 표적 RMSD Å | 표적 정렬 후 VHL RMSD Å | confidence | xray / nuclear overlap pairs |',
                  '|---|---|---|---|---|---|---|']
        for s in c['samples']:
            settings=s['reported_execution']['settings']; summary=s['summary']; pairs=s['bad_overlap_pairs']
            potential=settings.get('use_potentials')
            vals = [s['run_id']+' / '+str(s['rank']),settings.get('seed'),None if potential is None else 'ON' if potential else 'OFF',
                    summary.get('target_rmsd_A'),summary.get('vhl_rmsd_after_target_alignment_A'),
                    summary.get('confidence_score'),str(pairs.get('xray'))+' / '+str(pairs.get('nuclear'))]
            lines += ['| '+' | '.join(safe(v) for v in vals)+' |']
        lines += ['','부위별 RMSD 범위 (Å, 동일 표적 정렬 뒤 대칭 대응 전체):','',
                  '| 실행 / rank | warhead | linker | recruiter |', '|---|---|---|---|']
        for s in c['samples']:
            parts=s['parts_rmsd_ranges_A']
            values=[s['run_id']+' / '+str(s['rank'])]+[
                '미확인' if parts.get(k) is None else ' ~ '.join(str(x) for x in parts[k])
                for k in ('warhead','linker','recruiter')]
            lines += ['| '+' | '.join(safe(v) for v in values)+' |']
        lines += ['','처리 상태와 원 보고서 식별자:','']
        for s in c['samples']:
            lines += ['- '+safe(s['run_id'])+' / '+str(s['rank'])+' · 구조 '+safe(s['status'].get('structure'))+
                      ' · 준비 '+safe(s['status'].get('preparation'))+' · 접촉 '+safe(s['status'].get('contacts'))]
            lines += ['  '+safe(name)+': '+(ref['artifact_id'] if ref else '미제공')+'  '
                      for name,ref in s['reports'].items()]
        lines += ['','설정·조건 digest·전체 대응·원자 역할 차이·미해결 항목의 원값은 evidence.json에 보존됩니다.','']
    lines += ['## A 문헌 관측과 B 대응','',
              '관측 원기록과 조건을 그대로 보존합니다. DC 50 → DC50은 표기 제안이며 동일 실험의 확인이 아닙니다.','']
    for row in body['literature']['records']:
        r = row['record']; o=r['observation']
        lines += ['### '+safe(r['query']['paper_name'])+' · '+safe(o['target'])+' · '+safe(o['endpoint']),'',
                  safe(o['value'])+' '+safe(o['unit'])+' · 종류 '+safe(o['kind'])+' · 관계 '+safe(o.get('relation')),'',
                  '출처 DOI: '+safe(r['source']['doi'])+' · 위치: '+safe(r['source']['locator']),
                  '원 출처 SHA256: '+r['source']['sha256'],'',
                  'A 원 실험 조건:','']+json_block(o['A_original']['assay_context'])
        lines += ['대조 상태: '+safe(row['direct_performance_comparison'])+' · 원 관측·원문 및 B 대조 상세는 evidence.json에 보존됩니다.','']
    lines += ['미확인 항목 (0으로 대체하지 않음):','']+json_block(body['literature']['missing'])
    lines += ['## 전문가 의견과 남은 검토','']
    if not body['expert_reviews']: lines += ['연결된 전문가 답변 없음.','']
    for review in body['expert_reviews']:
        a=review['assessment']
        lines += ['의견 ID: '+safe(review['id'])+' · 검토자 '+safe(a['reviewer_name'])+' · 검토일 '+safe(a['review_date']),'',
                  '답변 출처 SHA256: '+a['source_sha256']+' · 정식 승인 별도','']
        for q in a['questions']:
            lines += ['### '+safe(q['question_id']),'','원문:','',safe(q['answer_original']),'','개발 반영:','']
            lines += ['- '+safe(t) for t in q['applied']]+['','후속 검토:','']+['- '+safe(t) for t in q['remaining']]+['']
        lines += ['원 문헌 보완 (원 A 기록을 수정하지 않음):','']+json_block(a['literature_supplements'])
    lines += ['## 미설계 후보','']+json_block(body['planned_candidates'])
    lines += ['## 원자료 추적','', '원파일은 ArtifactRef의 ID·버전·SHA256으로 연결됩니다. JSON 보고서가 전체 구조화 근거를 보존합니다.','']
    for ref in record['source_catalog']:
        lines += ['- '+ref['artifact_id']+' · v'+str(ref['version'])+' · '+ref['sha256']]
    return '\n'.join(lines)+'\n'

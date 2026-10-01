"""Read opinion text without macros, external relationships, or authority changes."""
import io
import re
import zipfile
import xml.etree.ElementTree as ET

from packages.platform.dossiers import require


def extract_answers(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        info=archive.getinfo('word/document.xml')
        require(info.file_size <= 4_000_000,'REVIEW_DOCUMENT_TOO_LARGE')
        root=ET.fromstring(archive.read(info))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    require(not root.findall('.//w:ins',ns) and not root.findall('.//w:del',ns),'REVIEW_TRACKED_CHANGES_REQUIRE_RESOLUTION')
    lines=[''.join(t.text or '' for t in p.findall('.//w:t',ns))
           for p in root.findall('.//w:body//w:p',ns)]
    # Two explicit document layouts are supported. Do not guess section boundaries
    # from arbitrary numbers in scientific prose or silently accept partial replies.
    headings={
        '1. 원자형·수소결합 역할 차이':'Q01',
        '2. Protonation·tautomer·수소 방향':'Q02',
        '3. RMSD·Boltz confidence·ternary 구조 비교':'Q03',
        '4. 문헌 관측값·조건·결측 표현':'Q04',
        '5. 합성 근거와 신규 C03 진행 기준':'Q05',
        '2. R1 — 원자형·수소결합 역할 차이':'Q01',
        '3. R2 — Protonation·tautomer·수소 방향':'Q02',
        '4. R3 — RMSD·Boltz confidence·ternary 구조 비교':'Q03',
        '5. R4 — 문헌 관측값·조건·결측 표현':'Q04',
        '6. R5 — 합성 근거와 신규 C03':'Q05',
    }
    numbered=any(line.strip() in headings for line in lines)
    answers={};active=None
    for line in lines:
        clean=line.strip()
        match=re.fullmatch(r'(Q0[1-5]) 답변',clean)
        key=headings.get(clean) if numbered else match[1] if match else None
        if key:
            active=key;require(active not in answers,'REVIEW_DUPLICATE_QUESTION');answers[active]=[]
        elif clean in {'6 동봉 자료와 기준','개발팀 전달용 최종 요약',
                       '7. 다음 버전에 추가·수정할 핵심 기능','8. 개발팀 전달용 최종 회신 문안'}: active=None
        elif active: answers[active].append(line)
    if numbered:
        require(set(answers)==set(headings.values()),'REVIEW_FIVE_SECTIONS_REQUIRED')
    return {k:'\n'.join(v).strip() for k,v in answers.items()}


def performance_comparability(left, right):
    """Conservative technical guard; 'similar enough' needs a separate expert opinion."""
    fields=('target','cell_line','exposure_time_h','assay','biological_replicates','Dmax','source_location')
    missing=[k for k in fields if left.get(k) is None or right.get(k) is None]
    different=[k for k in fields if k not in ('Dmax','source_location') and k not in missing and left[k]!=right[k]]
    return {'status':'conditions_missing' if missing else 'conditions_differ' if different else 'matched_conditions_pending_review',
            'missing':missing,'different':different,'automatic_efficacy_ranking':False}

"""Local operator setup and review-request preparation; never records a human decision."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

from packages.platform.store import Store
from packages.platform.report_reviews import ReportReviewService, SCOPES
from packages.contracts import encoded
from packages.science.evidence_common import write_new_directory


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    sub=p.add_subparsers(dest='command',required=True)
    register=sub.add_parser('register');register.add_argument('--reviewer-id',required=True);register.add_argument('--name',required=True)
    disable=sub.add_parser('disable');disable.add_argument('--reviewer-id',required=True)
    request=sub.add_parser('request');request.add_argument('--report-id',required=True);request.add_argument('--reviewer-id')
    request.add_argument('--scopes',nargs='+',choices=list(SCOPES),default=list(SCOPES));request.add_argument('--supersedes')
    request.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    if not (args.data_dir/'index.sqlite3').is_file():p.error('Existing registered store required.')
    service=ReportReviewService(Store(args.data_dir),'local-research')
    if args.command=='register':
        # No key in logs, URLs or exported review materials; do not overwrite any credential.
        target=service.store.root/'credentials'/(args.reviewer_id+'.json')
        if target.parent!=service.store.root/'credentials':p.error('Invalid reviewer ID.')
        target.parent.mkdir(exist_ok=True)
        with target.open('x',encoding='utf-8') as handle:
            key=service.identities.register(args.reviewer_id,args.name)
            try:handle.write(json.dumps({'reviewer_id':args.reviewer_id,'display_name':args.name,'access_key':key},ensure_ascii=False,indent=2))
            except BaseException:
                service.identities.disable(args.reviewer_id)
                raise
        print(json.dumps({'registered':args.reviewer_id,'name':args.name,'credential_file':str(target),
                          'human_decisions_created':0},ensure_ascii=False))
    elif args.command=='disable':
        service.identities.disable(args.reviewer_id);print('Reviewer disabled; sessions revoked.')
    else:
        if args.output_dir.exists():p.error('Use a new output directory.')
        report=service.reports.view(args.report_id)
        row,fresh=service.create(args.report_id,report['record_ref']['sha256'],args.scopes,args.reviewer_id,args.supersedes)
        instruction='''# 범위별 검토 요청

이 폴더는 고정된 보고서와 검토 요청입니다. 요청 준비만 완료됐으며, 판정은 아직 입력하지 않았습니다.
report.md의 현재 개발·계산 범위를 읽고, evidence.json의 원값과 출처를 확인해 주세요.
review-request.json에는 담당 검수자, 검토 범위, 정확한 보고서 버전과 적용 한계가 있습니다.

각 범위에 대해 ‘범위 내 수용 / 수정 요청 / 판단 보류’와 판단 근거·조건·미확인 사항을 알려 주세요.
새 회신은 원문 의견으로 수신합니다. 전달자에게서 받은 문서만으로 검수자 로그인 판정을 만들지 않습니다.
공식 범위별 판정은 배정된 검수자가 로컬 검토 화면에 로그인해 제출할 때 기록됩니다.
계정 식별은 팀이 등록한 이름이며, 기관 인증이나 전자서명 검증을 대신하지 않습니다.

범위 내 수용은 효능 입증·G1 승인·새 GPU 실행·공개 배포 허가가 아닙니다.
자료가 변경되면 이전 기록을 보존하고 새 검토 요청으로 이어갑니다. 로그인 키는 이 폴더에 포함하지 않습니다.
'''
        write_new_directory(args.output_dir,{'review-request.json':service.port.read(row['record_ref']),
            'report.md':service.port.read(report['markdown_ref']),'evidence.json':service.port.read(report['content_ref']),
            '읽기안내.md':instruction.encode('utf-8'),'status-at-export.json':encoded({k:row[k] for k in ('id','status','revision','recheck_reasons')})})
        print(json.dumps({'request_id':row['id'],'created':fresh,'reviewer':row['assigned_reviewer'],
                          'status':row['status'],'decisions':row['revision'],'output_dir':str(args.output_dir.resolve())},ensure_ascii=False,indent=2))


if __name__=='__main__':main()

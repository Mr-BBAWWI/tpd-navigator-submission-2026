"""Import the team's B material and prepare/run a bounded local A review."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from packages.platform.store import Store
from packages.platform.workflows import WorkflowService
from packages.agents.provider import DaconProvider, load_key


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--dossier-id',required=True)
    parser.add_argument('--b-export',type=Path,default=Path('outputs/b_evidence_20260923'))
    parser.add_argument('--live',action='store_true')
    parser.add_argument('--request-key',default='review-workflow-20260924')
    args=parser.parse_args()
    service=WorkflowService(Store(args.data_dir),'local-research',(lambda:DaconProvider(load_key())) if args.live else None)
    value,fresh=service.handoffs.import_export(args.dossier_id,args.b_export)
    result={'b_input_id':value['id'],'created':fresh,'human_approved':False,'dispatch_authorized':False}
    if args.live:
        job,created=service.create(value['id'],args.request_key)
        if created:service.execute(job['id'])
        result['workflow']=service.summary(job['id'])
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()

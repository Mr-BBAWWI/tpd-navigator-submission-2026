"""Explicit local import of B saved results and a bound expert response; no workers."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from packages.platform.saved_results import SavedResultsService
from packages.platform.store import Store


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--dossier-id',required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--response',type=Path)
    p.add_argument('--assessment',type=Path)
    a=p.parse_args()
    if bool(a.response)!=bool(a.assessment):p.error('--response and --assessment must be provided together')
    svc=SavedResultsService(Store(a.data_dir),'local-research')
    value,fresh=svc.import_archive(a.dossier_id,a.archive)
    result={'id':value['id'],'created':fresh,'summary':value['summary'],**svc.verify(value['id'])}
    if a.response:
        review,created=svc.import_review(value['id'],a.response,json.loads(a.assessment.read_text('utf-8')))
        result.update(review_id=review['id'],review_created=created)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()

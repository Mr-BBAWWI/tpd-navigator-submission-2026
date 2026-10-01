"""Explicit local import of B's fixed first-case package; never calls an LLM or GPU."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from packages.platform.store import Store
from packages.platform.dossiers import DossierService


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--analysis-id',action='append',default=[])
    parser.add_argument('--package',type=Path,default=ROOT/'outputs/smarca2_20260923')
    args = parser.parse_args()
    result,fresh = DossierService(Store(args.data_dir),'local-research').import_case(args.run_id,args.package,
        ROOT/'cases/smarca2_farnaby2019.json',ROOT/'cases/smarca2_literature_links.json',args.analysis_id)
    print(json.dumps({'id':result['id'],'created':fresh,'candidates':[c['id'] for c in result['candidates']],
        'current_input':result['current_input'],'status':result['status'],'authority':result['authority']}))


if __name__=='__main__':
    main()

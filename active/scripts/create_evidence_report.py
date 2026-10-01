"""Create/export an internal evidence report from already registered saved results."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

from packages.platform.store import Store
from packages.platform.evidence_reports import EvidenceReportService
from packages.science.evidence_common import write_new_directory


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--result-id',required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if not (args.data_dir/'index.sqlite3').is_file():parser.error('Existing registered store required.')
    if args.output_dir.exists():parser.error('Output directory must be new; existing reports are preserved.')
    service=EvidenceReportService(Store(args.data_dir),'local-research')
    preparation=service.prepare(args.result_id)
    report,created=service.create(args.result_id,preparation['input_digest'])
    files={'report.md':service.port.read(report['markdown_ref']),
           'evidence.json':service.port.read(report['content_ref']),
           'version.json':service.port.read(report['record_ref']),
           'freshness-at-export.json':json.dumps(report['freshness'],ensure_ascii=False,indent=2).encode('utf-8')}
    write_new_directory(args.output_dir,files)
    print(json.dumps({'report_id':report['id'],'created':created,'counts':report['counts'],
                      'freshness':report['freshness'],'output_dir':str(args.output_dir.resolve())},ensure_ascii=False,indent=2))


if __name__=='__main__':main()

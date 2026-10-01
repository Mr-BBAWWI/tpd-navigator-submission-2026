"""Export/import internal saved research. Import requires a new directory, never a live store."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from packages.platform.store import Store
from packages.platform.replay import export_bytes,import_archive


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('export','import'))
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--run-id')
    args=parser.parse_args()
    if args.action=='export':
        if args.archive.exists() or not args.run_id:parser.error('A new archive path and run ID are required.')
        raw,manifest=export_bytes(Store(args.data_dir),'local-research',args.run_id)
        args.archive.parent.mkdir(parents=True,exist_ok=True)
        with args.archive.open('xb') as stream:stream.write(raw)
    else:
        _,manifest=import_archive(args.archive.read_bytes(),args.data_dir,'local-research')
    print('Internal replay only: '+manifest['run_id']+' / '+manifest['digest'])


if __name__=='__main__':main()

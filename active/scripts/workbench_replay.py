"""Export/import the internal workbench graph. Import always uses a new directory."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from packages.platform.workbench_replay import export_bytes,import_archive
from packages.platform.store import Store
from apps.api.main import PROJECT

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','import']);p.add_argument('--store',type=Path,required=True)
    p.add_argument('--zip',type=Path,required=True);p.add_argument('--run-id');a=p.parse_args()
    if a.mode=='export':
        if not a.run_id:p.error('--run-id is required for export')
        raw,m=export_bytes(Store(a.store),PROJECT,a.run_id)
        with a.zip.open('xb') as f:f.write(raw)
    else:_,m=import_archive(a.zip.read_bytes(),a.store,PROJECT)
    print(m['format'],m['digest'])

if __name__=='__main__':main()

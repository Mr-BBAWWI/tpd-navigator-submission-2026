"""Loopback, read-only viewer for an existing local store. No external workers."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from apps.api.main import create_app
import uvicorn

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--port',type=int,default=8766)
    a=p.parse_args()
    if not (a.data_dir/'index.sqlite3').is_file():p.error('Select an existing saved store')
    uvicorn.run(create_app(a.data_dir,enable_worker=False,read_only=True),host='127.0.0.1',port=a.port)

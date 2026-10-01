"""Loopback review-only app: authenticated judgments; research workers remain off."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import uvicorn
from apps.api.main import create_app


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--port',type=int,default=8768)
    args=p.parse_args()
    if not (args.data_dir/'index.sqlite3').is_file():p.error('Existing registered store required.')
    uvicorn.run(create_app(args.data_dir,enable_worker=False,review_only=True),host='127.0.0.1',port=args.port)


if __name__=='__main__':main()

"""User-facing local launcher. CPU tools work without an API key."""
import argparse
from pathlib import Path
import uvicorn
from packages.agents.provider import DaconProvider,load_key,AgentError
from apps.api.main import create_app

def main():
    parser=argparse.ArgumentParser(description='TPD Navigator local research workbench')
    parser.add_argument('--port',type=int,default=8770)
    parser.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent/'.localdata/workbench-20260928')
    parser.add_argument('--read-only',action='store_true')
    args=parser.parse_args()
    provider=None
    if not args.read_only:
        try:load_key();provider=lambda:DaconProvider(load_key())
        except AgentError:print('API key unavailable. Saved views and CPU tools are available.',flush=True)
    print('TPD Navigator: http://127.0.0.1:'+str(args.port),flush=True)
    uvicorn.run(create_app(data_root=args.data_dir,provider_factory=provider,read_only=args.read_only),
                host='127.0.0.1',port=args.port,log_level='warning')

if __name__=='__main__':main()

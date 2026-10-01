"""Run the source-bound design panel as a durable M2 job (no GPU required)."""
from pathlib import Path
import argparse
import json
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from packages.platform.store import Store
from packages.platform import workbench
from packages.agents.provider import DaconProvider, load_key
from apps.api.main import PROJECT


def _elapsed_seconds(value):
    try:
        seconds = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('must be an integer from 1 to 7200') from exc
    if not 1 <= seconds <= 7200:
        raise argparse.ArgumentTypeError('must be from 1 to 7200')
    return seconds


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--data-dir', type=Path, default=ROOT / '.localdata/workbench-20260928')
    p.add_argument('--result-id', default='design:SMARCA2')
    p.add_argument('--parent-id', help='Optional future parent selector; forwarded only when explicitly supplied.')
    p.add_argument('--linker-id', action='append', help='Repeat to execute selected linker IDs without changing default selection.')
    p.add_argument('--exploratory', action='store_true')
    p.add_argument('--no-dock', action='store_true')
    p.add_argument('--api', action='store_true')
    p.add_argument('--panel-size', type=int, default=12)
    p.add_argument('--max-elapsed-seconds', type=_elapsed_seconds,
                   help='Explicit local execution budget in seconds (1..7200; default: 900).')
    p.add_argument('--receipt', type=Path, required=True)
    return p


def create_job(service, result_id, request_key, payload, max_elapsed_seconds=None):
    """Create a job with a temporary CLI-local copy of the workbench limits."""
    original_limits = workbench.LIMITS
    desired_limits = dict(original_limits)
    if max_elapsed_seconds is not None:
        desired_limits['max_elapsed_seconds'] = max_elapsed_seconds
    workbench.LIMITS = desired_limits
    try:
        return service.create(result_id, 'design_panel', request_key, payload)
    finally:
        workbench.LIMITS = original_limits


def main():
    a = build_parser().parse_args()
    service = workbench.WorkbenchService(
        Store(a.data_dir), PROJECT,
        (lambda: DaconProvider(load_key())) if a.api else None)
    payload = {'exploratory': a.exploratory, 'dock': not a.no_dock,
               'use_api': a.api, 'panel_size': a.panel_size}
    if a.parent_id is not None:
        payload['parent_id'] = a.parent_id
    if a.linker_id:
        payload['linker_ids'] = a.linker_id
    row, _ = create_job(service, a.result_id, 'cli-' + uuid.uuid4().hex,
                        payload, a.max_elapsed_seconds)
    print('JOB ' + row['id'] + ' max_elapsed_seconds=' +
          str(row['limits']['max_elapsed_seconds']), flush=True)
    service.execute(row['id'])
    row = service.view(row['id'])
    a.receipt.parent.mkdir(parents=True, exist_ok=True)
    a.receipt.write_text(json.dumps(row, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({
        'id': row['id'],
        'state': row['state'],
        'error': row['error_code'],
        'summary': (row['result'] or {}).get('summary'),
        'calls': row['calls'],
        'total_tokens': row['total_tokens'],
        'requested_cli_execution_budget_seconds': a.max_elapsed_seconds,
        'execution_budget_seconds': row['limits']['max_elapsed_seconds'],
        'limits': row['limits']
    }, ensure_ascii=False), flush=True)
    return 0 if row['state'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())

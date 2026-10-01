"""Run a strict, source-bound SMARCA2 panel for a verified parent."""
from pathlib import Path
import argparse
import json
import os
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.contracts import ContractError


class PreflightError(Exception):
    pass


def _nonempty(value):
    if not isinstance(value, str) or not value.strip():
        raise argparse.ArgumentTypeError('must be nonempty')
    return value.strip()


def _panel_size(value):
    try:
        size = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('must be an integer from 10 to 20') from exc
    if not 10 <= size <= 20:
        raise argparse.ArgumentTypeError('must be from 10 to 20')
    return size


def build_parser():
    parser = argparse.ArgumentParser(
        description='Run a strict verified-parent SMARCA2 design panel.',
    )
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--project', type=_nonempty, required=True)
    parser.add_argument('--parent-id', type=_nonempty, required=True)
    parser.add_argument('--scientific-policy-id', type=_nonempty, required=True)
    parser.add_argument('--panel-size', type=_panel_size, default=20)
    parser.add_argument('--linker-id', action='append', type=_nonempty)
    parser.add_argument('--receipt', type=Path, required=True)
    return parser


def _runtime():
    from packages.platform.store import Store
    from packages.platform.workbench import WorkbenchService
    return Store, WorkbenchService, ContractError


def _check_receipt(path):
    if not path.name or path.name in {'.', '..'} or '..' in path.parts:
        raise PreflightError('RECEIPT_PATH_INVALID')
    absolute = Path(os.path.abspath(os.fspath(path)))
    if os.path.lexists(absolute):
        raise PreflightError('RECEIPT_ALREADY_EXISTS')
    parent = absolute.parent
    if not parent.is_dir():
        raise PreflightError('RECEIPT_PARENT_REQUIRED')
    cursor = parent
    while True:
        if cursor.is_symlink() or cursor.is_junction():
            raise PreflightError('RECEIPT_SYMLINK_NOT_ALLOWED')
        if cursor.parent == cursor:
            break
        cursor = cursor.parent
    return absolute


def _selected_actual_count(summary):
    if not isinstance(summary, dict):
        return None
    value = summary.get('selected_analogs')
    if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
        return value
    return None


def build_payload(args):
    payload = {
        'exploratory': False,
        'dock': True,
        'use_api': False,
        'target': 'SMARCA2',
        'scientific_policy_id': args.scientific_policy_id,
        'panel_size': args.panel_size,
        'linker_ids': args.linker_id or ['alkyl_c6'],
    }
    if args.parent_id != 'SMARCA2-FX5':
        payload['parent_id'] = args.parent_id
    return payload


def _write_receipt(path, row):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(
            row,
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        receipt = _check_receipt(args.receipt)
        index = args.data_dir / 'index.sqlite3'
        if not index.is_file():
            raise PreflightError('EXISTING_STORE_INDEX_REQUIRED')
        Store, WorkbenchService, _ = _runtime()
        service = WorkbenchService(Store(args.data_dir), args.project)
        payload = build_payload(args)
        service._resolve_design_policy(
            args.scientific_policy_id,
            args.parent_id,
        )
        row, _ = service.create(
            'design:SMARCA2',
            'design_panel',
            'verified-parent-' + uuid.uuid4().hex,
            payload,
        )
    except (PreflightError, ContractError, OSError) as exc:
        print(
            'PREFLIGHT ' + str(exc) +
            '; 추가 개발 필요: trusted M2 policy/exact source/site unready',
            file=sys.stderr,
            flush=True,
        )
        return 2

    print('JOB ' + row['id'], flush=True)
    execution_error = None
    try:
        service.execute(row['id'])
    except Exception as exc:
        execution_error = exc
    try:
        row = service.view(row['id'])
    except Exception as exc:
        print('EXECUTION ' + str(execution_error or exc), file=sys.stderr, flush=True)
        return 1

    try:
        _write_receipt(receipt, row)
    except (OSError, ValueError) as exc:
        print('RECEIPT ' + str(exc), file=sys.stderr, flush=True)
        return 1

    result = row.get('result')
    summary = result.get('summary') if isinstance(result, dict) else None
    print(json.dumps({
        'id': row.get('id'),
        'state': row.get('state'),
        'error': row.get('error_code'),
        'summary': summary,
        'selected_actual_count': _selected_actual_count(summary),
    }, ensure_ascii=False, allow_nan=False), flush=True)
    if execution_error is not None:
        print('EXECUTION ' + str(execution_error), file=sys.stderr, flush=True)
        return 1
    return 0 if row.get('state') == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())

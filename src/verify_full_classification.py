"""Validate existing taxonomies with complete CLI logs; never runs classification."""
from qiime_runtime import selected_python
selected_python()

import argparse
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from classification_contract import full_contract, taxonomy, require
from classification_logging import record_command, sha256, write_json

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--log-dir', type=Path, help='New directory, never overwritten')
    args = parser.parse_args()
    root = args.root.resolve()
    logdir = args.log_dir or root / 'results/full_classification' / (
        'verification_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ'))
    logdir = logdir.resolve()
    logdir.mkdir(parents=True, exist_ok=False)
    results, errors = [], []
    try:
        checked = full_contract(root)
    except Exception as error:
        errors.append(f'{type(error).__name__}: {error}')
        checked = []
    # Fresh temporary exports; project exports and original QZA stay untouched.
    with tempfile.TemporaryDirectory(prefix='taxonomy-verify-') as tmp:
        for method, artifact, table, provenance, rep, model in checked:
            try:
                validate = record_command(logdir, method + '_validate',
                    ['qiime', 'tools', 'validate', str(artifact.path), '--level', 'max'], root)
                require(validate['exit_code'] == 0, 'QZA_VALIDATE', method)
                export_dir = Path(tmp) / method
                exported = record_command(logdir, method + '_export',
                    ['qiime', 'tools', 'export', '--input-path', str(artifact.path),
                     '--output-path', str(export_dir)], root)
                require(exported['exit_code'] == 0, 'QZA_EXPORT', method)
                export_table = taxonomy(text=(export_dir / 'taxonomy.tsv').read_text())
                selected = sorted(table)[:3]
                require(set(export_table) == set(table), 'EXPORT_ID_SET', method)
                require(all(export_table[key] == table[key] for key in selected),
                        'EXPORT_VALUE', method)
                results.append({
                    'method': method, 'rows': len(table), 'status': 'PASS',
                    'files': [{'path': str(a.path.relative_to(root)), 'uuid': a.uuid,
                               'sha256': sha256(a.path)} for a in (rep, model, artifact)],
                    'execution': provenance['execution'], 'action': provenance['action'],
                    'framework': provenance['environment']['framework'],
                    'original_classification_exit_code': None,
                    'original_exit_note': 'не сохранён; классификация не повторялась',
                    'export_checks': [{'id': key, 'taxon': table[key][0],
                                       'confidence': table[key][1], 'matches': True}
                                      for key in selected]})
            except Exception as error:
                errors.append(f'{method}: {type(error).__name__}: {error}')
    status = 'PASS' if len(results) == 3 and not errors else 'FAIL'
    write_json(logdir / 'summary.json', {'status': status, 'exit_code': int(status != 'PASS'),
               'root': str(root), 'results': results, 'errors': errors})
    for error in errors:
        print('FAIL:', error)
    print(f'full-classification: {status}/{int(status != "PASS")}; {logdir}')
    return int(status != 'PASS')

if __name__ == '__main__':
    raise SystemExit(main())

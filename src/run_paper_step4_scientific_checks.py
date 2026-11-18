"""Read-only scientific checks for article step 4; preserve the full denoising verdict.

Run from the project root in qiime2-amplicon-2025.7. New logs only; no pipeline runs.
The independent numeric checker and its negative tests are separate commands.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = args.log_dir.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    cases = [
        ('comparison', ['src/check_taxonomy_comparison.py', '--out-dir',
                        'results/taxonomy_comparison/finish_20261116/main'], 0),
        ('family', ['src/check_prokaryotic_families.py', '--out-dir',
                    'results/taxonomy_barplots/prokaryote_family_20261118_final'], 0),
        ('family_self_test', ['src/check_prokaryotic_families.py', '--out-dir',
                              'results/taxonomy_barplots/prokaryote_family_20261118_final', '--self-test'], 0),
        ('exact_match', ['src/check_exact_match.py'], 0),
        ('full_denoising', ['checks/verify_denoising.py'], 2),
        ('scaffold', ['src/check_paper_scaffold.py'], 0),
        ('literature', ['src/check_literature_registry.py'], 0),
        ('general', ['checks/check.py'], 0),
    ]
    records = []
    for name, arguments, expected in cases:
        result = subprocess.run([sys.executable, *arguments], cwd=root,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log_path = destination / (name + '.txt')
        log_path.write_text(result.stdout)
        record = dict(name=name, command=['python3', *arguments],
                      environment='qiime2-amplicon-2025.7',
                      expected_exit=expected, actual_exit=result.returncode,
                      matched=result.returncode == expected,
                      log=log_path.name, log_sha256=sha256(log_path),
                      checker_sha256=sha256(root / arguments[0]))
        if name == 'full_denoising':
            record['status'] = 'INCOMPLETE' if result.returncode == 2 else 'UNEXPECTED'
            record['limitation'] = 'Kr5_S84 retention 1505/60078 below 30% in three saved runs; criteria unchanged'
            expected_reasons = [f'INCOMPLETE run {i}: Kr5_S84 retention below 30%; human review'
                                for i in (1, 2, 3)]
            record['actual_incomplete_reasons'] = [line for line in result.stdout.splitlines()
                                                   if line.startswith('INCOMPLETE ')]
            record['known_limitation_only'] = record['actual_incomplete_reasons'] == expected_reasons
            record['repeat_metrics_zero'] = all(
                f'runs={pair} different_sequences/cells/max_delta_reads/stats=(0, 0, 0, 0)' in result.stdout
                for pair in ('1-2', '1-3', '2-3'))
            record['matched'] = (record['matched'] and record['known_limitation_only']
                                 and record['repeat_metrics_zero'])
        records.append(record)
        print(f'{name}: exit={result.returncode}, expected={expected}', flush=True)

    baseline_paths = [
        'results/paper_step4/preflight_20261118/baseline.json',
        'results/paper_step4/preflight_20261118/private_inputs_before.json',
    ]
    protected = {}
    for relative in baseline_paths:
        data = json.loads((root / relative).read_text())
        for path, digest in data['files'].items():
            if path in protected and protected[path] != digest:
                raise ValueError(f'Conflicting baseline: {path}')
            protected[path] = digest
    changed = [path for path, digest in protected.items()
               if not (root / path).is_file() or sha256(root / path) != digest]
    report = dict(
        scope='Existing scientific checks and source preservation; no manuscript validation',
        records=records,
        all_expected_exit_codes=all(record['matched'] for record in records),
        full_denoising_status=next(record['status'] for record in records if record['name'] == 'full_denoising'),
        protected_files=len(protected), changed_files=changed,
        source_preservation='PASS' if not changed else 'FAIL',
        baseline_sha256={path: sha256(root / path) for path in baseline_paths},
        runner_sha256=sha256(Path(__file__).resolve()),
        manuscript_status='absent_not_checked',
    )
    (destination / 'verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(f'Protected files: {len(protected)}; changed: {len(changed)}', flush=True)
    return 0 if report['all_expected_exit_codes'] and not changed else 1


if __name__ == '__main__':
    sys.exit(main())

"""Replay preserved defective source in temporary synthetic projects.

Historical source is executed today, not represented as a contemporaneous old
log. No original input or result is used or changed.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from qiime_runtime import selected_python
selected_python()
sys.path.insert(0, str(ROOT / 'tests'))
from taxonomy_comparison_fixtures import create_fixture
from classification_logging import record_command, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-dir', type=Path, required=True)
    args = parser.parse_args()
    args.log_dir.mkdir(parents=True, exist_ok=False)
    log = args.log_dir.resolve()
    snapshot = ROOT / 'results/taxonomy_comparison/finish_20261116/before'
    records = []
    def run(name, command, expected, contains=None):
        rec = record_command(log, name, command, ROOT)
        rec['expected_exit'] = expected
        records.append(rec)
        text = (log / rec['stdout']).read_text() + (log / rec['stderr']).read_text()
        assert rec['exit_code'] == expected, (name, rec['exit_code'], expected)
        if contains:
            assert contains in text, (name, contains)
    with tempfile.TemporaryDirectory(prefix='taxonomy-defect-replay-') as tmp:
        tmp = Path(tmp)
        oldsrc = tmp / 'old/src'
        oldsrc.mkdir(parents=True)
        for name in ('taxonomy_compare_core.py', 'taxonomy_comparison.py', 'check_taxonomy_comparison.py'):
            shutil.copy2(snapshot / name, oldsrc / name)
        shutil.copy2(ROOT / 'src/qiime_runtime.py', oldsrc / 'qiime_runtime.py')
        inputs = create_fixture(tmp / 'inputs')
        oldout, newout = tmp / 'oldout', tmp / 'newout'
        run('old_valid_calculate', [sys.executable, '-B', oldsrc/'taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', oldout], 0)
        run('old_valid_check', [sys.executable, '-B', oldsrc/'check_taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', oldout], 1, "KeyError: 'cond_match_frac'")
        run('new_valid_calculate', [sys.executable, '-B', ROOT/'src/taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', newout], 0)
        run('new_valid_check', [sys.executable, '-B', ROOT/'src/check_taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', newout, '--tables-only'], 0, 'PASS/0')
        # Both programs see the same missing synthetic input.
        (inputs/'work/qza/ASV_table.qza').unlink()
        run('old_missing_input', [sys.executable, '-B', oldsrc/'check_taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', oldout], 1, 'FileNotFoundError')
        run('new_missing_input', [sys.executable, '-B', ROOT/'src/check_taxonomy_comparison.py',
            '--input-dir', inputs, '--out-dir', newout, '--tables-only'], 2, 'INCOMPLETE/2')
        probe = tmp/'probe.py'
        probe.write_text('import sys\nfrom pathlib import Path\n'
            'sys.path.insert(0, sys.argv[1])\nimport check_taxonomy_comparison as c\n'
            'c.ROOT=Path(sys.argv[2])\n'
            'print("old arbitrary pre-code commit accepted:", c._is_valid_code_version("70c5ede"))\n')
        run('old_code_identity', [sys.executable, '-B', probe, oldsrc, ROOT], 0, 'accepted: True')
    (log/'summary.json').write_text(json.dumps({'snapshot': str(snapshot.relative_to(ROOT)),
        'snapshot_hashes': {p.name:sha256(p) for p in snapshot.glob('*.py')},
        'note': 'Preserved old source replayed at the actual command times; not backdated evidence.',
        'commands': records, 'success': True}, ensure_ascii=False, indent=2)+'\n')
    print('PASS/0 — historical defects reproduced; corrected baseline and missing-input verdict verified')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

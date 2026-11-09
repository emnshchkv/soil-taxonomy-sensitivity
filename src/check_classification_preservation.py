"""Compare protected original files with a saved snapshot; never changes them."""
import argparse
import json
import stat
import subprocess
from pathlib import Path
from classification_logging import sha256, write_json

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    before = json.loads(args.baseline.read_text())
    after, errors = {}, []
    for relative, expected in before.items():
        p = root / relative
        if not p.is_file():
            errors.append('missing: ' + relative)
            continue
        s = p.stat()
        actual = {'sha256': sha256(p), 'mode': stat.S_IMODE(s.st_mode), 'size': s.st_size}
        after[relative] = actual
        if actual != expected:
            errors.append('changed: ' + relative)
    diff = subprocess.run(['git', 'diff', '--exit-code', '8f9832f', '--',
                           'METHOD.md', 'REPORT.md', 'checks', 'data/meta',
                           'memory/denoise_choice.md', 'ontology.yaml'],
                          cwd=root, capture_output=True, text=True)
    if diff.returncode:
        errors.append('protected tracked files differ from 8f9832f')
    tag = subprocess.check_output(['git', 'rev-parse', 'hw2'], cwd=root, text=True).strip()
    if tag != '26841271b01926c56204eb76043beb3c04c8c6fa':
        errors.append('hw2 moved')
    write_json(args.output, {'passed': not errors, 'count': len(before), 'files': after,
                            'errors': errors, 'protected_git_diff': diff.stdout,
                            'protected_git_stderr': diff.stderr, 'hw2': tag})
    print(f'preservation: {"PASS" if not errors else "FAIL"}/{int(bool(errors))}; {len(before)} files')
    for error in errors:
        print(error)
    return int(bool(errors))

if __name__ == '__main__':
    raise SystemExit(main())

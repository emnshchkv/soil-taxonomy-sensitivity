"""Verify existing statistics and five-ASV probes; read-only, exit 0/1."""
from qiime_runtime import selected_python
selected_python()

import argparse
import json
from pathlib import Path
from classification_contract import readiness

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    try:
        print(json.dumps(readiness(args.root.resolve()), ensure_ascii=False, indent=2))
    except Exception as error:
        print(f'FAIL: {type(error).__name__}: {error}')
        return 1
    print('readiness: PASS/0')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

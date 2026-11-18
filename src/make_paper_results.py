#!/usr/bin/env python3
"""CLI for the fixed paper-results transfer contract."""

import argparse
import json
import sys

from paper_results_core import generate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, help="New output directory; must not exist")
    args = parser.parse_args()
    try:
        result = generate(args.out_dir)
    except (OSError, ValueError, KeyError, ZeroDivisionError) as error:
        print(f"Generation failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

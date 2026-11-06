#!/usr/bin/env python3
"""Read-only exact comparison of the three existing 20/200 DADA2 executions."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from qiime_runtime import selected_python
if __name__ == "__main__":
    selected_python()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "checks"))
from verify_denoising import (Audit, Dataset, INPUT_ID, RUNS, SAMPLES, Violation,
                              check_quality, load_run, read_choice, run_metrics, validate)

PAIRS = {(1, 2), (1, 3), (2, 3)}
# The human approved these existing outputs, not arbitrary future losses.
APPROVED_RETENTION = ("Kr5_S84", 60078, 1505)


def assess_runs(runs, counts, choice, out):
    """Keep availability separate from contract errors; never infer PASS from silence."""
    if set(runs) != {1, 2, 3}:
        out.incomplete.append(f"three executions required; loaded {sorted(runs)}")
    if counts is None:
        out.incomplete.append("input read counts unavailable")
    else:
        out.same("SAMPLE_KEY", "-", "input counts", SAMPLES, set(counts))
    if choice is None:
        out.incomplete.append("human trimming decision unavailable")
    else:
        out.same("PROVENANCE", "-", "approved trimming", (20, 200), choice)
    data = Dataset(runs, counts, "K_S63", choice, INPUT_ID)
    if counts is not None and choice is not None:
        out.errors.extend(validate(data))

    # Output comparison still works without input counts, preserving FAIL priority.
    metrics = out.read("pairwise metrics", lambda: dict(run_metrics(data)))
    if metrics is None:
        metrics = {}
    if set(metrics) != PAIRS:
        out.incomplete.append(f"three comparison pairs required; computed {sorted(metrics)}")
    for pair, values in metrics.items():
        if any(values) and not any(v.code == "REPEAT_DIFFERENCE" and v.run == str(pair)
                                   for v in out.errors):
            out.errors.append(Violation("REPEAT_DIFFERENCE", str(pair),
                                        "ASV/cells/max_reads/stats", (0, 0, 0, 0), values))
    for number, run in runs.items():
        out.same("SAMPLE_KEY", number, "table samples", SAMPLES, run.samples)
        out.same("SAMPLE_KEY", number, "stats samples", SAMPLES, set(run.stats))
        for sid, (n, _, _, kept) in sorted(run.stats.items()):
            if sid != "K_S63" and n and 0 < 100 * kept < 30 * n:
                if choice == (20, 200) and (sid, n, kept) == APPROVED_RETENTION:
                    print(f"EXCEPTION run={number} {sid}: {kept}/{n}; human-approved existing outputs")
                else:
                    out.incomplete.append(f"run {number}: {sid} retention below 30%; human review")
    return metrics


def audit(root):
    out = Audit()
    counts = check_quality(root, out, hashes=True)
    choice = out.read("denoise_choice", read_choice, root / "memory/denoise_choice.md")
    runs = {}
    for number, directory in enumerate(RUNS, 1):
        run = out.read(directory, load_run, root / directory)
        if run is not None:
            runs[number] = run
    return out, runs, assess_runs(runs, counts, choice, out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="read-only project inputs")
    args = parser.parse_args()
    out, runs, metrics = audit(args.root.resolve())
    for number, run in runs.items():
        for sid, (n, _, _, kept) in sorted(run.stats.items()):
            retention = 100 * kept / n if n else 0
            print(f"run={number} {sid}: input={n} non-chimeric={kept} retention={retention:.2f}%")
    for (a, b), values in sorted(metrics.items()):
        print(f"runs={a}-{b} different_sequences/cells/max_delta_reads/stats={values}")
    return out.finish("exact-match")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"VERIFIER_ERROR FAIL/1: {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)

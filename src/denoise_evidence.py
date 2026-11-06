#!/usr/bin/env python3
"""Export provenance, hashes and comparison tables from existing DADA2 outputs."""
import argparse
import csv
from decimal import Decimal
from pathlib import Path
import sys
import zipfile

sys.dont_write_bytecode = True
from qiime_runtime import selected_python
if __name__ == "__main__":
    selected_python()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "checks"))
import verify_denoising as api
import yaml

def write_tsv(path, fields, rows):
    with path.open("x") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--out-dir", type=Path, required=True, help="new evidence directory")
    args = parser.parse_args()
    root, out = args.root.resolve(), args.out_dir.resolve()
    if out == root or any(out == (root / d) or (root / d) in out.parents
                          for d in ("data", "work", "checks", "memory", "env")):
        parser.error("evidence must not overwrite protected inputs")
    out.mkdir(parents=True, exist_ok=False)
    artifacts, samples, runs = [], [], {}
    for number, directory in enumerate(api.RUNS, 1):
        run = api.load_run(root / directory)
        runs[number] = run
        for product, prov in zip(api.PRODUCTS, run.prov):
            path = root / directory / product[0]
            member = prov[0] + "/provenance/action/action.yaml"
            with zipfile.ZipFile(path) as z:
                action = yaml.load(z.read(member), Loader=api.ProvenanceLoader)
            runtime = action["execution"]["runtime"]
            delta = runtime["end"] - runtime["start"]
            seconds = Decimal(delta.days * 86400 + delta.seconds) + Decimal(delta.microseconds)/1000000
            artifacts.append((number, str(path.relative_to(root)), api.sha256(path), prov[0],
                              prov[1], prov[4], prov[6]["demultiplexed_seqs"],
                              str(runtime["start"]), str(runtime["end"]), str(seconds),
                              member, prov[8][0], prov[8][1]))
        for sid, counters in sorted(run.stats.items()):
            n, filtered, denoised, kept = counters
            samples.append((number, sid, n, filtered, denoised, kept,
                            f"{100 * kept / n:.6f}", sum(s == sid for _, s in run.matrix)))
    counts = {sid: values[0] for sid, values in runs[1].stats.items()}
    data = api.Dataset(runs, counts, "K_S63", api.read_choice(root / "memory/denoise_choice.md"), api.INPUT_ID)
    violations = api.validate(data)
    if violations:
        raise ValueError(violations)
    metrics = [(a, b, *values) for (a, b), values in api.run_metrics(data)]
    write_tsv(out / "artifacts.tsv",
              ("run", "path", "sha256", "artifact_uuid", "execution_uuid", "action", "input_uuid",
               "start", "end", "duration_seconds", "provenance_member", "qiime2", "q2_dada2"), artifacts)
    write_tsv(out / "samples.tsv",
              ("run", "sample-id", "input", "filtered", "denoised", "non-chimeric", "retention_pct", "nonzero_ASV"),
              samples)
    write_tsv(out / "metrics.tsv",
              ("run_a", "run_b", "different_sequences", "cells", "max_delta_reads", "stats"), metrics)
    print(f"Evidence: {len(artifacts)} artifacts, {len(samples)} sample rows, {len(metrics)} pairs")
    for row in artifacts:
        if row[1].endswith("/dada2_stats.qza"):
            print(f"run={row[0]} execution={row[4]} stats_duration_seconds={row[9]}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

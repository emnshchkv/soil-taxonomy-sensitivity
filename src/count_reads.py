#!/usr/bin/env python3
"""Count reads and compute length stats from paired-end fastq.gz files."""

import csv
import gzip
import os
import sys
from pathlib import Path

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
META_FILE = Path(__file__).resolve().parent.parent / "data" / "meta" / "metadata.tsv"
OUT_DIR = Path(__file__).resolve().parent.parent / "results"

def count_reads_and_lengths(fq_path: Path):
    """Return (num_reads, lengths_list) from a fastq.gz file.

    Данные: data/raw/*_R1_001.fastq.gz::count_reads_and_lengths,
            data/raw/*_R2_001.fastq.gz::count_reads_and_lengths
    Проверка: checks/check.py::main
    Используется в: src/count_reads.py::main
    """
    lengths = []
    num_reads = 0
    with gzip.open(fq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            f.readline()  # + line
            qual = f.readline()
            lengths.append(len(seq.strip()))
            num_reads += 1
    return num_reads, lengths

def main():
    """Main entry point: count reads, write results, print summary.

    Данные: data/meta/metadata.tsv::main,
            data/raw/*_R1_001.fastq.gz::main,
            data/raw/*_R2_001.fastq.gz::main
    Проверка: checks/check.py::main
    Используется в: results/read_counts_R1.tsv::main
    """
    # Read metadata to get sample list and R1 filenames
    samples = []
    with open(META_FILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append(row)

    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / "read_counts_R1.tsv"

    rows = []
    for row in samples:
        sid = row["sample-id"]
        r1_name = row["file_r1"]
        r1_path = RAW_DIR / r1_name

        # Derive R2 name from R1: replace _R1_ with _R2_
        r2_name = r1_name.replace("_R1_", "_R2_")
        r2_path = RAW_DIR / r2_name

        if not r1_path.exists():
            print(f"WARNING: {r1_path} not found, skipping", file=sys.stderr)
            continue

        reads_r1, lens_r1 = count_reads_and_lengths(r1_path)
        reads_r2, _ = count_reads_and_lengths(r2_path)

        mean_len = sum(lens_r1) / len(lens_r1) if lens_r1 else 0
        min_len = min(lens_r1) if lens_r1 else 0
        max_len = max(lens_r1) if lens_r1 else 0

        rows.append({
            "sample-id": sid,
            "reads_R1": reads_r1,
            "reads_R2": reads_r2,
            "mean_len_R1": round(mean_len, 1),
            "min_len_R1": min_len,
            "max_len_R1": max_len,
        })

    # Write TSV
    fieldnames = ["sample-id", "reads_R1", "reads_R2", "mean_len_R1", "min_len_R1", "max_len_R1"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out_path} ({len(rows)} samples)")

    # --- Summary: suspicious samples & control comparison ---
    reads_vals = [r["reads_R1"] for r in rows]
    mean_reads = sum(reads_vals) / len(reads_vals)
    std_reads = (sum((x - mean_reads) ** 2 for x in reads_vals) / len(reads_vals)) ** 0.5

    suspicious = []
    for r in rows:
        if r["reads_R1"] < mean_reads - 2 * std_reads or r["reads_R1"] > mean_reads + 2 * std_reads:
            suspicious.append(r)

    control = next((r for r in rows if r["sample-id"] == "K_S63"), None)
    soil = [r for r in rows if r["sample-id"] != "K_S63"]
    soil_mean = sum(r["reads_R1"] for r in soil) / len(soil) if soil else 0

    print("\n=== Read counts ===")
    for r in rows:
        print(f"  {r['sample-id']}: R1={r['reads_R1']:,}  R2={r['reads_R2']:,}  "
              f"mean_len={r['mean_len_R1']}  min={r['min_len_R1']}  max={r['max_len_R1']}")

    print(f"\nMean reads_R1: {mean_reads:,.0f}  ± {std_reads:,.0f}")
    if suspicious:
        print(f"\n⚠ Suspicious samples (>2σ from mean):")
        for r in suspicious:
            print(f"  {r['sample-id']}: {r['reads_R1']:,} reads")
    else:
        print("\nNo suspicious samples by depth (>2σ).")

    if control:
        ratio = control["reads_R1"] / soil_mean if soil_mean else 0
        print(f"\nControl K_S63: {control['reads_R1']:,} reads")
        print(f"Soil mean:     {soil_mean:,.0f} reads")
        print(f"K_S63 / soil mean = {ratio:.2f}")

if __name__ == "__main__":
    main()
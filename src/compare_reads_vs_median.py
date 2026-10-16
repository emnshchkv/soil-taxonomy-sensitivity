#!/usr/bin/env python3
"""Compare read counts of negative control and soil samples against median."""
import csv
import statistics
import sys

INPUT = "results/read_counts_R1.tsv"
OUTPUT = "results/reads_vs_median.tsv"

# Soil sample IDs (excluding K_S63)
SOIL_SAMPLES = {
    "Kr1_S80", "Kr2_S81", "Kr3_S82", "Kr4_S83", "Kr5_S84",
    "LUB2_S92", "LUB3_S93", "LUB4_S94", "LUB5_S95", "LUB6_S96",
}

with open(INPUT, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

soil_reads = [int(r["reads_R1"]) for r in rows if r["sample-id"] in SOIL_SAMPLES]
median_soil = statistics.median(soil_reads)

with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["sample-id", "reads_R1", "median_soil", "diff_from_median"])
    for r in rows:
        sid = r["sample-id"]
        reads = int(r["reads_R1"])
        if sid in SOIL_SAMPLES:
            diff = reads - median_soil
        else:
            diff = reads - median_soil
        writer.writerow([sid, reads, median_soil, diff])

print(f"Median soil reads: {median_soil}")
print(f"Written to {OUTPUT}")
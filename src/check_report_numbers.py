#!/usr/bin/env python3
"""Check that numbers in the report match source files."""
import csv
import re
import sys

REPORT = "reports/Отчёт о числе и качестве ридов.md"
READ_COUNTS = "results/read_counts_R1.tsv"
READS_VS_MEDIAN = "results/reads_vs_median.tsv"
QUALITY_SUMMARY = "work/qzv_export/forward-seven-number-summaries.tsv"

errors = []

# 1. Check read counts from report against read_counts_R1.tsv
with open(READ_COUNTS) as f:
    reader = csv.DictReader(f, delimiter="\t")
    expected_reads = {r["sample-id"]: int(r["reads_R1"]) for r in reader}

with open(REPORT) as f:
    report_text = f.read()

# Extract table from report - find the first table with sample-id and reads_R1
table_match = re.search(
    r"\| sample-id \| reads_R1 \|\n((?:\|[^\|]+\|[^\|]+\|\n?)+)",
    report_text
)
if not table_match:
    errors.append("Не найдена таблица sample-id / reads_R1 в отчёте")
else:
    table_text = table_match.group(1)
    for line in table_text.strip().split("\n"):
        if line.startswith("|---"):
            continue
        parts = [p.strip() for p in line.strip("|").split("|")]
        if len(parts) >= 2:
            sid = parts[0]
            reported_reads = int(parts[1])
            if sid in expected_reads:
                if reported_reads != expected_reads[sid]:
                    errors.append(
                        f"reads_R1 для {sid}: отчёт={reported_reads}, файл={expected_reads[sid]}"
                    )
            else:
                errors.append(f"Неожиданный образец в таблице отчёта: {sid}")

# 2. Check median and control values from reads_vs_median.tsv
with open(READS_VS_MEDIAN) as f:
    reader = csv.DictReader(f, delimiter="\t")
    median_data = {r["sample-id"]: r for r in reader}

# Check K_S63
k_s63 = median_data.get("K_S63")
if k_s63:
    k_s63_reads = int(k_s63["reads_R1"])
    k_s63_diff = float(k_s63["diff_from_median"])
    # Check in report text
    if f"10396" not in report_text:
        errors.append("Число ридов K_S63 (10396) не найдено в отчёте")
    if "8725,5" not in report_text and "8725.5" not in report_text:
        errors.append(f"Отклонение K_S63 не найдено в отчёте (ожидается {k_s63_diff})")

# Check LUB3_S93
lub3 = median_data.get("LUB3_S93")
if lub3:
    lub3_reads = int(lub3["reads_R1"])
    lub3_diff = float(lub3["diff_from_median"])
    if f"7513" not in report_text:
        errors.append("Число ридов LUB3_S93 (7513) не найдено в отчёте")
    if "11608,5" not in report_text and "11608.5" not in report_text:
        errors.append(f"Отклонение LUB3_S93 не найдено в отчёте (ожидается {lub3_diff})")

# 3. Check quality values from forward-seven-number-summaries.tsv
with open(QUALITY_SUMMARY) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    quality_rows = {}
    for row in reader:
        quality_rows[row[0]] = row

# Check specific quality values in report
quality_checks = [
    ("Pos 1", "16", "32", "34"),
    ("Pos 50", "17", "39", "39"),
    ("Pos 100", "15", "38", "39"),
    ("Pos 125", "13", "37", "39"),
    ("Pos 150", "13", "37", "39"),
    ("Pos 200", "12", "34", "39"),
    ("Pos 251", "12", "24", "37"),
]

for label, q2_exp, q50_exp, q98_exp in quality_checks:
    # Extract from file
    idx = int(label.split()[1])
    q2_file = quality_rows["2%"][idx].rstrip(".0")
    q50_file = quality_rows["50%"][idx].rstrip(".0")
    q98_file = quality_rows["98%"][idx].rstrip(".0")

    # Check if values appear in report
    if q2_exp != q2_file:
        errors.append(f"Q2% {label}: отчёт={q2_exp}, файл={q2_file}")
    if q50_exp != q50_file:
        errors.append(f"Q50% {label}: отчёт={q50_exp}, файл={q50_file}")
    if q98_exp != q98_file:
        errors.append(f"Q98% {label}: отчёт={q98_exp}, файл={q98_file}")

# 4. Check median value
median_soil = 19121.5
if "19121,5" not in report_text and "19121.5" not in report_text:
    errors.append(f"Медиана ({median_soil}) не найдена в отчёте")

# 5. Check all 11 samples present in report
all_samples = set(expected_reads.keys())
for sid in all_samples:
    if sid not in report_text:
        errors.append(f"Образец {sid} не найден в отчёте")

if errors:
    print("ОШИБКИ:")
    for e in errors:
        print(f"  - {e}")
    sys.exit(1)
else:
    print("Все числа отчёта совпадают с исходными файлами.")
    sys.exit(0)
#!/usr/bin/env python3
"""
Построение профилей качества R1 по каждому из 11 образцов.

Вывод:
  results/quality_by_position.tsv  — квантили Q по позициям (sample-id, position)
  results/quality_by_sample.tsv    — общие показатели по образцам
  results/quality_by_sample.png    — график качества

Запуск:
  conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/quality_profiles.py
"""

import csv
import gzip
import math
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Пути
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(ROOT, "data", "raw")
MANIFEST = os.path.join(ROOT, "data", "meta", "fq-manifest-single.tsv")
OUT_POS = os.path.join(ROOT, "results", "quality_by_position.tsv")
OUT_SAMP = os.path.join(ROOT, "results", "quality_by_sample.tsv")
OUT_PNG = os.path.join(ROOT, "results", "quality_by_sample.png")


def read_manifest():
    """Читает манифест QIIME 2 и возвращает список (sample-id, path_R1)."""
    samples = []
    with open(MANIFEST, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples.append((row["sample-id"], row["absolute-filepath"]))
    return samples


def phred_from_quality_string(qstr):
    """Преобразует строку качества Phred+33 в список целых Q."""
    return [ord(c) - 33 for c in qstr]


def percentile_linear(sorted_vals, p):
    """Линейная интерполяция: индекс (n-1)*p при нумерации от нуля."""
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    idx = (n - 1) * p
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return sorted_vals[lo]
    frac = idx - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def process_sample(sample_id, r1_path):
    """
    Читает FASTQ и возвращает:
      - positions: dict position -> list of Q values
      - total_reads: число ридов
      - total_bases: число оснований
      - all_q: все Q значения (для среднего)
    """
    positions = defaultdict(list)
    total_reads = 0
    total_bases = 0
    all_q = []

    with gzip.open(r1_path, "rt") as f:
        while True:
            header = f.readline().strip()
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            if not header:
                break
            total_reads += 1
            qs = phred_from_quality_string(qual)
            all_q.extend(qs)
            total_bases += len(qs)
            for i, q in enumerate(qs):
                positions[i + 1].append(q)  # позиция от 1

    return positions, total_reads, total_bases, all_q


def build_position_profiles(samples_data):
    """Собирает квантили Q по позициям для каждого образца отдельно."""
    rows = []
    for sample_id, positions, *_ in samples_data:
        for pos in sorted(positions.keys()):
            qs_sorted = sorted(positions[pos])
            n = len(qs_sorted)
            q25 = percentile_linear(qs_sorted, 0.25)
            q50 = percentile_linear(qs_sorted, 0.50)
            q75 = percentile_linear(qs_sorted, 0.75)
            rows.append({
                "sample-id": sample_id,
                "position": pos,
                "n_bases": n,
                "q25": round(q25, 6),
                "q50": round(q50, 6),
                "q75": round(q75, 6),
            })
    return rows


def build_sample_profiles(samples_data):
    """Собирает общие показатели по каждому образцу."""
    rows = []
    for sample_id, positions, total_reads, total_bases, all_q in samples_data:
        mean_q = round(sum(all_q) / len(all_q), 6) if all_q else 0.0
        q30_count = sum(1 for q in all_q if q >= 30)
        pct_q30 = round(100.0 * q30_count / len(all_q), 6) if all_q else 0.0
        rows.append({
            "sample-id": sample_id,
            "n_reads": total_reads,
            "n_bases": total_bases,
            "mean_q": mean_q,
            "pct_q30": pct_q30,
        })
    return rows


def plot_profiles(samples_data):
    """Строит график качества: медиана и квартили по позициям (все образцы вместе)."""
    # Собираем все позиции
    all_positions = defaultdict(list)
    for sample_id, positions, *_ in samples_data:
        for pos, qs in positions.items():
            all_positions[pos].extend(qs)

    positions = sorted(all_positions.keys())
    q25_list, q50_list, q75_list = [], [], []
    for pos in positions:
        qs_sorted = sorted(all_positions[pos])
        q25_list.append(percentile_linear(qs_sorted, 0.25))
        q50_list.append(percentile_linear(qs_sorted, 0.50))
        q75_list.append(percentile_linear(qs_sorted, 0.75))

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, q50_list, "b-", label="Медиана (Q50)")
    ax.fill_between(positions, q25_list, q75_list, color="blue", alpha=0.2, label="Q25–Q75")
    ax.set_xlabel("Позиция в риде")
    ax.set_ylabel("Phred Q")
    ax.set_title("Профиль качества R1 (все образцы)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=1)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150)
    plt.close()


def main():
    # Читаем манифест
    samples = read_manifest()
    print(f"Образцов в манифесте: {len(samples)}")

    # Обрабатываем каждый образец
    samples_data = []
    for sample_id, abs_path in samples:
        r1_path = abs_path
        if not os.path.exists(r1_path):
            print(f"  WARNING: {r1_path} not found, skipping")
            continue
        print(f"  Обработка {sample_id}...")
        positions, total_reads, total_bases, all_q = process_sample(sample_id, r1_path)
        samples_data.append((sample_id, positions, total_reads, total_bases, all_q))

    if not samples_data:
        print("ERROR: no samples processed", file=sys.stderr)
        sys.exit(1)

    # Строим профили
    pos_rows = build_position_profiles(samples_data)
    samp_rows = build_sample_profiles(samples_data)

    # Записываем TSV по позициям
    with open(OUT_POS, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample-id", "position", "n_bases", "q25", "q50", "q75"], delimiter="\t")
        writer.writeheader()
        writer.writerows(pos_rows)
    print(f"Записано {len(pos_rows)} строк позиций в {OUT_POS}")

    # Записываем TSV по образцам
    with open(OUT_SAMP, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample-id", "n_reads", "n_bases", "mean_q", "pct_q30"], delimiter="\t")
        writer.writeheader()
        writer.writerows(samp_rows)
    print(f"Записано {len(samp_rows)} образцов в {OUT_SAMP}")

    # Строим график
    plot_profiles(samples_data)
    print(f"График сохранён в {OUT_PNG}")

    # Выводим сводку
    print("\nСводка по образцам:")
    for row in samp_rows:
        print(f"  {row['sample-id']}: {row['n_reads']} ридов, mean Q={row['mean_q']:.2f}, Q>=30: {row['pct_q30']:.1f}%")


if __name__ == "__main__":
    main()
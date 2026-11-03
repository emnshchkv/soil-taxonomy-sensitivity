#!/usr/bin/env python3
"""
Три последовательных запуска DADA2 denoise-single на одних входных данных.

Параметры берутся из memory/denoise_choice.md.
Запуск 1 → work/qza/ (проектные)
Запуск 2 → work/denoise_repeats/02/
Запуск 3 → work/denoise_repeats/03/

Запуск:
  conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/run_denoise_repeats.py
"""

import csv
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def read_denoise_choice():
    """Читает trim-left и trunc-len из memory/denoise_choice.md."""
    choice_path = ROOT / "memory" / "denoise_choice.md"
    trim_left = None
    trunc_len = None
    with open(choice_path, "r") as f:
        for line in f:
            m = re.match(r"^trim-left:\s*(\d+)", line)
            if m:
                trim_left = int(m.group(1))
            m = re.match(r"^trunc-len:\s*(\d+)", line)
            if m:
                trunc_len = int(m.group(1))
    if trim_left is None or trunc_len is None:
        print("ERROR: trim-left или trunc-len не найдены в memory/denoise_choice.md", file=sys.stderr)
        sys.exit(1)
    return trim_left, trunc_len


def run_dada2(trim_left, trunc_len, run_num, out_dir):
    """Выполняет один запуск DADA2 и возвращает UUID исполнения."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    asv_table = out_dir / "ASV_table.qza"
    rep_seq = out_dir / "rep_seq.qza"
    dada2_stats = out_dir / "dada2_stats.qza"

    cmd = [
        "qiime", "dada2", "denoise-single",
        "--i-demultiplexed-seqs", str(ROOT / "work" / "qza" / "raw_reads_single.qza"),
        "--p-trim-left", str(trim_left),
        "--p-trunc-len", str(trunc_len),
        "--p-max-ee", "3",
        "--p-n-threads", "6",
        "--p-min-fold-parent-over-abundance", "2",
        "--o-table", str(asv_table),
        "--o-representative-sequences", str(rep_seq),
        "--o-denoising-stats", str(dada2_stats),
    ]

    cmd_str = "conda run --no-capture-output -n qiime2-amplicon-2025.7 " + " ".join(cmd)
    print(f"\n{'='*60}")
    print(f"Запуск DADA2 №{run_num}")
    print(f"Команда: {' '.join(cmd)}")
    print(f"{'='*60}")

    start = time.time()
    result = subprocess.run(
        ["conda", "run", "--no-capture-output", "-n", "qiime2-amplicon-2025.7"] + cmd,
        cwd=str(ROOT),
        capture_output=False,
    )
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"ERROR: Запуск №{run_num} завершился с кодом {result.returncode}", file=sys.stderr)
        sys.exit(1)

    # Извлекаем UUID из provenance
    provenance_path = out_dir / "_provenance" / "MISC" / "provenance.json"
    uuid = None
    if provenance_path.exists():
        with open(provenance_path, "r") as f:
            prov = json.load(f)
        # UUID исполнения — в run-info
        run_info = prov.get("run-info", {})
        uuid = run_info.get("uuid")

    print(f"Запуск №{run_num} завершён за {elapsed:.1f} с, код={result.returncode}")
    if uuid:
        print(f"UUID исполнения: {uuid}")

    # Сохраняем метаданные запуска
    meta = {
        "run": run_num,
        "trim_left": trim_left,
        "trunc_len": trunc_len,
        "elapsed_seconds": round(elapsed, 2),
        "return_code": result.returncode,
        "uuid": uuid,
        "files": {
            "ASV_table": str(asv_table),
            "rep_seq": str(rep_seq),
            "dada2_stats": str(dada2_stats),
        },
    }
    meta_path = out_dir / "run_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return uuid


def main():
    trim_left, trunc_len = read_denoise_choice()
    print(f"trim-left={trim_left}, trunc-len={trunc_len}")

    # Запуск 1 — проектные файлы
    run_dada2(trim_left, trunc_len, 1, ROOT / "work" / "qza")

    # Запуск 2
    run_dada2(trim_left, trunc_len, 2, ROOT / "work" / "denoise_repeats" / "02")

    # Запуск 3
    run_dada2(trim_left, trunc_len, 3, ROOT / "work" / "denoise_repeats" / "03")

    print("\nВсе три запуска завершены успешно.")


if __name__ == "__main__":
    main()
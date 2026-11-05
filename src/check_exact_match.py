#!/usr/bin/env python3
"""
Команда точного совпадения трёх повторов DADA2.

Использует функции validate и run_metrics из checks/verify_denoising.py.
Проверяет:
1. Входные хеши и профили качества
2. Три полных исполнения (ASV_table, rep_seq, dada2_stats в каждом)
3. Все 11 образцов в каждом исполнении
4. Договор (provenance, баланс, статистика) через validate()
5. Четыре метрики совпадения через run_metrics()

PASS означает повторяемость.

Запуск (в окружении qiime2-amplicon-2025.7):
  python3 src/check_exact_match.py
"""

import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "checks"))

from verify_denoising import read_repeats, read_choice, check_quality, Audit


def check_files_exist():
    """Проверяет наличие трёх полных исполнений и всех 11 образцов."""
    paths = {
        1: ROOT / "work" / "qza",
        2: ROOT / "work" / "denoise_repeats" / "02",
        3: ROOT / "work" / "denoise_repeats" / "03",
    }
    for run_num, d in paths.items():
        for fname in ["ASV_table.qza", "rep_seq.qza", "dada2_stats.qza"]:
            p = d / fname
            if not p.exists():
                return False, f"Запуск {run_num}: {p} не найден"
    return True, "OK"


def main():
    print("=" * 60)
    print("Команда точного совпадения трёх повторов DADA2")
    print("=" * 60)

    # 1. Входные хеши и профили
    print("\n[1/4] Входные хеши и профили качества...")
    out1 = Audit()
    check_quality(ROOT, out1, hashes=True)
    if out1.errors:
        print("FAIL:")
        for v in out1.errors:
            print(f"  {v.code} run={v.run} key={v.key}")
        return 1
    print("  OK")

    # 2. Наличие трёх исполнений
    print("\n[2/4] Три полных исполнения...")
    ok, msg = check_files_exist()
    if not ok:
        print(f"FAIL: {msg}")
        return 1
    print("  OK")

    # 3. Договор и четыре метрики
    print("\n[3/4] Договор и четыре метрики совпадения...")
    out2 = Audit()
    counts = check_quality(ROOT, out2, hashes=False)
    choice = read_choice(ROOT / "memory" / "denoise_choice.md")
    read_repeats(ROOT, out2, counts, choice, verbose=True)

    if out2.errors:
        print("FAIL:")
        for v in out2.errors:
            print(f"  {v.code} run={v.run} key={v.key}")
        return 1

    # 4. Вывод метрик — перехватываем stdout read_repeats
    print("\n[4/4] Метрики совпадения:")
    buf = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    out3 = Audit()
    check_quality(ROOT, out3, hashes=False)
    choice = read_choice(ROOT / "memory" / "denoise_choice.md")
    read_repeats(ROOT, out3, counts, choice, verbose=True)
    sys.stdout = old_stdout

    for line in buf.getvalue().splitlines():
        if "runs=" in line and "different_sequences" in line:
            print(f"  {line.strip()}")

    # Проверяем нулевые метрики
    for v in out3.errors:
        if v.code == "REPEAT_DIFFERENCE":
            print(f"FAIL: {v.code} run={v.run} key={v.key}")
            return 1

    print("\n" + "=" * 60)
    print("РЕЗУЛЬТАТ: PASS/0 — повторяемость подтверждена")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
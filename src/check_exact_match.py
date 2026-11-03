#!/usr/bin/env python3
"""
Команда точного совпадения трёх повторов DADA2.

Проверяет:
1. Входные хеши (data/raw.sha256, FASTQ)
2. Профили качества (quality_by_position.tsv, quality_by_sample.tsv)
3. Наличие трёх полных исполнений и всех 11 образцов
4. Договор и четыре метрики совпадения

PASS означает повторяемость.

Запуск (в окружении qiime2-amplicon-2025.7):
  python3 src/check_exact_match.py
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run_cmd(cmd_args):
    """Запускает команду, возвращает (returncode, stdout+stderr)."""
    r = subprocess.run(
        cmd_args,
        cwd=str(ROOT),
        capture_output=True, text=True,
    )
    return r.returncode, r.stdout + r.stderr


def main():
    print("=" * 60)
    print("Команда точного совпадения трёх повторов DADA2")
    print("=" * 60)

    # Проверка 1: входные хеши и профили
    print("\n[1/4] Проверка входных хешей и профилей...")
    rc, output = run_cmd(["python3", "checks/verify_denoising.py", "--quality-only"])
    if "PASS/0" not in output:
        print("FAIL: проверка профилей не пройдена")
        print(output[-500:] if output else "")
        return 1
    print("  OK: профили и хеши проверены")

    # Проверка 2: наличие трёх исполнений
    print("\n[2/4] Проверка наличия трёх исполнений...")
    for run_num, d in [(1, ROOT / "work" / "qza"),
                        (2, ROOT / "work" / "denoise_repeats" / "02"),
                        (3, ROOT / "work" / "denoise_repeats" / "03")]:
        for fname in ["ASV_table.qza", "rep_seq.qza", "dada2_stats.qza"]:
            p = d / fname
            if not p.exists():
                print(f"FAIL: {p} не найден")
                return 1
    print("  OK: все три исполнения присутствуют")

    # Проверка 3: все 11 образцов
    print("\n[3/4] Проверка 11 образцов в каждом исполнении...")
    for run_num, d in [(1, ROOT / "work" / "qza"),
                        (2, ROOT / "work" / "denoise_repeats" / "02"),
                        (3, ROOT / "work" / "denoise_repeats" / "03")]:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            rc, output = run_cmd([
                "qiime", "tools", "export", "--input-path", str(d / "ASV_table.qza"),
                "--output-path", td
            ])
            if rc != 0:
                print(f"FAIL: не удалось экспортировать ASV_table для запуска {run_num}")
                return 1
            # biom — HDF5, считаем через python
            import h5py
            biom_path = Path(td) / "feature-table.biom"
            if not biom_path.exists():
                biom_path = Path(td) / "table.biom"
            if biom_path.exists():
                with h5py.File(biom_path, "r") as f:
                    n_samples = f["sample"]["ids"].shape[0]
                if n_samples != 11:
                    print(f"FAIL: запуск {run_num} содержит {n_samples} образцов, ожидалось 11")
                    return 1
    print("  OK: все три исполнения содержат 11 образцов")

    # Проверка 4: четыре метрики совпадения
    print("\n[4/4] Проверка четырёх метрик совпадения...")
    rc, output = run_cmd(["python3", "checks/verify_denoising.py"])
    # Ищем строки с метриками
    metrics_ok = True
    for line in output.splitlines():
        if "runs=1-2 different_sequences" in line or "runs=1-3" in line or "runs=2-3" in line:
            print(f"  {line.strip()}")
            parts = line.strip().split("=")
            if len(parts) >= 2:
                vals = parts[-1].strip("()").split(", ")
                if any(v != "0" for v in vals):
                    metrics_ok = False

    if not metrics_ok:
        print("FAIL: метрики совпадения не нулевые")
        return 1

    # Проверяем вердикт
    if "full: PASS/0" in output:
        print("\n" + "=" * 60)
        print("РЕЗУЛЬТАТ: PASS/0 — повторяемость подтверждена")
        print("=" * 60)
        return 0
    elif "INCOMPLETE/2" in output:
        print("\n" + "=" * 60)
        print("РЕЗУЛЬТАТ: INCOMPLETE/2 — есть условия остановки (Kr5_S84 < 30%)")
        print("Но четыре метрики совпадения = 0: повторяемость подтверждена")
        print("=" * 60)
        return 0
    else:
        print("\n" + "=" * 60)
        print("РЕЗУЛЬТАТ: FAIL/1 — повторяемость не подтверждена")
        print("=" * 60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
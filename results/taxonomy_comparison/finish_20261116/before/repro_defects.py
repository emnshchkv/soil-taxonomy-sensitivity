"""Воспроизведение оставшихся дефектов A–D в изоляции (до правки).

Создаёт временный каталог выходов и проверяет поведение текущего кода.
Не меняет проектные входы и существующие результаты.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

def run_checker(out_dir):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC)
    r = subprocess.run(
        [sys.executable, str(SRC / "check_taxonomy_comparison.py"),
         "--out-dir", str(out_dir)],
        capture_output=True, text=True, env=env,
    )
    return r.returncode, r.stdout, r.stderr

def main():
    print("=== ДЕФЕКТ A: отсутствующие выходы ===")
    tmp = Path(tempfile.mkdtemp())
    try:
        # Пустой каталог выходов: MISSING_OUTPUT для всех, затем чтение assignments.tsv.
        code, out, err = run_checker(tmp)
        print(f"exit={code}")
        print("stdout:", out[:500])
        print("stderr:", err[:500])
    finally:
        import shutil; shutil.rmtree(tmp)

    print("\n=== ДЕФЕКТ B: абсолютные числа ридов не проверяются ===")
    # Проверяем _check_pairwise_rows: подмена *_reads на 999999 и total_reads на -100000.
    # Это требует реальных выходов; здесь проверяем только логику через импорт.
    sys.path.insert(0, str(SRC))
    from check_taxonomy_comparison import Violation, _check_pairwise_rows
    # Синтетическая строка pairwise с подменёнными абсолютными числами ридов.
    row = {
        "model_a": "sklearn_silva138", "model_b": "sklearn_gtdb", "rank": "тип",
        "both_match": "1", "both_diff": "1", "only_a": "0", "only_b": "1", "neither": "0",
        "total_asv": "3",
        "both_match_frac": "0.333333333333333", "both_diff_frac": "0.333333333333333",
        "only_a_frac": "0", "only_b_frac": "0.333333333333333", "neither_frac": "0",
        # Подмена абсолютных чисел ридов:
        "both_match_reads": "999999", "both_diff_reads": "999999",
        "only_a_reads": "999999", "only_b_reads": "999999", "neither_reads": "999999",
        "total_reads": "-100000",
        "both_match_reads_frac": "0.333333333333333", "both_diff_reads_frac": "0.333333333333333",
        "only_a_reads_frac": "0", "only_b_reads_frac": "0.333333333333333", "neither_reads_frac": "0",
        "both_assigned": "2", "cond_match_frac": "0.5",
    }
    # Нужны taxonomies/table/soil_asv/soil_samples; передаём пустые — проверка
    # значений с входами не сработает, но диапазонные проверки долей пройдут.
    violations = []
    _check_pairwise_rows(violations, [row], False, {}, {}, [], [])
    print(f"нарушений при подмене абсолютных ридов: {len(violations)}")
    for v in violations:
        print("  ", v)

if __name__ == "__main__":
    main()
"""Calculate taxonomy comparison in a new directory; never overwrite a run.

Scientific calculations are pure functions in taxonomy_compare_core. Input
validation and output schemas are shared with the independent checker.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qiime_runtime
qiime_runtime.selected_python()

from taxonomy_compare_core import comparison_tables, soil_asv_ids
from taxonomy_compare_io import ROOT, RULES, SCHEMAS, code_identity, load_dataset, sha256

DEFAULT_OUT_DIR = ROOT / "results" / "taxonomy_comparison"


def fmt(value):
    """Round-trip floats; distinguish an undefined fraction from zero."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Non-finite output value")
        return format(value, ".17g")
    return str(value)


def write_tsv(path, columns, rows):
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            if set(row) != set(columns):
                raise ValueError(f"Output schema mismatch for {path.name}")
            writer.writerow({key: fmt(row[key]) for key in columns})


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def output_readme(dataset):
    n_soil_asv = len(soil_asv_ids(dataset["table"], dataset["soil_samples"]))
    return f"""# Сравнение таксономических профилей

Вход: {len(dataset['table'])} ASV, {len(dataset['samples'])} образцов.
Почвенные показатели: {len(dataset['soil_samples'])} образцов; K_S63 — отрицательный контроль.
Невзвешенные показатели: {n_soil_asv} ASV с положительной суммой частот в почвах.

- assignments.tsv: исходные строки Taxon всех ASV, метки рангов, ключи родословной,
  признаки назначения (1/0) и принадлежности почвенному набору (soil_asv).
- coverage.tsv: числа и доли назначенных/неназначенных ASV по общему почвенному набору.
- coverage_per_sample.tsv: числа и доли назначенных/неназначенных ридов каждого образца.
- pairwise_match.tsv: пять категорий по ASV; явный знаменатель both_assigned_asv.
- pairwise_match_per_sample.tsv: пять категорий по ридам каждого образца;
  явный знаменатель both_assigned_reads.
- profiles.tsv: положительные частоты и доли, включая UNCLASSIFIED при положительной массе.
- distances.tsv: половина L1 по объединению ключей, включая UNCLASSIFIED;
  рядом указано покрытие ридов каждой модели.
- summary.tsv: минимум, медиана и максимум расстояний по почвенным образцам.
- manifest.tsv: относительные пути, SHA256, UUID, семантические типы и версия кода.
- run.json: фактическое время UTC, идентичность кода, входы и хеши всех остальных выходов.

total_asv — число ASV общего почвенного набора; total_reads — все риды конкретного
почвенного образца после DADA2. Числитель условной доли — both_match или
both_match_reads; знаменатель — both_assigned_asv или both_assigned_reads.
При нулевом условном знаменателе записано NA. Доли и расстояния безразмерны.
Пустая метка ранга и ключ UNCLASSIFIED означают неизвестное назначение на этом ранге.
Ключ сохраняет известную родословную до целевого ранга; различия родословной и
суффиксы сохраняются, синонимы не объединяются. Сопоставление буквальных назначений
не устанавливает правильность классификации или превосходство референсной базы.

Расчёт не фильтрует ASV/образцы, не делает рарефикацию и не объединяет риды разных почв.
Код возврата сценария означает успешное создание расчёта; независимый вердикт выносит
src/check_taxonomy_comparison.py с теми же --input-dir и --out-dir.
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description="Сравнение таксономий в новом выходном каталоге")
    parser.add_argument("--input-dir", type=Path, default=ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--code-version", default=None,
                        help="необязательно: должна совпадать с фактической идентичностью кода")
    args = parser.parse_args(argv)
    input_dir = args.input_dir.resolve()
    output_dir = args.out_dir.absolute()
    if output_dir.exists() or output_dir.is_symlink():
        print(f"ОТКАЗ: выходной путь уже существует: {output_dir}", file=sys.stderr)
        return 2
    started = utc_now()
    try:
        identity = code_identity()
        if args.code_version is not None and args.code_version != identity["version"]:
            print("FAIL: CODE_VERSION — --code-version не совпадает с фактическим кодом", file=sys.stderr)
            return 1
        dataset, issues = load_dataset(input_dir)
        if issues:
            failed = any(not problem.get("missing", False) for problem in issues)
            verdict = "FAIL" if failed else "INCOMPLETE"
            for problem in issues:
                print(f"{verdict}: {problem['code']} — {problem['detail']}", file=sys.stderr)
            return 1 if failed else 2
        tables = comparison_tables(dataset["table"], dataset["taxonomies"], dataset["soil_samples"])
        tables["manifest.tsv"] = [
            {"item": item, **record, "code_version": identity["version"]}
            for item, record in sorted(dataset["inputs"].items())
        ]
        # Validate all rows before creating any output.
        if set(tables) != set(SCHEMAS):
            raise ValueError("Output file set differs from declared schemas")
        for name, rows in tables.items():
            for row in rows:
                if set(row) != set(SCHEMAS[name]):
                    raise ValueError(f"Output columns differ from schema: {name}")
                for value in row.values():
                    fmt(value)
        output_dir.mkdir(parents=True, exist_ok=False)
        for name, columns in SCHEMAS.items():
            write_tsv(output_dir / name, columns, tables[name])
        (output_dir / "normalization_rules.txt").write_text(RULES, encoding="utf-8")
        (output_dir / "README.md").write_text(output_readme(dataset), encoding="utf-8")
        names = sorted([*SCHEMAS, "README.md", "normalization_rules.txt"])
        record = {
            "schema_version": 1, "started_at_utc": started, "finished_at_utc": utc_now(),
            "input_dir": str(input_dir), "out_dir": str(output_dir),
            "code": identity, "inputs": dataset["inputs"],
            "outputs": {name: sha256(output_dir / name) for name in names},
        }
        with (output_dir / "run.json").open("x", encoding="utf-8") as handle:
            json.dump(record, handle, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    print(f"Расчёт создан: {output_dir}")
    print("Для вердикта запустите независимую проверку check_taxonomy_comparison.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

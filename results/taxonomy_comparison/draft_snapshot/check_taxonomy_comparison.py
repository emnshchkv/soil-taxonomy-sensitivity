"""Команда проверки сравнения таксономических профилей.

Контролирует точные наборы ID, отсутствие дубликатов и потерянных ASV,
десять почвенных образцов, целочисленность/неотрицательность частот,
сохранение сумм, полноту рангов и пар, допустимые диапазоны метрик
и отсутствие неявных нулей вместо пропусков.

Коды: 0 — PASS, 1 — FAIL, 2 — INCOMPLETE.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import biom

from taxonomy_compare_core import MAIN_RANKS, MODELS, PAIRS, SEPARATORS, parse_taxonomy_rows

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "taxonomy_comparison"

TOL = 1e-12

TAXONOMY_FILES = {
    "sklearn_silva138": "work/qza/taxonomy_sklearn_silva138.qza",
    "sklearn_gtdb": "work/qza/taxonomy_sklearn_gtdb.qza",
    "sklearn_gg2": "work/qza/taxonomy_sklearn_gg2.qza",
}


class Violation:
    def __init__(self, code: str, detail: str):
        self.code = code
        self.detail = detail

    def __str__(self) -> str:
        return f"{self.code}: {self.detail}"


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def check_required_files(violations: List[Violation]) -> bool:
    required = [
        "assignments.tsv", "normalization_rules.txt", "coverage.tsv",
        "pairwise_match.tsv", "profiles.tsv", "distances.tsv",
        "summary.tsv", "manifest.tsv", "README.md",
    ]
    ok = True
    for name in required:
        if not (OUT_DIR / name).exists():
            violations.append(Violation("MISSING_OUTPUT", name))
            ok = False
    return ok


def check_inputs(violations: List[Violation]) -> Tuple[Dict[str, Dict[str, int]], List[str], List[str], Dict[str, str]]:
    """Проверяет входы и возвращает (table, feature_ids, sample_ids, metadata)."""
    import h5py
    import io
    with open(ROOT / "work/qza/ASV_table.qza", "rb") as fh:
        data = fh.read()
    with h5py.File(io.BytesIO(data), "r") as f:
        t = biom.Table.from_hdf5(f)
    table: Dict[str, Dict[str, int]] = {}
    feature_ids = list(t.ids(axis="observation"))
    sample_ids = list(t.ids(axis="sample"))
    for fid, counts in zip(feature_ids, t.iter_data(axis="observation")):
        row = {sid: int(c) for sid, c in zip(sample_ids, counts)}
        table[fid] = row

    # Целочисленность/неотрицательность частот.
    for fid, row in table.items():
        for sid, c in row.items():
            if c < 0:
                violations.append(Violation("NEGATIVE_COUNT", f"{fid} {sid} {c}"))

    # Метаданные.
    metadata: Dict[str, str] = {}
    with open(ROOT / "data/meta/metadata.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            metadata[r["sample-id"]] = r["sample_type"]

    soil = [s for s in sample_ids if metadata.get(s) == "soil"]
    if len(soil) != 10:
        violations.append(Violation("SOIL_SAMPLES", f"ожидалось 10, получено {len(soil)}"))

    # Таксономии: наборы ID совпадают, нет дубликатов.
    for model in MODELS:
        rows: List[Tuple[str, str]] = []
        with open(ROOT / TAXONOMY_FILES[model]) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                rows.append((r["Feature ID"], r["Taxon"]))
        ids = [x[0] for x in rows]
        if len(set(ids)) != len(ids):
            violations.append(Violation("DUP_TAXON_ID", model))
        if set(ids) != set(feature_ids):
            violations.append(Violation("TAXON_ID_MISMATCH", model))
        # Разбор не должен падать.
        try:
            parse_taxonomy_rows(rows, SEPARATORS[model])
        except ValueError as e:
            violations.append(Violation("TAXON_PARSE", f"{model}: {e}"))

    return table, feature_ids, sample_ids, metadata


def check_assignments(violations: List[Violation], feature_ids: List[str]) -> None:
    rows = read_tsv(OUT_DIR / "assignments.tsv")
    ids = [r["feature_id"] for r in rows]
    if len(ids) != len(set(ids)):
        violations.append(Violation("DUP_ASSIGNMENT_ID", "assignments.tsv"))
    if set(ids) != set(feature_ids):
        violations.append(Violation("ASSIGNMENT_ID_MISMATCH", "assignments.tsv"))
    # Полнота рангов.
    for r in rows:
        for model in MODELS:
            for rank in MAIN_RANKS:
                col = f"{model}_{rank}"
                if col not in r:
                    violations.append(Violation("MISSING_RANK_COL", col))


def check_coverage(violations: List[Violation]) -> None:
    rows = read_tsv(OUT_DIR / "coverage.tsv")
    seen = set()
    for r in rows:
        key = (r["model"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_COVERAGE", str(key)))
        seen.add(key)
        if r["model"] not in MODELS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_COVERAGE_KEY", str(key)))
        # Доли в [0,1], суммы долей = 1.
        a = float(r["assigned_asv_frac"]); u = float(r["unassigned_asv_frac"])
        if not (0 - TOL <= a <= 1 + TOL and 0 - TOL <= u <= 1 + TOL):
            violations.append(Violation("COVERAGE_FRAC_RANGE", str(key)))
        if abs(a + u - 1.0) > TOL:
            violations.append(Violation("COVERAGE_FRAC_SUM", str(key)))
        ar = float(r["assigned_reads_frac"]); ur = float(r["unassigned_reads_frac"])
        if not (0 - TOL <= ar <= 1 + TOL and 0 - TOL <= ur <= 1 + TOL):
            violations.append(Violation("COVERAGE_READS_FRAC_RANGE", str(key)))
        if abs(ar + ur - 1.0) > TOL:
            violations.append(Violation("COVERAGE_READS_FRAC_SUM", str(key)))
        # Сумма ASV и ридов согласована.
        if int(r["assigned_asv"]) + int(r["unassigned_asv"]) != int(r["total_asv"]):
            violations.append(Violation("COVERAGE_ASV_SUM", str(key)))
        if int(r["assigned_reads"]) + int(r["unassigned_reads"]) != int(r["total_reads"]):
            violations.append(Violation("COVERAGE_READS_SUM", str(key)))


def check_pairwise(violations: List[Violation]) -> None:
    rows = read_tsv(OUT_DIR / "pairwise_match.tsv")
    seen = set()
    for r in rows:
        key = (r["model_a"], r["model_b"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_PAIRWISE", str(key)))
        seen.add(key)
        if (r["model_a"], r["model_b"]) not in PAIRS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_PAIRWISE_KEY", str(key)))
        # Категории покрывают весь набор.
        n = (int(r["both_match"]) + int(r["both_diff"]) + int(r["only_a"])
             + int(r["only_b"]) + int(r["neither"]))
        if n != int(r["total_asv"]):
            violations.append(Violation("PAIRWISE_COVER", str(key)))
        # Доли в [0,1], сумма = 1.
        fracs = [float(r[f"{c}_frac"]) for c in
                 ("both_match", "both_diff", "only_a", "only_b", "neither")]
        if any(not (0 - TOL <= x <= 1 + TOL) for x in fracs):
            violations.append(Violation("PAIRWISE_FRAC_RANGE", str(key)))
        if abs(sum(fracs) - 1.0) > TOL:
            violations.append(Violation("PAIRWISE_FRAC_SUM", str(key)))
        # Взвешенные доли.
        rfracs = [float(r[f"{c}_reads_frac"]) for c in
                  ("both_match", "both_diff", "only_a", "only_b", "neither")]
        if any(not (0 - TOL <= x <= 1 + TOL) for x in rfracs):
            violations.append(Violation("PAIRWISE_READS_FRAC_RANGE", str(key)))
        if abs(sum(rfracs) - 1.0) > TOL:
            violations.append(Violation("PAIRWISE_READS_FRAC_SUM", str(key)))
        # Условная доля: NA при нулевом знаменателе, иначе в [0,1].
        both_assigned = int(r["both_assigned"])
        cond = r["cond_match_frac"]
        if both_assigned == 0:
            if cond != "NA":
                violations.append(Violation("COND_NOT_NA", str(key)))
        else:
            if cond == "NA":
                violations.append(Violation("COND_NA_WITH_DENOM", str(key)))
            else:
                cv = float(cond)
                if not (0 - TOL <= cv <= 1 + TOL):
                    violations.append(Violation("COND_RANGE", str(key)))
                if abs(cv - int(r["both_match"]) / both_assigned) > TOL:
                    violations.append(Violation("COND_VALUE", str(key)))


def check_profiles(violations: List[Violation], table: Dict[str, Dict[str, int]], soil_samples: List[str]) -> None:
    rows = read_tsv(OUT_DIR / "profiles.tsv")
    # Сумма частот по (model, sample, rank) = total_reads образца.
    sums: Dict[Tuple[str, str, str], int] = {}
    totals: Dict[Tuple[str, str, str], int] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        sums[key] = sums.get(key, 0) + int(r["count"])
        totals[key] = int(r["total_reads"])
        if r["model"] not in MODELS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_PROFILE_KEY", str(key)))
        if r["sample"] not in soil_samples:
            violations.append(Violation("BAD_PROFILE_SAMPLE", str(key)))
        # Доли в [0,1].
        f = float(r["frac"])
        if not (0 - TOL <= f <= 1 + TOL):
            violations.append(Violation("PROFILE_FRAC_RANGE", str(key)))
    # Полнота: каждый (model, sample, rank) должен присутствовать.
    for model in MODELS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                key = (model, sample, rank)
                if key not in sums:
                    violations.append(Violation("MISSING_PROFILE", str(key)))
                    continue
                # Сумма частот = total_reads образца (сохранение массы).
                if sums[key] != totals[key]:
                    violations.append(Violation("PROFILE_MASS", str(key)))
                # Сумма долей = 1.
                # (пересчитываем по строкам)
    # Сумма долей = 1 по каждому ключу.
    frac_sums: Dict[Tuple[str, str, str], float] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        frac_sums[key] = frac_sums.get(key, 0.0) + float(r["frac"])
    for key, s in frac_sums.items():
        if abs(s - 1.0) > TOL:
            violations.append(Violation("PROFILE_FRAC_SUM", str(key)))


def check_distances(violations: List[Violation], soil_samples: List[str]) -> None:
    rows = read_tsv(OUT_DIR / "distances.tsv")
    seen = set()
    for r in rows:
        key = (r["model_a"], r["model_b"], r["sample"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_DISTANCE", str(key)))
        seen.add(key)
        if (r["model_a"], r["model_b"]) not in PAIRS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_DISTANCE_KEY", str(key)))
        if r["sample"] not in soil_samples:
            violations.append(Violation("BAD_DISTANCE_SAMPLE", str(key)))
        d = float(r["distance"])
        if not (0 - TOL <= d <= 1 + TOL):
            violations.append(Violation("DISTANCE_RANGE", str(key)))
    # Полнота пар/рангов/образцов.
    for a, b in PAIRS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                if (a, b, sample, rank) not in seen:
                    violations.append(Violation("MISSING_DISTANCE", f"{a} {b} {sample} {rank}"))


def check_summary(violations: List[Violation]) -> None:
    rows = read_tsv(OUT_DIR / "summary.tsv")
    seen = set()
    for r in rows:
        key = (r["model_a"], r["model_b"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_SUMMARY", str(key)))
        seen.add(key)
        if (r["model_a"], r["model_b"]) not in PAIRS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_SUMMARY_KEY", str(key)))
        if int(r["n_samples"]) != 10:
            violations.append(Violation("SUMMARY_N", str(key)))
        lo, med, hi = float(r["min"]), float(r["median"]), float(r["max"])
        if not (0 - TOL <= lo <= med <= hi <= 1 + TOL):
            violations.append(Violation("SUMMARY_ORDER", str(key)))


def check_manifest(violations: List[Violation]) -> None:
    rows = read_tsv(OUT_DIR / "manifest.tsv")
    items = {r["item"] for r in rows}
    for model in MODELS:
        if model not in items:
            violations.append(Violation("MISSING_MANIFEST", model))
    if "ASV_table" not in items or "metadata" not in items:
        violations.append(Violation("MISSING_MANIFEST", "inputs"))


def main() -> int:
    violations: List[Violation] = []

    if not check_required_files(violations):
        print("INCOMPLETE/2 — отсутствуют обязательные выходы")
        for v in violations:
            print(f"  {v}")
        return 2

    table, feature_ids, sample_ids, metadata = check_inputs(violations)
    soil_samples = [s for s in sample_ids if metadata.get(s) == "soil"]

    check_assignments(violations, feature_ids)
    check_coverage(violations)
    check_pairwise(violations)
    check_profiles(violations, table, soil_samples)
    check_distances(violations, soil_samples)
    check_summary(violations)
    check_manifest(violations)

    if violations:
        print("FAIL/1")
        for v in violations:
            print(f"  {v}")
        return 1

    print("PASS/0 — сравнение таксономических профилей корректно")
    return 0


if __name__ == "__main__":
    sys.exit(main())
"""Сценарий расчёта сравнения таксономических профилей.

Читает таблицу частот, метаданные и три таксономии, считает метрики
и пишет выходы в results/taxonomy_comparison/. Чтение файлов отделено
от чистых вычислительных функций (src/taxonomy_compare_core.py).
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import biom

from taxonomy_compare_core import (
    MAIN_RANKS,
    MODELS,
    PAIRS,
    ParsedTaxon,
    coverage,
    parse_taxonomy_rows,
    pairwise_match,
    profile,
    profile_distance,
    soil_asv_ids,
)

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "taxonomy_comparison"

# Разделители рангов по моделям (фактический формат).
SEPARATORS = {
    "sklearn_silva138": ";",
    "sklearn_gtdb": ";",
    "sklearn_gg2": "; ",
}

TAXONOMY_FILES = {
    "sklearn_silva138": "work/qza/taxonomy_sklearn_silva138.qza",
    "sklearn_gtdb": "work/qza/taxonomy_sklearn_gtdb.qza",
    "sklearn_gg2": "work/qza/taxonomy_sklearn_gg2.qza",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def read_biom(path: Path) -> Tuple[Dict[str, Dict[str, int]], List[str], List[str]]:
    """Читает BIOM таблицу (HDF5). Возвращает (table, feature_ids, sample_ids)."""
    import h5py
    import io
    with open(path, "rb") as fh:
        data = fh.read()
    with h5py.File(io.BytesIO(data), "r") as f:
        t = biom.Table.from_hdf5(f)
    table: Dict[str, Dict[str, int]] = {}
    for fid, counts in zip(t.ids(axis="observation"), t.iter_data(axis="observation")):
        row = {sid: int(c) for sid, c in zip(t.ids(axis="sample"), counts)}
        table[fid] = row
    return table, list(t.ids(axis="observation")), list(t.ids(axis="sample"))


def read_metadata(path: Path) -> Dict[str, str]:
    """Читает metadata.tsv. Возвращает {sample_id: sample_type}."""
    out: Dict[str, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[row["sample-id"]] = row["sample_type"]
    return out


def read_taxonomy(path: Path, separator: str) -> Dict[str, ParsedTaxon]:
    """Читает taxonomy.tsv (Feature ID, Taxon). Возвращает словарь."""
    rows: List[Tuple[str, str]] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append((row["Feature ID"], row["Taxon"]))
    return parse_taxonomy_rows(rows, separator)


def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Чтение входов ---
    table, feature_ids, sample_ids = read_biom(ROOT / "work/qza/ASV_table.qza")
    metadata = read_metadata(ROOT / "data/meta/metadata.tsv")

    soil_samples = [s for s in sample_ids if metadata.get(s) == "soil"]
    if len(soil_samples) != 10:
        raise SystemExit(f"Ожидалось 10 почвенных образцов, получено {len(soil_samples)}")

    taxonomies: Dict[str, Dict[str, ParsedTaxon]] = {}
    for model in MODELS:
        taxonomies[model] = read_taxonomy(
            ROOT / TAXONOMY_FILES[model], SEPARATORS[model]
        )

    # Проверка: наборы ID совпадают с таблицей.
    for model, tax in taxonomies.items():
        if set(tax) != set(feature_ids):
            raise SystemExit(f"Набор ID {model} не совпадает с таблицей ASV")

    # --- Почвенный набор ASV (ненулевая сумма в почвенных образцах) ---
    soil_asv = soil_asv_ids(table, soil_samples)
    soil_asv_set = set(soil_asv)

    # --- Выходы ---

    # 1. Назначения по ASV и рангам (все 3374 ID).
    assign_rows: List[List[object]] = []
    for fid in feature_ids:
        row = [fid]
        for model in MODELS:
            pt = taxonomies[model][fid]
            row.append(pt.raw)
            for rank in MAIN_RANKS:
                label = pt.ranks.get(rank)
                row.append(label if label is not None else "")
        assign_rows.append(row)
    assign_header = ["feature_id"] + [
        f"{model}_taxon" for model in MODELS
    ] + [
        f"{model}_{rank}" for model in MODELS for rank in MAIN_RANKS
    ]
    write_tsv(OUT_DIR / "assignments.tsv", assign_header, assign_rows)

    # 2. Правила нормализации/сопоставления — статический текст.
    rules_text = (
        "Нормализация только техническая: пробелы, разделители, обозначения рангов.\n"
        "Значимые суффиксы GTDB, названия кластеров и различия родословной сохраняются.\n"
        "Служебные метки (Unassigned, Unknown_Family, пустые суффиксы) не считаются известными таксонами.\n"
        "Синонимы и переименования не объединяются без проверяемого соответствия источнику.\n"
        "Основной анализ — буквальные нормализованные назначения с учётом родословной.\n"
    )
    (OUT_DIR / "normalization_rules.txt").write_text(rules_text)

    # 3. Покрытие классификацией.
    cov = coverage(taxonomies, soil_asv, table, soil_samples)
    cov_rows: List[List[object]] = []
    for model in MODELS:
        for rank in MAIN_RANKS:
            c = cov[model][rank]
            cov_rows.append([
                model, rank,
                c["assigned_asv"], c["unassigned_asv"], c["total_asv"],
                f"{c['assigned_asv_frac']:.12f}", f"{c['unassigned_asv_frac']:.12f}",
                c["assigned_reads"], c["unassigned_reads"], c["total_reads"],
                f"{c['assigned_reads_frac']:.12f}", f"{c['unassigned_reads_frac']:.12f}",
            ])
    write_tsv(
        OUT_DIR / "coverage.tsv",
        ["model", "rank", "assigned_asv", "unassigned_asv", "total_asv",
         "assigned_asv_frac", "unassigned_asv_frac",
         "assigned_reads", "unassigned_reads", "total_reads",
         "assigned_reads_frac", "unassigned_reads_frac"],
        cov_rows,
    )

    # 4. Пять категорий попарного совпадения.
    match_rows: List[List[object]] = []
    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            m = pairwise_match(taxonomies[a], taxonomies[b], soil_asv, table, soil_samples, rank)
            cond = "NA" if m["cond_match_frac"] is None else f"{m['cond_match_frac']:.12f}"
            match_rows.append([
                a, b, rank,
                m["both_match"], m["both_diff"], m["only_a"], m["only_b"], m["neither"],
                m["total_asv"],
                f"{m['both_match_frac']:.12f}", f"{m['both_diff_frac']:.12f}",
                f"{m['only_a_frac']:.12f}", f"{m['only_b_frac']:.12f}",
                f"{m['neither_frac']:.12f}",
                m["both_match_reads"], m["both_diff_reads"], m["only_a_reads"],
                m["only_b_reads"], m["neither_reads"], m["total_reads"],
                f"{m['both_match_reads_frac']:.12f}", f"{m['both_diff_reads_frac']:.12f}",
                f"{m['only_a_reads_frac']:.12f}", f"{m['only_b_reads_frac']:.12f}",
                f"{m['neither_reads_frac']:.12f}",
                m["both_assigned"], cond,
            ])
    write_tsv(
        OUT_DIR / "pairwise_match.tsv",
        ["model_a", "model_b", "rank",
         "both_match", "both_diff", "only_a", "only_b", "neither", "total_asv",
         "both_match_frac", "both_diff_frac", "only_a_frac", "only_b_frac", "neither_frac",
         "both_match_reads", "both_diff_reads", "only_a_reads", "only_b_reads",
         "neither_reads", "total_reads",
         "both_match_reads_frac", "both_diff_reads_frac", "only_a_reads_frac",
         "only_b_reads_frac", "neither_reads_frac",
         "both_assigned", "cond_match_frac"],
        match_rows,
    )

    # 5. Частоты и доли профилей по образцу/рангу/модели.
    prof_rows: List[List[object]] = []
    for model in MODELS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                p, total = profile(taxonomies[model], table, sample, rank)
                for key in sorted(p):
                    prof_rows.append([
                        model, sample, rank, key, p[key], total,
                        f"{p[key] / total:.12f}" if total else "0.000000000000",
                    ])
    write_tsv(
        OUT_DIR / "profiles.tsv",
        ["model", "sample", "rank", "key", "count", "total_reads", "frac"],
        prof_rows,
    )

    # 6. Расстояния профилей по образцу/рангу/паре.
    dist_rows: List[List[object]] = []
    for a, b in PAIRS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                pa, _ = profile(taxonomies[a], table, sample, rank)
                pb, _ = profile(taxonomies[b], table, sample, rank)
                d = profile_distance(pa, pb)
                dist_rows.append([a, b, sample, rank, f"{d:.12f}"])
    write_tsv(
        OUT_DIR / "distances.tsv",
        ["model_a", "model_b", "sample", "rank", "distance"],
        dist_rows,
    )

    # 7. Сводная таблица: медиана и диапазон расстояний по десяти почвам.
    summary_rows: List[List[object]] = []
    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            vals = [
                float(r[4]) for r in dist_rows
                if r[0] == a and r[1] == b and r[3] == rank
            ]
            vals_sorted = sorted(vals)
            n = len(vals_sorted)
            med = (
                vals_sorted[n // 2]
                if n % 2 == 1
                else (vals_sorted[n // 2 - 1] + vals_sorted[n // 2]) / 2
            )
            summary_rows.append([
                a, b, rank, n,
                f"{min(vals):.12f}", f"{med:.12f}", f"{max(vals):.12f}",
            ])
    write_tsv(
        OUT_DIR / "summary.tsv",
        ["model_a", "model_b", "rank", "n_samples", "min", "median", "max"],
        summary_rows,
    )

    # 8. Манифест входов с UUID/SHA256 и версией кода.
    manifest_rows: List[List[object]] = []
    for model in MODELS:
        p = ROOT / TAXONOMY_FILES[model]
        manifest_rows.append([model, str(p), sha256(p)])
    manifest_rows.append(["ASV_table", str(ROOT / "work/qza/ASV_table.qza"), sha256(ROOT / "work/qza/ASV_table.qza")])
    manifest_rows.append(["metadata", str(ROOT / "data/meta/metadata.tsv"), sha256(ROOT / "data/meta/metadata.tsv")])
    write_tsv(
        OUT_DIR / "manifest.tsv",
        ["item", "path", "sha256"],
        manifest_rows,
    )

    # 9. README выходного каталога.
    readme = f"""# results/taxonomy_comparison — сравнение таксономических профилей

Этап 4. Вход: 3374 ASV, 11 образцов; почвенные профили — 10 образцов (sample_type=soil).
K_S63 (negative_control) в почвенные профили не входит. Kr5_S84 сохранён (1505/60078, 2.51%).

Почвенный набор ASV для невзвешенных показателей: {len(soil_asv)} ASV с ненулевой суммой
частот в десяти почвенных образцах (знаменатель записан в каждом файле как total_asv).

## Файлы
- assignments.tsv — назначения по всем {len(feature_ids)} ASV и рангам (тип, класс, порядок, семейство, род).
- normalization_rules.txt — правила нормализации и сопоставления.
- coverage.tsv — покрытие классификацией: доли ASV и доли ридов (не смешиваются).
- pairwise_match.tsv — пять категорий попарного совпадения назначений (невзвешенные и взвешенные доли).
- profiles.tsv — частоты и доли профилей по образцу/рангу/модели; UNCLASSIFIED — неназначенная масса.
- distances.tsv — попарное расстояние профилей D = 0.5*sum|p-q| по образцу/рангу/паре.
- summary.tsv — медиана и диапазон расстояний по десяти почвам.
- manifest.tsv — SHA256 входов.

## Знаменатели и единицы
- total_asv — число ASV в почвенном наборе (знаменатель невзвешенных долей).
- total_reads — сумма ридов по почвенным образцам (знаменатель взвешенных долей).
- total_reads в profiles.tsv — общее число ридов образца после DADA2.
- distance — доля, диапазон [0,1].
- cond_match_frac — условная доля совпадений среди назначенных обеими; NA при нулевом знаменателе.
"""
    (OUT_DIR / "README.md").write_text(readme)

    print(f"Расчёт завершён. Выходы в {OUT_DIR}")
    print(f"Почвенных образцов: {len(soil_samples)}; почвенных ASV: {len(soil_asv)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
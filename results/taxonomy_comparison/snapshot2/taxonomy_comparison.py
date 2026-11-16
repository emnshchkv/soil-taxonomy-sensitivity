"""Сценарий расчёта сравнения таксономических профилей.

Читает таблицу частот, метаданные и три таксономии, считает метрики
и пишет выходы в results/taxonomy_comparison/. Чтение файлов отделено
от чистых вычислительных функций (src/taxonomy_compare_core.py).

Исправления 2026-11-11 (по поручению):
- QZA читаются по ZIP-структуре (data/feature-table.biom, data/taxonomy.tsv,
  data/dna-sequences.fasta) и metadata.yaml, а не прямым открытием как HDF5/TSV;
- manifest.tsv содержит UUID и версию кода, включает rep_seq.qza и VERSIONS.md;
- взвешенные показатели считаются отдельно для каждого почвенного образца;
- скрипт принимает отдельные пути входов/выходов и отказывается затирать
  существующий расчёт.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
import os
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

# Выбор Python до сторонних импортов.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import qiime_runtime  # noqa: E402
qiime_runtime.selected_python()  # noqa: E402

import biom  # noqa: E402

from taxonomy_compare_core import (  # noqa: E402
    MAIN_RANKS,
    MODELS,
    PAIRS,
    SEPARATORS,
    ParsedTaxon,
    coverage,
    coverage_per_sample,
    parse_taxonomy_rows,
    pairwise_match,
    pairwise_match_per_sample,
    profile,
    profile_distance,
    soil_asv_ids,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT_DIR = ROOT / "results" / "taxonomy_comparison"

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


def qza_uuid(path: Path) -> str:
    """UUID из metadata.yaml внутри QZA (ZIP)."""
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith("metadata.yaml")]
        # Корневой metadata.yaml — тот, что не в provenance/.
        root_meta = [n for n in names if "/provenance/" not in n]
        if not root_meta:
            raise ValueError(f"Нет корневого metadata.yaml в {path}")
        meta = z.read(root_meta[0]).decode("utf-8")
    for line in meta.splitlines():
        if line.startswith("uuid:"):
            return line.split(":", 1)[1].strip()
    raise ValueError(f"Нет uuid в metadata.yaml {path}")


def qza_semantic_type(path: Path) -> str:
    """Семантический тип из корневого metadata.yaml внутри QZA (ZIP)."""
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith("metadata.yaml")]
        root_meta = [n for n in names if "/provenance/" not in n]
        if not root_meta:
            raise ValueError(f"Нет корневого metadata.yaml в {path}")
        meta = z.read(root_meta[0]).decode("utf-8")
    for line in meta.splitlines():
        if line.startswith("type:"):
            return line.split(":", 1)[1].strip()
    raise ValueError(f"Нет type в metadata.yaml {path}")


def qza_archive_root(path: Path) -> str:
    """Корень архива QZA (общий префикс путей внутри ZIP)."""
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
    if not names:
        raise ValueError(f"Пустой архив {path}")
    # Корень — общий префикс до первого '/'.
    roots = {n.split("/")[0] for n in names}
    if len(roots) != 1:
        raise ValueError(f"Неоднозначный корень архива {path}: {roots}")
    return roots.pop()


def read_biom_from_qza(path: Path) -> Tuple[Dict[str, Dict[str, int]], List[str], List[str]]:
    """Читает BIOM таблицу из data/feature-table.biom внутри QZA (ZIP)."""
    import tempfile
    with zipfile.ZipFile(path) as z:
        biom_names = [n for n in z.namelist() if n.endswith("data/feature-table.biom")]
        if not biom_names:
            raise ValueError(f"Нет data/feature-table.biom в {path}")
        data = z.read(biom_names[0])
    # biom_open не принимает BytesIO; пишем во временный файл.
    with tempfile.NamedTemporaryFile(suffix=".biom", delete=False) as tf:
        tf.write(data)
        tmp_path = tf.name
    try:
        with biom.util.biom_open(tmp_path) as f:
            t = biom.Table.from_hdf5(f)
    finally:
        os.unlink(tmp_path)
    table: Dict[str, Dict[str, int]] = {}
    feature_ids = list(t.ids(axis="observation"))
    sample_ids = list(t.ids(axis="sample"))
    for fid, counts in zip(feature_ids, t.iter_data(axis="observation")):
        row = {}
        for sid, c in zip(sample_ids, counts):
            # Проверка конечности/неотрицательности/целочисленности ДО int(c).
            if isinstance(c, float) and not math.isfinite(c):
                raise ValueError(f"Не конечная частота {fid} {sid}: {c}")
            if c < 0:
                raise ValueError(f"Отрицательная частота {fid} {sid}: {c}")
            if not float(c).is_integer():
                raise ValueError(f"Дробная частота {fid} {sid}: {c}")
            row[sid] = int(c)
        table[fid] = row
    return table, feature_ids, sample_ids


def read_taxonomy_from_qza(path: Path, separator: str) -> Dict[str, ParsedTaxon]:
    """Читает taxonomy.tsv из data/taxonomy.tsv внутри QZA (ZIP)."""
    with zipfile.ZipFile(path) as z:
        tax_names = [n for n in z.namelist() if n.endswith("data/taxonomy.tsv")]
        if not tax_names:
            raise ValueError(f"Нет data/taxonomy.tsv в {path}")
        text = z.read(tax_names[0]).decode("utf-8")
    rows: List[Tuple[str, str]] = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        rows.append((row["Feature ID"], row["Taxon"]))
    return parse_taxonomy_rows(rows, separator)


def read_repseq_ids(path: Path) -> List[str]:
    """Читает ID последовательностей из data/dna-sequences.fasta внутри QZA (ZIP)."""
    with zipfile.ZipFile(path) as z:
        fasta_names = [n for n in z.namelist() if n.endswith("data/dna-sequences.fasta")]
        if not fasta_names:
            raise ValueError(f"Нет data/dna-sequences.fasta в {path}")
        text = z.read(fasta_names[0]).decode("utf-8")
    ids = []
    for line in text.splitlines():
        if line.startswith(">"):
            ids.append(line[1:].strip())
    return ids


def read_metadata(path: Path) -> Dict[str, str]:
    """Читает metadata.tsv. Возвращает {sample_id: sample_type}."""
    out: Dict[str, str] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[row["sample-id"]] = row["sample_type"]
    return out


def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def fmt(x: float) -> str:
    """Формат с достаточной точностью для допуска 1e-12."""
    return f"{x:.15g}"


def _detect_code_version() -> str:
    """Версия кода: текущий коммит git или SHA256 снимка кода.

    Если рабочий каталог чист и есть коммит — возвращает короткий хеш коммита.
    Иначе — SHA256 объединённого содержимого трёх скриптов (снимок кода).
    """
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=ROOT,
        )
        if r.returncode == 0:
            commit = r.stdout.strip()
            # Проверка чистоты рабочего каталога для этих файлов.
            st = subprocess.run(
                ["git", "status", "--porcelain", "src/taxonomy_comparison.py",
                 "src/taxonomy_compare_core.py", "src/check_taxonomy_comparison.py"],
                capture_output=True, text=True, cwd=ROOT,
            )
            if st.returncode == 0 and not st.stdout.strip():
                return commit
    except Exception:
        pass
    # Снимок кода: SHA256 объединённых скриптов.
    h = hashlib.sha256()
    for name in ["taxonomy_compare_core.py", "taxonomy_comparison.py",
                 "check_taxonomy_comparison.py"]:
        with open(ROOT / "src" / name, "rb") as f:
            h.update(f.read())
    return f"snapshot-{h.hexdigest()[:12]}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Расчёт сравнения таксономических профилей")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--force", action="store_true",
                    help="перезаписать существующий расчёт (по умолчанию отказ)")
    ap.add_argument("--code-version", default=None,
                    help="версия кода для manifest (коммит или SHA256 снимка)")
    args = ap.parse_args()

    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()) and not args.force:
        print(f"ОТКАЗ: каталог {out_dir} не пуст. Используйте --force для перезаписи.", file=sys.stderr)
        return 2
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Чтение входов ---
    table, feature_ids, sample_ids = read_biom_from_qza(ROOT / "work/qza/ASV_table.qza")
    metadata = read_metadata(ROOT / "data/meta/metadata.tsv")

    soil_samples = [s for s in sample_ids if metadata.get(s) == "soil"]
    if len(soil_samples) != 10:
        raise SystemExit(f"Ожидалось 10 почвенных образцов, получено {len(soil_samples)}")

    taxonomies: Dict[str, Dict[str, ParsedTaxon]] = {}
    for model in MODELS:
        taxonomies[model] = read_taxonomy_from_qza(
            ROOT / TAXONOMY_FILES[model], SEPARATORS[model]
        )

    # Проверка семантических типов QZA и корня архива.
    expected_types = {
        "ASV_table": "FeatureTable[Frequency]",
        "rep_seq": "FeatureData[Sequence]",
    }
    for model in MODELS:
        expected_types[model] = "FeatureData[Taxonomy]"
    for item, rel in [
        ("ASV_table", "work/qza/ASV_table.qza"),
        ("rep_seq", "work/qza/rep_seq.qza"),
    ] + [(m, TAXONOMY_FILES[m]) for m in MODELS]:
        p = ROOT / rel
        st = qza_semantic_type(p)
        if st != expected_types[item]:
            raise SystemExit(f"Семантический тип {item}: {st!r} != {expected_types[item]!r}")
        qza_archive_root(p)  # проверка однозначного корня архива

    # Проверка: наборы ID совпадают с таблицей.
    for model, tax in taxonomies.items():
        if set(tax) != set(feature_ids):
            raise SystemExit(f"Набор ID {model} не совпадает с таблицей ASV")

    # rep_seq ID совпадают с таблицей.
    repseq_ids = read_repseq_ids(ROOT / "work/qza/rep_seq.qza")
    if set(repseq_ids) != set(feature_ids):
        raise SystemExit("Набор ID rep_seq.qza не совпадает с таблицей ASV")

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
    # Заголовок соответствует порядку значений в строке: для каждой модели
    # сначала исходная строка Taxon, затем пять рангов.
    assign_header = ["feature_id"]
    for model in MODELS:
        assign_header.append(f"{model}_taxon")
        assign_header.extend(f"{model}_{rank}" for rank in MAIN_RANKS)
    write_tsv(out_dir / "assignments.tsv", assign_header, assign_rows)

    # 2. Правила нормализации/сопоставления — статический текст.
    rules_text = (
        "Нормализация только техническая: пробелы, разделители, обозначения рангов.\n"
        "Значимые суффиксы GTDB, названия кластеров и различия родословной сохраняются.\n"
        "Служебные метки (Unassigned, Unknown_Family, пустые суффиксы) не считаются известными таксонами.\n"
        "Синонимы и переименования не объединяются без проверяемого соответствия источнику.\n"
        "Основной анализ — буквальные нормализованные назначения с учётом родословной.\n"
        "Ключ родословной — от домена до целевого ранга; род под разными семействами — разные ключи.\n"
        "Дублированный ранг (p__P1;p__P2) — диагностика, а не молчаливый выбор.\n"
    )
    (out_dir / "normalization_rules.txt").write_text(rules_text)

    # 3. Покрытие классификацией.
    cov = coverage(taxonomies, soil_asv, table, soil_samples)
    cov_rows: List[List[object]] = []
    for model in MODELS:
        for rank in MAIN_RANKS:
            c = cov[model][rank]
            cov_rows.append([
                model, rank,
                c["assigned_asv"], c["unassigned_asv"], c["total_asv"],
                fmt(c["assigned_asv_frac"]), fmt(c["unassigned_asv_frac"]),
                c["assigned_reads"], c["unassigned_reads"], c["total_reads"],
                fmt(c["assigned_reads_frac"]), fmt(c["unassigned_reads_frac"]),
            ])
    write_tsv(
        out_dir / "coverage.tsv",
        ["model", "rank", "assigned_asv", "unassigned_asv", "total_asv",
         "assigned_asv_frac", "unassigned_asv_frac",
         "assigned_reads", "unassigned_reads", "total_reads",
         "assigned_reads_frac", "unassigned_reads_frac"],
        cov_rows,
    )

    # 3b. Покрытие классификацией по каждому почвенному образцу.
    cov_sample_rows: List[List[object]] = []
    for sample in soil_samples:
        cov_s = coverage_per_sample(taxonomies, soil_asv, table, sample)
        for model in MODELS:
            for rank in MAIN_RANKS:
                c = cov_s[model][rank]
                cov_sample_rows.append([
                    model, sample, rank,
                    c["assigned_asv"], c["unassigned_asv"], c["total_asv"],
                    fmt(c["assigned_asv_frac"]), fmt(c["unassigned_asv_frac"]),
                    c["assigned_reads"], c["unassigned_reads"], c["total_reads"],
                    fmt(c["assigned_reads_frac"]), fmt(c["unassigned_reads_frac"]),
                ])
    write_tsv(
        out_dir / "coverage_per_sample.tsv",
        ["model", "sample", "rank", "assigned_asv", "unassigned_asv", "total_asv",
         "assigned_asv_frac", "unassigned_asv_frac",
         "assigned_reads", "unassigned_reads", "total_reads",
         "assigned_reads_frac", "unassigned_reads_frac"],
        cov_sample_rows,
    )

    # 4. Пять категорий попарного совпадения (по почвенному набору ASV).
    match_rows: List[List[object]] = []
    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            m = pairwise_match(taxonomies[a], taxonomies[b], soil_asv, table, soil_samples, rank)
            cond = "NA" if m["cond_match_frac"] is None else fmt(m["cond_match_frac"])
            match_rows.append([
                a, b, rank,
                m["both_match"], m["both_diff"], m["only_a"], m["only_b"], m["neither"],
                m["total_asv"],
                fmt(m["both_match_frac"]), fmt(m["both_diff_frac"]),
                fmt(m["only_a_frac"]), fmt(m["only_b_frac"]),
                fmt(m["neither_frac"]),
                int(m["both_match_reads"]), int(m["both_diff_reads"]), int(m["only_a_reads"]),
                int(m["only_b_reads"]), int(m["neither_reads"]), int(m["total_reads"]),
                fmt(m["both_match_reads_frac"]), fmt(m["both_diff_reads_frac"]),
                fmt(m["only_a_reads_frac"]), fmt(m["only_b_reads_frac"]),
                fmt(m["neither_reads_frac"]),
                m["both_assigned"], cond,
            ])
    write_tsv(
        out_dir / "pairwise_match.tsv",
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

    # 4b. Попарное совпадение по каждому почвенному образцу (взвешенные по образцу).
    match_sample_rows: List[List[object]] = []
    for a, b in PAIRS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                m = pairwise_match_per_sample(
                    taxonomies[a], taxonomies[b], soil_asv, table, sample, rank
                )
                cond = "NA" if m["cond_match_frac"] is None else fmt(m["cond_match_frac"])
                match_sample_rows.append([
                    a, b, sample, rank,
                    m["both_match"], m["both_diff"], m["only_a"], m["only_b"], m["neither"],
                    m["total_asv"],
                    fmt(m["both_match_frac"]), fmt(m["both_diff_frac"]),
                    fmt(m["only_a_frac"]), fmt(m["only_b_frac"]),
                    fmt(m["neither_frac"]),
                    int(m["both_match_reads"]), int(m["both_diff_reads"]), int(m["only_a_reads"]),
                    int(m["only_b_reads"]), int(m["neither_reads"]), int(m["total_reads"]),
                    fmt(m["both_match_reads_frac"]), fmt(m["both_diff_reads_frac"]),
                    fmt(m["only_a_reads_frac"]), fmt(m["only_b_reads_frac"]),
                    fmt(m["neither_reads_frac"]),
                    m["both_assigned"], cond,
                ])
    write_tsv(
        out_dir / "pairwise_match_per_sample.tsv",
        ["model_a", "model_b", "sample", "rank",
         "both_match", "both_diff", "only_a", "only_b", "neither", "total_asv",
         "both_match_frac", "both_diff_frac", "only_a_frac", "only_b_frac", "neither_frac",
         "both_match_reads", "both_diff_reads", "only_a_reads", "only_b_reads",
         "neither_reads", "total_reads",
         "both_match_reads_frac", "both_diff_reads_frac", "only_a_reads_frac",
         "only_b_reads_frac", "neither_reads_frac",
         "both_assigned", "cond_match_frac"],
        match_sample_rows,
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
                        fmt(p[key] / total) if total else "0",
                    ])
    write_tsv(
        out_dir / "profiles.tsv",
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
                dist_rows.append([a, b, sample, rank, fmt(d)])
    write_tsv(
        out_dir / "distances.tsv",
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
                fmt(min(vals)), fmt(med), fmt(max(vals)),
            ])
    write_tsv(
        out_dir / "summary.tsv",
        ["model_a", "model_b", "rank", "n_samples", "min", "median", "max"],
        summary_rows,
    )

    # 8. Манифест входов с UUID/SHA256 и версией кода.
    code_version = args.code_version or _detect_code_version()
    manifest_rows: List[List[object]] = []
    for model in MODELS:
        p = ROOT / TAXONOMY_FILES[model]
        manifest_rows.append([model, str(p), sha256(p), qza_uuid(p), code_version])
    for item, rel in [
        ("ASV_table", "work/qza/ASV_table.qza"),
        ("rep_seq", "work/qza/rep_seq.qza"),
    ]:
        p = ROOT / rel
        manifest_rows.append([item, str(p), sha256(p), qza_uuid(p), code_version])
    for item, rel in [
        ("metadata", "data/meta/metadata.tsv"),
        ("VERSIONS", "data/refs/VERSIONS.md"),
    ]:
        p = ROOT / rel
        manifest_rows.append([item, str(p), sha256(p), "", code_version])
    write_tsv(
        out_dir / "manifest.tsv",
        ["item", "path", "sha256", "uuid", "code_version"],
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
- coverage.tsv — покрытие классификацией по почвенному набору ASV: доли ASV и доли ридов (не смешиваются).
- coverage_per_sample.tsv — покрытие классификацией по каждому почвенному образцу (доли ридов образца).
- pairwise_match.tsv — пять категорий попарного совпадения назначений по почвенному набору ASV
  (невзвешенные доли по ASV и взвешенные по сумме ридов почвенных образцов).
- pairwise_match_per_sample.tsv — те же категории по каждому почвенному образцу
  (взвешенные по частотам ридов конкретного образца).
- profiles.tsv — частоты и доли профилей по образцу/рангу/модели; UNCLASSIFIED — неназначенная масса.
- distances.tsv — попарное расстояние профилей D = 0.5*sum|p-q| по образцу/рангу/паре.
- summary.tsv — медиана и диапазон расстояний по десяти почвам.
- manifest.tsv — SHA256, UUID и версия кода для всех входов.

## Знаменатели и единицы
- total_asv — число ASV в почвенном наборе (знаменатель невзвешенных долей).
- total_reads в pairwise_match.tsv — сумма ридов по почвенным образцам (знаменатель взвешенных долей).
- total_reads в pairwise_match_per_sample.tsv — число ридов конкретного образца.
- total_reads в profiles.tsv — общее число ридов образца после DADA2.
- distance — доля, диапазон [0,1].
- cond_match_frac — условная доля совпадений среди назначенных обеими; NA при нулевом знаменателе.
- Ключи профилей — ключи родословной от домена до целевого ранга; UNCLASSIFIED — неназначенная масса.
"""
    (out_dir / "README.md").write_text(readme)

    print(f"Расчёт завершён. Выходы в {out_dir}")
    print(f"Почвенных образцов: {len(soil_samples)}; почвенных ASV: {len(soil_asv)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
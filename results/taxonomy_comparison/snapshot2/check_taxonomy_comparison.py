"""Команда проверки сравнения таксономических профилей.

Контролирует точные наборы ID, отсутствие дубликатов и потерянных ASV,
десять почвенных образцов, целочисленность/неотрицательность частот,
сохранение сумм, полноту рангов и пар, допустимые диапазоны метрик
и отсутствие неявных нулей вместо пропусков.

Коды: 0 — PASS, 1 — FAIL, 2 — INCOMPLETE.

Исправления 2026-11-11 (по поручению):
- импорт SEPARATORS из taxonomy_compare_core (был отсутствующий импорт);
- выбор Python через src/qiime_runtime.py до сторонних импортов;
- QZA читаются по ZIP-структуре, а не прямым открытием как HDF5/TSV;
- частоты проверяются на конечность/неотрицательность/целочисленность ДО int(c);
- дубликаты выявляются до преобразования в словари;
- проверяются ключи родословной, взвешенные показатели по образцам,
  SHA256/UUID с настоящими входами, отсутствие неявных нулей.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
import os
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
    lineage_key,
    parse_taxonomy_rows,
    pairwise_match,
    pairwise_match_per_sample,
    profile,
    profile_distance,
    soil_asv_ids,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT_DIR = ROOT / "results" / "taxonomy_comparison"

TOL = 1e-12


def _is_valid_code_version(v: str) -> bool:
    """Версия кода — фактический коммит git или снимок кода (snapshot-<sha>)."""
    if v.startswith("snapshot-") and len(v) == len("snapshot-") + 12:
        return True
    # Проверка, что это существующий коммит в репозитории.
    try:
        import subprocess
        r = subprocess.run(
            ["git", "cat-file", "-e", v + "^{commit}"],
            capture_output=True, cwd=ROOT,
        )
        return r.returncode == 0
    except Exception:
        return False

TAXONOMY_FILES = {
    "sklearn_silva138": "work/qza/taxonomy_sklearn_silva138.qza",
    "sklearn_gtdb": "work/qza/taxonomy_sklearn_gtdb.qza",
    "sklearn_gg2": "work/qza/taxonomy_sklearn_gg2.qza",
}

INPUT_FILES = {
    "ASV_table": "work/qza/ASV_table.qza",
    "rep_seq": "work/qza/rep_seq.qza",
    "metadata": "data/meta/metadata.tsv",
    "VERSIONS": "data/refs/VERSIONS.md",
}
for m in MODELS:
    INPUT_FILES[m] = TAXONOMY_FILES[m]


class Violation:
    def __init__(self, code: str, detail: str):
        self.code = code
        self.detail = detail

    def __str__(self) -> str:
        return f"{self.code}: {self.detail}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def qza_uuid(path: Path) -> str:
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith("metadata.yaml")]
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
    roots = {n.split("/")[0] for n in names}
    if len(roots) != 1:
        raise ValueError(f"Неоднозначный корень архива {path}: {roots}")
    return roots.pop()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def check_required_files(out_dir: Path, violations: List[Violation]) -> None:
    required = [
        "assignments.tsv", "normalization_rules.txt", "coverage.tsv",
        "coverage_per_sample.tsv",
        "pairwise_match.tsv", "pairwise_match_per_sample.tsv",
        "profiles.tsv", "distances.tsv",
        "summary.tsv", "manifest.tsv", "README.md",
    ]
    for name in required:
        if not (out_dir / name).exists():
            violations.append(Violation("MISSING_OUTPUT", name))


def read_biom_from_qza(path: Path) -> Tuple[Dict[str, Dict[str, int]], List[str], List[str]]:
    """Читает BIOM таблицу из data/feature-table.biom внутри QZA (ZIP)."""
    import tempfile
    with zipfile.ZipFile(path) as z:
        biom_names = [n for n in z.namelist() if n.endswith("data/feature-table.biom")]
        if not biom_names:
            raise ValueError(f"Нет data/feature-table.biom в {path}")
        data = z.read(biom_names[0])
    with tempfile.NamedTemporaryFile(suffix=".biom", delete=False) as tf:
        tf.write(data)
        tmp_path = tf.name
    try:
        with biom.util.biom_open(tmp_path) as f:
            t = biom.Table.from_hdf5(f)
    finally:
        os.unlink(tmp_path)
    feature_ids = list(t.ids(axis="observation"))
    sample_ids = list(t.ids(axis="sample"))
    table: Dict[str, Dict[str, int]] = {}
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


def check_inputs(violations: List[Violation]) -> Tuple[Dict[str, Dict[str, int]], List[str], List[str], Dict[str, str], Dict[str, Dict[str, ParsedTaxon]]]:
    """Проверяет входы и возвращает (table, feature_ids, sample_ids, metadata, taxonomies)."""
    try:
        table, feature_ids, sample_ids = read_biom_from_qza(ROOT / "work/qza/ASV_table.qza")
    except ValueError as e:
        violations.append(Violation("BAD_FREQUENCY", str(e)))
        raise

    # Дубликаты ID в таблице.
    if len(set(feature_ids)) != len(feature_ids):
        violations.append(Violation("DUP_TABLE_ID", "ASV_table"))
    if len(set(sample_ids)) != len(sample_ids):
        violations.append(Violation("DUP_SAMPLE_ID", "ASV_table"))

    # Метаданные.
    metadata: Dict[str, str] = {}
    with open(ROOT / "data/meta/metadata.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            metadata[r["sample-id"]] = r["sample_type"]

    soil = [s for s in sample_ids if metadata.get(s) == "soil"]
    if len(soil) != 10:
        violations.append(Violation("SOIL_SAMPLES", f"ожидалось 10, получено {len(soil)}"))
    # Все 11 образцов в метаданных.
    if set(sample_ids) != set(metadata):
        violations.append(Violation("SAMPLE_METADATA_MISMATCH",
                                    f"таблица {len(sample_ids)}, метаданные {len(metadata)}"))

    # rep_seq ID.
    try:
        with zipfile.ZipFile(ROOT / "work/qza/rep_seq.qza") as z:
            fasta_names = [n for n in z.namelist() if n.endswith("data/dna-sequences.fasta")]
            if not fasta_names:
                raise ValueError("Нет data/dna-sequences.fasta")
            text = z.read(fasta_names[0]).decode("utf-8")
        repseq_ids = [l[1:].strip() for l in text.splitlines() if l.startswith(">")]
        if len(set(repseq_ids)) != len(repseq_ids):
            violations.append(Violation("DUP_REPSEQ_ID", "rep_seq"))
        if set(repseq_ids) != set(feature_ids):
            violations.append(Violation("REPSEQ_ID_MISMATCH", "rep_seq"))
    except Exception as e:
        violations.append(Violation("REPSEQ_READ", str(e)))

    # Таксономии: наборы ID совпадают, нет дубликатов, разбор не падает.
    taxonomies: Dict[str, Dict[str, ParsedTaxon]] = {}
    for model in MODELS:
        try:
            with zipfile.ZipFile(ROOT / TAXONOMY_FILES[model]) as z:
                tax_names = [n for n in z.namelist() if n.endswith("data/taxonomy.tsv")]
                if not tax_names:
                    raise ValueError("Нет data/taxonomy.tsv")
                text = z.read(tax_names[0]).decode("utf-8")
            rows: List[Tuple[str, str]] = []
            for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
                rows.append((r["Feature ID"], r["Taxon"]))
            ids = [x[0] for x in rows]
            if len(set(ids)) != len(ids):
                violations.append(Violation("DUP_TAXON_ID", model))
            if set(ids) != set(feature_ids):
                violations.append(Violation("TAXON_ID_MISMATCH", model))
            try:
                taxonomies[model] = parse_taxonomy_rows(rows, SEPARATORS[model])
            except ValueError as e:
                violations.append(Violation("TAXON_PARSE", f"{model}: {e}"))
        except Exception as e:
            violations.append(Violation("TAXON_READ", f"{model}: {e}"))

    # Семантические типы QZA и корень архива.
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
        try:
            st = qza_semantic_type(p)
            if st != expected_types[item]:
                violations.append(Violation("QZA_TYPE", f"{item}: {st!r} != {expected_types[item]!r}"))
            qza_archive_root(p)
        except Exception as e:
            violations.append(Violation("QZA_STRUCTURE", f"{item}: {e}"))

    return table, feature_ids, sample_ids, metadata, taxonomies


def check_assignments(violations: List[Violation], feature_ids: List[str],
                      taxonomies: Dict[str, Dict[str, ParsedTaxon]], out_dir: Path) -> None:
    rows = read_tsv(out_dir / "assignments.tsv")
    ids = [r["feature_id"] for r in rows]
    if len(ids) != len(set(ids)):
        violations.append(Violation("DUP_ASSIGNMENT_ID", "assignments.tsv"))
    if set(ids) != set(feature_ids):
        violations.append(Violation("ASSIGNMENT_ID_MISMATCH", "assignments.tsv"))
    # Полнота рангов и соответствие заголовка порядку значений.
    for r in rows:
        for model in MODELS:
            for rank in MAIN_RANKS:
                col = f"{model}_{rank}"
                if col not in r:
                    violations.append(Violation("MISSING_RANK_COL", col))
    # Сверка значений с входами: каждая колонка {model}_{rank} должна совпадать
    # с разобранным рангом из исходной таксономии; {model}_taxon — с исходной строкой.
    by_id = {r["feature_id"]: r for r in rows}
    for fid in feature_ids:
        r = by_id.get(fid)
        if r is None:
            continue
        for model in MODELS:
            pt = taxonomies[model][fid]
            if r.get(f"{model}_taxon") != pt.raw:
                violations.append(Violation(
                    "ASSIGNMENT_TAXON", f"{fid} {model}: {r.get(f'{model}_taxon')!r} != {pt.raw!r}"))
            for rank in MAIN_RANKS:
                expected = pt.ranks.get(rank)
                expected_s = expected if expected is not None else ""
                got = r.get(f"{model}_{rank}", "")
                if got != expected_s:
                    violations.append(Violation(
                        "ASSIGNMENT_RANK", f"{fid} {model} {rank}: {got!r} != {expected_s!r}"))


def check_coverage(violations: List[Violation], out_dir: Path,
                   taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                   table: Dict[str, Dict[str, int]], soil_asv: List[str],
                   soil_samples: List[str]) -> None:
    rows = read_tsv(out_dir / "coverage.tsv")
    seen = set()
    for r in rows:
        key = (r["model"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_COVERAGE", str(key)))
        seen.add(key)
        if r["model"] not in MODELS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_COVERAGE_KEY", str(key)))
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
        if int(r["assigned_asv"]) + int(r["unassigned_asv"]) != int(r["total_asv"]):
            violations.append(Violation("COVERAGE_ASV_SUM", str(key)))
        if int(r["assigned_reads"]) + int(r["unassigned_reads"]) != int(r["total_reads"]):
            violations.append(Violation("COVERAGE_READS_SUM", str(key)))
        # Сверка с входами: пересчёт покрытия по почвенному набору ASV.
        model, rank = r["model"], r["rank"]
        exp = coverage({model: taxonomies[model]}, soil_asv, table, soil_samples)[model][rank]
        if int(r["assigned_asv"]) != exp["assigned_asv"] or int(r["unassigned_asv"]) != exp["unassigned_asv"]:
            violations.append(Violation("COVERAGE_ASV_VALUE", str(key)))
        if int(r["assigned_reads"]) != exp["assigned_reads"] or int(r["unassigned_reads"]) != exp["unassigned_reads"]:
            violations.append(Violation("COVERAGE_READS_VALUE", str(key)))
        if abs(a - exp["assigned_asv_frac"]) > TOL or abs(u - exp["unassigned_asv_frac"]) > TOL:
            violations.append(Violation("COVERAGE_ASV_FRAC_VALUE", str(key)))
        if abs(ar - exp["assigned_reads_frac"]) > TOL or abs(ur - exp["unassigned_reads_frac"]) > TOL:
            violations.append(Violation("COVERAGE_READS_FRAC_VALUE", str(key)))
    # Полнота: все (model, rank).
    for model in MODELS:
        for rank in MAIN_RANKS:
            if (model, rank) not in seen:
                violations.append(Violation("MISSING_COVERAGE", f"{model} {rank}"))


def check_coverage_per_sample(violations: List[Violation], out_dir: Path,
                              taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                              table: Dict[str, Dict[str, int]], soil_asv: List[str],
                              soil_samples: List[str]) -> None:
    rows = read_tsv(out_dir / "coverage_per_sample.tsv")
    seen = set()
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_COVERAGE_SAMPLE", str(key)))
        seen.add(key)
        if r["model"] not in MODELS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_COVERAGE_SAMPLE_KEY", str(key)))
        if r["sample"] not in soil_samples:
            violations.append(Violation("BAD_COVERAGE_SAMPLE", str(key)))
        a = float(r["assigned_asv_frac"]); u = float(r["unassigned_asv_frac"])
        if abs(a + u - 1.0) > TOL:
            violations.append(Violation("COVERAGE_SAMPLE_FRAC_SUM", str(key)))
        ar = float(r["assigned_reads_frac"]); ur = float(r["unassigned_reads_frac"])
        if abs(ar + ur - 1.0) > TOL:
            violations.append(Violation("COVERAGE_SAMPLE_READS_FRAC_SUM", str(key)))
        if int(r["assigned_asv"]) + int(r["unassigned_asv"]) != int(r["total_asv"]):
            violations.append(Violation("COVERAGE_SAMPLE_ASV_SUM", str(key)))
        if int(r["assigned_reads"]) + int(r["unassigned_reads"]) != int(r["total_reads"]):
            violations.append(Violation("COVERAGE_SAMPLE_READS_SUM", str(key)))
        # Сверка с входами.
        model, sample, rank = r["model"], r["sample"], r["rank"]
        exp = coverage_per_sample({model: taxonomies[model]}, soil_asv, table, sample)[model][rank]
        if int(r["assigned_asv"]) != exp["assigned_asv"] or int(r["unassigned_asv"]) != exp["unassigned_asv"]:
            violations.append(Violation("COVERAGE_SAMPLE_ASV_VALUE", str(key)))
        if int(r["assigned_reads"]) != exp["assigned_reads"] or int(r["unassigned_reads"]) != exp["unassigned_reads"]:
            violations.append(Violation("COVERAGE_SAMPLE_READS_VALUE", str(key)))
        if abs(ar - exp["assigned_reads_frac"]) > TOL or abs(ur - exp["unassigned_reads_frac"]) > TOL:
            violations.append(Violation("COVERAGE_SAMPLE_READS_FRAC_VALUE", str(key)))
    for model in MODELS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                if (model, sample, rank) not in seen:
                    violations.append(Violation("MISSING_COVERAGE_SAMPLE", f"{model} {sample} {rank}"))


def _check_pairwise_rows(violations: List[Violation], rows: List[Dict[str, str]],
                         per_sample: bool,
                         taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                         table: Dict[str, Dict[str, int]], soil_asv: List[str],
                         soil_samples: List[str]) -> None:
    seen = set()
    for r in rows:
        if per_sample:
            key = (r["model_a"], r["model_b"], r["sample"], r["rank"])
        else:
            key = (r["model_a"], r["model_b"], r["rank"])
        if key in seen:
            violations.append(Violation("DUP_PAIRWISE", str(key)))
        seen.add(key)
        if (r["model_a"], r["model_b"]) not in PAIRS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_PAIRWISE_KEY", str(key)))
        n = (int(r["both_match"]) + int(r["both_diff"]) + int(r["only_a"])
             + int(r["only_b"]) + int(r["neither"]))
        if n != int(r["total_asv"]):
            violations.append(Violation("PAIRWISE_COVER", str(key)))
        fracs = [float(r[f"{c}_frac"]) for c in
                 ("both_match", "both_diff", "only_a", "only_b", "neither")]
        if any(not (0 - TOL <= x <= 1 + TOL) for x in fracs):
            violations.append(Violation("PAIRWISE_FRAC_RANGE", str(key)))
        if abs(sum(fracs) - 1.0) > TOL:
            violations.append(Violation("PAIRWISE_FRAC_SUM", str(key)))
        rfracs = [float(r[f"{c}_reads_frac"]) for c in
                  ("both_match", "both_diff", "only_a", "only_b", "neither")]
        if any(not (0 - TOL <= x <= 1 + TOL) for x in rfracs):
            violations.append(Violation("PAIRWISE_READS_FRAC_RANGE", str(key)))
        if abs(sum(rfracs) - 1.0) > TOL:
            violations.append(Violation("PAIRWISE_READS_FRAC_SUM", str(key)))
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
        # Сверка с входами: пересчёт категорий по почвенному набору ASV.
        a, b, rank = r["model_a"], r["model_b"], r["rank"]
        if per_sample:
            sample = r["sample"]
            exp = pairwise_match_per_sample(
                taxonomies[a], taxonomies[b], soil_asv, table, sample, rank)
        else:
            exp = pairwise_match(taxonomies[a], taxonomies[b], soil_asv, table, soil_samples, rank)
        if int(r["both_match"]) != exp["both_match"] or int(r["both_diff"]) != exp["both_diff"] \
           or int(r["only_a"]) != exp["only_a"] or int(r["only_b"]) != exp["only_b"] \
           or int(r["neither"]) != exp["neither"]:
            violations.append(Violation("PAIRWISE_COUNTS_VALUE", str(key)))
        if abs(float(r["both_match_frac"]) - exp["both_match_frac"]) > TOL \
           or abs(float(r["both_diff_frac"]) - exp["both_diff_frac"]) > TOL \
           or abs(float(r["only_a_frac"]) - exp["only_a_frac"]) > TOL \
           or abs(float(r["only_b_frac"]) - exp["only_b_frac"]) > TOL \
           or abs(float(r["neither_frac"]) - exp["neither_frac"]) > TOL:
            violations.append(Violation("PAIRWISE_FRAC_VALUE", str(key)))
        if abs(float(r["both_match_reads_frac"]) - exp["both_match_reads_frac"]) > TOL \
           or abs(float(r["both_diff_reads_frac"]) - exp["both_diff_reads_frac"]) > TOL \
           or abs(float(r["only_a_reads_frac"]) - exp["only_a_reads_frac"]) > TOL \
           or abs(float(r["only_b_reads_frac"]) - exp["only_b_reads_frac"]) > TOL \
           or abs(float(r["neither_reads_frac"]) - exp["neither_reads_frac"]) > TOL:
            violations.append(Violation("PAIRWISE_READS_FRAC_VALUE", str(key)))
        if int(r["both_assigned"]) != exp["both_assigned"]:
            violations.append(Violation("PAIRWISE_BOTH_ASSIGNED_VALUE", str(key)))
    return seen


def check_pairwise(violations: List[Violation], out_dir: Path,
                   taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                   table: Dict[str, Dict[str, int]], soil_asv: List[str],
                   soil_samples: List[str]) -> None:
    rows = read_tsv(out_dir / "pairwise_match.tsv")
    seen = _check_pairwise_rows(violations, rows, False,
                                taxonomies, table, soil_asv, soil_samples)
    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            if (a, b, rank) not in seen:
                violations.append(Violation("MISSING_PAIRWISE", f"{a} {b} {rank}"))


def check_pairwise_per_sample(violations: List[Violation], out_dir: Path,
                              soil_samples: List[str],
                              taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                              table: Dict[str, Dict[str, int]], soil_asv: List[str]) -> None:
    rows = read_tsv(out_dir / "pairwise_match_per_sample.tsv")
    seen = _check_pairwise_rows(violations, rows, True,
                                taxonomies, table, soil_asv, soil_samples)
    for a, b in PAIRS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                if (a, b, sample, rank) not in seen:
                    violations.append(Violation("MISSING_PAIRWISE_SAMPLE",
                                                f"{a} {b} {sample} {rank}"))


def check_profiles(violations: List[Violation], table: Dict[str, Dict[str, int]],
                   soil_samples: List[str], out_dir: Path,
                   taxonomies: Dict[str, Dict[str, ParsedTaxon]]) -> None:
    rows = read_tsv(out_dir / "profiles.tsv")
    sums: Dict[Tuple[str, str, str], int] = {}
    totals: Dict[Tuple[str, str, str], int] = {}
    keys_seen: Dict[Tuple[str, str, str], set] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        sums[key] = sums.get(key, 0) + int(r["count"])
        totals[key] = int(r["total_reads"])
        keys_seen.setdefault(key, set()).add(r["key"])
        if r["model"] not in MODELS or r["rank"] not in MAIN_RANKS:
            violations.append(Violation("BAD_PROFILE_KEY", str(key)))
        if r["sample"] not in soil_samples:
            violations.append(Violation("BAD_PROFILE_SAMPLE", str(key)))
        f = float(r["frac"])
        if not (0 - TOL <= f <= 1 + TOL):
            violations.append(Violation("PROFILE_FRAC_RANGE", str(key)))
    # Полнота и сохранение массы.
    for model in MODELS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                key = (model, sample, rank)
                if key not in sums:
                    violations.append(Violation("MISSING_PROFILE", str(key)))
                    continue
                if sums[key] != totals[key]:
                    violations.append(Violation("PROFILE_MASS", str(key)))
    # Дубликаты ключей внутри (model, sample, rank).
    seen_keys: Dict[Tuple[str, str, str], set] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        if r["key"] in seen_keys.get(key, set()):
            violations.append(Violation("DUP_PROFILE_KEY", str(key) + " " + r["key"]))
        seen_keys.setdefault(key, set()).add(r["key"])
    # Сумма долей = 1 по каждому ключу.
    frac_sums: Dict[Tuple[str, str, str], float] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        frac_sums[key] = frac_sums.get(key, 0.0) + float(r["frac"])
    for key, s in frac_sums.items():
        if abs(s - 1.0) > TOL:
            violations.append(Violation("PROFILE_FRAC_SUM", str(key)))
    # Сверка с входами: каждый ключ и частота профиля — из исходной таблицы и таксономии.
    by_key: Dict[Tuple[str, str, str], Dict[str, Dict[str, str]]] = {}
    for r in rows:
        key = (r["model"], r["sample"], r["rank"])
        by_key.setdefault(key, {})[r["key"]] = r
    for model in MODELS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                key = (model, sample, rank)
                exp, exp_total = profile(taxonomies[model], table, sample, rank)
                if exp_total != totals.get(key):
                    violations.append(Violation("PROFILE_TOTAL_VALUE", str(key)))
                exp_keys = set(exp)
                got_keys = set(by_key.get(key, {}))
                if exp_keys != got_keys:
                    violations.append(Violation("PROFILE_KEYS_VALUE", str(key)))
                for k, cnt in exp.items():
                    r = by_key.get(key, {}).get(k)
                    if r is None:
                        continue
                    if int(r["count"]) != cnt:
                        violations.append(Violation("PROFILE_COUNT_VALUE", f"{key} {k}"))
                    exp_frac = cnt / exp_total if exp_total else 0.0
                    if abs(float(r["frac"]) - exp_frac) > TOL:
                        violations.append(Violation("PROFILE_FRAC_VALUE", f"{key} {k}"))


def check_distances(violations: List[Violation], soil_samples: List[str], out_dir: Path,
                    taxonomies: Dict[str, Dict[str, ParsedTaxon]],
                    table: Dict[str, Dict[str, int]]) -> None:
    rows = read_tsv(out_dir / "distances.tsv")
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
        # Сверка с входами: пересчёт расстояния из профилей.
        a, b, sample, rank = r["model_a"], r["model_b"], r["sample"], r["rank"]
        pa, _ = profile(taxonomies[a], table, sample, rank)
        pb, _ = profile(taxonomies[b], table, sample, rank)
        exp_d = profile_distance(pa, pb)
        if abs(d - exp_d) > TOL:
            violations.append(Violation("DISTANCE_VALUE", str(key)))
    for a, b in PAIRS:
        for sample in soil_samples:
            for rank in MAIN_RANKS:
                if (a, b, sample, rank) not in seen:
                    violations.append(Violation("MISSING_DISTANCE", f"{a} {b} {sample} {rank}"))


def check_summary(violations: List[Violation], out_dir: Path) -> None:
    rows = read_tsv(out_dir / "summary.tsv")
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
    # Сверка сводки с индивидуальными расстояниями.
    dist_rows = read_tsv(out_dir / "distances.tsv")
    dist_by_key: Dict[Tuple[str, str, str], List[float]] = {}
    for r in dist_rows:
        key = (r["model_a"], r["model_b"], r["rank"])
        dist_by_key.setdefault(key, []).append(float(r["distance"]))
    for r in rows:
        key = (r["model_a"], r["model_b"], r["rank"])
        vals = sorted(dist_by_key.get(key, []))
        if len(vals) != int(r["n_samples"]):
            violations.append(Violation("SUMMARY_N_VALUE", str(key)))
            continue
        if abs(float(r["min"]) - vals[0]) > TOL or abs(float(r["max"]) - vals[-1]) > TOL:
            violations.append(Violation("SUMMARY_MINMAX_VALUE", str(key)))
        n = len(vals)
        med = vals[n // 2] if n % 2 == 1 else (vals[n // 2 - 1] + vals[n // 2]) / 2
        if abs(float(r["median"]) - med) > TOL:
            violations.append(Violation("SUMMARY_MEDIAN_VALUE", str(key)))
    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            if (a, b, rank) not in seen:
                violations.append(Violation("MISSING_SUMMARY", f"{a} {b} {rank}"))


def check_manifest(violations: List[Violation], out_dir: Path) -> None:
    rows = read_tsv(out_dir / "manifest.tsv")
    items = {r["item"]: r for r in rows}
    for item in INPUT_FILES:
        if item not in items:
            violations.append(Violation("MISSING_MANIFEST", item))
            continue
        r = items[item]
        # SHA256 совпадает с настоящим входом.
        actual = sha256(ROOT / INPUT_FILES[item])
        if r["sha256"] != actual:
            violations.append(Violation("MANIFEST_SHA256", f"{item}: {r['sha256']} != {actual}"))
        # UUID для QZA.
        if item in MODELS or item in ("ASV_table", "rep_seq"):
            if not r["uuid"]:
                violations.append(Violation("MANIFEST_NO_UUID", item))
            else:
                try:
                    actual_uuid = qza_uuid(ROOT / INPUT_FILES[item])
                    if r["uuid"] != actual_uuid:
                        violations.append(Violation("MANIFEST_UUID", f"{item}: {r['uuid']} != {actual_uuid}"))
                except Exception as e:
                    violations.append(Violation("MANIFEST_UUID_READ", f"{item}: {e}"))
        # Версия кода: должна быть фактическим коммитом или снимком кода.
        if not r["code_version"] or r["code_version"] == "не указано":
            violations.append(Violation("MANIFEST_NO_CODE_VERSION", item))
        elif not _is_valid_code_version(r["code_version"]):
            violations.append(Violation("MANIFEST_BAD_CODE_VERSION", f"{item}: {r['code_version']}"))


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Проверка сравнения таксономических профилей")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = ap.parse_args()
    out_dir = args.out_dir

    violations: List[Violation] = []

    check_required_files(out_dir, violations)

    try:
        table, feature_ids, sample_ids, metadata, taxonomies = check_inputs(violations)
    except ValueError:
        # Ошибка чтения/частот — аварийный выход, не штатный FAIL.
        print("ERROR — ошибка чтения входов или частот")
        for v in violations:
            print(f"  {v}")
        return 2
    soil_samples = [s for s in sample_ids if metadata.get(s) == "soil"]
    soil_asv = soil_asv_ids(table, soil_samples)

    check_assignments(violations, feature_ids, taxonomies, out_dir)
    check_coverage(violations, out_dir, taxonomies, table, soil_asv, soil_samples)
    check_coverage_per_sample(violations, out_dir, taxonomies, table, soil_asv, soil_samples)
    check_pairwise(violations, out_dir, taxonomies, table, soil_asv, soil_samples)
    check_pairwise_per_sample(violations, out_dir, soil_samples, taxonomies, table, soil_asv)
    check_profiles(violations, table, soil_samples, out_dir, taxonomies)
    check_distances(violations, soil_samples, out_dir, taxonomies, table)
    check_summary(violations, out_dir)
    check_manifest(violations, out_dir)

    if violations:
        # FAIL/1 имеет приоритет, если установлено нарушение данных.
        # INCOMPLETE/2 — только если нарушений нет, но не хватает обязательных выходов.
        data_violations = [v for v in violations if v.code != "MISSING_OUTPUT"]
        if data_violations:
            print("FAIL/1")
            for v in violations:
                print(f"  {v}")
            return 1
        print("INCOMPLETE/2 — отсутствуют обязательные выходы")
        for v in violations:
            print(f"  {v}")
        return 2

    print("PASS/0 — сравнение таксономических профилей корректно")
    return 0


if __name__ == "__main__":
    sys.exit(main())
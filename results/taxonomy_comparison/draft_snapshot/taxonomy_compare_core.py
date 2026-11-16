"""Чистые вычислительные функции сравнения таксономических профилей.

Чтение файлов отделено от вычислений: этот модуль не знает о файлах,
только о структурах данных (таблица частот, метаданные, таксономии).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

# Основные ранги в порядке глубины. Вид в основной анализ не входит.
MAIN_RANKS: Tuple[str, ...] = ("тип", "класс", "порядок", "семейство", "род")

# Обозначения рангов в строках Taxon каждой модели.
# «тип» = филум (p__); домен (d__) не входит в основные ранги.
RANK_PREFIXES: Dict[str, str] = {
    "тип": "p__",
    "класс": "c__",
    "порядок": "o__",
    "семейство": "f__",
    "род": "g__",
    "вид": "s__",
}

# Префиксы, которые распознаются, но не входят в основные ранги (домен).
NON_MAIN_PREFIXES = {"d__"}

# Служебные метки, не считающиеся известными таксонами.
SERVICE_LABELS = {"Unassigned", "Unknown_Family"}

# Модели и их пары.
MODELS: Tuple[str, ...] = ("sklearn_silva138", "sklearn_gtdb", "sklearn_gg2")
PAIRS: Tuple[Tuple[str, str], ...] = (
    ("sklearn_silva138", "sklearn_gtdb"),
    ("sklearn_silva138", "sklearn_gg2"),
    ("sklearn_gtdb", "sklearn_gg2"),
)


@dataclass
class ParsedTaxon:
    """Разобранная строка Taxon одной модели."""

    raw: str
    # Ранг -> нормализованная метка (без префикса ранга) или None, если ранга нет.
    ranks: Dict[str, Optional[str]]
    # True, если назначение известно (не служебная метка и не пустое).
    assigned: bool


def normalize_label(label: str) -> str:
    """Техническая нормализация метки: убрать пробелы по краям.

    Значимые суффиксы GTDB, названия кластеров и различия родословной
    сохраняются. Это только техническая нормализация, не объединение синонимов.
    """
    return label.strip()


def parse_taxon(raw: str, separator: str) -> ParsedTaxon:
    """Разобрать строку Taxon по фактическому формату модели.

    separator — разделитель рангов (';' или '; '). Ранги определяются по
    префиксам, а не по числу непустых элементов после удаления пропусков.
    """
    if raw == "Unassigned":
        return ParsedTaxon(raw=raw, ranks={r: None for r in MAIN_RANKS}, assigned=False)

    # Разделитель '; ' (Greengenes2) или ';' (SILVA, GTDB). Нормализуем к ';'.
    parts = [p.strip() for p in raw.split(";")]
    ranks: Dict[str, Optional[str]] = {r: None for r in MAIN_RANKS}
    assigned = True
    for part in parts:
        if not part:
            continue
        matched = False
        for rank, prefix in RANK_PREFIXES.items():
            if part.startswith(prefix):
                label = part[len(prefix):].strip()
                if rank in MAIN_RANKS:
                    if label == "" or label in SERVICE_LABELS:
                        ranks[rank] = None
                    else:
                        ranks[rank] = normalize_label(label)
                matched = True
                break
        if not matched:
            # Домен (d__) распознаётся, но не входит в основные ранги.
            if any(part.startswith(p) for p in NON_MAIN_PREFIXES):
                continue
            # Неожиданный формат — явная диагностика.
            raise ValueError(f"Неожиданный фрагмент Taxon: {part!r} в строке {raw!r}")

    # Если ни один ранг не назначен (например, только d__Bacteria без пустых
    # суффиксов — это назначение на уровне домена, но домен не входит в MAIN_RANKS).
    # Назначение считается известным, если есть хотя бы один непустой основной ранг.
    if not any(v is not None for v in ranks.values()):
        assigned = False
    return ParsedTaxon(raw=raw, ranks=ranks, assigned=assigned)


def parse_taxonomy_rows(
    rows: Sequence[Tuple[str, str]], separator: str
) -> Dict[str, ParsedTaxon]:
    """Разобрать строки таксономии (Feature ID, Taxon) в словарь.

    rows — последовательность пар (feature_id, taxon_string).
    Возвращает словарь feature_id -> ParsedTaxon.
    """
    result: Dict[str, ParsedTaxon] = {}
    for fid, taxon in rows:
        if fid in result:
            raise ValueError(f"Дубликат Feature ID в таксономии: {fid}")
        result[fid] = parse_taxon(taxon, separator)
    return result


def soil_asv_ids(
    table: Dict[str, Dict[str, int]], soil_samples: Sequence[str]
) -> List[str]:
    """ASV с ненулевой суммой частот в почвенных образцах.

    table — feature_id -> {sample_id -> count}.
    Возвращает список feature_id в порядке появления.
    """
    result = []
    for fid in table:
        total = sum(table[fid].get(s, 0) for s in soil_samples)
        if total > 0:
            result.append(fid)
    return result


def coverage(
    taxonomies: Dict[str, Dict[str, ParsedTaxon]],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    soil_samples: Sequence[str],
) -> Dict[str, Dict[str, object]]:
    """Покрытие классификацией: доля назначенных/неназначенных ASV и ридов.

    Возвращает {model: {rank: {assigned_asv, unassigned_asv, assigned_reads,
    unassigned_reads, ...}}}.
    """
    out: Dict[str, Dict[str, object]] = {}
    for model, tax in taxonomies.items():
        by_rank: Dict[str, Dict[str, object]] = {}
        for rank in MAIN_RANKS:
            assigned_asv = 0
            unassigned_asv = 0
            assigned_reads = 0
            unassigned_reads = 0
            for fid in asv_ids:
                pt = tax[fid]
                is_assigned = pt.assigned and pt.ranks.get(rank) is not None
                reads = sum(table[fid].get(s, 0) for s in soil_samples)
                if is_assigned:
                    assigned_asv += 1
                    assigned_reads += reads
                else:
                    unassigned_asv += 1
                    unassigned_reads += reads
            total_asv = assigned_asv + unassigned_asv
            total_reads = assigned_reads + unassigned_reads
            by_rank[rank] = {
                "assigned_asv": assigned_asv,
                "unassigned_asv": unassigned_asv,
                "total_asv": total_asv,
                "assigned_asv_frac": assigned_asv / total_asv if total_asv else 0.0,
                "unassigned_asv_frac": unassigned_asv / total_asv if total_asv else 0.0,
                "assigned_reads": assigned_reads,
                "unassigned_reads": unassigned_reads,
                "total_reads": total_reads,
                "assigned_reads_frac": assigned_reads / total_reads if total_reads else 0.0,
                "unassigned_reads_frac": unassigned_reads / total_reads if total_reads else 0.0,
            }
        out[model] = by_rank
    return out


def _label_for(pt: ParsedTaxon, rank: str) -> Optional[str]:
    """Метка ранга, если назначение известно, иначе None."""
    if not pt.assigned:
        return None
    return pt.ranks.get(rank)


def pairwise_match(
    tax_a: Dict[str, ParsedTaxon],
    tax_b: Dict[str, ParsedTaxon],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    soil_samples: Sequence[str],
    rank: str,
) -> Dict[str, object]:
    """Пять категорий попарного совпадения назначений на одном ранге.

    Категории: both_match, both_diff, only_a, only_b, neither.
    Возвращает невзвешенные доли по ASV, взвешенные частотами ридов,
    условную долю совпадений среди назначенных обеими и размер подгруппы.
    """
    both_match = both_diff = only_a = only_b = neither = 0
    both_match_reads = both_diff_reads = only_a_reads = only_b_reads = neither_reads = 0
    for fid in asv_ids:
        la = _label_for(tax_a[fid], rank)
        lb = _label_for(tax_b[fid], rank)
        reads = sum(table[fid].get(s, 0) for s in soil_samples)
        if la is not None and lb is not None:
            if la == lb:
                both_match += 1
                both_match_reads += reads
            else:
                both_diff += 1
                both_diff_reads += reads
        elif la is not None:
            only_a += 1
            only_a_reads += reads
        elif lb is not None:
            only_b += 1
            only_b_reads += reads
        else:
            neither += 1
            neither_reads += reads

    total_asv = len(asv_ids)
    total_reads = sum(
        sum(table[fid].get(s, 0) for s in soil_samples) for fid in asv_ids
    )
    both_assigned = both_match + both_diff
    cond_match = both_match / both_assigned if both_assigned else None

    return {
        "rank": rank,
        "both_match": both_match,
        "both_diff": both_diff,
        "only_a": only_a,
        "only_b": only_b,
        "neither": neither,
        "total_asv": total_asv,
        "both_match_frac": both_match / total_asv if total_asv else 0.0,
        "both_diff_frac": both_diff / total_asv if total_asv else 0.0,
        "only_a_frac": only_a / total_asv if total_asv else 0.0,
        "only_b_frac": only_b / total_asv if total_asv else 0.0,
        "neither_frac": neither / total_asv if total_asv else 0.0,
        "both_match_reads": both_match_reads,
        "both_diff_reads": both_diff_reads,
        "only_a_reads": only_a_reads,
        "only_b_reads": only_b_reads,
        "neither_reads": neither_reads,
        "total_reads": total_reads,
        "both_match_reads_frac": both_match_reads / total_reads if total_reads else 0.0,
        "both_diff_reads_frac": both_diff_reads / total_reads if total_reads else 0.0,
        "only_a_reads_frac": only_a_reads / total_reads if total_reads else 0.0,
        "only_b_reads_frac": only_b_reads / total_reads if total_reads else 0.0,
        "neither_reads_frac": neither_reads / total_reads if total_reads else 0.0,
        "both_assigned": both_assigned,
        "cond_match_frac": cond_match,  # None -> NA
    }


def profile(
    tax: Dict[str, ParsedTaxon],
    table: Dict[str, Dict[str, int]],
    sample: str,
    rank: str,
) -> Tuple[Dict[str, int], int]:
    """Профиль образца на ранге: {ключ: частота}, total_reads.

    Неназначенная масса сохраняется отдельной категорией UNCLASSIFIED.
    Ключи — нормализованные метки ранга; UNCLASSIFIED — для неназначенных.
    Возвращает (профиль, total_reads_образца).
    """
    prof: Dict[str, int] = {}
    total = 0
    for fid, counts in table.items():
        reads = counts.get(sample, 0)
        if reads == 0:
            continue
        total += reads
        pt = tax[fid]
        label = _label_for(pt, rank)
        key = label if label is not None else "UNCLASSIFIED"
        prof[key] = prof.get(key, 0) + reads
    return prof, total


def profile_distance(p: Dict[str, int], q: Dict[str, int]) -> float:
    """Расстояние профилей D = 0.5 * sum_t |p_t - q_t| по объединению ключей.

    p, q — профили {ключ: частота}. Нормируются на общее число ридов образца.
    Возвращает долю в [0, 1].
    """
    total_p = sum(p.values())
    total_q = sum(q.values())
    if total_p == 0 and total_q == 0:
        return 0.0
    keys = set(p) | set(q)
    d = 0.0
    for k in keys:
        fp = p.get(k, 0) / total_p if total_p else 0.0
        fq = q.get(k, 0) / total_q if total_q else 0.0
        d += abs(fp - fq)
    return 0.5 * d
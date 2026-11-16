"""Чистые вычислительные функции сравнения таксономических профилей.

Чтение файлов отделено от вычислений: этот модуль не знает о файлах,
только о структурах данных (таблица частот, метаданные, таксономии).

Исправления 2026-11-11 (по поручению):
- ключ родословной от домена до целевого ранга используется во всех
  расчётах и выходах (род G под разными семействами — разные ключи);
- дублированный ранг (p__P1;p__P2) — явная диагностика, а не молчаливый
  выбор последнего;
- взвешенные показатели считаются отдельно для каждого почвенного образца,
  невзвешенные — по общему почвенному набору ASV.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import fsum
from statistics import median
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

# Разделители рангов по моделям (фактический формат).
SEPARATORS: Dict[str, str] = {
    "sklearn_silva138": ";",
    "sklearn_gtdb": ";",
    "sklearn_gg2": "; ",
}


@dataclass
class ParsedTaxon:
    """Разобранная строка Taxon одной модели."""

    raw: str
    # Ранг -> нормализованная метка (без префикса ранга) или None, если ранга нет.
    ranks: Dict[str, Optional[str]]
    # True, если назначение известно (не служебная метка и не пустое).
    assigned: bool
    # Домен (d__) как есть, или None.
    domain: Optional[str] = None


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
    Дублированный ранг (например p__P1;p__P2) — явная диагностика.
    """
    if raw.strip() == "Unassigned":
        return ParsedTaxon(raw=raw, ranks={r: None for r in MAIN_RANKS}, assigned=False)

    parts = [p.strip() for p in raw.split(";")]
    ranks: Dict[str, Optional[str]] = {r: None for r in MAIN_RANKS}
    domain: Optional[str] = None
    seen_prefixes: Dict[str, str] = {}
    for part in parts:
        if not part:
            continue
        matched = False
        for rank, prefix in RANK_PREFIXES.items():
            if part.startswith(prefix):
                label = part[len(prefix):].strip()
                # Дублированный ранг — диагностика.
                if rank in seen_prefixes:
                    raise ValueError(
                        f"Дублированный ранг {prefix} в строке {raw!r}: "
                        f"{seen_prefixes[rank]!r} и {label!r}"
                    )
                seen_prefixes[rank] = label
                if rank in MAIN_RANKS:
                    if label == "" or label in SERVICE_LABELS:
                        ranks[rank] = None
                    else:
                        ranks[rank] = normalize_label(label)
                matched = True
                break
        if not matched:
            if any(part.startswith(p) for p in NON_MAIN_PREFIXES):
                # Домен (d__). Дублированный домен — диагностика.
                if "domain" in seen_prefixes:
                    raise ValueError(f"Дублированный домен в строке {raw!r}")
                label = part[3:].strip()
                seen_prefixes["domain"] = label
                domain = None if label in SERVICE_LABELS else (label or None)
                continue
            raise ValueError(f"Неожиданный фрагмент Taxon: {part!r} в строке {raw!r}")

    # Назначение считается известным, если есть хотя бы один непустой основной ранг.
    assigned = any(v is not None for v in ranks.values())
    return ParsedTaxon(raw=raw, ranks=ranks, assigned=assigned, domain=domain)


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


def lineage_key(pt: ParsedTaxon, rank: str) -> Optional[str]:
    """Ключ родословной от домена до целевого ранга (включительно).

    Возвращает None, если назначение на целевом ранге неизвестно
    (нет известного таксона на целевом ранге). Ключ учитывает всю
    родословную: род G под разными семействами — разные ключи.
    """
    if rank not in MAIN_RANKS:
        raise ValueError(f"Неизвестный основной ранг: {rank!r}")
    if not pt.assigned:
        return None
    label = pt.ranks.get(rank)
    if label is None:
        return None
    parts = []
    if pt.domain:
        parts.append(f"d__{pt.domain}")
    for r in MAIN_RANKS:
        if r == rank:
            parts.append(f"{RANK_PREFIXES[r]}{label}")
            break
        v = pt.ranks.get(r)
        if v is not None:
            parts.append(f"{RANK_PREFIXES[r]}{v}")
    return ";".join(parts)


def _label_for(pt: ParsedTaxon, rank: str) -> Optional[str]:
    """Метка ранга, если назначение известно, иначе None."""
    if not pt.assigned:
        return None
    return pt.ranks.get(rank)


def coverage(
    taxonomies: Dict[str, Dict[str, ParsedTaxon]],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    soil_samples: Sequence[str],
) -> Dict[str, Dict[str, object]]:
    """Покрытие классификацией: доля назначенных/неназначенных ASV и ридов.

    Возвращает {model: {rank: {...}}}. Доли ASV — по почвенному набору ASV;
    доли ридов — по сумме ридов почвенных образцов.
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
                is_assigned = lineage_key(pt, rank) is not None
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


def coverage_per_sample(
    taxonomies: Dict[str, Dict[str, ParsedTaxon]],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    sample: str,
) -> Dict[str, Dict[str, object]]:
    """Покрытие классификацией для одного образца.

    Возвращает {model: {rank: {...}}}. Доли ASV — по почвенному набору ASV
    (одинаковы для всех образцов); доли ридов — по частотам конкретного образца.
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
                is_assigned = lineage_key(pt, rank) is not None
                reads = table[fid].get(sample, 0)
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


def _pairwise_counts(
    tax_a: Dict[str, ParsedTaxon],
    tax_b: Dict[str, ParsedTaxon],
    asv_ids: Sequence[str],
    weights: Dict[str, float],
    rank: str,
) -> Dict[str, object]:
    """Пять категорий попарного совпадения на одном ранге с весами.

    weights — {feature_id: вес} (1 для невзвешенных, частота ридов для взвешенных).
    Возвращает счётчики и доли по категориям.
    """
    both_match = both_diff = only_a = only_b = neither = 0
    both_match_n = both_diff_n = only_a_n = only_b_n = neither_n = 0
    total_w = 0
    for fid in asv_ids:
        ka = lineage_key(tax_a[fid], rank)
        kb = lineage_key(tax_b[fid], rank)
        w = weights.get(fid, 0)
        total_w += w
        if ka is not None and kb is not None:
            if ka == kb:
                both_match += w
                both_match_n += 1
            else:
                both_diff += w
                both_diff_n += 1
        elif ka is not None:
            only_a += w
            only_a_n += 1
        elif kb is not None:
            only_b += w
            only_b_n += 1
        else:
            neither += w
            neither_n += 1

    both_assigned_w = both_match + both_diff
    both_assigned_n = both_match_n + both_diff_n
    # Условные доли совпадений среди назначенных обеими — раздельно по ASV и по ридам,
    # каждая со своим числителем и знаменателем. При нулевом знаменателе — None (NA).
    cond_match_asv = both_match_n / both_assigned_n if both_assigned_n else None
    cond_match_reads = both_match / both_assigned_w if both_assigned_w else None

    return {
        "both_match": both_match_n, "both_diff": both_diff_n,
        "only_a": only_a_n, "only_b": only_b_n, "neither": neither_n,
        "total_asv": len(asv_ids),
        "both_match_w": both_match, "both_diff_w": both_diff,
        "only_a_w": only_a, "only_b_w": only_b, "neither_w": neither,
        "total_w": total_w,
        "both_assigned_n": both_assigned_n,
        "both_assigned_w": both_assigned_w,
        "cond_match_asv_frac": cond_match_asv,
        "cond_match_reads_frac": cond_match_reads,
    }


def pairwise_match(
    tax_a: Dict[str, ParsedTaxon],
    tax_b: Dict[str, ParsedTaxon],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    soil_samples: Sequence[str],
    rank: str,
) -> Dict[str, object]:
    """Пять категорий попарного совпадения назначений на одном ранге.

    Невзвешенные доли — по почвенному набору ASV (вес 1).
    Взвешенные доли — по сумме ридов почвенных образцов (вес = частота).
    Возвращает невзвешенные и взвешенные доли, условную долю совпадений
    среди назначенных обеими и размер подгруппы.
    """
    # Невзвешенные.
    unweighted = _pairwise_counts(
        tax_a, tax_b, asv_ids, {fid: 1 for fid in asv_ids}, rank
    )
    # Взвешенные по сумме ридов почвенных образцов.
    weights = {
        fid: sum(table[fid].get(s, 0) for s in soil_samples)
        for fid in asv_ids
    }
    weighted = _pairwise_counts(tax_a, tax_b, asv_ids, weights, rank)

    total_asv = unweighted["total_asv"]
    total_w = weighted["total_w"]
    return {
        "rank": rank,
        "both_match": unweighted["both_match"],
        "both_diff": unweighted["both_diff"],
        "only_a": unweighted["only_a"],
        "only_b": unweighted["only_b"],
        "neither": unweighted["neither"],
        "total_asv": total_asv,
        "both_match_frac": unweighted["both_match"] / total_asv if total_asv else 0.0,
        "both_diff_frac": unweighted["both_diff"] / total_asv if total_asv else 0.0,
        "only_a_frac": unweighted["only_a"] / total_asv if total_asv else 0.0,
        "only_b_frac": unweighted["only_b"] / total_asv if total_asv else 0.0,
        "neither_frac": unweighted["neither"] / total_asv if total_asv else 0.0,
        "both_match_reads": weighted["both_match_w"],
        "both_diff_reads": weighted["both_diff_w"],
        "only_a_reads": weighted["only_a_w"],
        "only_b_reads": weighted["only_b_w"],
        "neither_reads": weighted["neither_w"],
        "total_reads": total_w,
        "both_match_reads_frac": weighted["both_match_w"] / total_w if total_w else 0.0,
        "both_diff_reads_frac": weighted["both_diff_w"] / total_w if total_w else 0.0,
        "only_a_reads_frac": weighted["only_a_w"] / total_w if total_w else 0.0,
        "only_b_reads_frac": weighted["only_b_w"] / total_w if total_w else 0.0,
        "neither_reads_frac": weighted["neither_w"] / total_w if total_w else 0.0,
        "both_assigned": weighted["both_assigned_n"],
        "both_assigned_reads": weighted["both_assigned_w"],
        "cond_match_asv_frac": unweighted["cond_match_asv_frac"],  # None -> NA
        "cond_match_reads_frac": weighted["cond_match_reads_frac"],  # None -> NA
    }


def pairwise_match_per_sample(
    tax_a: Dict[str, ParsedTaxon],
    tax_b: Dict[str, ParsedTaxon],
    asv_ids: Sequence[str],
    table: Dict[str, Dict[str, int]],
    sample: str,
    rank: str,
) -> Dict[str, object]:
    """Пять категорий попарного совпадения для одного образца.

    Взвешенные доли — по частотам ридов конкретного образца.
    Возвращает невзвешенные (по ASV) и взвешенные (по ридам образца) доли.
    """
    unweighted = _pairwise_counts(
        tax_a, tax_b, asv_ids, {fid: 1 for fid in asv_ids}, rank
    )
    weights = {fid: table[fid].get(sample, 0) for fid in asv_ids}
    weighted = _pairwise_counts(tax_a, tax_b, asv_ids, weights, rank)

    total_asv = unweighted["total_asv"]
    total_w = weighted["total_w"]
    return {
        "rank": rank,
        "both_match": unweighted["both_match"],
        "both_diff": unweighted["both_diff"],
        "only_a": unweighted["only_a"],
        "only_b": unweighted["only_b"],
        "neither": unweighted["neither"],
        "total_asv": total_asv,
        "both_match_frac": unweighted["both_match"] / total_asv if total_asv else 0.0,
        "both_diff_frac": unweighted["both_diff"] / total_asv if total_asv else 0.0,
        "only_a_frac": unweighted["only_a"] / total_asv if total_asv else 0.0,
        "only_b_frac": unweighted["only_b"] / total_asv if total_asv else 0.0,
        "neither_frac": unweighted["neither"] / total_asv if total_asv else 0.0,
        "both_match_reads": weighted["both_match_w"],
        "both_diff_reads": weighted["both_diff_w"],
        "only_a_reads": weighted["only_a_w"],
        "only_b_reads": weighted["only_b_w"],
        "neither_reads": weighted["neither_w"],
        "total_reads": total_w,
        "both_match_reads_frac": weighted["both_match_w"] / total_w if total_w else 0.0,
        "both_diff_reads_frac": weighted["both_diff_w"] / total_w if total_w else 0.0,
        "only_a_reads_frac": weighted["only_a_w"] / total_w if total_w else 0.0,
        "only_b_reads_frac": weighted["only_b_w"] / total_w if total_w else 0.0,
        "neither_reads_frac": weighted["neither_w"] / total_w if total_w else 0.0,
        "both_assigned": weighted["both_assigned_n"],
        "both_assigned_reads": weighted["both_assigned_w"],
        "cond_match_asv_frac": unweighted["cond_match_asv_frac"],  # None -> NA
        "cond_match_reads_frac": weighted["cond_match_reads_frac"],  # None -> NA
    }


def profile(
    tax: Dict[str, ParsedTaxon],
    table: Dict[str, Dict[str, int]],
    sample: str,
    rank: str,
) -> Tuple[Dict[str, int], int]:
    """Профиль образца на ранге: {ключ родословной: частота}, total_reads.

    Неназначенная масса сохраняется отдельной категорией UNCLASSIFIED.
    Ключи — ключи родословной до целевого ранга; UNCLASSIFIED — для неназначенных.
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
        key = lineage_key(pt, rank)
        if key is None:
            key = "UNCLASSIFIED"
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
    keys = sorted(set(p) | set(q))
    deltas = []
    for k in keys:
        fp = p.get(k, 0) / total_p if total_p else 0.0
        fq = q.get(k, 0) / total_q if total_q else 0.0
        deltas.append(abs(fp - fq))
    return 0.5 * fsum(deltas)


def comparison_tables(table, raw_taxonomies, soil_samples):
    """Build scientific output rows from validated in-memory inputs only.

    Unweighted tables describe the union of ASVs present in soils. Weighted
    tables describe one soil sample at a time; no pooled read metric is emitted.
    All counts remain integers, and zero conditional denominators yield None.
    """
    categories = ("both_match", "both_diff", "only_a", "only_b", "neither")
    taxonomies = {
        m: parse_taxonomy_rows(list(raw_taxonomies[m].items()), SEPARATORS[m])
        for m in MODELS
    }
    samples = sorted(soil_samples)
    features = sorted(table)
    soil_asv = sorted(soil_asv_ids(table, samples))
    if not samples or not soil_asv:
        raise ValueError("Пустой почвенный набор образцов или ASV")
    soil_set = set(soil_asv)
    out = {name: [] for name in (
        "assignments.tsv", "coverage.tsv", "coverage_per_sample.tsv",
        "pairwise_match.tsv", "pairwise_match_per_sample.tsv", "profiles.tsv",
        "distances.tsv", "summary.tsv")}
    for fid in features:
        row = {"feature_id": fid, "soil_asv": int(fid in soil_set)}
        for model in MODELS:
            pt = taxonomies[model][fid]
            row[f"{model}_taxon"] = pt.raw
            for rank in MAIN_RANKS:
                key = lineage_key(pt, rank)
                row[f"{model}_{rank}"] = pt.ranks[rank] or ""
                row[f"{model}_{rank}_key"] = key or "UNCLASSIFIED"
                row[f"{model}_{rank}_assigned"] = int(key is not None)
        out["assignments.tsv"].append(row)

    for model in MODELS:
        for rank in MAIN_RANKS:
            known = sum(lineage_key(taxonomies[model][f], rank) is not None for f in soil_asv)
            n = len(soil_asv)
            out["coverage.tsv"].append({
                "model": model, "rank": rank, "assigned_asv": known,
                "unassigned_asv": n - known, "total_asv": n,
                "assigned_asv_frac": known / n, "unassigned_asv_frac": (n - known) / n,
            })

    coverages = {}
    for sample in samples:
        coverages[sample] = coverage_per_sample(taxonomies, soil_asv, table, sample)
        for model in MODELS:
            for rank in MAIN_RANKS:
                c = coverages[sample][model][rank]
                if not c["total_reads"]:
                    raise ValueError(f"Почвенный образец {sample}: нулевой знаменатель")
                row = {"model": model, "sample": sample, "rank": rank}
                row.update({k: c[k] for k in (
                    "assigned_reads", "unassigned_reads", "total_reads",
                    "assigned_reads_frac", "unassigned_reads_frac")})
                out["coverage_per_sample.tsv"].append(row)

    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            # No read weights in this ASV-only output.
            m = _pairwise_counts(taxonomies[a], taxonomies[b], soil_asv, {}, rank)
            n = m["total_asv"]
            row = {"model_a": a, "model_b": b, "rank": rank, "total_asv": n,
                   "both_assigned_asv": m["both_assigned_n"],
                   "cond_match_asv_frac": m["cond_match_asv_frac"]}
            row.update({c: m[c] for c in categories})
            row.update({c + "_frac": m[c] / n for c in categories})
            out["pairwise_match.tsv"].append(row)
        for sample in samples:
            for rank in MAIN_RANKS:
                m = pairwise_match_per_sample(taxonomies[a], taxonomies[b], soil_asv,
                                              table, sample, rank)
                row = {"model_a": a, "model_b": b, "sample": sample, "rank": rank,
                       "total_reads": m["total_reads"],
                       "both_assigned_reads": m["both_assigned_reads"],
                       "cond_match_reads_frac": m["cond_match_reads_frac"]}
                row.update({c + "_reads": m[c + "_reads"] for c in categories})
                row.update({c + "_reads_frac": m[c + "_reads_frac"] for c in categories})
                out["pairwise_match_per_sample.tsv"].append(row)

    profiles = {}
    for model in MODELS:
        for sample in samples:
            for rank in MAIN_RANKS:
                p, total = profile(taxonomies[model], table, sample, rank)
                profiles[model, sample, rank] = p
                for key, count in sorted(p.items()):
                    if count > 0:
                        out["profiles.tsv"].append({
                            "model": model, "sample": sample, "rank": rank,
                            "key": key, "count": count, "total_reads": total,
                            "frac": count / total,
                        })

    for a, b in PAIRS:
        for rank in MAIN_RANKS:
            distances = []
            for sample in samples:
                distance = profile_distance(profiles[a, sample, rank], profiles[b, sample, rank])
                distances.append(distance)
                out["distances.tsv"].append({
                    "model_a": a, "model_b": b, "sample": sample, "rank": rank,
                    "distance": distance,
                    "assigned_reads_frac_a": coverages[sample][a][rank]["assigned_reads_frac"],
                    "assigned_reads_frac_b": coverages[sample][b][rank]["assigned_reads_frac"],
                })
            out["summary.tsv"].append({
                "model_a": a, "model_b": b, "rank": rank, "n_samples": len(samples),
                "min": min(distances), "median": median(distances), "max": max(distances),
            })
    return out

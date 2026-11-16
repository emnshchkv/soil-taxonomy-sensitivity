"""Синтетические примеры и отрицательные случаи для сравнения таксономий.

Ожидаемые ответы вычислены вручную, не той же функцией, что проверяется.
Все внесённые ошибки — только во временных копиях.

Исправления 2026-11-11 (по поручению):
- ключ родословной от домена до целевого ранга;
- дублированный ранг — диагностика;
- дробные/отрицательные/неконечные частоты — диагностика;
- интеграционные испытания читателя QZA и CLI-команды.
"""

from __future__ import annotations

import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from qiime_runtime import selected_python
selected_python()

from taxonomy_compare_core import (  # noqa: E402
    MAIN_RANKS,
    ParsedTaxon,
    lineage_key,
    parse_taxon,
    parse_taxonomy_rows,
    pairwise_match,
    pairwise_match_per_sample,
    profile,
    profile_distance,
    soil_asv_ids,
)

ROOT = Path(__file__).resolve().parent.parent


def _mk_tax(assignments, domain=None):
    """assignments: {fid: {rank: label or None}}."""
    out = {}
    for fid, ranks in assignments.items():
        out[fid] = ParsedTaxon(
            raw="", ranks=ranks,
            assigned=any(v is not None for v in ranks.values()),
            domain=domain,
        )
    return out


def _mk_table():
    # 3 ASV, 2 почвенных образца.
    return {
        "a1": {"s1": 10, "s2": 0},
        "a2": {"s1": 5, "s2": 5},
        "a3": {"s1": 0, "s2": 10},
    }


# --- Разбор строк Taxon ---


def test_parse_silva():
    p = parse_taxon(
        "d__Bacteria;p__Acidobacteriota;c__Holophagae;o__Subgroup_7;f__Subgroup_7;g__Subgroup_7;s__",
        ";",
    )
    assert p.assigned
    assert p.ranks["тип"] == "Acidobacteriota"
    assert p.ranks["род"] == "Subgroup_7"
    assert p.domain == "Bacteria"


def test_parse_gtdb_no_species():
    p = parse_taxon(
        "d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus",
        ";",
    )
    assert p.assigned
    assert p.ranks["род"] == "Staphylococcus"


def test_parse_gg2_space_separator():
    p = parse_taxon(
        "d__Bacteria; p__Acidobacteriota; c__Thermoanaerobaculia; o__UBA5066; f__; g__; s__",
        "; ",
    )
    assert p.assigned
    assert p.ranks["тип"] == "Acidobacteriota"
    assert p.ranks["порядок"] == "UBA5066"
    assert p.ranks["семейство"] is None  # пустой суффикс
    assert p.ranks["род"] is None


def test_parse_unassigned():
    p = parse_taxon("Unassigned", ";")
    assert not p.assigned
    assert all(v is None for v in p.ranks.values())


def test_parse_domain_only():
    # Только домен, без основных рангов — не известный таксон на основных рангах.
    p = parse_taxon("d__Bacteria", ";")
    assert not p.assigned


def test_parse_unknown_family():
    p = parse_taxon(
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Gammaproteobacteria_Incertae_Sedis;f__Unknown_Family;g__Acidibacter;s__",
        ";",
    )
    assert p.assigned
    assert p.ranks["семейство"] is None  # Unknown_Family — служебная метка
    assert p.ranks["род"] == "Acidibacter"


def test_parse_unexpected_format():
    try:
        parse_taxon("d__Bacteria;x__Foo", ";")
        assert False, "должна быть ошибка"
    except ValueError:
        pass


def test_parse_duplicate_id():
    rows = [
        ("a", "d__Bacteria;p__X;c__Y;o__Z;f__W;g__V"),
        ("a", "d__Bacteria;p__X;c__Y;o__Z;f__W;g__V"),
    ]
    try:
        parse_taxonomy_rows(rows, ";")
        assert False, "должна быть ошибка на дубликате"
    except ValueError:
        pass


def test_parse_duplicate_rank():
    # Дублированный ранг p__P1;p__P2 — диагностика, а не молчаливый выбор.
    try:
        parse_taxon("d__Bacteria;p__P1;p__P2;c__C;o__O;f__F;g__G", ";")
        assert False, "должна быть ошибка на дублированном ранге"
    except ValueError:
        pass


def test_parse_duplicate_domain():
    try:
        parse_taxon("d__Bacteria;d__Archaea;p__X;c__Y;o__Z;f__W;g__V", ";")
        assert False, "должна быть ошибка на дублированном домене"
    except ValueError:
        pass


# --- Ключ родословной ---


def test_lineage_key_distinguishes_family():
    # Род G под разными семействами F1/F2 — разные ключи.
    a = _mk_tax({"a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F1", "род": "G"}}, domain="Bacteria")
    b = _mk_tax({"a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F2", "род": "G"}}, domain="Bacteria")
    ka = lineage_key(a["a1"], "род")
    kb = lineage_key(b["a1"], "род")
    assert ka != kb
    assert ka == "d__Bacteria;p__B;c__C;o__O;f__F1;g__G"
    assert kb == "d__Bacteria;p__B;c__C;o__O;f__F2;g__G"


def test_lineage_key_unassigned_rank():
    a = _mk_tax({"a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None}})
    assert lineage_key(a["a1"], "род") is None


def test_lineage_key_truncates_at_rank():
    # Ключ на ранге «тип» обрезается на типе.
    a = _mk_tax({"a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G"}}, domain="Bacteria")
    assert lineage_key(a["a1"], "тип") == "d__Bacteria;p__B"


# --- Профили и расстояние ---


def test_soil_asv_ids():
    table = _mk_table()
    soil = soil_asv_ids(table, ["s1", "s2"])
    assert soil == ["a1", "a2", "a3"]


def test_profile_mass_and_unclassified():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    p, total = profile(tax, table, "s1", "род")
    assert total == 15  # 10 + 5
    assert p == {"d__Bacteria;p__B;c__C;o__O;f__F;g__G1": 10, "UNCLASSIFIED": 5}
    assert sum(p.values()) == total


def test_profile_distance_zero_identical():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s1", "род")
    assert profile_distance(p1, p2) == 0.0


def test_profile_distance_hand_computed():
    # s1: G1=10, UNCLASSIFIED=5 (total 15). s2: G1=5, G2=10 (total 15).
    # Объединение ключей: G1, G2, UNCLASSIFIED.
    # p_s1: G1=10/15=2/3, G2=0, UNCLASSIFIED=1/3.
    # p_s2: G1=5/15=1/3, G2=10/15=2/3, UNCLASSIFIED=0.
    # D = 0.5*(|2/3-1/3| + |0-2/3| + |1/3-0|) = 0.5*(1/3+2/3+1/3) = 0.5*(4/3) = 2/3.
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s2", "род")
    d = profile_distance(p1, p2)
    assert abs(d - 2.0 / 3.0) < 1e-12, f"ожидалось 2/3, получено {d}"


def test_profile_distance_symmetry_and_bounds():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s2", "род")
    d12 = profile_distance(p1, p2)
    d21 = profile_distance(p2, p1)
    assert abs(d12 - d21) < 1e-12
    assert 0.0 <= d12 <= 1.0


def test_profile_distance_order_independent():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    }, domain="Bacteria")
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s2", "род")
    d1 = profile_distance(p1, p2)
    table2 = {k: table[k] for k in ["a3", "a1", "a2"]}
    p1b, _ = profile(tax, table2, "s1", "род")
    p2b, _ = profile(tax, table2, "s2", "род")
    d2 = profile_distance(p1b, p2b)
    assert abs(d1 - d2) < 1e-12


def test_profile_distance_renamed_ids():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    }, domain="Bacteria")
    p1, _ = profile(tax, table, "s1", "род")
    rename = {"a1": "x1", "a2": "x2", "a3": "x3"}
    table2 = {rename[k]: v for k, v in table.items()}
    tax2 = {rename[k]: v for k, v in tax.items()}
    p1b, _ = profile(tax2, table2, "s1", "род")
    assert p1 == p1b


def test_fully_unassigned_rank():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
    })
    p, total = profile(tax, table, "s1", "род")
    assert p == {"UNCLASSIFIED": 15}
    assert total == 15


def test_zero_denominator_conditional():
    table = _mk_table()
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
    })
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
    })
    m = pairwise_match(tax_a, tax_b, ["a1", "a2", "a3"], table, ["s1", "s2"], "род")
    assert m["both_assigned"] == 0
    assert m["cond_match_asv_frac"] is None  # NA
    assert m["cond_match_reads_frac"] is None  # NA


def test_pairwise_categories_cover():
    table = _mk_table()
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
    }, domain="Bacteria")
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G9"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    }, domain="Bacteria")
    m = pairwise_match(tax_a, tax_b, ["a1", "a2", "a3"], table, ["s1", "s2"], "род")
    assert m["both_match"] == 1
    assert m["both_diff"] == 1
    assert m["only_a"] == 0
    assert m["only_b"] == 1
    assert m["neither"] == 0
    assert m["total_asv"] == 3
    assert m["both_assigned"] == 2
    assert abs(m["cond_match_asv_frac"] - 0.5) < 1e-12
    # Условная доля по ридам: a1=10 (match), a2=5 (diff) в s1, a2=5 (diff) в s2.
    # Назначены обеими: a1, a2. Риды назначенных обеими = 10+5+5 = 20; совпали a1=10.
    assert abs(m["cond_match_reads_frac"] - 0.5) < 1e-12


def test_pairwise_family_distinguishes():
    # Род G под разными семействами — НЕ совпадение.
    table = _mk_table()
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F1", "род": "G"},
    }, domain="Bacteria")
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F2", "род": "G"},
    }, domain="Bacteria")
    m = pairwise_match(tax_a, tax_b, ["a1"], table, ["s1", "s2"], "род")
    assert m["both_match"] == 0
    assert m["both_diff"] == 1


def test_conditional_fractions_separate_weights():
    # Два ASV с частотами 9 и 1; первый совпал.
    # Доля ASV = 1/2 = 0.5; доля ридов = 9/10 = 0.9. Знаменатели раздельные.
    table = {
        "a1": {"s1": 9, "s2": 0},
        "a2": {"s1": 1, "s2": 0},
    }
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G9"},
    }, domain="Bacteria")
    m = pairwise_match(tax_a, tax_b, ["a1", "a2"], table, ["s1"], "род")
    assert m["both_match"] == 1
    assert m["both_diff"] == 1
    assert m["both_assigned"] == 2
    assert abs(m["cond_match_asv_frac"] - 0.5) < 1e-12
    assert abs(m["cond_match_reads_frac"] - 0.9) < 1e-12


def test_pairwise_per_sample_weights():
    # Разные образцы с неодинаковыми весами.
    table = _mk_table()
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    }, domain="Bacteria")
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G9"},
    }, domain="Bacteria")
    # s1: a1=10 (match), a2=5 (diff). s2: a1=0, a2=5 (diff).
    m1 = pairwise_match_per_sample(tax_a, tax_b, ["a1", "a2"], table, "s1", "род")
    m2 = pairwise_match_per_sample(tax_a, tax_b, ["a1", "a2"], table, "s2", "род")
    assert m1["both_match_reads"] == 10
    assert m1["both_diff_reads"] == 5
    assert m2["both_match_reads"] == 0
    assert m2["both_diff_reads"] == 5
    assert abs(m1["both_match_reads_frac"] - 10 / 15) < 1e-12
    assert abs(m2["both_match_reads_frac"] - 0.0) < 1e-12


def test_mass_preservation():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    }, domain="Bacteria")
    for sample in ["s1", "s2"]:
        p, total = profile(tax, table, sample, "род")
        assert sum(p.values()) == total


def test_duplicate_empty_domain():
    for raw in ("d__;d__Bacteria;p__P", "d__;d__;p__P"):
        try:
            parse_taxon(raw, ";")
        except ValueError:
            continue
        raise AssertionError("Repeated empty domain must be rejected")


def test_missing_rank_does_not_shift():
    value = parse_taxon("d__Bacteria;p__P;o__O;f__Unknown_Family;g__G", ";")
    assert value.ranks["класс"] is None
    assert value.ranks["семейство"] is None
    assert lineage_key(value, "род") == "d__Bacteria;p__P;o__O;g__G"


def test_asv_and_read_zero_denominators():
    table = {"a": {"s1": 9, "s2": 0}, "b": {"s1": 1, "s2": 0},
             "c": {"s1": 0, "s2": 5}}
    a = {fid: parse_taxon(raw, ";") for fid, raw in
         {"a": "g__G", "b": "g__F", "c": "Unassigned"}.items()}
    b = {fid: parse_taxon(raw, ";") for fid, raw in
         {"a": "g__G", "b": "g__X", "c": "Unassigned"}.items()}
    result = pairwise_match_per_sample(a, b, list(table), table, "s2", "род")
    assert result["both_assigned"] == 2
    assert result["both_assigned_reads"] == 0
    assert result["cond_match_asv_frac"] == 0.5
    assert result["cond_match_reads_frac"] is None


def test_five_categories_manual():
    labels_a = ["g__G", "g__F", "g__A", "Unassigned", "Unassigned"]
    labels_b = ["g__G", "g__X", "Unassigned", "g__B", "Unassigned"]
    ids = list("abcde")
    table = {k: {"soil": w} for k, w in zip(ids, [9, 1, 2, 3, 4])}
    a = {k: parse_taxon(v, ";") for k, v in zip(ids, labels_a)}
    b = {k: parse_taxon(v, ";") for k, v in zip(ids, labels_b)}
    result = pairwise_match(a, b, ids, table, ["soil"], "род")
    for category, weight in zip(("both_match", "both_diff", "only_a", "only_b", "neither"), [9, 1, 2, 3, 4]):
        assert result[category] == 1
        assert result[category + "_reads"] == weight
    assert result["total_asv"] == 5
    assert result["total_reads"] == 19
    assert result["cond_match_asv_frac"] == 0.5
    assert result["cond_match_reads_frac"] == 0.9


def test_properties_at_multiple_sizes():
    # Five properties at four sizes: conservation, axis permutation, consistent
    # ID renaming, distance symmetry/bounds, and invariance to count scaling.
    for n in (2, 3, 5, 8):
        ids = [f"a{i}" for i in range(n)]
        table = {fid: {"s1": i + 1, "s2": n - i} for i, fid in enumerate(ids)}
        tax = {fid: parse_taxon("Unassigned" if i % 3 == 0 else f"d__Bacteria;g__G{i % 2}", ";")
               for i, fid in enumerate(ids)}
        p, total = profile(tax, table, "s1", "род")
        q, total_q = profile(tax, table, "s2", "род")
        assert sum(p.values()) == total == n * (n + 1) // 2
        assert sum(q.values()) == total_q == total
        for freq in p.values():
            assert isinstance(freq, int) and freq >= 0
        permuted = {fid: dict(reversed(list(table[fid].items()))) for fid in reversed(ids)}
        assert profile(dict(reversed(list(tax.items()))), permuted, "s1", "род") == (p, total)
        renamed = {"z" + fid: counts for fid, counts in table.items()}
        assert profile({"z" + fid: value for fid, value in tax.items()}, renamed, "s1", "род") == (p, total)
        d = profile_distance(p, q)
        assert 0 <= d <= 1 and abs(d - profile_distance(q, p)) < 1e-12
        assert profile_distance(p, p) == 0
        scaled = {fid: {sid: 7 * count for sid, count in row.items()} for fid, row in table.items()}
        ps, _ = profile(tax, scaled, "s1", "род")
        qs, _ = profile(tax, scaled, "s2", "род")
        assert abs(d - profile_distance(ps, qs)) < 1e-12


# CLI tests use only disposable generated inputs. Every output mutation also
# updates its recorded SHA, so detection requires content validation, not just
# noticing that a file changed.
import argparse
import hashlib
import json
import unittest
from datetime import datetime, timezone

from classification_logging import record_command
from taxonomy_compare_io import CODE_FILES, SCHEMAS, code_identity, sha256
from taxonomy_comparison_fixtures import (FEATURES, MODELS, SAMPLES, VERSION,
    alter_member, create_fixture, matrix, write_table)

LOG_ROOT = None
COMMAND_RECORDS = []


def read_rows(path):
    with Path(path).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_rows(path, headers, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ComparisonCLI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="taxonomy-comparison-tests-")
        cls.base = Path(cls.temp.name)
        cls.inputs = create_fixture(cls.base / "inputs")
        cls.outputs = cls.base / "outputs"
        cls.input_before = {str(p.relative_to(cls.inputs)): sha256(p)
                            for p in cls.inputs.rglob("*") if p.is_file()}
        cls.call("baseline_calculation", "taxonomy_comparison.py", cls.inputs, cls.outputs, 0)
        cls.call("baseline_tables", "check_taxonomy_comparison.py", cls.inputs, cls.outputs, 0, "--tables-only")

    @classmethod
    def tearDownClass(cls):
        after = {str(p.relative_to(cls.inputs)): sha256(p)
                 for p in cls.inputs.rglob("*") if p.is_file()}
        if after != cls.input_before:
            raise AssertionError("Synthetic source inputs changed during tests")
        cls.temp.cleanup()

    @classmethod
    def call(cls, name, script, inputs, outputs, expected, *extras):
        command = [sys.executable, "-B", str(ROOT / "src" / script),
                   "--input-dir", str(inputs), "--out-dir", str(outputs), *extras]
        logs = LOG_ROOT or cls.base / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        name = f"{len(COMMAND_RECORDS):03d}_{name}"
        result = record_command(logs, name, command, ROOT)
        result["expected_exit"] = expected
        result["metadata_file"] = f"{name}.command.json"
        result["name"] = name
        result["stdout_sha256"] = sha256(logs / result["stdout"])
        result["stderr_sha256"] = sha256(logs / result["stderr"])
        (logs / f"{name}.command.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
        COMMAND_RECORDS.append(result)
        stdout = (logs / result["stdout"]).read_text()
        stderr = (logs / result["stderr"]).read_text()
        if result["exit_code"] != expected:
            raise AssertionError(f"{name}: expected {expected}, got {result['exit_code']}\n{stdout}\n{stderr}")
        if "Traceback (most recent call last)" in stdout + stderr:
            raise AssertionError(f"{name}: uncaught exception\n{stdout}\n{stderr}")
        return stdout + stderr

    def setUp(self):
        self.case = self.base / self.id().rsplit(".", 1)[-1]
        self.case.mkdir()
        self.ins = self.case / "inputs"
        self.out = self.case / "outputs"
        shutil.copytree(self.inputs, self.ins)
        shutil.copytree(self.outputs, self.out)

    def check(self, expected=1, *extras):
        output = self.call(self.id().rsplit(".", 1)[-1], "check_taxonomy_comparison.py",
                           self.ins, self.out, expected, "--tables-only", *extras)
        self.assertNotIn("CODE_CONTENT:", output,
                         "Source changed during the suite; cannot use unrelated failure as evidence")
        return output

    def rehash_output(self, name):
        path = self.out / "run.json"
        run = json.loads(path.read_text())
        run["outputs"][name] = sha256(self.out / name)
        path.write_text(json.dumps(run, ensure_ascii=False, indent=2) + "\n")

    def change_cell(self, name, column, value):
        headers, rows = read_rows(self.out / name)
        rows[0][column] = value
        write_rows(self.out / name, headers, rows)
        self.rehash_output(name)

    def test_expected_answer_from_manual_fixture(self):
        _, rows = read_rows(self.out / "pairwise_match.tsv")
        row = next(r for r in rows if r["model_a"] == MODELS[0] and r["model_b"] == MODELS[1] and r["rank"] == "род")
        for category in ("both_match", "both_diff", "only_a", "only_b"):
            self.assertEqual(int(row[category]), 1)
        self.assertEqual(int(row["neither"]), 3369)
        self.assertEqual(int(row["total_asv"]), 3373)
        self.assertEqual(int(row["both_assigned_asv"]), 2)
        self.assertEqual(float(row["cond_match_asv_frac"]), 0.5)
        _, rows = read_rows(self.out / "pairwise_match_per_sample.tsv")
        rows = [r for r in rows if r["model_a"] == MODELS[0] and r["model_b"] == MODELS[1] and r["rank"] == "род"]
        ordinary = next(r for r in rows if r["sample"] == "Kr1_S80")
        zero = next(r for r in rows if r["sample"] == "Kr2_S81")
        self.assertEqual(int(ordinary["total_reads"]), 3387)
        self.assertEqual(int(ordinary["both_assigned_reads"]), 10)
        self.assertEqual(float(ordinary["cond_match_reads_frac"]), 0.9)
        self.assertEqual(int(zero["both_assigned_reads"]), 0)
        self.assertEqual(zero["cond_match_reads_frac"], "NA")
        self.check(0)

    def test_complete_tables_without_reproduction_are_incomplete(self):
        self.call("no_acceptance", "check_taxonomy_comparison.py", self.ins, self.out, 2)

    def test_refuses_overwrite(self):
        before = {p.name: sha256(p) for p in self.out.iterdir() if p.is_file()}
        self.call("overwrite_refused", "taxonomy_comparison.py", self.ins, self.out, 2)
        self.assertEqual(before, {p.name: sha256(p) for p in self.out.iterdir() if p.is_file()})

    def test_missing_input(self):
        (self.ins / "work/qza/ASV_table.qza").unlink()
        self.check(2)

    def test_missing_summary_dependency(self):
        (self.out / "distances.tsv").unlink()
        self.check(2)

    def test_missing_all_outputs(self):
        shutil.rmtree(self.out)
        self.check(2)

    def test_bad_input_outranks_missing_output(self):
        (self.out / "distances.tsv").unlink()
        path = self.ins / "data/meta/metadata.tsv"
        path.write_text(path.read_text().replace("K_S63\tnegative_control", "K_S63\tsoil"))
        self.assertIn("SAMPLE_ROLE", self.check(1))

    def test_duplicate_metadata(self):
        path = self.ins / "data/meta/metadata.tsv"
        path.write_text(path.read_text() + "Kr1_S80\tsoil\n")
        self.assertIn("DUPLICATE", self.check())

    def test_wrong_metadata_ids(self):
        path = self.ins / "data/meta/metadata.tsv"
        path.write_text(path.read_text().replace("Kr1_S80", "unexpected"))
        self.assertIn("SAMPLE_SET", self.check())

    def test_wrong_role_without_missing_files(self):
        path = self.ins / "data/meta/metadata.tsv"
        path.write_text(path.read_text().replace("K_S63\tnegative_control", "K_S63\tsoil"))
        self.assertIn("SAMPLE_ROLE", self.check())

    def test_wrong_table_samples(self):
        samples = SAMPLES.copy()
        samples[1] = "unexpected"
        write_table(self.ins, samples=samples)
        self.assertIn("SAMPLE_SET", self.check())

    def test_wrong_qza_type(self):
        write_table(self.ins, kind="FeatureData[Taxonomy]")
        self.check()

    def test_wrong_qza_uuid(self):
        def replace_uuid(contents):
            lines = contents.decode().splitlines()
            return "\n".join("uuid: 00000000-0000-0000-0000-000000000001" if line.startswith("uuid:")
                             else line for line in lines).encode() + b"\n"
        alter_member(self.ins / "work/qza/ASV_table.qza", "metadata.yaml", replace_uuid)
        self.assertIn("QZA_UUID", self.check())

    def test_missing_asv_in_table(self):
        write_table(self.ins, values=matrix()[:-1], features=FEATURES[:-1])
        self.assertIn("ASV_COUNT", self.check())

    def test_duplicate_taxonomy_ids(self):
        alter_member(self.ins / f"work/qza/taxonomy_{MODELS[0]}.qza", "data/taxonomy.tsv",
                     lambda b: b + b.splitlines(keepends=True)[1])
        self.assertIn("DUPLICATE", self.check())

    def test_duplicate_fasta_ids(self):
        alter_member(self.ins / "work/qza/rep_seq.qza", "data/dna-sequences.fasta",
                     lambda b: b + b"\n>synthetic_0000\n" + b"A" * 180 + b"\n")
        self.assertIn("DUPLICATE", self.check())

    def test_wrong_taxonomy_ids(self):
        alter_member(self.ins / f"work/qza/taxonomy_{MODELS[0]}.qza", "data/taxonomy.tsv",
                     lambda b: b.replace(b"synthetic_0000", b"foreign_0000", 1))
        self.assertIn("ASV_SET", self.check())

    def test_bad_taxonomy_prefix(self):
        alter_member(self.ins / f"work/qza/taxonomy_{MODELS[0]}.qza", "data/taxonomy.tsv",
                     lambda b: b.replace(b"p__P", b"x__P", 1))
        self.check()

    def test_empty_output(self):
        headers, _ = read_rows(self.out / "coverage.tsv")
        write_rows(self.out / "coverage.tsv", headers, [])
        self.rehash_output("coverage.tsv")
        self.check()

    def test_empty_output_file(self):
        (self.out / "coverage.tsv").write_bytes(b"")
        self.rehash_output("coverage.tsv")
        self.assertIn("OUTPUT_HEADER", self.check())

    def test_extra_tsv_cell(self):
        path = self.out / "coverage.tsv"
        lines = path.read_text().splitlines()
        lines[1] += "\textra"
        path.write_text("\n".join(lines) + "\n")
        self.rehash_output(path.name)
        self.assertIn("OUTPUT_ROW", self.check())

    def test_consistent_scaling_keeps_fractions_but_is_wrong(self):
        path = self.out / "pairwise_match_per_sample.tsv"
        headers, rows = read_rows(path)
        for row in rows:
            for column in headers:
                if column.endswith("_reads"):
                    row[column] = str(2 * int(row[column]))
        write_rows(path, headers, rows)
        self.rehash_output(path.name)
        self.assertIn("OUTPUT_VALUE", self.check())

    def test_balanced_profile_transfer_is_wrong(self):
        path = self.out / "profiles.tsv"
        headers, rows = read_rows(path)
        group = tuple(rows[0][key] for key in ("model", "sample", "rank"))
        same = [row for row in rows if tuple(row[key] for key in ("model", "sample", "rank")) == group]
        donor = next(row for row in same if int(row["count"]) > 1)
        receiver = next(row for row in same if row is not donor)
        before = sum(int(row["count"]) for row in same)
        for row, delta in ((donor, -1), (receiver, 1)):
            row["count"] = str(int(row["count"]) + delta)
            row["frac"] = format(int(row["count"]) / int(row["total_reads"]), ".17g")
        self.assertEqual(before, sum(int(row["count"]) for row in same))
        write_rows(path, headers, rows)
        self.rehash_output(path.name)
        self.assertIn("OUTPUT_VALUE", self.check())

    def test_duplicate_output_key(self):
        headers, rows = read_rows(self.out / "profiles.tsv")
        rows.append(rows[0])
        write_rows(self.out / "profiles.tsv", headers, rows)
        self.rehash_output("profiles.tsv")
        self.check()

    def test_mutation_fail_then_restoration_pass(self):
        path = self.out / "pairwise_match_per_sample.tsv"
        original = path.read_bytes()
        self.change_cell(path.name, "both_match_reads", "999")
        self.assertIn("OUTPUT_VALUE", self.check(1))
        path.write_bytes(original)
        self.rehash_output(path.name)
        self.check(0)

    def test_missing_output_row(self):
        headers, rows = read_rows(self.out / "distances.tsv")
        write_rows(self.out / "distances.tsv", headers, rows[:-1])
        self.rehash_output("distances.tsv")
        self.check()

    def test_bad_snapshot_identity(self):
        path = self.out / "run.json"
        value = json.loads(path.read_text())
        value["code"]["version"] = "snapshot-" + "f" * 64
        path.write_text(json.dumps(value))
        self.check()

    def test_older_commit_without_source_is_invalid(self):
        old = subprocess.run(["git", "rev-parse", "70c5ede^{commit}"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout.strip()
        path = self.out / "run.json"
        value = json.loads(path.read_text())
        value["code"].update(version=old, commit=old, clean=True)
        path.write_text(json.dumps(value))
        self.assertIn("CODE_COMMIT", self.check())

    def test_valid_row_reordering(self):
        for name in SCHEMAS:
            headers, rows = read_rows(self.out / name)
            write_rows(self.out / name, headers, list(reversed(rows)))
            self.rehash_output(name)
        self.check(0)


def _domain_case(value):
    def test(self):
        values = matrix()
        values[0, 1] = value
        write_table(self.ins, values=values)
        self.assertIn("COUNT_DOMAIN", self.check())
    return test


for _name, _value in (("negative", -1), ("fraction", 0.5), ("nan", float("nan")), ("infinity", float("inf"))):
    setattr(ComparisonCLI, "test_frequency_" + _name, _domain_case(_value))


def _output_case(filename, column, value):
    def test(self):
        self.change_cell(filename, column, value)
        self.assertIn("OUTPUT_VALUE", self.check())
    return test


for _name, _filename, _column, _value in (
    ("raw_taxon", "assignments.tsv", "sklearn_silva138_taxon", "Unassigned"),
    ("rank_label", "assignments.tsv", "sklearn_silva138_род", "Other"),
    ("lineage", "assignments.tsv", "sklearn_silva138_род_key", "g__Other"),
    ("status", "assignments.tsv", "sklearn_silva138_род_assigned", "0"),
    ("soil_flag", "assignments.tsv", "soil_asv", "0"),
    ("coverage", "coverage.tsv", "assigned_asv", "999"),
    ("coverage_reads", "coverage_per_sample.tsv", "assigned_reads", "999"),
    ("asv_conditional", "pairwise_match.tsv", "cond_match_asv_frac", "0.123"),
    ("read_conditional", "pairwise_match_per_sample.tsv", "cond_match_reads_frac", "0.123"),
    ("read_denominator", "pairwise_match_per_sample.tsv", "both_assigned_reads", "999"),
    ("both_match_reads", "pairwise_match_per_sample.tsv", "both_match_reads", "999"),
    ("both_diff_reads", "pairwise_match_per_sample.tsv", "both_diff_reads", "999"),
    ("only_a_reads", "pairwise_match_per_sample.tsv", "only_a_reads", "999"),
    ("only_b_reads", "pairwise_match_per_sample.tsv", "only_b_reads", "999"),
    ("neither_reads", "pairwise_match_per_sample.tsv", "neither_reads", "999"),
    ("total_reads", "pairwise_match_per_sample.tsv", "total_reads", "999"),
    ("total_asv", "pairwise_match.tsv", "total_asv", "999"),
    ("frequency", "profiles.tsv", "count", "0.5"),
    ("profile_fraction", "profiles.tsv", "frac", "NaN"),
    ("distance", "distances.tsv", "distance", "0.123"),
    ("distance_coverage", "distances.tsv", "assigned_reads_frac_a", "0.123"),
    ("summary", "summary.tsv", "median", "0.123"),
    ("manifest_hash", "manifest.tsv", "sha256", "0" * 64),
    ("manifest_uuid", "manifest.tsv", "uuid", "00000000-0000-0000-0000-000000000001"),
):
    setattr(ComparisonCLI, "test_mutated_" + _name, _output_case(_filename, _column, _value))


def run_all(log_dir=None):
    global LOG_ROOT
    if log_dir is not None:
        LOG_ROOT = Path(log_dir)
        LOG_ROOT.mkdir(parents=True, exist_ok=False)
    start = datetime.now(timezone.utc).isoformat()
    suite = unittest.TestSuite()
    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            suite.addTest(unittest.FunctionTestCase(test))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(ComparisonCLI))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    summary = {"schema": "taxonomy-comparison-tests-v1", "fixture_version": VERSION,
               "start_utc": start, "end_utc": datetime.now(timezone.utc).isoformat(),
               "tests_run": result.testsRun, "failures": len(result.failures),
               "errors": len(result.errors), "skipped": len(result.skipped), "success": result.wasSuccessful(),
               "code": code_identity(),
               "sources": {name: sha256(ROOT / name) for name in CODE_FILES},
               "commands": COMMAND_RECORDS,
               "details": [{"test": str(test), "traceback": error}
                           for test, error in result.failures + result.errors]}
    if LOG_ROOT:
        (LOG_ROOT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(f"{'PASS/0' if result.wasSuccessful() else 'FAIL/1'} — {result.testsRun} tests")
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-dir", type=Path, help="New directory for complete command evidence")
    args = parser.parse_args()
    sys.exit(run_all(args.log_dir))



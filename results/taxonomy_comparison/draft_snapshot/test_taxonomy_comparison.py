"""Синтетические примеры и отрицательные случаи для сравнения таксономий.

Ожидаемые ответы вычислены вручную, не той же функцией, что проверяется.
Все внесённые ошибки — только во временных копиях.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from taxonomy_compare_core import (  # noqa: E402
    MAIN_RANKS,
    ParsedTaxon,
    parse_taxon,
    parse_taxonomy_rows,
    profile,
    profile_distance,
    soil_asv_ids,
)

# --- Разбор строк Taxon ---


def test_parse_silva():
    p = parse_taxon(
        "d__Bacteria;p__Acidobacteriota;c__Holophagae;o__Subgroup_7;f__Subgroup_7;g__Subgroup_7;s__",
        ";",
    )
    assert p.assigned
    assert p.ranks["тип"] == "Acidobacteriota"
    assert p.ranks["род"] == "Subgroup_7"


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
    # Неожиданный фрагмент — явная диагностика.
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


# --- Профили и расстояние ---


def _mk_table():
    # 3 ASV, 2 почвенных образца.
    return {
        "a1": {"s1": 10, "s2": 0},
        "a2": {"s1": 5, "s2": 5},
        "a3": {"s1": 0, "s2": 10},
    }


def _mk_tax(assignments):
    # assignments: {fid: {rank: label or None}}
    out = {}
    for fid, ranks in assignments.items():
        out[fid] = ParsedTaxon(raw="", ranks=ranks, assigned=any(v is not None for v in ranks.values()))
    return out


def test_soil_asv_ids():
    table = _mk_table()
    soil = soil_asv_ids(table, ["s1", "s2"])
    assert soil == ["a1", "a2", "a3"]


def test_profile_mass_and_unclassified():
    table = _mk_table()
    # a1 и a3 назначены на роде, a2 — нет.
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    })
    p, total = profile(tax, table, "s1", "род")
    assert total == 15  # 10 + 5
    assert p == {"G1": 10, "UNCLASSIFIED": 5}
    # Сумма частот = total.
    assert sum(p.values()) == total


def test_profile_distance_zero_identical():
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
    })
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
    })
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
    })
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s2", "род")
    d12 = profile_distance(p1, p2)
    d21 = profile_distance(p2, p1)
    assert abs(d12 - d21) < 1e-12
    assert 0.0 <= d12 <= 1.0


def test_profile_distance_order_independent():
    # Перестановка строк таблицы не меняет расстояние.
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    })
    p1, _ = profile(tax, table, "s1", "род")
    p2, _ = profile(tax, table, "s2", "род")
    d1 = profile_distance(p1, p2)
    # Переставленная таблица.
    table2 = {k: table[k] for k in ["a3", "a1", "a2"]}
    p1b, _ = profile(tax, table2, "s1", "род")
    p2b, _ = profile(tax, table2, "s2", "род")
    d2 = profile_distance(p1b, p2b)
    assert abs(d1 - d2) < 1e-12


def test_profile_distance_renamed_ids():
    # Согласованное переименование ASV ID не меняет профиль.
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    })
    p1, _ = profile(tax, table, "s1", "род")
    # Переименованные ID.
    rename = {"a1": "x1", "a2": "x2", "a3": "x3"}
    table2 = {rename[k]: v for k, v in table.items()}
    tax2 = {rename[k]: v for k, v in tax.items()}
    p1b, _ = profile(tax2, table2, "s1", "род")
    assert p1 == p1b


def test_fully_unassigned_rank():
    # Полностью неназначенный ранг: вся масса в UNCLASSIFIED.
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
    # Нулевой знаменатель условного совпадения -> NA (None), не 0.
    # Обе модели полностью неназначены на роде.
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
    from taxonomy_compare_core import pairwise_match
    m = pairwise_match(tax_a, tax_b, ["a1", "a2", "a3"], table, ["s1", "s2"], "род")
    assert m["both_assigned"] == 0
    assert m["cond_match_frac"] is None  # NA


def test_pairwise_categories_cover():
    # Пять категорий покрывают весь набор.
    table = _mk_table()
    tax_a = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": None},
    })
    tax_b = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G9"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    })
    from taxonomy_compare_core import pairwise_match
    m = pairwise_match(tax_a, tax_b, ["a1", "a2", "a3"], table, ["s1", "s2"], "род")
    # a1: both_match; a2: both_diff; a3: only_b.
    assert m["both_match"] == 1
    assert m["both_diff"] == 1
    assert m["only_a"] == 0
    assert m["only_b"] == 1
    assert m["neither"] == 0
    assert m["total_asv"] == 3
    assert m["both_assigned"] == 2
    assert abs(m["cond_match_frac"] - 0.5) < 1e-12


def test_mass_preservation():
    # Сохранение массы: сумма частот профиля = total_reads образца.
    table = _mk_table()
    tax = _mk_tax({
        "a1": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G1"},
        "a2": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G2"},
        "a3": {"тип": "B", "класс": "C", "порядок": "O", "семейство": "F", "род": "G3"},
    })
    for sample in ["s1", "s2"]:
        p, total = profile(tax, table, sample, "род")
        assert sum(p.values()) == total


def run_all():
    tests = [
        test_parse_silva, test_parse_gtdb_no_species, test_parse_gg2_space_separator,
        test_parse_unassigned, test_parse_domain_only, test_parse_unknown_family,
        test_parse_unexpected_format, test_parse_duplicate_id,
        test_soil_asv_ids, test_profile_mass_and_unclassified,
        test_profile_distance_zero_identical, test_profile_distance_hand_computed,
        test_profile_distance_symmetry_and_bounds, test_profile_distance_order_independent,
        test_profile_distance_renamed_ids, test_fully_unassigned_rank,
        test_zero_denominator_conditional, test_pairwise_categories_cover,
        test_mass_preservation,
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  ok: {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL: {t.__name__}: {e}")
        except Exception as e:
            failed += 1
            print(f"  ERROR: {t.__name__}: {e!r}")
    if failed:
        print(f"FAIL/1 — {failed} тестов не прошли")
        return 1
    print("PASS/0 — все синтетические примеры прошли")
    return 0


if __name__ == "__main__":
    sys.exit(run_all())
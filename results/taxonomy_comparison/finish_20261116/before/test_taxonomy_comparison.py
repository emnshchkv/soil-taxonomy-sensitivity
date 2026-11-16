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


# --- Интеграционные испытания читателя QZA и CLI ---


def _make_qza(tmp: Path, name: str, uuid: str, data_files: dict):
    """Создаёт минимальный QZA (ZIP) с metadata.yaml и data-файлами."""
    qza = tmp / name
    with zipfile.ZipFile(qza, "w") as z:
        z.writestr(f"{uuid}/metadata.yaml", f"uuid: {uuid}\ntype: test\nformat: test\n")
        for rel, content in data_files.items():
            z.writestr(f"{uuid}/data/{rel}", content)
    return qza


def test_qza_reader_biom():
    """Интеграционный тест читателя BIOM из QZA (ZIP)."""
    import biom
    from taxonomy_comparison import read_biom_from_qza
    tmp = Path(tempfile.mkdtemp())
    try:
        # Создаём BIOM в памяти (плотная матрица).
        t = biom.Table([[1, 2], [3, 4]], ["a1", "a2"], ["s1", "s2"],
                       input_is_dense=True)
        biom_path = tmp / "table.biom"
        with biom.util.biom_open(str(biom_path), "w") as f:
            t.to_hdf5(f, "test")
        qza = _make_qza(tmp, "table.qza", "u1",
                        {"feature-table.biom": biom_path.read_bytes()})
        table, fids, sids = read_biom_from_qza(qza)
        assert fids == ["a1", "a2"]
        assert sids == ["s1", "s2"]
        assert table["a1"] == {"s1": 1, "s2": 2}
        assert table["a2"] == {"s1": 3, "s2": 4}
    finally:
        shutil.rmtree(tmp)


def test_qza_reader_taxonomy():
    from taxonomy_comparison import read_taxonomy_from_qza
    tmp = Path(tempfile.mkdtemp())
    try:
        content = "Feature ID\tTaxon\tConfidence\na1\td__Bacteria;p__X;c__Y;o__Z;f__W;g__V\t0.9\n"
        qza = _make_qza(tmp, "tax.qza", "u2", {"taxonomy.tsv": content})
        tax = read_taxonomy_from_qza(qza, ";")
        assert set(tax) == {"a1"}
        assert tax["a1"].ranks["род"] == "V"
    finally:
        shutil.rmtree(tmp)


def test_cli_refuses_overwrite():
    """CLI отказывается затирать существующий расчёт без --force."""
    tmp = Path(tempfile.mkdtemp())
    try:
        (tmp / "existing.txt").write_text("x")
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")
        r = subprocess.run(
            [sys.executable, str(ROOT / "src/taxonomy_comparison.py"),
             "--out-dir", str(tmp)],
            capture_output=True, text=True, env=env,
        )
        assert r.returncode == 2, f"ожидался отказ (2), получено {r.returncode}"
        assert "ОТКАЗ" in r.stderr
    finally:
        shutil.rmtree(tmp)


def test_cli_accepts_force():
    """CLI с --force перезаписывает существующий расчёт."""
    tmp = Path(tempfile.mkdtemp())
    try:
        (tmp / "existing.txt").write_text("x")
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")
        r = subprocess.run(
            [sys.executable, str(ROOT / "src/taxonomy_comparison.py"),
             "--out-dir", str(tmp), "--force"],
            capture_output=True, text=True, env=env,
        )
        # Может упасть на реальных входах, но не должен отказать из-за непустого каталога.
        assert "ОТКАЗ" not in r.stderr
    finally:
        shutil.rmtree(tmp)


def run_all():
    tests = [
        test_parse_silva, test_parse_gtdb_no_species, test_parse_gg2_space_separator,
        test_parse_unassigned, test_parse_domain_only, test_parse_unknown_family,
        test_parse_unexpected_format, test_parse_duplicate_id,
        test_parse_duplicate_rank, test_parse_duplicate_domain,
        test_lineage_key_distinguishes_family, test_lineage_key_unassigned_rank,
        test_lineage_key_truncates_at_rank,
        test_soil_asv_ids, test_profile_mass_and_unclassified,
        test_profile_distance_zero_identical, test_profile_distance_hand_computed,
        test_profile_distance_symmetry_and_bounds, test_profile_distance_order_independent,
        test_profile_distance_renamed_ids, test_fully_unassigned_rank,
        test_zero_denominator_conditional, test_pairwise_categories_cover,
        test_pairwise_family_distinguishes, test_pairwise_per_sample_weights,
        test_conditional_fractions_separate_weights,
        test_mass_preservation,
        test_qza_reader_biom, test_qza_reader_taxonomy,
        test_cli_refuses_overwrite, test_cli_accepts_force,
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
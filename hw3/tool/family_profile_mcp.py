"""family_profile_mcp — MCP-инструмент «family-profile» для проекта бактериальных сообществ.

Один инструмент: `family_profile` читает готовую таблицу
`results/paper_step4/20261118/family_sample_counts.tsv` и возвращает для одной пробы
счётчики ридов по семействам для трёх моделей таксономии (sklearn_silva138,
sklearn_gtdb, sklearn_gg2). Нужен, когда требуется назвать числа к рисунку статьи
без ручного поиска строки в TSV и без риска перепутать модель или пробу.

Источник данных только читается: файл открывается на чтение, ничего не пишется.

Проверка без агента:  python3 hw3/tool/family_profile_mcp.py --selftest
Вызов из командной строки: python3 hw3/tool/family_profile_mcp.py --sample Kr2_S81
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

try:                       # пакет mcp 1.x
    from mcp.server.fastmcp import FastMCP
except ImportError:
    try:                   # пакет mcp 2.x
        from mcp.server.mcpserver import MCPServer as FastMCP
    except ImportError:    # без пакета mcp работают режимы --selftest и --sample
        class FastMCP:     # type: ignore[no-redef]
            def __init__(self, name: str) -> None:
                self.name = name

            def tool(self):
                return lambda f: f

            def run(self) -> None:
                sys.exit("нет пакета mcp: запускайте интерпретатором из окружения учебного агента")

if sys.platform == "win32":  # протокол MCP — UTF-8
    for _stream in (sys.stdin, sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
        except Exception:
            pass

mcp = FastMCP("family-profile")

# относительный путь к таблице от корня проекта
TSV_REL = Path("results") / "paper_step4" / "20261118" / "family_sample_counts.tsv"
MODELS: tuple[str, ...] = ("sklearn_silva138", "sklearn_gtdb", "sklearn_gg2")

INT_FIELDS: tuple[str, ...] = ("post_dada2_reads", "selected_reads", "unclassified_family_reads")
PCT_FIELD = "unclassified_pct_of_selected"


def default_tsv_path() -> Path:
    """Путь к таблице: переменная окружения, затем корень проекта от файла модуля, затем cwd."""
    env = os.environ.get("FAMILY_PROFILE_TSV")
    if env:
        return Path(env).expanduser()
    here = Path(__file__).resolve()
    candidates = [here.parents[2] / TSV_REL, Path.cwd() / TSV_REL]
    for cand in candidates:
        if cand.is_file():
            return cand
    return candidates[0]


def _read_rows(tsv_path: Path) -> list[dict[str, str]]:
    """Прочитать TSV с разделителем табуляции; ошибки чтения — исключением."""
    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in ("model", "sample", *INT_FIELDS, PCT_FIELD) if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"в таблице нет столбцов: {', '.join(missing)}")
        return [dict(row) for row in reader]


def _as_int(row: dict[str, str], field: str, where: str) -> int:
    raw = (row.get(field) or "").strip()
    try:
        return int(raw)
    except ValueError:
        raise ValueError(f"{where}: поле {field} не целое число: {raw!r}") from None


def profile_at(sample_id: str, tsv_path: Path) -> dict:
    """Собрать ответ инструмента по явному пути к таблице (используется и в самопроверке)."""
    if not isinstance(sample_id, str) or not sample_id.strip():
        return {"error": "не задан sample_id: передайте идентификатор пробы, например Kr2_S81"}
    sample_id = sample_id.strip()

    if not tsv_path.is_file():
        return {"error": f"нет таблицы {tsv_path} — проверьте путь к results/paper_step4/20261118/family_sample_counts.tsv"}

    try:
        rows = _read_rows(tsv_path)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        return {"error": f"таблица {tsv_path} не читается: {exc}"}

    selected = [r for r in rows if (r.get("sample") or "").strip() == sample_id]
    if not selected:
        known = sorted({(r.get("sample") or "").strip() for r in rows})
        return {"error": f"неизвестная проба {sample_id!r}: в таблице нет строк с таким sample",
                "known_samples": known}

    records: list[dict] = []
    for model in MODELS:
        hits = [r for r in selected if (r.get("model") or "").strip() == model]
        if not hits:
            return {"error": f"для пробы {sample_id} нет строки модели {model}"}
        if len(hits) > 1:
            return {"error": f"для пробы {sample_id} и модели {model} найдено строк: {len(hits)}, ожидалась одна"}
        row = hits[0]
        where = f"{sample_id}/{model}"
        try:
            counts = {field: _as_int(row, field, where) for field in INT_FIELDS}
        except ValueError as exc:
            return {"error": str(exc)}

        pct = (row.get(PCT_FIELD) or "").strip()
        if counts["selected_reads"] == 0:
            # при нуле отобранных ридов процент неопределим
            pct = "NA"
        elif pct in ("", "NA"):
            return {"error": f"{where}: при {counts['selected_reads']} отобранных ридах процент не задан"}

        records.append({
            "model": model,
            "post_dada2_reads": counts["post_dada2_reads"],
            "selected_reads": counts["selected_reads"],
            "unclassified_family_reads": counts["unclassified_family_reads"],
            "unclassified_pct_of_selected": pct,
        })

    return {
        "sample_id": sample_id,
        "source": TSV_REL.as_posix(),
        "records": records,
    }


@mcp.tool()
def family_profile(sample_id: str) -> dict:
    """Данные к рисунку статьи: счётчики ридов по семействам для одной пробы почвы или контроля.

    Вызывать, когда нужно назвать числа к рисунку или таблице статьи по семействам:
    сколько ридов осталось после DADA2, сколько отобрано для профиля семейств, сколько
    ридов осталось без назначенного семейства и какой это процент от отобранных, —
    сразу для трёх моделей таксономии (sklearn_silva138, sklearn_gtdb, sklearn_gg2).

    sample_id — идентификатор пробы ровно как в столбце sample таблицы, например
    "Kr2_S81", "LUB2_S92", "K_S63" (отрицательный контроль). Совпадение точное,
    регистр важен; список допустимых значений — в поле known_samples ответа на ошибку.

    Возвращает словарь: sample_id, source (относительный путь к исходному TSV) и records —
    три записи по одной на модель с полями model, post_dada2_reads, selected_reads,
    unclassified_family_reads (целые числа) и unclassified_pct_of_selected (строка из TSV
    с сохранённой точностью; "NA", если отобранных ридов ноль и процент неопределим).
    При неизвестной пробе, отсутствии таблицы, дубликате или пропуске строки модели
    возвращает словарь с полем error. Только чтение, файлы проекта не меняются.
    """
    return profile_at(sample_id, default_tsv_path())


# ---------------------------------------------------------------- самопроверка

def _check(cond: bool, message: str) -> None:
    if not cond:
        raise SystemExit(f"FAIL: {message}")


def _triple(record: dict) -> tuple[int, int, str]:
    return (record["selected_reads"], record["unclassified_family_reads"],
            record["unclassified_pct_of_selected"])


def selftest() -> None:
    """Проверка на известных значениях таблицы и на подменённых таблицах во временном каталоге."""
    import tempfile

    tsv = default_tsv_path()
    _check(tsv.is_file(), f"нет таблицы {tsv}")

    kr2 = profile_at("Kr2_S81", tsv)
    _check("error" not in kr2, f"Kr2_S81: {kr2.get('error')}")
    _check(kr2["sample_id"] == "Kr2_S81", "Kr2_S81: неверный sample_id")
    _check(kr2["source"] == TSV_REL.as_posix(), f"Kr2_S81: неверный source {kr2['source']}")
    _check([r["model"] for r in kr2["records"]] == list(MODELS), "Kr2_S81: не тот набор моделей")
    for rec in kr2["records"]:
        _check(rec["post_dada2_reads"] == 8904, f"Kr2_S81/{rec['model']}: post_dada2_reads {rec['post_dada2_reads']} != 8904")
        _check(isinstance(rec["post_dada2_reads"], int), "post_dada2_reads не int")
        _check(isinstance(rec["selected_reads"], int), "selected_reads не int")
        _check(isinstance(rec["unclassified_family_reads"], int), "unclassified_family_reads не int")
        _check(isinstance(rec["unclassified_pct_of_selected"], str), "процент должен быть строкой")
    expected = {
        "sklearn_silva138": (0, 0, "NA"),
        "sklearn_gtdb": (8013, 8013, "100"),
        "sklearn_gg2": (8098, 8098, "100"),
    }
    for rec in kr2["records"]:
        _check(_triple(rec) == expected[rec["model"]],
               f"Kr2_S81/{rec['model']}: {_triple(rec)} != {expected[rec['model']]}")

    lub2 = profile_at("LUB2_S92", tsv)
    _check("error" not in lub2, f"LUB2_S92: {lub2.get('error')}")
    silva = next(r for r in lub2["records"] if r["model"] == "sklearn_silva138")
    _check(silva["selected_reads"] == 25197, f"LUB2_S92/SILVA: selected_reads {silva['selected_reads']} != 25197")
    _check(silva["unclassified_family_reads"] == 2182,
           f"LUB2_S92/SILVA: unclassified_family_reads {silva['unclassified_family_reads']} != 2182")
    _check(silva["unclassified_pct_of_selected"].startswith("8.65976108266857165535579632495932055403421042187562011350558"),
           f"LUB2_S92/SILVA: процент потерял точность: {silva['unclassified_pct_of_selected']}")

    unknown = profile_at("NO_SUCH_SAMPLE_42", tsv)
    _check("error" in unknown, "неизвестная проба не дала ошибку")
    _check("NO_SUCH_SAMPLE_42" in unknown["error"], "в ошибке нет идентификатора пробы")

    # подменённые таблицы: отсутствие файла, дубликат, пропуск строки модели
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        missing = profile_at("Kr2_S81", tmpdir / "absent.tsv")
        _check("error" in missing and "absent.tsv" in missing["error"], "отсутствие таблицы не дало ошибку")

        header = ["scope", "model", "sample", "sample_type", "post_dada2_reads", "selected_reads",
                  "excluded_reads", "selected_bacteria_reads", "selected_archaea_reads",
                  "named_family_reads", "unclassified_family_reads", "selected_pct_of_post_dada2",
                  "unclassified_pct_of_selected", "zero_selected"]

        def write(name: str, rows: list[list[str]]) -> Path:
            path = tmpdir / name
            with path.open("w", encoding="utf-8", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(header)
                writer.writerows(rows)
            return path

        def row(model: str, sample: str = "X1") -> list[str]:
            return ["family_prokaryote_selected", model, sample, "soil", "100", "90", "10",
                    "90", "0", "80", "10", "90", "11.1111111111111111111111111111111111111111111111111111111111", "0"]

        dup = write("dup.tsv", [row(m) for m in MODELS] + [row("sklearn_gtdb")])
        res = profile_at("X1", dup)
        _check("error" in res and "найдено строк" in res["error"], f"дубликат не дал ошибку: {res}")

        gap = write("gap.tsv", [row("sklearn_silva138"), row("sklearn_gg2")])
        res = profile_at("X1", gap)
        _check("error" in res and "sklearn_gtdb" in res["error"], f"пропуск модели не дал ошибку: {res}")

        ok = write("ok.tsv", [row(m) for m in MODELS])
        res = profile_at("X1", ok)
        _check("error" not in res, f"корректная таблица дала ошибку: {res}")
        _check(_triple(res["records"][0]) == (90, 10, "11.1111111111111111111111111111111111111111111111111111111111"),
               f"корректная таблица: {_triple(res['records'][0])}")

        # ноль отобранных ридов: процент неопределим, даже если в TSV стоит число
        def zero_row(model: str) -> list[str]:
            return ["family_prokaryote_selected", model, "X2", "soil", "100", "0", "100",
                    "0", "0", "0", "0", "0", "0", "1"]

        zero = write("zero.tsv", [zero_row(m) for m in MODELS])
        res = profile_at("X2", zero)
        _check("error" not in res, f"нулевая проба дала ошибку: {res}")
        for rec in res["records"]:
            _check(_triple(rec) == (0, 0, "NA"), f"нулевая проба/{rec['model']}: {_triple(rec)} != (0, 0, 'NA')")

    print("family_profile: проверки пройдены")
    print("  Kr2_S81:", json.dumps(kr2, ensure_ascii=False))
    print("  LUB2_S92/SILVA:", json.dumps(silva, ensure_ascii=False))
    print("  неизвестная проба:", json.dumps(unknown, ensure_ascii=False))


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    elif "--sample" in sys.argv:        # тот же инструмент из командной строки: --sample <id>
        print(json.dumps(family_profile(sys.argv[sys.argv.index("--sample") + 1]), ensure_ascii=False, indent=1))
    else:
        mcp.run()

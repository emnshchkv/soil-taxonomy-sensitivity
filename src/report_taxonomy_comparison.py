"""Render or check a report from completed, saved comparison evidence only."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BASE = ROOT / "results/taxonomy_comparison/finish_20261116"
DEFAULT_OUTPUT = ROOT / "reports/Сравнение таксономических профилей.md"
MODELS = ("sklearn_silva138", "sklearn_gtdb", "sklearn_gg2")
LABELS = dict(zip(MODELS, ("SILVA", "GTDB", "Greengenes2")))
RANKS = ("тип", "класс", "порядок", "семейство", "род")
PAIRS = ((MODELS[0], MODELS[1]), (MODELS[0], MODELS[2]), (MODELS[1], MODELS[2]))


def rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wiki(path):
    return "[[" + str(Path(path).resolve().relative_to(ROOT)) + "]]"


def table(headers, records):
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    lines.extend("| " + " | ".join(map(str, row)) + " |" for row in records)
    return "\n".join(lines)


def fraction(value):
    return "NA" if value is None else f"{value:.6f}"


def spread(values, percent=False):
    values = [float(v) for v in values if v not in (None, "NA")]
    if not values:
        return "NA"
    scale = 100 if percent else 1
    form = ".2f" if percent else ".6f"
    output = [format(scale * v, form) for v in (statistics.median(values), min(values), max(values))]
    return f"{output[0]} [{output[1]}; {output[2]}]" + ("%" if percent else "")


def render(base=DEFAULT_BASE):
    base = Path(base).resolve()
    main = base / "main"
    run = read_json(main / "run.json")
    acceptance = read_json(main / "acceptance.json")
    tests = read_json(base / "tests_final/summary.json")
    checkout = read_json(base / "clean_checkout.json")
    reproduction = read_json(base / "reproduction_comparison.json")
    protected = read_json(base / "protected_after.json")
    final_check = read_json(base / "execution_logs/full_acceptance.command.json")
    if final_check["exit_code"] != 0 or not tests["success"] or tests["failures"] or tests["errors"]:
        raise ValueError("Full acceptance and all synthetic tests must pass before rendering")
    if reproduction["violations"] or protected["violations"]:
        raise ValueError("Reproduction/protected inputs contain violations")
    if checkout["status_before"] or checkout["status_after"] or not checkout["readonly_inputs"]:
        raise ValueError("Reproduction was not from a clean checkout with read-only inputs")
    if not run["code"]["clean"] or run["code"]["commit"] != acceptance["code_commit"]:
        raise ValueError("Comparison source was not committed")
    for name, sha in run["outputs"].items():
        if digest(main / name) != sha:
            raise ValueError(f"Changed comparison output: {name}")
    inputs_before = read_json(base / "inputs_before.json")
    if inputs_before != read_json(base / "inputs_after.json"):
        raise ValueError("Input identities changed")
    metadata_path = ROOT / run["inputs"]["metadata"]["path"]
    if digest(metadata_path) != run["inputs"]["metadata"]["sha256"]:
        raise ValueError("Metadata no longer match the comparison input")
    metadata = rows(metadata_path)
    soil_samples = sorted(row["sample-id"] for row in metadata if row["sample_type"] == "soil")
    assignments = rows(main / "assignments.tsv")
    total_asv = len(assignments)
    soil_asv = sum(int(row["soil_asv"]) for row in assignments)
    coverage = rows(main / "coverage.tsv")
    weighted = rows(main / "coverage_per_sample.tsv")
    distances = rows(main / "summary.tsv")
    sample_distances = rows(main / "distances.tsv")
    pairs = rows(main / "pairwise_match.tsv")
    weighted_pairs = rows(main / "pairwise_match_per_sample.tsv")
    denoise_path = ROOT / "results/exact_match_reproduction/provenance/samples.tsv"
    kr5 = next(row for row in rows(denoise_path) if row["run"] == "1" and row["sample-id"] == "Kr5_S84")
    retained, incoming = int(kr5["non-chimeric"]), int(kr5["input"])
    retention = 100 * retained / incoming
    source = run["code"]["commit"]
    protocol = acceptance["protocol_commit"]
    date = datetime.fromisoformat(final_check["end_utc"]).date().isoformat()

    coverage_rows = []
    for model in MODELS:
        for rank in RANKS:
            item = next(r for r in coverage if r["model"] == model and r["rank"] == rank)
            sample_rows = [r for r in weighted if r["model"] == model and r["rank"] == rank]
            if sorted(r["sample"] for r in sample_rows) != soil_samples:
                raise ValueError("Coverage sample set differs")
            coverage_rows.append([LABELS[model], rank, f"{item['assigned_asv']} / {item['total_asv']}",
                                  f"{100 * float(item['assigned_asv_frac']):.2f}%",
                                  spread([r["assigned_reads_frac"] for r in sample_rows], True)])

    distance_rows = []
    conditional_rows = []
    unknown_zero_groups = []
    for a, b in PAIRS:
        label = LABELS[a] + "–" + LABELS[b]
        zero_samples = sorted(r["sample"] for r in sample_distances
                              if r["model_a"] == a and r["model_b"] == b and r["rank"] == "род"
                              and float(r["distance"]) == 0 and float(r["assigned_reads_frac_a"]) == 0
                              and float(r["assigned_reads_frac_b"]) == 0)
        if zero_samples:
            unknown_zero_groups.append(label + ": " + ", ".join(zero_samples))
        for rank in RANKS:
            item = next(r for r in distances if r["model_a"] == a and r["model_b"] == b and r["rank"] == rank)
            distance_rows.append([label, rank, item["n_samples"],
                                  fraction(float(item["median"])), fraction(float(item["min"])),
                                  fraction(float(item["max"]))])
        item = next(r for r in pairs if r["model_a"] == a and r["model_b"] == b and r["rank"] == "род")
        per_sample = [r for r in weighted_pairs if r["model_a"] == a and r["model_b"] == b and r["rank"] == "род"]
        values = [r["cond_match_reads_frac"] for r in per_sample]
        asv_fraction = None if item["cond_match_asv_frac"] == "NA" else float(item["cond_match_asv_frac"])
        conditional_rows.append([label, f"{item['both_match']} / {item['both_assigned_asv']}",
                                 fraction(asv_fraction), spread(values),
                                 f"{sum(v != 'NA' for v in values)} / {len(per_sample)}"])

    restore_commands = [c for c in tests["commands"]
                        if c["name"].split("_", 1)[-1] == "test_mutation_fail_then_restoration_pass"]
    if [c["exit_code"] for c in restore_commands] != [1, 0]:
        raise ValueError("Missing recorded mutation FAIL -> repair PASS")
    failure_links = " → ".join(f"`{'FAIL/1' if c['exit_code'] == 1 else 'PASS/0'}` " +
        wiki(base / "tests_final" / c["stdout"]) for c in restore_commands)
    command_labels = {"tests": "Синтетические испытания", "main": "Основной расчёт",
                      "tables_main": "Проверка таблиц", "reproduction": "Расчёт из чистого коммита",
                      "tables_reproduction": "Проверка воспроизведённых таблиц",
                      "project_check": "Общая проверка проекта"}
    timing_rows = []
    calculation_seconds = 0.0
    for key, label in command_labels.items():
        record_path = main / acceptance["commands"][key]
        command = read_json(record_path)
        if command["exit_code"] != 0:
            raise ValueError(f"Failed required command: {key}")
        start, end = (datetime.fromisoformat(command[x]) for x in ("start_utc", "end_utc"))
        if key in {"main", "reproduction"}:
            calculation_seconds += (end-start).total_seconds()
        timing_rows.append([label, f"{(end-start).total_seconds():.3f}", command["exit_code"], wiki(record_path)])

    output_rows = [[name, count] for name, count in sorted(reproduction["rows"].items())]
    base_rel = str(base.relative_to(ROOT))
    text = f'''---
type: Отчёт
---
# Сравнение таксономических профилей

Дата завершённой проверки: {date} (UTC). Этап 4 выполнен: рассчитаны описательные
метрики трёх существующих таксономий, таблицы проверены независимым пересчётом и
воспроизведены из чистой локальной копии коммита. Полный вердикт — **PASS/0**.
Это техническое принятие расчёта; биологическая правильность назначений не установлена.

## Методика и фиксированный набор

Предварительная методика: `e6dec7b`; уточнение фактических входов и родословных:
`790e2de`; уточнение окончательной схемы до нового расчёта: `{protocol}` (`e7ffbeb`).
Источник вычислений и тестов: `{source}`. История и даты прежних коммитов сохранены.
Правила: [[METHOD_COMPARISON.md]], основные §§0–7 и последующие датированные дополнения.

Исходные данные содержат **{total_asv} ASV** и **{len(metadata)} образцов**.
В почвенный анализ включены **{len(soil_samples)} образцов**: {', '.join(soil_samples)}.
Общий почвенный знаменатель — **{soil_asv} ASV** с положительной суммой частот
в почвах; ещё **{total_asv-soil_asv} ASV** встречаются только вне почвенного набора.
Таблица назначений сохраняет все {total_asv} ID. K_S63 исключён из почвенных профилей.

Сравниваются SILVA 138.1, GTDB r220 и Greengenes2 2024.09 при одном алгоритме
`classify-sklearn`, ранее выполненном с `confidence=0.7`. Использованы результаты
первого DADA2 20/200. DADA2 и классификация в этом этапе повторно не запускались.
Рарефакция, новые фильтры, удаление Kr5 и объединение ридов разных почв не применялись.

Нормализация техническая: префиксы рангов, разделители и пробелы. Значимые суффиксы
и родословные сохраняются; синонимы не объединяются. Ключ охватывает известных
предков до целевого ранга. Неназначенная масса остаётся `UNCLASSIFIED`.
Основные ранги: тип, класс, порядок, семейство, род; вид не анализируется.

## Покрытие классификацией

Доля ASV использует общий почвенный знаменатель. Доля ридов рассчитывается для
каждого образца отдельно, относительно всех его ридов после DADA2. В последнем
столбце приведены медиана [минимум; максимум] по {len(soil_samples)} почвам.

{table(['Модель', 'Ранг', 'Назначено ASV / всего', 'Доля ASV', 'Доля ридов: медиана [мин; макс]'], coverage_rows)}

Источники: {wiki(main / 'coverage.tsv')} и {wiki(main / 'coverage_per_sample.tsv')}.

## Расстояния профилей

`D = 0.5 × Σ |p − q|` по объединению ключей, включая `UNCLASSIFIED`; единицы —
доля, диапазон [0,1]. Это расстояние в пространстве выбранных названий и
родословных, а не оценка биологической ошибки модели. Покрытие приведено выше;
общая неназначенная масса не считается подтверждением биологического согласия.

{table(['Пара', 'Ранг', 'Образцов', 'Медиана D', 'Минимум', 'Максимум'], distance_rows)}

На ранге рода нулевое расстояние при нулевом покрытии обеих моделей наблюдается
для следующих пар и образцов: {'; '.join(unknown_zero_groups) or 'не наблюдается'}.
В этих случаях обе массы полностью `UNCLASSIFIED`; ноль не означает совпадения
известных биологических таксонов.

Источники: {wiki(main / 'distances.tsv')}, {wiki(main / 'summary.tsv')}.
Профили по каждому образцу и рангу: {wiki(main / 'profiles.tsv')}.
Частоты целочисленные, общая масса каждого образца сохраняется.

## Условное совпадение назначений рода

Условная доля относится только к назначениям, известным обеим моделям.
Совместное отсутствие назначения не считается совпадением. Доли ASV и доли
ридов имеют разные знаменатели; при нулевом знаменателе записывается `NA`.
Сводка по ридам вычислена только по определённым значениям, их число указано явно.

{table(['Пара', 'Совпавших / совместно назначенных ASV', 'Условная доля ASV', 'Доля ридов: медиана [мин; макс]', 'Определено образцов / всего'], conditional_rows)}

Полные пять категорий (`both_match`, `both_diff`, `only_a`, `only_b`, `neither`)
и знаменатели для всех рангов сохранены в {wiki(main / 'pairwise_match.tsv')} и
{wiki(main / 'pairwise_match_per_sample.tsv')}. Исходные строки, выделенные ранги,
ключи и статусы назначения: {wiki(main / 'assignments.tsv')}.

## Проверки, ошибки и журналы

Синтетические испытания: **{tests['tests_run']} тестов**, **{len(tests['commands'])} запусков CLI**,
ошибок {tests['errors']}, провалов {tests['failures']}, пропусков {tests.get('skipped', 0)}.
Все фактические коды совпали с заранее указанными ожидаемыми. Генератор и ручные
ответы: [[tests/fixtures/taxonomy-comparison/README.md]]. Проверяются различие
0.5 по ASV и 0.9 по ридам, отдельный нулевой знаменатель ридов, пять категорий,
родословные, сохранение массы, перестановки, переименование ID и масштабирование.
Отрицательные случаи включают пропуски, дубликаты, неверные ID/роли/UUID,
отрицательные, дробные и неконечные частоты, подмену чисел и версии кода.

В сохранённой демонстрации изменён счётчик `both_match_reads` выходной таблицы,
затем исходное значение восстановлено: {failure_links}.
SHA выходного файла при мутации пересчитан: ошибку обнаруживает проверка значений,
а не одно несовпадение хеша. Дополнительно испытан перенос одного рида между
ключами профиля с сохранением общей массы. Все мутации выполнялись во временных копиях.

Полные команды, stdout/stderr, времена и коды сохранены в
{wiki(base / 'tests_final/summary.json')} и {wiki(base / 'execution_logs')}.
Исторические попытки не заменены новыми: {wiki(base / 'tests_development/summary.json')}
фиксирует изменение исходного кода во время первого испытания, последующий полный
прогон {wiki(base / 'tests_development2/summary.json')} прошёл; дополнительные
случаи проверены в {wiki(base / 'tests_added/summary.json')}.
Воспроизведение ошибок прежнего снимка: {wiki(base / 'historical_replay')}.

{table(['Команда', 'Длительность, с', 'Код', 'Запись команды'], timing_rows)}

Длительности относятся к текущим процессам сравнения и проверок; это не время
первичного DADA2 или классификации. Сумма времени основного и воспроизводящего
расчётов: **{calculation_seconds:.3f} с**, без тестов и проверок.
Расход токенов недоступен и не оценивался.

## Воспроизведение и сохранность

Второй расчёт выполнен из чистого локального клона коммита `{checkout['commit']}`.
Частные входы скопированы только для чтения. До и после сверены SHA256, UUID и
набор входов; состояние кода до и после чистое. Это локальное воспроизведение
тем же исполнителем, не проверка другим человеком. Реальные FASTQ/QZA в Git не включаются.

Все следующие таблицы совпали по ключам и значениям с абсолютным допуском
{reproduction['tolerance']:g} для долей; идентификаторы и счётчики сравниваются точно.

{table(['Таблица', 'Строк в каждом расчёте'], output_rows)}

Свидетельства: {wiki(base / 'clean_checkout.json')},
{wiki(base / 'reproduction_comparison.json')}, {wiki(main / 'acceptance.json')}.
Защищённых файлов проверено **{protected['checked']}**; изменений SHA256, размеров
или прав не найдено. Исходные данные, классификаторы, старые результаты и проверки
сохранены: {wiki(base / 'protected_after.json')}.

Проверка существующего результата из корня проекта:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_taxonomy_comparison.py --out-dir {base_rel}/main
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/report_taxonomy_comparison.py --check
```

Первой команде нужны локальные частные входы. Публичные синтетические испытания
работают без реальных данных:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_taxonomy_comparison.py --log-dir /tmp/taxonomy-tests-new
```

Каталог журналов должен быть новым. `--tables-only` проверяет только таблицы,
а полный режим дополнительно требует тесты, чистое воспроизведение и сохранность.
Общая проверка `checks/check.py` также завершилась кодом 0, но имеет свой договор.

## Ограничения и границы вывода

Kr5_S84 сохранён: **{retained} / {incoming} ридов = {retention:.2f}%** после DADA2
({wiki(denoise_path)}). Низкая сохранность не исправлена; причины потерь и пересечения
отрицательного контроля с почвенными ASV по согласованному объёму не исследовались.
Полный верификатор DADA2 сохраняет своё `INCOMPLETE/2` из-за порога 30%; это
не отменяется техническим PASS текущего сравнения.

Сравниваются три обученные модели/базы при одном алгоритме. Эталонной таксономии
и положительного контроля нет; лучшая модель не выбирается. Сопоставление
синонимов между базами не подтверждено, поэтому различия могут быть номенклатурными.
Праймеры, регион обучения и полнота референсного гена документально не установлены.
Результаты относятся к данным десяти почв при фиксированной обработке и не
обобщаются на все почвы Москвы. Статистическая гипотеза, p-value, уровень
значимости и порог биологического согласия на этом описательном этапе не задавались.
Экологическая интерпретация и исследовательский вывод этапа 5 здесь не выполнялись.

## Авторство

Первичные DADA2 и классификации, а также первоначальные наброски сравнения
подготовлены предыдущими исполнителями проекта. Доработка кода и проверки,
текущие расчёты сравнения, локальное воспроизведение и этот отчёт выполнены
ИИ-ассистентом по поручению пользователя. Параметры обработки и ограничения
задавал пользователь. Проверка не представляется независимой работой другого человека.

## Связи
- тип: [[Отчёт]]
- основанНа: [[Результаты сравнения таксономических профилей]]
- основанНа: [[Доказательства воспроизведения сравнения]]
- выполняет: [[Оформить проверенные результаты сравнения]]
- выполняет: [[Оформить отчёт и завершить этап сравнения]]
- выполняет: [[Завершить исправления и этап 4 сравнения таксономий]]
'''
    return text


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    try:
        expected = render(args.base_dir)
        if args.check:
            if not args.output.is_file() or args.output.read_text(encoding="utf-8") != expected:
                raise ValueError("Report differs from deterministic rendering of saved evidence")
            print("PASS/0 — report numbers and text match saved comparison evidence")
        else:
            with args.output.open("x", encoding="utf-8") as stream:
                stream.write(expected)
            print(f"Created: {args.output}")
        return 0
    except (OSError, ValueError, KeyError, StopIteration) as exc:
        print(f"FAIL/1 — {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

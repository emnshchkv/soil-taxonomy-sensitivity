# Выполнение сравнения 16 ноября 2026 года

Этап 4 завершён. Версия кода: `b813ac47ffee43a8f9eb40d1b317b96e9e31cfe3`. Предварительная методика `e6dec7b`, дополнения `790e2de` и `e7ffbeb` записаны до основного расчёта. Итоговая метка `taxonomy-comparison-v1` включает отчёт и свидетельства; метки созданы локально, публикация требует отдельного поручения.

- `main/`: 3374 назначения, показатели по 3309 почвенным ASV и десяти почвам; полный вердикт PASS/0.
- `reproduction/`: чистое локальное воспроизведение того же коммита. Все девять TSV совпали; детали в `reproduction_comparison.json`.
- `tests_final/`: 89 тестов, 60 вызовов CLI, PASS/0. Все входы синтетические, включая случаи FAIL → исправление → PASS.
- `historical_replay/`: прежний код из `before/` повторён в изоляции; сохранены настоящие коды и traceback прежних дефектов.
- `execution_logs/`: полные команды, времена, stdout/stderr и exit-коды основного цикла.
- `protected_before.json`, `protected_start.json`, `protected_after.json`: 85 защищённых файлов не изменились.
- `inputs_before.json`, `inputs_after.json`, `clean_checkout.json`: хеши, UUID и чистота копии; частные входы копии закрыты от записи.
- `tests_development/`, `tests_development2/`, `tests_added/`, `replay_development/`: промежуточные попытки, включая обнаруженное изменение исходника при первом прогоне. Они не подменяют финальные журналы.

## Проверить сохранённый результат

Из корня проекта с исходными частными входами:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_taxonomy_comparison.py --out-dir results/taxonomy_comparison/finish_20261116/main
```

`--tables-only` ограничивает проверку таблицами, не завершает этап. Общий `checks/check.py` относится к подготовке/классификации и имеет отдельный PASS. Старый критерий 30% сохранности DADA2 не отменён.

## Синтетическое воспроизведение без реальных входов

В чистой копии репозитория на метке `taxonomy-comparison-v1` с существующим окружением:

```sh
git checkout taxonomy-comparison-v1
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_taxonomy_comparison.py --log-dir /tmp/taxonomy-tests-new
```

Каталог журналов должен отсутствовать. Это доступная внешнему проверяющему синтетическая проверка; в ней не нужны FASTQ, частные QZA или классификаторы.

## Повторить именно записанный реальный расчёт

Для полного повторения есть частные входы из `main/manifest.tsv`. Создать чистую локальную копию репозитория, перейти на записанный в `main/run.json` коммит `b813ac47ffee43a8f9eb40d1b317b96e9e31cfe3`, скопировать эти семь входов с проверкой SHA256, закрыть копии от записи и выполнить расчёт в новый каталог. Точные уже исполненные команды и пути сохранены в `execution_logs/clone.command.json`, `checkout.command.json`, `reproduction.command.json` и `tables_reproduction.command.json`. Источником была локальная копия, не другой человек.

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/taxonomy_comparison.py --input-dir /path/to/private-input-copy --out-dir /tmp/taxonomy-reproduction-new
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_taxonomy_comparison.py --input-dir /path/to/private-input-copy --out-dir /tmp/taxonomy-reproduction-new --tables-only
```

Таблицы сравниваются по ключам и значениям (допуск 1e-12), а путь и время запусков не являются биологическими различиями. Новая версия Git может иметь тот же код, но другую запись `code_version`; для точного воспроизведения всех девяти TSV используется указанный исходный коммит.

Форматы и знаменатели — `main/README.md`, `main/normalization_rules.txt` и METHOD_COMPARISON.md §8.5. Итоговый отчёт — `reports/Сравнение таксономических профилей.md`.

# Проверка существующих результатов классификации — 2026-11-09

Научные расчёты DADA2 и classify-sklearn не запускались. Проверялись существующие QZA. Исходный код договора и регрессий зафиксирован в d73deae.

## Результаты

| Проверка | Результат | Доказательство |
|---|---|---|
| Готовность, статистика 11 образцов, пять последовательностей и три пробы | PASS/0 | initial/readiness.command.json и stdout/stderr |
| Регрессии договора | PASS/0, 29 ожидаемых исходов | regressions/summary.json |
| Две функции старого кода на дубликатах | старый ложный PASS/0; новый FAIL/1 | regressions/legacy_duplicate_*.stdout.txt и duplicate_stats/duplicate_probe.* |
| Первый штатный validate в ограниченной среде | FAIL/1 до чтения артефакта | full_verification/ и full_verification.command.json |
| Повтор штатной проверки с доступом к сведениям CPU | PASS/0; шесть команд validate/export с кодом 0 | full_verification_retry/summary.json и *.command.json |
| Общая проверка проекта | PASS/0 | project_check.command.json и stdout/stderr |
| Сохранность исходных файлов | PASS/0; 83 файла, SHA256/размеры/права совпадают | protected_before.json, protected_after.json, preservation.* |

Каждая запись *.command.json содержит точный argv, рабочий каталог, время UTC и настоящий exit_code; stdout и stderr записаны отдельно целиком, включая предупреждения. Пустой stderr-файл означает отсутствие stderr, а не потерю журнала. initial/ и regressions/ уже зафиксированы в d73deae.

Первая попытка full_verification не дошла до проверки QZA: импорт QIIME 2 вызвал ошибку psutil.cpu_count() == None в ограниченной среде. Повтор выполнен тем же кодом с доступом к системной информации, в новый каталог; окружение и входы не изменялись. Журнал неудачи сохранён.

## Воспроизведение

Нужны частные локальные QZA и существующее окружение qiime2-amplicon-2025.7. Из корня проекта:

~~~sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/verify_classifier_readiness.py
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_classification_contract.py --log-dir /tmp/classification-regressions-new
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/verify_full_classification.py
conda run --no-capture-output -n qiime2-amplicon-2025.7 /bin/sh -c 'exec "$CONDA_PREFIX/bin/python3" -B checks/check.py'
~~~

Каталог регрессий должен быть новым. Полная проверка по умолчанию создаёт новый каталог results/full_classification/verification_*; выводит его путь. Экспорты выполняются во временный каталог и сверяются с архивом, рабочие экспорты не перезаписываются. Скрипты verify_* сами выбирают Python conda-окружения; checks/check.py остаётся неизменным и вызывается явно через интерпретатор окружения.

Регрессии создают временные копии входов, изменяют только их и сравнивают SHA256 оригиналов до/после. Два legacy-случая исполняют отдельные функции из 8f9832f; для них нужна история Git. Они показывают старый дефект, а не запуск всего старого верификатора.

## Происхождение и ограничения

- model_provenance.json содержит адресные выдержки provenance трёх моделей, версии и SHA256; генератор src/classifier_provenance.py, команда в model_provenance.command.json.
- environment.* фиксирует Python и версии установленных пакетов.
- full_verification_retry/summary.json содержит фактические UUID, SHA256, параметры, времена и UUID исполнения из provenance, а также три проверенных ASV для каждой модели.
- Исходные журналы *_run.log в родительском каталоге относятся к классификации 2026-11-08. Их команды и exit-коды не были сохранены. Проверочные exit-коды 2026-11-09 нельзя выдавать за коды тех запусков.
- Прежние case*.log, verify_*log и results/classifier_readiness/verify_full_log.txt — исторические попытки, не полное доказательство новой реализации. Они сохранены без переписывания.
- Эта локальная проверка требует закрытых данных и не является независимой проверкой другим человеком. Учебный синтетический комплект DADA2 и метка hw2 сохранены отдельно.
- PASS означает выполнение технического договора; он не подтверждает биологическую правильность таксономии или соответствие региона обучения ампликону.

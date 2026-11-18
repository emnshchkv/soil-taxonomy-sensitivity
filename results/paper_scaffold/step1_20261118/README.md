# Проверки адаптации каркаса и глоссария

Шаг 1, 2026-11-18. Базовый коммит до работы — 4b4c216; согласование структуры — d51e34a. Это свидетельства подготовки каркаса, не результаты нового научного эксперимента и не готовность статьи к подаче.

| Свидетельство | Содержание |
|---|---|
| [baseline.json](baseline.json) | SHA256 1198 существовавших файлов: отслеживаемые научные файлы, checks/, ontology.yaml, заметки классов, .obsidian и два локальных шаблона. |
| [preservation.json](preservation.json) | Все файлы снимка сохранены, manuscript не содержит файлов; PASS. |
| [scaffold_check.json](scaffold_check.json), [stdout](scaffold_check.stdout.txt), [stderr](scaffold_check.stderr.txt) | Проверка полного каркаса и глоссария, код 0. |
| [project_check.json](project_check.json), [stdout](project_check.stdout.txt), [stderr](project_check.stderr.txt) | Прежняя общая проверка проекта, код 0; включает неотслеживаемые научные входы. |
| [tests_check.json](tests_check.json) | Сводка фактического запуска семи изолированных синтетических тестов подагентом; код 0. Полный stdout этого запуска отдельно не сохранялся. |
| [final_documents.json](final_documents.json) | Итоговая структурная проверка, ссылки новых документов, закрытие задач, сохранность файлов и прежних записей журнала. |

Воспроизведение из корня, без научных расчётов:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_paper_scaffold.py
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_paper_scaffold.py
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/check.py
```

Смысл определений сверялся чтением кода и таблиц; формальный checker проверяет структуру и ссылки. При расширении постановки и литературы после нового поручения условия структурной проверки шага 1 следует пересмотреть: сейчас она намеренно требует `state: deferred` для этих заготовок. Настоящий комплект статьи и его будущая приёмка остаются за пределами этой задачи.

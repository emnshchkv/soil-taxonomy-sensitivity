# Шаг 4 подготовки статьи — проверенные результаты

2026-11-18. Комплект подготовлен по practice.md с адаптацией к существующим src/ и results/. Научные обработки не повторялись; рукопись не создавалась и не восстанавливалась (`absent_not_checked`).

## Что использовать дальше

- [Выбор материалов](../../research/results_selection.md): две компактные таблицы для основного текста, четыре существующих PDF для приложения, допустимые наблюдения с источниками и ограничениями.
- [Контракт](../../research/paper_results_contract.md): зафиксирован до генерации; разделяет исходное сравнение всех ридов после DADA2 и модельный отбор Bacteria + Archaea.
- [Описание числового комплекта](20261118/README.md), [карта 37 чисел](20261118/numbers.tsv), [макросы LaTeX](20261118/numbers.tex).
- [Покрытие](20261118/comparison_coverage.tsv), [расстояния и условное совпадение](20261118/comparison_pairs.tsv), [семейства по образцам](20261118/family_sample_counts.tsv): 15, 15 и 33 строки. Контроль явно помечен и не входит в почвенные сводки.
- [Основная сравнительная таблица TeX](20261118/comparison_main.tex), [модельные N и доли неопределённых семейств TeX](20261118/family_main.tex), [реестр четырёх PDF](20261118/figure_sources.tsv).

Это готовые данные и фрагменты для будущего текста. Финальная русская подпись, вёрстка в шаблоне, все числа и цитаты фактической рукописи проверяются на шаге 5. Компактные фрагменты проверены по содержимому; локальная сборка LaTeX не выполнялась.

## Свидетельства

| Проверка | Результат и файл |
|---|---|
| Выбор наблюдений и источников | PASS/0: пять групп утверждений, десять ссылок; [selection_validation.json](preflight_20261118/selection_validation.json). |
| Два запуска генератора | PASS: все десять файлов, включая manifest.json, совпали побайтно; [reproduction.json](validation_20261118/reproduction.json). |
| Независимая сверка нового комплекта | PASS/0, 1520 сверок; [numeric_verification.json](validation_20261118/numeric_verification.json). |
| Изолированные отрицательные испытания | PASS/0, 16 CLI-кейсов: исходный PASS, 14 ожидаемых FAIL, восстановленный PASS; [numeric_tests.json](validation_20261118/numeric_tests.json). Временные выходные хеши обновлялись: обнаружены содержательные ошибки, а не только несовпадение SHA. |
| Исходное сравнение, полный режим | PASS/0; [comparison.txt](validation_20261118/scientific/comparison.txt). |
| Окончательные семейства и их self-test | PASS/0; [family.txt](validation_20261118/scientific/family.txt), [family_self_test.txt](validation_20261118/scientific/family_self_test.txt). |
| Точное совпадение трёх DADA2 | PASS/0 при ранее согласованном исключении Kr5; [exact_match.txt](validation_20261118/scientific/exact_match.txt). |
| Полный исходный критерий DADA2 | **INCOMPLETE/2**: Kr5_S84, 1505/60078 ридов, около 2.51% <30% во всех трёх повторах. Четыре метрики каждой пары нулевые; других причин INCOMPLETE нет. [full_denoising.txt](validation_20261118/scientific/full_denoising.txt). |
| Каркас, литература, общая checks/check.py | Все exit 0; [scaffold.txt](validation_20261118/scientific/scaffold.txt), [literature.txt](validation_20261118/scientific/literature.txt), [general.txt](validation_20261118/scientific/general.txt). |
| Сохранность прежних файлов | 1209 файлов без изменений: 1167 ранее отслеживаемых и 42 дополнительных входа; [verification.json](validation_20261118/scientific/verification.json). |

Подмены охватывают значение макроса, модель/источник, два вида знаменателей, NA→0, перенос одного рида при прежнем N, число/роль/состав почв, модельные заголовки TeX, ложный статус рукописи и пропуск строки. Испытания меняют только временную копию; 33 файла числового комплекта и его источников/генератора также проверены до и после тестов.

PASS проверки переноса относится к соответствию сохранённым таблицам. Он не устанавливает истинную таксономию, превосходство модели или экологические причины. INCOMPLETE полного денойзинга остаётся отдельным ограничением; критерии не ослаблялись.

## Воспроизведение

Команды из корня проекта, Python окружения qiime2-amplicon-2025.7:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/make_paper_results.py --out-dir results/paper_step4/NEW
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_paper_results.py --out-dir results/paper_step4/20261118
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_paper_results.py --out-dir results/paper_step4/20261118
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/run_paper_step4_scientific_checks.py --log-dir results/paper_step4/validation_NEW
```

Каталоги NEW и validation_NEW должны отсутствовать. Последняя команда запускает восемь существующих проверок, сохраняет фактические коды и контрольные суммы, отдельно проверяет, что причины INCOMPLETE ограничены известным Kr5. Её exit 0 означает совпадение с указанными исходами и сохранность файлов; полный DADA2 при этом остаётся INCOMPLETE/2.

Числовой checker не импортирует генератор: ожидаемые селекторы и арифметика заданы независимо. Генератор использует рациональную арифметику, checker — Decimal; проверяются реальные значения/подписи TeX, состав строк и SHA256. Источники и код генератора — [manifest.json](20261118/manifest.json); команды новых проверок сохранены рядом с JSON-результатами.

---
state: fixed_for_step4
date: 2026-11-18
---
# Контракт числового комплекта шага 4

Зафиксирован до генерации нового комплекта. Это договор о переносе существующих результатов, а не изменение METHOD_COMPARISON.md и не предварительная регистрация уже просмотренных результатов. Выбор материалов — [[research/results_selection.md]]. Рукопись отсутствует и получает статус `absent_not_checked`.

## Область и источники

1. `comparison_all_post_dada2`: `results/taxonomy_comparison/finish_20261116/main/`. Пять рангов в порядке тип, класс, порядок, семейство, род; модели `sklearn_silva138`, `sklearn_gtdb`, `sklearn_gg2`; пары SILVA–GTDB, SILVA–GG2, GTDB–GG2. Все десять почвенных образцов; контроль исключён. Знаменатель покрытия ридов и профилей — все риды образца после DADA2. Покрытие ASV — все ASV, присутствующие в почвах. Сохранены исходная буквальная нормализация родословной и категория UNCLASSIFIED, включая прежнюю обработку uncultured.
2. `family_prokaryote_selected`: `results/taxonomy_barplots/prokaryote_family_20261118_final/`. Отбор Bacteria + Archaea отдельно каждой моделью; исключены явно помеченные органеллы, Eukaryota и неизвестный домен. Конечное точное f__uncultured без учёта регистра означает UNCLASSIFIED; ASV при этом не удаляется. Знаменатель семейственных долей — selected_reads данной модели и образца. При N=0 доля NA. Порог 2% из манифеста относится только к отображению семейств, не к удалению ридов.
3. Метаданные — `data/meta/metadata.tsv`; параметры обрезки — `memory/denoise_choice.md`; статистика первого сохранённого исполнения — `results/exact_match_reproduction/provenance/samples.tsv`, run=1. Количество всех/почвенных ASV — assignments.tsv; количество моделей — уникальные модели coverage.tsv. Порог отображения — parameters.min_percent семейственного манифеста.

Контроль K_S63 сохраняется в подробной семейственной таблице как negative_control. В почвенные сводки и основные таблицы он не входит. Данные/код прежних научных расчётов не изменяются; новые выходы — только `results/paper_step4/20261118/`. Снимок защищённых файлов до работы — [[results/paper_step4/preflight_20261118/baseline.json]].

## Выходы и расчёты

Генератор `src/make_paper_results.py --out-dir <новый каталог>` отказывает при существующем каталоге. Вспомогательные функции допустимы только в новых файлах src/. Повторный запуск идёт в отдельный временный каталог и даёт побайтно те же выходы; абсолютные рабочие пути и время генерации в детерминированные файлы не включаются.

| Файл | Ключ, строки | Обязательное содержание |
|---|---|---|
| comparison_coverage.tsv | model, rank; 15 строк | scope; assigned_asv, total_asv, assigned_asv_pct; read_coverage_min_pct, read_coverage_median_pct, read_coverage_max_pct; n_samples=10. Доли пересчитываются из целых assigned/total, сводка по отдельным почвам. |
| comparison_pairs.tsv | model_a, model_b, rank; 15 строк | scope; distance_min, distance_median, distance_max, n_samples; both_match, both_diff, only_a, only_b, neither, total_asv, both_assigned_asv; conditional_asv_pct; conditional_read_min_pct, conditional_read_median_pct, conditional_read_max_pct, conditional_defined_samples. Условный знаменатель — both_assigned, нулевой знаменатель даёт NA и исключается только из условной сводки. |
| family_sample_counts.tsv | model, sample; 33 строки | scope, sample_type; post_dada2_reads, selected_reads, excluded_reads, selected_bacteria_reads, selected_archaea_reads; named_family_reads, unclassified_family_reads; selected_pct_of_post_dada2, unclassified_pct_of_selected, zero_selected. Полная сумма family_profiles.count равна N; UNCLASSIFIED выбирается по полному ключу. |
| figure_sources.tsv | figure_id; 4 строки | scope, model (или all), path, sha256, source_table, denominator, display_threshold_pct. Регистрируются PDF трёх моделей и обзора в figures/taxonomy/family; исходники не копируются и не перерисовываются. |
| numbers.tsv | id; 37 строк | id, scope, source_file, row_key, column, operation, numerator, denominator, n, unit, raw_value, digits, display_value. Составной источник/ключ допустим как однозначная JSON-запись. Значения выводятся из источников, не вводятся вручную. |
| numbers.tex | имена из numbers.tsv | 37 команд newcommand, только ASCII-буквы в имени; фактические значения совпадают с display_value. |
| comparison_main.tex | 5 строк рангов | Медиана покрытия ридов (%) каждой модели и медиана D каждой пары; явный порядок заголовков SILVA, GTDB, GG2; исходная область сравнения. Диапазоны доступны в TSV. |
| family_main.tex | 10 строк почв | sample, всего после DADA2, затем N и U% каждой модели; U% — UNCLASSIFIED/N; N=0 → NA. |
| manifest.json | один объект | версия схемы, manuscript_status=absent_not_checked, относительные пути и SHA256 всех использованных источников/кода генератора и выходов; никаких утверждений о проверке рукописи. |

Показатели и селекторы карты должны быть достаточны, чтобы независимо восстановить число. Для всех таблиц источник/ключ/правило также фиксируются этим контрактом и в README комплекта. Регистр рисунков дополняет связь с исходным family_barplot_data.tsv, масками, манифестом и сохранённой проверкой семейств.

37 макросов:

- ProjectSamplesAll, ProjectSamplesSoil, ProjectSamplesControl, ProjectAsvAll, ProjectAsvSoil, ProjectTrimLeft, ProjectTruncLen, ProjectAsvLength, ProjectKrFiveInput, ProjectKrFiveRetained, ProjectKrFiveRetentionPct, ProjectModels, SelectedDisplayThresholdPct.
- Для Silva, Gtdb, GgTwo: Cmp<модель>FamilyAsvPct, Cmp<модель>FamilyReadMedianPct, Selected<модель>ZeroSoilSamples.
- Для SilvaGtdb, SilvaGgTwo, GtdbGgTwo: Cmp<пара>FamilyDistanceMedian, Cmp<пара>GenusDistanceMedian, Cmp<пара>GenusConditionalAsvPct, Cmp<пара>GenusConditionalReadMedianPct, Cmp<пара>GenusDefinedSamples.

В числовых TSV сохраняется точность вычисления, минимум 12 значащих цифр для нецелых. Для отображения используются Decimal и ROUND_HALF_EVEN: проценты — два знака, D — шесть, целые — без десятичной части. Проценты сначала вычисляются из целых числителей/знаменателей; медиана для чётного числа — среднее двух центральных значений. NA остаётся строкой NA.

## Независимая проверка

`src/check_paper_results.py --out-dir <каталог>` не импортирует генератор и не доверяет селекторам/операциям из проверяемой карты: ожидаемые источники, состав строк и расчёты заданы независимо. Проверяются значения всех TSV, происхождение 37 макросов, реальные значения и подписи обоих TeX-фрагментов, исходные/выходные хеши, роли образцов, уникальность ключей и область знаменателей. Коды: PASS/0, FAIL/1; отсутствующие обязательные файлы также FAIL, поскольку комплект заявлен как завершённый.

Испытания выполняются в изолированных копиях. Верный комплект принимается; подмена макроса, модели/источника, знаменателя, NA на 0, перенос одного рида, изменение состава почв и перестановка подписей должны обнаруживаться даже после обновления временных выходных хешей. Восстановленный комплект снова проходит. Сохраняются команды, результаты и проверка повторной генерации. Источники в проекте не изменяются.

Дополнительно запускаются существующие полные проверки сравнения, семейственных профилей и их self-test, check_exact_match.py, verify_denoising.py, checks/check.py, проверки каркаса и литературы. Ожидаемый INCOMPLETE/2 полного денойзинга из-за Kr5 (1505/60078 при 20/200) записывается отдельно, без ослабления исходного критерия. Успех переноса чисел не является доказательством биологической правильности назначений.

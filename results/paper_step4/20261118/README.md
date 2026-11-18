# Числовой комплект шага 4

Дата договора: 2026-11-18. Контракт: research/paper_results_contract.md.
Рукопись: absent_not_checked. Здесь перенесены сохранённые результаты; DADA2,
классификация, научные сравнения и семейственные рисунки не выполнялись заново.

comparison_coverage.tsv: 15 строк, ключ model/rank. Целые assigned_asv/total_asv
из comparison main/coverage.tsv; read_coverage — проценты assigned_reads/total_reads
из coverage_per_sample.tsv отдельно по десяти почвам, затем min/median/max.
comparison_pairs.tsv: 15 строк, ключ model_a/model_b/rank. Пять категорий и общий
ASV-знаменатель из pairwise_match.tsv. Условный процент — both_match/both_assigned.
По ридам источник pairwise_match_per_sample.tsv; нулевой both_assigned_reads
исключается только из условной сводки; число определённых почв показано отдельно.
Расстояния — сохранённый distance из distances.tsv, сводка по десяти почвам.
Область comparison_all_post_dada2 использует все риды после DADA2 и прежнюю
буквальную обработку таксонов; uncultured не переопределяется.

family_sample_counts.tsv: 33 строки, ключ model/sample; sample_type берётся из
metadata.tsv и сверяется с исходниками. Источник filter_counts.tsv, а UNCLASSIFIED
выбирается по точному key из family_profiles.tsv. post_dada2_reads = total_reads
filter_counts; selected_reads = N. named_family_reads = N − UNCLASSIFIED.
selected_pct_of_post_dada2 = 100*N/post_dada2_reads; unclassified_pct_of_selected
= 100*UNCLASSIFIED/N. При N=0 последняя доля NA, не ноль. Контроль K_S63 остаётся
в подробной таблице, но не входит в почвенные сводки или основной TeX-фрагмент.
Область family_prokaryote_selected — модельный отбор Bacteria + Archaea с
исключением органелл/Eukaryota/неизвестного домена; точное f__uncultured здесь
считается неопределённым семейством. Это назначения моделей, не биологическая истина.

figure_sources.tsv регистрирует четыре существующих PDF figures/taxonomy/family;
файлы не копируются. Source_table — family_barplot_data.tsv. Знаменатель каждой
панели — selected_reads(model,sample), порог отображения — parameters.min_percent
prokaryotic_family_manifest.json. В manifest.json также зарегистрированы маска,
полные профили, сохранённая проверка и манифест переноса рисунков.

numbers.tsv содержит 37 скаляров и их полное происхождение. row_key — JSON
селектор: строка означает точное равенство, массив sample — заданный набор строк
в указанном порядке. column — поле либо JSON-массив полей числителя/знаменателя.
operation: count_rows, count_distinct, value, subtract, ratio_pct,
median_ratio_pct или median_decimal. numerator/denominator у ratio_pct — целые
значения; у median_ratio_pct — массивы целых в порядке sample; у subtract — два
операнда. n — число исходных выбранных строк (для subtract — два параметра).
raw_value: точная рациональная арифметика с десятичной сериализацией precision=60;
для печати Decimal ROUND_HALF_EVEN: проценты 2 знака, D 6, целые 0. NA сохраняется.
Медиана чётного числа значений — среднее двух центральных.

numbers.tex определяет 37 ASCII-команд newcommand. comparison_main.tex:
пять рангов, медианы покрытия и D; family_main.tex: десять почв, объём до отбора,
N и процент UNCLASSIFIED для каждой модели. Это простые tabular-фрагменты,
не рукопись и не проверка окончательной вёрстки. Исходные русские ранги TSV
соответствуют Phylum/Class/Order/Family/Genus в TeX.

manifest.json фиксирует SHA256 всех использованных источников, двух новых
файлов генератора и девяти выходов. Сам manifest не включает собственный хеш.
Пути относительны корню проекта или каталогу комплекта; время и абсолютные
рабочие пути не сохраняются, поэтому повторный запуск побайтно воспроизводим.

Запуск из корня проекта, выходной каталог должен отсутствовать:
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/make_paper_results.py --out-dir results/paper_step4/NEW

Генерация сама не является независимой проверкой. Приёмка выполняется отдельно
src/check_paper_results.py; отсутствие рукописи остаётся absent_not_checked.

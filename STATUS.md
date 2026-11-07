# Состояние проекта

Выполнена [[Завершить сводку денойзинга и проверить готовность классификаторов]] (AGENT_CLASSIFIER_READINESS_PROMPT.txt).

## Сводка денойзинга
- Создана недостающая `work/qzv/dada2_stats.qzv` штатной командой `qiime metadata tabulate` (SampleData[DADA2Stats] → Metadata) и проверена через `qiime tools validate` — валидна (level=max).
- Статистика первого запуска (run=1) из `results/exact_match_reproduction/provenance/samples.tsv` сверена с QZA `dada2_stats.qza`: все 11 образцов, input, filtered, denoised, non-chimeric совпадают точно. Kr5_S84: 1505/60078 = 2.51% — согласовано.

## Проверка классификаторов
- Три файла — тип `TaxonomicClassifier`, SHA256 совпадают с `data/refs/VERSIONS.md`, обучены в QIIME 2 2024.5.0 (feature-classifier 2024.5.0) на scikit-learn 1.4.2 (= установленному).
- Источники из provenance: SILVA 138.1 NR99, GTDB r220 SSU reps, Greengenes 2024.09 backbone (полная длина). Регион обучения не установлен.
- Пробная классификация: 5 ASV (первые по лексикографической сортировке ID из rep_seq.qza), classify-sklearn confidence=0.7 n_jobs=1, все три exit 0, ровно 5 ASV без дубликатов, Confidence∈[0,1], UUID корректны.

## Результаты
- `results/classifier_readiness/` — команды, выводы, коды, выбранные 5 ID, таблица проверки моделей.
- `work/classifier_probe/` — артефакты пробной классификации (не в Git).
- `reports/Готовность к классификации.md` — отчёт.

## checks/check.py
Код 1, причина — отсутствуют три полные таксономии `work/qza/taxonomy_sklearn_*.qza` (полная классификация не выполнялась по инструкции). Новых нарушений нет.

## Дальше
Полная классификация тремя методами и сравнение методов — следующий этап, не выполнялся. DADA2 не запускался.
# Журнал команд, выводов и кодов завершения

Этап: завершение сводки денойзинга и проверка готовности классификаторов (AGENT_CLASSIFIER_READINESS_PROMPT.txt). Дата: 2026-11-07. Окружение: qiime2-amplicon-2025.7.

## 1. Сводка денойзинга

### 1.1 Создание dada2_stats.qzv
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime metadata tabulate --m-input-file work/qza/dada2_stats.qza --o-visualization work/qzv/dada2_stats.qzv
```
Вывод:
```
/Users/elizavetamenshchikova/miniconda3/envs/qiime2-amplicon-2025.7/lib/python3.10/site-packages/q2_demux/_summarize/_visualizer.py:11: UserWarning: pkg_resources is deprecated as an API. ...
Saved Visualization to: work/qzv/dada2_stats.qzv
```
Код завершения: 0

### 1.2 Валидация dada2_stats.qzv
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime tools validate work/qzv/dada2_stats.qzv
```
Вывод:
```
Result work/qzv/dada2_stats.qzv appears to be valid at level=max.
```
Код завершения: 0

### 1.3 Экспорт dada2_stats.qza для сверки
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime tools export --input-path work/qza/dada2_stats.qza --output-path work/qzv_export/dada2_stats
```
Вывод: `Exported work/qza/dada2_stats.qza as DADA2StatsDirFmt to directory work/qzv_export/dada2_stats`
Код завершения: 0

Сверка: все 11 образцов (input, filtered, denoised, non-chimeric) совпадают точно между samples.tsv (run=1) и QZA stats.tsv; сохранность совпадает (QZA 2 знака, samples.tsv 6). Kr5_S84: 1505/60078 = 2.51%.

## 2. Проверка классификаторов

### 2.1 SHA256
Команда: `shasum -a 256 data/refs/*.qza`
```
5c5900c10ad601b12cf0d2280d52cc0f7852627e89070595cc3899a6d690ee27  data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza
07aadcf7472d9cc6f853f6b4615348619f1a3eceb56c1fb1b6d8dbb20554765f  data/refs/gtdb_classifier_r220.qza
c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616  data/refs/silva-138-99-nb-classifier.qza
```
Все совпадают с VERSIONS.md.

### 2.2 UUID и тип (qiime tools peek)
```
silva-138-99-nb-classifier: UUID 70b4b5f4-8fce-40bd-b508-afacbc12a5ed, Type TaxonomicClassifier
gtdb_classifier_r220:      UUID 5d5461cc-6a51-434b-90ab-040f388e4221, Type TaxonomicClassifier
gg_2024.09...:             UUID 49ccfb0a-155d-404b-80b4-818d2aeb53b2, Type TaxonomicClassifier
```

### 2.3 Provenance (из action.yaml внутри QZA)
Все три: QIIME 2 2024.5.0, feature-classifier 2024.5.0, fit_classifier_naive_bayes, python 3.9.19, linux-x86_64, scikit-learn 1.4.2.
- silva: источник rescript get_silva_data SILVA 138.1 SSURef_NR99; filter_seqs_length_by_taxon min_lens 900/1200/1400 (полная длина).
- gtdb: импорт ar53_ssu_reps_r220.fna + bac120_ssu_reps_r220.fna (GTDB r220 SSU reps).
- gg: импорт 2024.9.backbone.fna (Greengenes 2024.09 backbone).

Установленный scikit-learn: 1.4.2 (совпадает).

## 3. Пробная классификация (5 ASV)

Выбор: первые пять по лексикографической сортировке Feature ID из rep_seq.qza (3374 ASV). Список — results/classifier_readiness/probe5_ids.md.

Подготовка:
```
qiime tools export --input-path work/qza/rep_seq.qza --output-path work/classifier_probe/rep_seq   (exit 0)
qiime tools import --input-path work/classifier_probe/probe5.fasta --output-path work/classifier_probe/probe5_seqs.qza --type 'FeatureData[Sequence]'   (exit 0)
```

### 3.1 silva138
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime feature-classifier classify-sklearn --i-classifier data/refs/silva-138-99-nb-classifier.qza --i-reads work/classifier_probe/probe5_seqs.qza --p-confidence 0.7 --p-n-jobs 1 --o-classification work/classifier_probe/taxonomy_probe_silva138.qza
```
Вывод: `Saved FeatureData[Taxonomy] to: work/classifier_probe/taxonomy_probe_silva138.qza`
Код завершения: 0

### 3.2 gtdb
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime feature-classifier classify-sklearn --i-classifier data/refs/gtdb_classifier_r220.qza --i-reads work/classifier_probe/probe5_seqs.qza --p-confidence 0.7 --p-n-jobs 1 --o-classification work/classifier_probe/taxonomy_probe_gtdb.qza
```
Вывод: `Saved FeatureData[Taxonomy] to: work/classifier_probe/taxonomy_probe_gtdb.qza`
Код завершения: 0

### 3.3 gg2
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 qiime feature-classifier classify-sklearn --i-classifier data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza --i-reads work/classifier_probe/probe5_seqs.qza --p-confidence 0.7 --p-n-jobs 1 --o-classification work/classifier_probe/taxonomy_probe_gg2.qza
```
Вывод: `Saved FeatureData[Taxonomy] to: work/classifier_probe/taxonomy_probe_gg2.qza`
Код завершения: 0

### 3.4 Проверка результатов
- Все три: ровно 5 ASV, без дубликатов, Confidence ∈ [0,1].
- Входной UUID reads f191f736-c066-42be-aa44-80e29282c3a8 (probe5_seqs.qza) — совпадает во всех трёх.
- UUID классификаторов в provenance: silva138=70b4b5f4..., gtdb=5d5461cc..., gg2=49ccfb0a... — совпадают с peek.
- Параметры: confidence=0.7, n_jobs=1 во всех трёх.

## 4. checks/check.py
Команда:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/check.py
```
Вывод (хвост):
```
ПРОВЕРКА НЕ ПРОЙДЕНА:
 - work/qza/taxonomy_sklearn_silva138.qza: файла нет
 - work/qza/taxonomy_sklearn_gtdb.qza: файла нет
 - work/qza/taxonomy_sklearn_gg2.qza: файла нет
```
Код завершения: 1
Причина: отсутствуют три полные таксономии (полная классификация не выполнялась по инструкции). Новых нарушений нет; разделы импорта, денойзинга и классификаторов (кроме полных таксономий) пройдены.
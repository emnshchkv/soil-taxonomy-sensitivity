> Актуальная проверка 2026-11-09: [repair_20261109/README.md](repair_20261109/README.md). Код договора — d73deae, новые полные журналы — 807bc3f. Текст ниже описывает предыдущую проверку 8f9832f; её краткий вывод не содержит stdout/stderr внутренних validate/export. Для актуального вердикта используйте новые журналы.

# Полная проверка трёх готовых таксономий

Дата: 2026-11-09. Окружение: qiime2-amplicon-2025.7, q2cli 2025.7.0.
Скрипт: `src/verify_full_classification.py` (воспроизводимый).
Полный вывод: `results/full_classification/verify_full_classification.log` (код 0).

## Команды

Проверка каждой модели выполнялась скриптом `src/verify_full_classification.py`,
который внутри вызывает:

```
qiime tools validate work/qza/taxonomy_sklearn_<метод>.qza
qiime tools export --input-path work/qza/taxonomy_sklearn_<метод>.qza --output-path work/qzv_export/full_classification/<метод>
```

Запуск скрипта:
```
conda run --no-capture-output -n qiime2-amplicon-2025.7 /bin/sh -c 'exec "$CONDA_PREFIX/bin/python3" -B src/verify_full_classification.py'
```
Код завершения: 0 (PASS).

## Итоги проверки по моделям

| Метод | validate | Уник. ID | Совпадение с rep_seq | Confidence∈[0,1] | Taxon непустой | action | reads UUID | classifier UUID | confidence | n_jobs | probe3 сверка |
|---|---|---|---|---|---|---|---|---|---|---|---|
| silva138 | PASS | 3374 | PASS | PASS | PASS | classify_sklearn | 978af6f8…5f3c | 70b4b5f4…a5ed | 0.7 | 1 | PASS |
| gtdb | PASS | 3374 | PASS | PASS | PASS | classify_sklearn | 978af6f8…5f3c | 5d5461cc…e4221 | 0.7 | 1 | PASS |
| gg2 | PASS | 3374 | PASS | PASS | PASS | classify_sklearn | 978af6f8…5f3c | 49ccfb0a…b53b2 | 0.7 | 1 | PASS |

- `reads` UUID = `978af6f8-91bb-4df2-9e19-b395076e5f3c` — UUID `work/qza/rep_seq.qza`.
- `classifier` UUID соответствует модели: silva138=`70b4b5f4-8fce-40bd-b508-afacbc12a5ed`,
  gtdb=`5d5461cc-6a51-434b-90ab-040f388e4221`, gg2=`49ccfb0a-155d-404b-80b4-818d2aeb53b2`.
- probe3 — первые три ASV по лексикографической сортировке ID:
  `0021cfa2357aba7fe9d1b92465676a8d`, `0027dadcdc187da02f8f422ae4690999`,
  `002d2685b94200ea0724b24286fec6ba`; их Taxon и Confidence сверены с `qiime tools export`.

## SHA256 и UUID файлов

| Файл | SHA256 | UUID |
|---|---|---|
| work/qza/rep_seq.qza | 2f81bcce47316cfac022c82a6b4f2e3206eb1a3fa817a81883aff2bc858cda86 | 978af6f8-91bb-4df2-9e19-b395076e5f3c |
| data/refs/silva-138-99-nb-classifier.qza | c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616 | 70b4b5f4-8fce-40bd-b508-afacbc12a5ed |
| data/refs/gtdb_classifier_r220.qza | 07aadcf7472d9cc6f853f6b4615348619f1a3eceb56c1fb1b6d8dbb20554765f | 5d5461cc-6a51-434b-90ab-040f388e4221 |
| data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza | 5c5900c10ad601b12cf0d2280d52cc0f7852627e89070595cc3899a6d690ee27 | 49ccfb0a-155d-404b-80b4-818d2aeb53b2 |
| work/qza/taxonomy_sklearn_silva138.qza | 66d7b18f8349c6adcd2e0d8a74f860e6674f5691e5bdbea8cd885cb0381505d6 | 1e17faa9-825e-4930-ba09-f02e38391ca4 |
| work/qza/taxonomy_sklearn_gtdb.qza | a542ba9b762e42112b7dc95a05de6b9870b4b8e22b01acd0847b22201090f73e | 64a8d4d3-2751-4a3d-a0b2-711c358793fc |
| work/qza/taxonomy_sklearn_gg2.qza | e972b452e45cfc991e69712ded06b3ef47d5915a85caecfbdc06016ab2439771 | fdc19c64-00f2-4fae-b8ac-387eadd6e74f |

## Времена из provenance (корневое действие classify_sklearn)

| Метод | start | end | duration |
|---|---|---|---|
| silva138 | 2026-11-08T05:57:37.780711+03:00 | 2026-11-08T06:15:27.353262+03:00 | 17 min 49 s |
| gtdb | 2026-11-08T06:36:32.691838+03:00 | 2026-11-08T07:26:28.642475+03:00 | 49 min 55 s |
| gg2 | 2026-11-08T08:16:41.890256+03:00 | 2026-11-08T08:32:59.983326+03:00 | 16 min 18 s |

## Реконструкция команд классификации

Команды полной классификации не сохранялись в журналах `*_run.log` (они содержат
только вывод QIIME 2 без команд и кодов). Команда восстановлена по provenance
(корневое действие `classify_sklearn`, параметры confidence=0.7, n_jobs=1,
вход rep_seq.qza, классификатор из data/refs/) и обозначена как **реконструированная**:

```
qiime feature-classifier classify-sklearn \
  --i-classifier data/refs/<модель>.qza \
  --i-reads work/qza/rep_seq.qza \
  --p-confidence 0.7 \
  --p-n-jobs 1 \
  --o-classification work/qza/taxonomy_sklearn_<метод>.qza
```

Исходный код первоначального процесса: **не сохранён** (журналы `*_run.log` не
содержат команд). Классификация не повторялась ради журнала.

## Экспорты

- `work/qzv_export/full_classification/silva138/taxonomy.tsv`
- `work/qzv_export/full_classification/gtdb/taxonomy.tsv`
- `work/qzv_export/full_classification/gg2/taxonomy.tsv`

## Ограничения

- Проверка подтверждает техническую корректность результатов (тип, набор ID,
  Confidence, provenance, UUID). Биологическая правильность таксономических
  назначений и соответствие региона классификатора региону данных (V4–V5)
  этой проверкой не устанавливаются.
- Kr5_S84: 1505/60078 = 2.51% сохранённых ридов — согласовано ранее; потери
  не исследуются (решение пользователя).
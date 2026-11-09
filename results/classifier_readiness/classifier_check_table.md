# Проверка готовых классификаторов и результатов

Актуализировано 2026-11-09. Это таблица технических проверок, без сравнения таксономических назначений.

Все три файла имеют тип TaxonomicClassifier. Обучение: fit_classifier_naive_bayes, QIIME 2 / feature-classifier 2024.5.0, scikit-learn 1.4.2. Текущее окружение: QIIME 2 / feature-classifier / q2cli 2025.7.0, scikit-learn 1.4.2.

| Метод | UUID классификатора | SHA256 | Наблюдаемые сведения | Проба / полный результат |
|---|---|---|---|---|
| sklearn_silva138 | 70b4b5f4-8fce-40bd-b508-afacbc12a5ed | c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616 | get_silva_data: SILVA 138.1, SSURef_NR99; filter_seqs_length_by_taxon 900/1200/1400 | PASS / PASS |
| sklearn_gtdb | 5d5461cc-6a51-434b-90ab-040f388e4221 | 07aadcf7472d9cc6f853f6b4615348619f1a3eceb56c1fb1b6d8dbb20554765f | Импорты bac120_ssu_reps_r220.fna и ar53_ssu_reps_r220.fna | PASS / PASS |
| sklearn_gg2 | 49ccfb0a-155d-404b-80b4-818d2aeb53b2 | 5c5900c10ad601b12cf0d2280d52cc0f7852627e89070595cc3899a6d690ee27 | Импорты 2024.9.backbone.fna и 2024.9.backbone.tax | PASS / PASS |

## Проверенные условия

Пять исходных проб: точный набор первых пяти ID из текущего rep_seq.qza, совпадение последовательностей с probe5_seqs.qza, отсутствие дубликатов, Confidence в [0,1], правильные типы, действие, параметры и фактические UUID.
Полные результаты: 3374 уникальных ASV в каждом, точное совпадение набора с rep_seq.qza, конечные Confidence в [0,1], непустой Taxon; confidence=0.7, n_jobs=1. Все validate/export — 0; три выбранных ASV на модель сверены с экспортом.
Числа 0 здесь — коды новых проверок, не восстановленные коды первоначальной классификации.

## Доказательства

- results/full_classification/repair_20261109/initial/readiness.stdout.txt
- results/full_classification/repair_20261109/full_verification_retry/summary.json
- results/full_classification/repair_20261109/model_provenance.json
- results/full_classification/repair_20261109/environment.stdout.txt
- results/full_classification/repair_20261109/project_check.stdout.txt

## Ограничения происхождения

Наблюдаемые версии/имена импортов перечислены выше; внешний источник получения готового QZA и его автор не установлены. Регион обучения, праймеры и полная длина гена не подтверждены. Фильтр минимальной длины SILVA не доказывает полноту гена, а слова full-length в имени Greengenes не заменяют проверку исходных последовательностей. Подготовка до импорта неизвестна. VERSIONS.md согласован с этими ограничениями. Прежние журналы сохранены как история.

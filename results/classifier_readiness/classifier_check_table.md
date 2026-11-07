# Проверка трёх классификаторов

Проверено 2026-11-07. Окружение: qiime2-amplicon-2025.7, q2cli 2025.7.0, scikit-learn 1.4.2.

| Метод | Файл | SHA256 | UUID | Тип артефакта | Сведения об обучении (provenance) | Результат пробы | Ограничения |
|---|---|---|---|---|---|---|---|
| sklearn_silva138 | data/refs/silva-138-99-nb-classifier.qza | c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616 | 70b4b5f4-8fce-40bd-b508-afacbc12a5ed | TaxonomicClassifier | QIIME 2 2024.5.0, feature-classifier 2024.5.0, fit_classifier_naive_bayes; источник: rescript get_silva_data SILVA 138.1 SSURef_NR99; фильтрация по длине 900/1200/1400 (полная длина); scikit-learn 1.4.2 | PASS (exit 0, 5 ASV, Confidence∈[0,1]) | Происхождение региона/праймеров из имени файла не выводилось; полная длина подтверждена фильтрацией по длине, регион не установлен |
| sklearn_gtdb | data/refs/gtdb_classifier_r220.qza | 07aadcf7472d9cc6f853f6b4615348619f1a3eceb56c1fb1b6d8dbb20554765f | 5d5461cc-6a51-434b-90ab-040f388e4221 | TaxonomicClassifier | QIIME 2 2024.5.0, feature-classifier 2024.5.0, fit_classifier_naive_bayes; источник: импорт ar53_ssu_reps_r220.fna + bac120_ssu_reps_r220.fna (GTDB r220 SSU reps, археи+бактерии); scikit-learn 1.4.2 | PASS (exit 0, 5 ASV, Confidence∈[0,1]) | Полная длина SSU reps; регион/праймеры не установлены |
| sklearn_gg2 | data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza | 5c5900c10ad601b12cf0d2280d52cc0f7852627e89070595cc3899a6d690ee27 | 49ccfb0a-155d-404b-80b4-818d2aeb53b2 | TaxonomicClassifier | QIIME 2 2024.5.0, feature-classifier 2024.5.0, fit_classifier_naive_bayes; источник: импорт 2024.9.backbone.fna (Greengenes 2024.09 backbone); scikit-learn 1.4.2 | PASS (exit 0, 5 ASV, Confidence∈[0,1]) | Полная длина (backbone); регион/праймеры не установлены |

## Сопоставление с VERSIONS.md

- SHA256 всех трёх файлов совпадают с data/refs/VERSIONS.md.
- Версия scikit-learn в каждом классификаторе (1.4.2) совпадает с установленной в окружении (1.4.2) — загрузка совместима.
- VERSIONS.md указывает источник «QIIME 2» и «full-length» для всех трёх. Provenance подтверждает: обучены в QIIME 2 2024.5.0 (feature-classifier 2024.5.0), полная длина. Конкретный регион/праймеры в provenance не указаны — «не установлено».

## Примечание о происхождении

«QIIME 2» в VERSIONS.md — это среда обучения, а не источник референсной базы. Источники баз установлены из provenance: SILVA 138.1 (NR99), GTDB r220 (SSU reps), Greengenes 2024.09 (backbone).
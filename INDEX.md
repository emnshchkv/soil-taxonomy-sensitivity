# INDEX.md — где что лежит

- `AGENTS.md` — правила работы агента в проекте.
- `STATUS.md` — где работа сейчас; читать в начале сессии.
- `diary.md` — журнал проекта; только дописывается.
- `cascade.sh` — <назначение уточнить>.
- `data/raw/` — риды `*_R1_001.fastq.gz` и `*_R2_001.fastq.gz` по 11 образцам; только чтение; в анализ идёт R1.
- `data/raw.sha256` — контрольные суммы `data/raw/`.
- `data/meta/metadata.tsv` — метаданные: `sample-id`, `sample_type` (`negative_control` или `soil`), `site_code`, `site`, `file_r1`.
- `data/meta/fq-manifest-single.tsv` — манифест QIIME 2 (`SingleEndFastqManifestPhred33V2`, только R1).
- `data/refs/silva-138-99-nb-classifier.qza` — классификатор метода `sklearn_silva138`.
- `data/refs/gtdb_classifier_r220.qza` — классификатор метода `sklearn_gtdb`.
- `data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza` — классификатор метода `sklearn_gg2`.
- `data/refs/VERSIONS.md` — источник, дата и sha256 классификаторов.
- `src/count_reads.py` — подсчёт ридов (шаг 1); пишет `results/read_counts_R1.tsv`.
- `checks/check.py` — критерий готовности этапа; только чтение.
- `env/environment.yml`, `env/explicit.txt` — описание окружения `qiime2-amplicon-2025.7`; только справка.
- `memory/denoise_choice.md` — `trim-left` и `trunc-len`, решение человека.
- `memory/L2.md` — сжатие основной статьи проекта; уровни сжатия `memory/L1.md`, `memory/L3.md` и по частям `memory/L1/part_01.md`–`part_11.md`.
- `memory/state.md` — прежняя сводка состояния, заменена `STATUS.md`; не читать.
- `work/qza/raw_reads_single.qza` — импорт ридов; `work/qzv/reads_summary.qzv` — сводка качества. Остальные артефакты QIIME 2 — тоже в `work/qza/` и `work/qzv/`; git их не отслеживает.
- `results/read_counts_R1.tsv` — число и длина ридов по образцам; остальные итоговые таблицы — тоже в `results/`.
- `paper/` — литература: `bokulich18.pdf`, `mcdonald23.pdf`, `parks21.pdf`, `robeson21.pdf`; `robeson21.txt` — текст статьи; `paper/parts/part_01.md`–`part_11.md` — статья по частям.

- src/verify_classifier_readiness.py — проверка актуальной статистики и трёх проб на пяти ASV; без расчётов.
- src/verify_full_classification.py — проверка готовых полных таксономий, validate/export и новые полные журналы.
- results/full_classification/repair_20261109/README.md — актуальные доказательства, воспроизведение, регрессии и ограничения исторических журналов.
- reports/Готовность к классификации.md — текущий отчёт о готовности и выполненной полной классификации.


## Этап 4 — сравнение готовых таксономий

- [[METHOD_COMPARISON.md]] — предварительная методика, дополнения и журнал попыток.
- [[src/taxonomy_comparison.py]] — расчёт; [[src/check_taxonomy_comparison.py]] — независимая проверка и полный вердикт.
- [[tests/test_taxonomy_comparison.py]] — синтетические испытания; [[tests/fixtures/taxonomy-comparison/README.md]] — ручной эталон.
- [[results/taxonomy_comparison/finish_20261116/README.md]] — входы, команды, таблицы и воспроизведение.
- [[reports/Сравнение таксономических профилей.md]] — численные результаты и ограничения.
- [[Завершить этап сравнения в отдельной ветке]] — текущая задача завершения.

## Подготовка статьи AI4Science-YSDA 2026

- [[AGENT_ARTICLE_PROMPT.txt]] — готовый промпт для написания и проверки методологической статьи по этапам 1–4.
- [[Проверка требований шаблона статьи AI4Science 2026]] — формат, структура и правила шаблона.
- [[Проверка готовности проекта к написанию статьи]] — актуальный аудит результатов и оставшихся вопросов подготовки.
- [[Проверка промпта агента для статьи]] — проверка полноты, источников и границ промпта.
- [[Подготовить промпт для написания статьи и проверить готовность проекта]] — выполненное поручение; статья ещё не написана.
- [[Проверка результата агента по написанию статьи]] — аудит запуска 17 ноября: сделанное, остановка по лимиту, препятствия и состояние научных результатов.

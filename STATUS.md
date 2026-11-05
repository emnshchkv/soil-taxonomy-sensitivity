# Состояние проекта

Выполнена [[Исправить и проверить верификатор DADA2]]. План — коммит 6d12744; исправленная реализация и синтетические испытания — 866ae21. --self-test PASS/0: все случаи METHOD.md §5 и семь свойств на четырёх размерах матрицы. Интеграционная проверка настоящих архивных QZA через штатный QIIME 2 и сверка архивных профилей по текущим FASTQ тоже PASS/0. Сохранены results/caught_error_before.txt (FAIL/1) и results/caught_error_after.txt (PASS/0).

Все 67 исходных защищённых файлов, включая checks/check.py, memory/denoise_choice.md и METHOD.md §§0–7, неизменны; состав data/work не изменён. checks только для чтения. src/qiime_runtime.py выбирает Python и R одного существующего окружения для процесса.

## Выполнено в текущей сессии

1. **Профили качества** — `src/quality_profiles.py`, `results/quality_by_position.tsv`, `results/quality_by_sample.tsv`, `results/quality_by_sample.png`. Проверка `--quality-only`: PASS/0.
2. **Решение человека** — `memory/denoise_choice.md`: trim-left=20, trunc-len=200, исключение Kr5_S84 < 30%.
3. **Три запуска DADA2** — `work/qza/` (запуск 1), `work/denoise_repeats/02/` (запуск 2), `work/denoise_repeats/03/` (запуск 3). Параметры: max-ee=3, n-threads=6, min-fold-parent-over-abundance=2.
4. **Проверка совпадения** — `src/check_exact_match.py` использует `validate` и `run_metrics` из верификатора. Все четыре метрики = 0 для трёх пар.
5. **Воспроизведение** — самотест PASS/0, CASE_REPRO_TEMP: FAIL/1 → PASS/0 после repair.
6. **METHOD.md §§8–10** — времена из provenance, токены неизвестно, повторяемость ≠ правильность.

## Текущий вердикт

- **Повторяемость**: PASS/0 (метрики 0,0,0,0)
- **Полный режим**: INCOMPLETE/2 (Kr5_S84 retention 2.51% < 30%, пользователь разрешил)
- **Классификация и сравнение методов**: не начинались

## Дальнейшая работа

1. Заполнить REPORT.md
2. Создать метку hw2
3. Классификация тремя методами (sklearn_silva138, sklearn_gtdb, sklearn_gg2)
4. Сравнение классификаторов
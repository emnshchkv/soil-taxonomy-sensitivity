# Проверка верификатора

Из корня проекта, в существующем окружении:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/verify_denoising.py --self-test
```

Самотест включает все случаи METHOD.md §5, известный ответ для квантилей, сохранённые позиционные и сводные TSV, приоритет FAIL над отсутствующими данными и семь свойств на четырёх размерах матрицы. Минимальные синтетические QZA описаны в `fixtures/denoise-clean/README.md`; они не заменяют проверку настоящего формата.

Сверка неизменности проектных входов, результатов и fixtures, включая права, плюс отдельные заполненные файлы-сторожа:

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 tests/test_verifier_integration.py
```

Во время этого теста другие процессы не должны записывать файлы в `results/`, `work/`, `memory/` и `data/`: изменение любого файла будет обнаружено.

Для независимой проверки настоящего читателя добавьте `--real-run /путь/к/архивному/набору` с тремя QZA одного запуска 20/200. Тест выполняет штатные QIIME 2 `Artifact.validate(level="max")`, `view` и `export_data` и сравнивает все последовательности, ячейки и счётчики. Он не запускает DADA2. Дополнительный `--quality-archive /путь/к/архиву.tar` читает только `results/quality_by_position.tsv` и `results/quality_by_sample.tsv` из локального архивного снимка; пересчёт выполняется по текущим FASTQ, повреждение и восстановление — во временной копии. Частные архивы не нужны публичному самотесту.

Воспроизведение внесённой ошибки (CASE — новое, ещё не существующее имя в временном каталоге):

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/denoise_examples.py create "$CASE"
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/verify_denoising.py --example-dir "$CASE"
# Ожидаются REPEAT_DIFFERENCE и код 1.
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/denoise_examples.py repair "$CASE"
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/verify_denoising.py --example-dir "$CASE"
# Ожидается PASS/0.
```

Журналы: `results/verifier_selftest.txt`, `results/verifier_integration.txt`, `results/caught_error_before.txt`, `results/caught_error_after.txt`. До/после сверки исходных файлов — `results/verifier_protected_*.json`. Неудачные интеграционные попытки сохранены отдельно: ограничение CPU в sandbox, параллельное создание журналов и смешение системного R с conda. Последнее устранено выбором PATH внутри процесса; пакеты не изменены.

Полный научный вердикт остаётся отдельной командой без аргументов: отсутствие новых профилей, решения и трёх реальных результатов даёт INCOMPLETE/2. Здесь проверяется реализация верификатора, а не повторяется научный эксперимент.

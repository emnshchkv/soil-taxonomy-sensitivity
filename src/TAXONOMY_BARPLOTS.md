# Код таксономических барплотов

Основной набор — семейства Bacteria + Archaea по SILVA, GTDB и Greengenes2. [Готовые картинки](../figures/taxonomy/family/README.md), [таблицы и правила](../results/taxonomy_barplots/prokaryote_family_20261118_final/README.md), [диагностика Kr](../results/taxonomy_barplots/kr_diagnostics_20261118/README.md).

| Файл | Назначение |
|---|---|
| [plot_prokaryotic_families.py](plot_prokaryotic_families.py) | Читает готовые QZA, отбирает бактерии и археи, рассчитывает семейства и вызывает рендерер |
| [render_prokaryotic_families.py](render_prokaryotic_families.py) | Рисует общий обзор и отдельные модели в PNG/PDF/SVG |
| [check_prokaryotic_families.py](check_prokaryotic_families.py) | Независимо проверяет отбор, числа, группы и контрольные суммы; поддерживает изолированный самотест |
| [audit_kr_taxonomy.py](audit_kr_taxonomy.py) | Прослеживает расхождения Kr до ASV и исходных назначений; проверяет прежнюю фильтрацию |

Скрипты используют соседние модули `taxonomy_compare_io.py`, `taxonomy_compare_core.py` и `classification_contract.py`. Запускать из корня репозитория в существующем окружении `qiime2-amplicon-2025.7`; нужны локальные входы `work/qza/` и `data/meta/`. Их точные версии и SHA256 записаны в манифесте результата. В Git исходные QZA не включены.

## Проверить готовые картинки и числа

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_prokaryotic_families.py --out-dir results/taxonomy_barplots/prokaryote_family_20261118_final
```

Дополнительный флаг `--self-test` проверяет обнаружение подмены в изолированной временной копии.

## Построить новый комплект

Каталог `family_new` должен отсутствовать. Классификация и DADA2 повторно не запускаются.

```sh
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/plot_prokaryotic_families.py --out-dir results/taxonomy_barplots/family_new --min-percent 2
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_prokaryotic_families.py --out-dir results/taxonomy_barplots/family_new
```

Рисунки появятся в `results/taxonomy_barplots/family_new/figures/`, таблицы и манифест — рядом. `--min-percent` регулирует визуальное объединение малых долей, а не отбор ASV. Каталог `figures/taxonomy/family/` содержит сохранённые копии окончательного запуска от 2026-11-18; новый расчёт их не перезаписывает.

## Предыдущий набор только с бактериями

[plot_taxonomy_profiles.py](plot_taxonomy_profiles.py), [render_taxonomy_barplots.py](render_taxonomy_barplots.py), [check_taxonomy_barplots.py](check_taxonomy_barplots.py) и [check_taxonomy_domain_audit.py](check_taxonomy_domain_audit.py) относятся к прежнему набору пяти рангов. Его [описание и команды](../results/taxonomy_barplots/BACTERIA_ONLY_20261118.md) сохранены отдельно; текущий запрос с археями выполняют скрипты из таблицы выше.

## ai4science doctor

```
=== AI4Science doctor ===
...
MCP research: OK (5 инструментов: openalex_search, openalex_work, openalex_links, wikidata_search, wikidata_sparql)
MCP family-profile: OK (1 инструментов: family_profile)
...
```

## Запрос пользователя

```
> Для пробы Kr2_S81 покажи, сколько ридов вошло в семейственный профиль у каждой из трёх моделей и у какой модели доля ридов без назначенного семейства неопределима. Объясни почему.
```

## Вызов инструмента агентом

```
Запрашиваю данные рисунка по пробе Kr2_S81.

  ────────────────────────────────────────
  ▸ family_profile family-profile-mcp
    sample_id: Kr2_S81
```

## Результат работы агента


#### Kr2_S81 — риды в семейственном профиле

Источник: `results/paper_step4/20261118/family_sample_counts.tsv` (строка `family_prokaryote_selected`), сверено с `results/taxonomy_barplots/prokaryote_family_20261118_final/filter_counts.tsv`. После DADA2 у пробы 8904 рида.

| Модель           | Вошло в профиль (N) | Без семейства |    Доля без семейства |
|------------------|---------------------|---------------|-----------------------|
| sklearn_silva138 |               **0** |             0 | **NA (неопределима)** |
| sklearn_gtdb     |            **8013** |          8013 |                 100 % |
| sklearn_gg2      |            **8098** |          8098 |                 100 % |
**Неопределима у SILVA.** Доля считается как `unclassified_family_reads / selected_reads`, а у SILVA `selected_reads = 0`: знаменатель нулевой, деление не определено, поэтому в таблице стоит `NA` (по контракту `research/paper_results_contract.md`: «При N=0 доля NA»). Это не пропуск данных и не ошибка агрегации.

**Почему у SILVA ноль отобранных.** Разбиение всех 8904 ридов по назначению SILVA (`kr_diagnostics_20261118/README.md`, `filter_counts.tsv`):

- 7780 — Eukaryota (исключаются);
- 1113 — неизвестный домен (исключаются);
- 11 — Bacteria, но все 11 несут явную органелльную метку (Mitochondria/Chloroplast), поэтому тоже исключаются;
- Archaea — 0.

Итого Bacteria + Archaea после отбора = 0. Отбор идёт отдельно по каждой модели, поэтому GTDB и Greengenes2 на тех же ридах сохраняют 8013 и 8098: у них эти ASV назначены как Bacteria (например, ASV `121fc77a…` — 1002 рида, который SILVA относит к Eukaryota/Didymellaceae, а GTDB и GG2 только к Bacteria). У GTDB/GG2 все отобранные риды имеют лишь `d__Bacteria` без семейства — отсюда 100 %, но знаменатель ненулевой, и доля определена.

Оговорка: это наблюдаемые назначения моделей, а не установленное биологическое происхождение ASV; «лучшая» модель здесь не выбирается.
___
---
type: Проверка
---
# Проверка переименования каталога статьи в paper

2026-11-18: по подтверждению пользователя каталог manuscript/ переименован в paper/. Статья — [[paper/ai4science_ysda2026/paper.pdf]], исходник — [[paper/ai4science_ysda2026/paper.tex]]. Сохранены все 27 файлов, включая служебные файлы сборки; PDF, TeX, BibTeX, стили и вычитка побайтно прежние. В README.md и claims.tsv обновлены только пути.

Рабочие ссылки, промпты, .gitignore, команды сборки/чтения и проверки используют paper/. Проверка литературы отделяет текущий каталог от прежнего source_base в Git; исторический манифест, журналы и исходные критерии завершённых задач сохранены. Каталог references/ по-прежнему содержит литературу.

Проверки src/check_article.py, src/check_literature_registry.py, src/check_paper_scaffold.py, checks/check.py и существующие tests/test_article.py, tests/test_literature_registry.py через conda qiime2-amplicon-2025.7 завершились с кодом 0.
Финальный проверяющий скрипт статьи выполнил 646 сверок, в том числе проверку PDF. Существующие испытания статьи и литературы проходят после изменения путей. 168 файлов прежних свидетельств, замороженных чисел, литературы, рисунков и защищённых проверок сохранены. Научные расчёты и сборка PDF не повторялись, подача не выполнялась.

Свидетельства: [[results/article_relocation/20261118/before.json]], [[results/article_relocation/20261118/verification.json]], [[results/article_relocation/20261118/article.json]], [[results/article_relocation/20261118/literature.json]].

## Связи
- тип: [[Проверка]]
- решает: [[Переименовать каталог рукописи в paper]]

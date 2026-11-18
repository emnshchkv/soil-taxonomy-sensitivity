---
type: Проверка
---
# Проверка исключения промптов статьи из Git

2026-11-18. По запросу пользователя из индекса Git исключены только AGENT_ARTICLE_PROMPT.txt и AGENT_ARTICLE_STEP5_PROMPT.txt. Для них добавлены отдельные правила в .gitignore; локальные файлы сохранены без изменения. results/, tasks/ и другие пути из отслеживания не исключались.

Проверка: git ls-files -- AGENT_ARTICLE_PROMPT.txt AGENT_ARTICLE_STEP5_PROMPT.txt — пустой вывод; git check-ignore -v показывает оба правила .gitignore. SHA256 до и после совпадают:

- AGENT_ARTICLE_PROMPT.txt: ac8c36859bc464ec1079aea76daeed6af91227f4dfb344b890b6f2ed2037581b.
- AGENT_ARTICLE_STEP5_PROMPT.txt: 493c98956c2b39b6f7036975cd154d5cf7f1c0be3501390959473bad2b1afa65.

conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 checks/check.py — код 0; q2cli version 2025.7.0. Старые коммиты не переписаны; push не выполнялся. Само добавление уже отслеживаемых файлов в .gitignore не удаляет их из Git; удаление из индекса меняет следующую версию, но не старые коммиты. Передача такого коммита на сервер — отдельный push.

## Связи
- тип: [[Проверка]]
- решает: [[Исключить промпты статьи из отслеживания Git]]

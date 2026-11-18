"""Read-only checks for the adapted paper scaffold, not for article acceptance.

Run through the project's qiime2-amplicon-2025.7 environment.  --layout-only
checks the first subtask without requiring the glossary.  No scientific
verifier, LaTeX command, network request, or filesystem write is performed.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


LAYOUT_FILES = (
    "AGENTS.md", "STATUS.md", "INDEX.md", "diary.md", "PAPER_PROJECT.md",
    "WRITING.md", "ontology.yaml", "research/problem.md", "lit/lit.md",
    "src/check_paper_scaffold.py",
)
LINK_DOCUMENTS = ("PAPER_PROJECT.md", "WRITING.md", "research/problem.md", "lit/lit.md")
SYMBOL_HEADER = ("Обозначение", "Значение", "Единицы", "В коде или данных")
TERM_HEADER = ("Термин", "Значение", "В коде или данных")
LINK_RE = re.compile(r"(?<!!)\[[^\]\n]+\]\(\s*(<[^>\n]+>|[^)\n]+)\s*\)")


def without_fences(text: str) -> str:
    """Ignore examples inside fenced code blocks while retaining line numbers."""
    lines = []
    fence = None
    for line in text.splitlines():
        match = re.match(r"^\s*(`{3,}|~{3,})", line)
        if match:
            marker = match.group(1)
            if fence is None:
                fence = marker
            elif marker[0] == fence[0] and len(marker) >= len(fence):
                fence = None
            lines.append("")
        else:
            lines.append(line if fence is None else "")
    return "\n".join(lines)


def local_links(text: str, document: Path) -> list[tuple[str, Path]]:
    links = []
    for match in LINK_RE.finditer(text):
        target = match.group(1).strip()
        if target.startswith("<") and target.endswith(">"):
            target = target[1:-1]
        else:
            target = re.split(r"\s+[\"']", target, maxsplit=1)[0]
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        links.append((target, (document.parent / unquote(parsed.path)).resolve()))
    return links


def table_cells(line: str) -> tuple[str, ...]:
    return tuple(cell.strip() for cell in re.split(r"(?<!\\)\|", line.strip().strip("|")))


def markdown_tables(text: str):
    lines = without_fences(text).splitlines()
    for index in range(len(lines) - 1):
        if not lines[index].lstrip().startswith("|"):
            continue
        header = table_cells(lines[index])
        separator = table_cells(lines[index + 1])
        if len(header) != len(separator) or not all(re.fullmatch(r":?-{3,}:?", cell) for cell in separator):
            continue
        rows = []
        for row_index in range(index + 2, len(lines)):
            if not lines[row_index].lstrip().startswith("|"):
                break
            rows.append((row_index + 1, table_cells(lines[row_index])))
        yield header, rows


def check_scaffold(root: Path, layout_only: bool = False) -> list[str]:
    errors: list[str] = []
    documents: dict[str, str] = {}

    def require(ok: bool, reason: str) -> None:
        if not ok:
            errors.append(reason)

    def read(relative: str) -> str:
        if relative in documents:
            return documents[relative]
        try:
            value = (root / relative).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            errors.append(f"{relative}: не удалось прочитать ({exc})")
            value = ""
        documents[relative] = value
        return value

    for relative in LAYOUT_FILES + (() if layout_only else ("GLOSSARY.md",)):
        require(bool(read(relative).strip()), f"{relative}: файл отсутствует или пуст")

    try:
        import yaml
    except ImportError:
        errors.append("PyYAML недоступен: используйте существующее окружение qiime2-amplicon-2025.7; пакеты автоматически не устанавливаются")
        return errors

    def yaml_mapping(text: str, label: str) -> dict:
        try:
            value = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            errors.append(f"{label}: неверный YAML ({exc})")
            return {}
        if not isinstance(value, dict):
            errors.append(f"{label}: ожидается YAML-словарь")
            return {}
        return value

    def frontmatter(text: str, label: str) -> dict:
        match = re.match(r"\A---\s*\n(.*?)\n---(?:\s*\n|$)", text, re.S)
        if not match:
            errors.append(f"{label}: нет YAML-полей в начале файла")
            return {}
        return yaml_mapping(match.group(1), label)

    for relative in ("research/problem.md", "lit/lit.md"):
        metadata = frontmatter(read(relative), relative)
        require(metadata.get("state") == "deferred", f"{relative}: ожидалось state: deferred; шаги постановки и литературы не входят в приёмку этого каркаса")

    ontology = yaml_mapping(read("ontology.yaml"), "ontology.yaml")
    classes = ontology.get("classes", {})
    if not isinstance(classes, dict):
        errors.append("ontology.yaml: classes должен быть словарём")
        classes = {}
    rules = ontology.get("rules", [])
    require(len(classes) >= 5, "ontology.yaml: требуется не меньше пяти классов")
    require(isinstance(rules, list) and len(rules) >= 3, "ontology.yaml: требуется не меньше трёх правил")
    for name, definition in classes.items():
        if not isinstance(name, str) or "/" in name or "\\" in name:
            errors.append(f"ontology.yaml: недопустимое имя класса {name!r}")
            continue
        require(isinstance(definition, dict), f"ontology.yaml: описание класса {name} должно быть словарём")
        if isinstance(definition, dict):
            parent = definition.get("parent")
            require(parent is None or parent in classes, f"ontology.yaml: неизвестный родитель {parent!r} у {name}")
        relative = f"ontology/{name}.md"
        note = read(relative)
        require(frontmatter(note, relative).get("kind") == "класс", f"{relative}: требуется kind: класс")
        require(bool(re.search(rf"^# {re.escape(name)}\s*$", note, re.M)), f"{relative}: заголовок должен называть класс {name}")
    properties = ontology.get("properties", {})
    if not isinstance(properties, dict):
        errors.append("ontology.yaml: properties должен быть словарём")
        properties = {}
    for name, definition in properties.items():
        if not isinstance(definition, dict):
            errors.append(f"ontology.yaml: отношение {name} должно быть словарём")
            continue
        for field in ("domain", "range"):
            value = definition.get(field)
            values = value if isinstance(value, list) else [value]
            require(bool(values) and all(isinstance(item, str) and item in classes for item in values), f"ontology.yaml: неизвестный или пустой {field} отношения {name}")

    for relative, keys in (
        (".obsidian/app.json", ("showUnsupportedFiles",)),
        (".obsidian/graph.json", ("showAttachments", "showOrphans", "showArrow")),
    ):
        try:
            settings = json.loads(read(relative))
        except json.JSONDecodeError as exc:
            errors.append(f"{relative}: неверный JSON ({exc})")
            continue
        require(isinstance(settings, dict), f"{relative}: ожидается JSON-объект")
        if isinstance(settings, dict):
            for key in keys:
                require(settings.get(key) is True, f"{relative}: {key} должен быть true")

    for relative in LINK_DOCUMENTS + (() if layout_only else ("GLOSSARY.md",)):
        for target, destination in local_links(without_fences(read(relative)), root / relative):
            require(destination.exists(), f"{relative}: не существует цель ссылки {target}")

    if not layout_only:
        glossary = read("GLOSSARY.md")
        tables = list(markdown_tables(glossary))
        for expected in (SYMBOL_HEADER, TERM_HEADER):
            matching = [(header, rows) for header, rows in tables if header == expected]
            require(bool(matching), f"GLOSSARY.md: нет таблицы с колонками {' | '.join(expected)}")
            for header, rows in matching:
                require(bool(rows), f"GLOSSARY.md: таблица {header[0]} пуста")
                seen = set()
                for line, cells in rows:
                    require(len(cells) == len(header) and all(cells), f"GLOSSARY.md:{line}: незаполненная или повреждённая строка таблицы")
                    if len(cells) != len(header):
                        continue
                    symbol = cells[0].strip("`$").strip()
                    require(symbol not in seen, f"GLOSSARY.md:{line}: повтор обозначения/термина {symbol}")
                    seen.add(symbol)
                    sources = local_links(cells[-1], root / "GLOSSARY.md")
                    require(any(path.is_file() for _, path in sources), f"GLOSSARY.md:{line}: колонка «В коде или данных» должна ссылаться на существующий локальный файл")
        headings = [line.lower() for line in glossary.splitlines() if line.startswith("#")]
        require(any(re.search(r"этап\s*4|сравнен", line) for line in headings), "GLOSSARY.md: нет раздела обозначений исходного сравнения")
        require(any("семейств" in line for line in headings), "GLOSSARY.md: нет отдельного раздела семейственных профилей")
        require(any(cells and cells[0].strip("`$").strip() == "NA" for _, rows in tables for _, cells in rows), "GLOSSARY.md: нет отдельной строки с определением NA")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1], help="корень проекта (по умолчанию — родитель src/)")
    parser.add_argument("--layout-only", action="store_true", help="проверить только размещение и онтологию, без GLOSSARY.md")
    args = parser.parse_args()
    errors = check_scaffold(args.root.resolve(), args.layout_only)
    if errors:
        print("FAIL/1 — каркас проекта статьи")
        for error in errors:
            print(f"- {error}")
        return 1
    scope = "размещение и онтология; глоссарий не проверялся" if args.layout_only else "размещение, онтология и структура глоссария"
    print(f"PASS/0 — каркас проекта статьи: {scope}")
    print("Это не приёмка научных результатов, содержания статьи или готовности к подаче; научный смысл определений проверяет автор.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

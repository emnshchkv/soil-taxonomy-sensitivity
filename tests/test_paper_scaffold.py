"""CLI contracts for the paper scaffold; all fixtures are synthetic and temporary.

Run with the existing qiime2-amplicon-2025.7 Python. No real scientific input
is copied, and the verifier is invoked from a working directory outside the
project to check its explicit --root contract.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER = Path(__file__).resolve().parents[1] / "src" / "check_paper_scaffold.py"
SYMBOL_ROW = "| p | Доля | безразмерная доля | [Код](src/example.py) |"
GLOSSARY = """# Синтетический глоссарий

## Этап 4 — сравнение

| Обозначение | Значение | Единицы | В коде или данных |
| --- | --- | --- | --- |
{symbol_row}

## Семейственные профили

| Обозначение | Значение | Единицы | В коде или данных |
| --- | --- | --- | --- |
| N | Знаменатель | риды | [Код](src/example.py) |

| Термин | Значение | В коде или данных |
| --- | --- | --- |
| NA | Доля не определена при нулевом знаменателе | [Код](src/example.py) |
""".format(symbol_row=SYMBOL_ROW)


class ScaffoldCLIContracts(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="paper-scaffold-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "project"
        self.root.mkdir()
        for name in ("AGENTS.md", "STATUS.md", "INDEX.md", "diary.md", "WRITING.md"):
            self.write(name, "# Синтетический файл\n")
        self.write("PAPER_PROJECT.md", "# Каркас\n\n[Пример](src/example.py)\n")
        self.write("src/example.py", "# Synthetic code reference; never executed.\n")
        self.write("src/check_paper_scaffold.py", "# Layout marker only; CLI uses the real checker.\n")
        for name in ("research/problem.md", "lit/lit.md"):
            self.write(name, "---\nstate: deferred\n---\n# Отложено\n")
        self.write(".obsidian/app.json", json.dumps({"showUnsupportedFiles": True}))
        self.write(".obsidian/graph.json", json.dumps({
            "showAttachments": True, "showOrphans": True, "showArrow": True,
        }))
        self.write("ontology.yaml", """classes:
  Образец: {parent: null, definition: проба}
  Данные: {parent: null, definition: набор}
  Метод: {parent: null, definition: способ}
  Результат: {parent: null, definition: выход}
  Проверка: {parent: null, definition: сверка}
properties:
  входитВ: {domain: Образец, range: Данные}
  проверенЧерез: {domain: Результат, range: Проверка}
rules:
  - {id: П1, text: образец входит в данные}
  - {id: П2, text: результат проверен}
  - {id: П3, text: метод описан}
""")
        for name in ("Образец", "Данные", "Метод", "Результат", "Проверка"):
            self.write(f"ontology/{name}.md", f"---\nkind: класс\n---\n# {name}\n")
        self.write("GLOSSARY.md", GLOSSARY)

    def write(self, relative: str, content: str):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def replace(self, relative: str, old: str, new: str):
        path = self.root / relative
        text = path.read_text(encoding="utf-8")
        self.assertIn(old, text, "Fault injection must target existing fixture text")
        path.write_text(text.replace(old, new), encoding="utf-8")

    def snapshot(self):
        return {
            str(path.relative_to(self.root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in self.root.rglob("*") if path.is_file()
        }

    def run_checker(self, *options: str):
        before = self.snapshot()
        result = subprocess.run(
            [sys.executable, str(CHECKER), "--root", str(self.root), *options],
            cwd=self.temporary.name, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(before, self.snapshot(), "Checker must not change the supplied project")
        return result

    def assert_failure(self, result, reason: str):
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("FAIL/1 — каркас проекта статьи", result.stdout)
        self.assertIn(reason, result.stdout)
        self.assertNotIn("Traceback", result.stderr)

    def test_clean_scaffold_passes(self):
        result = self.run_checker()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS/0 — каркас проекта статьи", result.stdout)
        self.assertIn("Это не приёмка научных результатов", result.stdout)

    def test_problem_draft_passes_without_claiming_content_acceptance(self):
        self.replace("research/problem.md", "state: deferred", "state: draft")
        result = self.run_checker()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Это не приёмка научных результатов, содержания статьи", result.stdout)

    def test_unknown_problem_state_is_rejected(self):
        self.replace("research/problem.md", "state: deferred", "state: accepted")
        self.assert_failure(self.run_checker(), "research/problem.md: ожидалось state: deferred или draft")

    def test_literature_draft_passes_without_claiming_content_acceptance(self):
        self.replace("lit/lit.md", "state: deferred", "state: draft")
        result = self.run_checker()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Это не приёмка научных результатов, содержания статьи", result.stdout)

    def test_unknown_literature_state_is_rejected(self):
        self.replace("lit/lit.md", "state: deferred", "state: accepted")
        self.assert_failure(self.run_checker(), "lit/lit.md: ожидалось state: deferred или draft")

    def test_layout_only_does_not_accept_missing_glossary_as_complete(self):
        (self.root / "GLOSSARY.md").unlink()
        layout = self.run_checker("--layout-only")
        self.assertEqual(layout.returncode, 0, layout.stdout + layout.stderr)
        self.assertIn("глоссарий не проверялся", layout.stdout)
        self.assert_failure(self.run_checker(), "GLOSSARY.md: файл отсутствует или пуст")

    def test_broken_local_link_is_rejected(self):
        self.replace("PAPER_PROJECT.md", "src/example.py", "src/missing.py")
        self.assert_failure(self.run_checker(), "PAPER_PROJECT.md: не существует цель ссылки src/missing.py")

    def test_duplicate_symbol_is_rejected(self):
        self.replace("GLOSSARY.md", SYMBOL_ROW, SYMBOL_ROW + "\n" + SYMBOL_ROW)
        self.assert_failure(self.run_checker(), "повтор обозначения/термина p")

    def test_empty_unit_is_rejected(self):
        self.replace("GLOSSARY.md", "| безразмерная доля |", "|  |")
        self.assert_failure(self.run_checker(), "незаполненная или повреждённая строка таблицы")

    def test_malformed_table_row_is_rejected(self):
        self.replace("GLOSSARY.md", "| p | Доля |", "| p | лишняя колонка | Доля |")
        self.assert_failure(self.run_checker(), "незаполненная или повреждённая строка таблицы")

    def test_unknown_ontology_range_is_rejected(self):
        self.replace("ontology.yaml", "range: Данные", "range: НеизвестныйКласс")
        self.assert_failure(self.run_checker(), "неизвестный или пустой range отношения входитВ")


if __name__ == "__main__":
    unittest.main(verbosity=2)

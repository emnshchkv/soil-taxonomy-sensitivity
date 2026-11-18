"""Synthetic CLI checks for class-note creation, synchronization and read-only QA."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'src' / 'build_ontology_notes.py'
SCHEMA = '''classes:
  Объект:
    parent: null
    definition: общее понятие
  Подтип:
    parent: Объект
    definition: частное понятие
  Цель:
    parent: null
    definition: цель отношения
properties:
  направленК: {domain: Подтип, range: Цель}
  относитсяК: {domain: [Объект, Подтип], range: Цель}
attributes:
  режим: {domain: Подтип, values: [известен, неизвестен]}
rules: []
'''


class ClassNoteCLI(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='ontology-notes-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'project'
        self.root.mkdir()
        (self.root / 'ontology.yaml').write_text(SCHEMA, encoding='utf-8')

    def snapshot(self):
        return {str(p.relative_to(self.root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in self.root.rglob('*') if p.is_file()}

    def run_cli(self, *args, expected=0):
        before = self.snapshot()
        result = subprocess.run([sys.executable, str(SCRIPT), '--root', str(self.root), *args],
                                cwd=self.tmp.name, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual((self.root / 'ontology.yaml').read_text(), SCHEMA)
        if '--check' in args:
            self.assertEqual(before, self.snapshot(), '--check must not write files')
        return result

    def test_default_creates_missing_notes(self):
        result = self.run_cli()
        self.assertIn('ontology/Подтип.md', result.stdout)
        self.assertEqual(len(list((self.root / 'ontology').glob('*.md'))), 3)
        self.assertIn('PASS/0', self.run_cli('--check').stdout)

    def test_default_preserves_existing_note(self):
        directory = self.root / 'ontology'
        directory.mkdir()
        note = directory / 'Подтип.md'
        note.write_text('Предыдущий текст\n')
        self.run_cli()
        self.assertEqual(note.read_text(), 'Предыдущий текст\n')
        self.assertTrue((directory / 'Объект.md').exists())

    def test_check_detects_missing_without_creating_directory(self):
        result = self.run_cli('--check', expected=1)
        self.assertIn('отсутствует: ontology/Подтип.md', result.stdout)
        self.assertFalse((self.root / 'ontology').exists())

    def test_check_detects_stale_note_without_repairing_it(self):
        self.run_cli()
        note = self.root / 'ontology' / 'Подтип.md'
        note.write_text(note.read_text() + 'Ручная добавка\n')
        result = self.run_cli('--check', expected=1)
        self.assertIn('устарела: ontology/Подтип.md', result.stdout)
        self.assertTrue(note.read_text().endswith('Ручная добавка\n'))

    def test_update_repairs_stale_parent_and_direct_domain(self):
        self.run_cli()
        note = self.root / 'ontology' / 'Подтип.md'
        note.write_text('Устаревшая заметка\n')
        result = self.run_cli('--update')
        self.assertIn('Updated: ontology/Подтип.md', result.stdout)
        text = note.read_text()
        self.assertIn('- подклассОт: [[Объект]]', text)
        self.assertIn('Отношения от этого класса: направленК, относитсяК.', text)
        self.assertIn('Атрибуты: режим (известен, неизвестен).', text)
        target = (self.root / 'ontology' / 'Цель.md').read_text()
        self.assertIn('Отношения от этого класса: нет.', target)
        self.assertNotIn('подклассОт:', target)
        self.assertIn('PASS/0', self.run_cli('--check').stdout)

    def test_update_is_idempotent_and_preserves_unrelated_note(self):
        self.run_cli('--update')
        unrelated = self.root / 'ontology' / 'Личная заметка.md'
        unrelated.write_text('Не входит в схему\n')
        before = self.snapshot()
        result = self.run_cli('--update')
        self.assertEqual(before, self.snapshot())
        self.assertIn('Created: none', result.stdout)
        self.assertIn('Updated: none', result.stdout)
        self.assertIn('PASS/0', self.run_cli('--check').stdout)


if __name__ == '__main__':
    unittest.main(verbosity=2)

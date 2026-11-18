"""Generate class notes from ontology.yaml, without changing the schema.

The default creates missing notes and preserves existing ones.  --update
synchronizes schema class notes; --check only compares their canonical text.
Notes whose names are absent from the schema are never removed or changed.
"""
from qiime_runtime import selected_python
selected_python()

import argparse
from pathlib import Path
import yaml


def canonical_notes(root):
    """Render all expected notes before making any filesystem changes."""
    schema = yaml.safe_load((root / 'ontology.yaml').read_text())
    if not isinstance(schema, dict) or not isinstance(schema.get('classes'), dict):
        raise ValueError('ontology.yaml: classes должен быть словарём')
    classes = schema['classes']
    properties = schema.get('properties', {})
    attributes_schema = schema.get('attributes', {})
    if not isinstance(properties, dict) or not isinstance(attributes_schema, dict):
        raise ValueError('ontology.yaml: properties и attributes должны быть словарями')
    rendered = {}
    for name, item in classes.items():
        if (not isinstance(name, str) or not name.strip()
                or name in ('.', '..') or any(c in name for c in '/\\:#|[]')):
            raise ValueError(f'ontology.yaml: недопустимое имя класса {name!r}')
        if not isinstance(item, dict) or not isinstance(item.get('definition'), str):
            raise ValueError(f'ontology.yaml: отсутствует определение класса {name}')
        parent_name = item.get('parent')
        if parent_name is not None and parent_name not in classes:
            raise ValueError(f'ontology.yaml: неизвестный родитель {parent_name!r}')
        relations, attributes = [], []
        for relation, rule in properties.items():
            domain = rule['domain'] if isinstance(rule['domain'], list) else [rule['domain']]
            if name in domain:
                relations.append(str(relation))
        for attribute, rule in attributes_schema.items():
            domain = rule['domain'] if isinstance(rule['domain'], list) else [rule['domain']]
            if name in domain:
                attributes.append(str(attribute) + ' (' + ', '.join(map(str, rule['values'])) + ')')
        parent = f'- подклассОт: [[{parent_name}]]\n' if parent_name else ''
        rendered[name + '.md'] = (
            f"---\nkind: класс\n---\n# {name}\n\n{item['definition']}\n\n"
            f"Отношения от этого класса: {', '.join(relations) or 'нет'}.\n"
            f"Атрибуты: {'; '.join(attributes) or 'нет'}.\n\n## Связи\n{parent}"
        )
    return rendered


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent,
                        help='корень проекта; по умолчанию родитель src/')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--update', action='store_true', help='обновить существующие заметки классов по схеме')
    mode.add_argument('--check', action='store_true', help='сравнить заметки со схемой без записи')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        rendered = canonical_notes(root)
        directory = root / 'ontology'
        missing, stale = [], []
        for filename, expected in rendered.items():
            path = directory / filename
            if path.is_symlink():
                raise ValueError(f'Ссылка вместо заметки класса: {path}')
            if not path.exists():
                missing.append(filename)
            elif path.read_bytes() != expected.encode('utf-8'):
                stale.append(filename)
        if directory.is_symlink():
            raise ValueError(f'Ссылка вместо каталога заметок: {directory}')
        if args.check:
            if missing or stale:
                print('FAIL/1 — заметки классов не соответствуют ontology.yaml')
                for name in missing:
                    print(f'- отсутствует: ontology/{name}')
                for name in stale:
                    print(f'- устарела: ontology/{name}')
                return 1
            print(f'PASS/0 — заметки классов соответствуют ontology.yaml ({len(rendered)})')
            return 0
        directory.mkdir(exist_ok=True)
        created, updated = [], []
        for filename in missing:
            with (directory / filename).open('x', encoding='utf-8', newline='\n') as stream:
                stream.write(rendered[filename])
            created.append('ontology/' + filename)
        if args.update:
            for filename in stale:
                with (directory / filename).open('w', encoding='utf-8', newline='\n') as stream:
                    stream.write(rendered[filename])
                updated.append('ontology/' + filename)
    except (OSError, UnicodeError, ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
        print(f'FAIL/1 — построение заметок классов: {exc}')
        return 1
    print('Created:', ', '.join(created) or 'none')
    if args.update:
        print('Updated:', ', '.join(updated) or 'none')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

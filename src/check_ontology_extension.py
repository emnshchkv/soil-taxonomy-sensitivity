"""Read-only validation of the article ontology extension, not the whole historical graph."""
from qiime_runtime import selected_python
selected_python()

import argparse
import copy
import json
from pathlib import Path
import re
import yaml


RULE_IDS = {'Т1', 'Т2', 'Т3', 'Т4', 'Т5'}
BASELINE = 'results/ontology_extension/20261118/ontology_before.yaml'


def as_list(value):
    return value if isinstance(value, list) else [value]


def ancestors(name, classes):
    chain = set()
    while name is not None:
        if name not in classes or name in chain:
            raise ValueError(f'Неизвестный класс или цикл родителей: {name}')
        chain.add(name)
        name = classes[name].get('parent')
    return chain


def evaluate(schema, nodes, rules):
    """Nodes contain resolved target classes, so the same checks work on synthetic input."""
    errors, counts = [], {}
    classes, properties = schema['classes'], schema['properties']
    isa = lambda child, parent: parent in ancestors(child, classes)
    for name, node in nodes.items():
        ancestors(node['type'], classes)
        for relation, target_type, target_name in node['links']:
            prop = properties.get(relation)
            if prop is None:
                errors.append(f'{name}: неизвестное отношение {relation}')
                continue
            if not any(isa(node['type'], c) for c in as_list(prop['domain'])):
                errors.append(f'{name}: недопустимый domain отношения {relation}')
            if not any(isa(target_type, c) for c in as_list(prop['range'])):
                errors.append(f'{name}: недопустимый range {relation} → {target_name} ({target_type})')
    for rule in rules:
        check, identifier = rule['check'], rule['id']
        relevant = [(name, node) for name, node in nodes.items() if isa(node['type'], check['class'])]
        counts[identifier] = len(relevant)
        for name, node in relevant:
            if not any(rel == check['property'] and isa(target_type, check['target'])
                       for rel, target_type, _ in node['links']):
                errors.append(f'{identifier}: {name} не имеет {check["property"]} → {check["target"]}')
    return errors, counts


def read_note(path):
    text = path.read_text(encoding='utf-8')
    match = re.match(r'\A---\s*\n(.*?)\n---(?:\s*\n|$)', text, re.S)
    metadata = yaml.safe_load(match[1]) if match else {}
    return (metadata or {}).get('type'), text


def verify(root):
    schema = yaml.safe_load((root / 'ontology.yaml').read_text(encoding='utf-8'))
    old = yaml.safe_load((root / BASELINE).read_text(encoding='utf-8'))
    classes, properties = schema['classes'], schema['properties']
    extension = schema['extensions']['article_preparation']
    errors = []
    if len(classes) < 5 or len(schema['rules']) < 3:
        errors.append('Схема должна содержать не меньше пяти классов и трёх правил')
    for name in classes:
        ancestors(name, classes)
    for name, prop in properties.items():
        for field in ('domain', 'range'):
            if not as_list(prop[field]) or any(c not in classes for c in as_list(prop[field])):
                errors.append(f'{name}: неизвестный или пустой {field}')
    for name, definition in old['classes'].items():
        if classes.get(name) != definition:
            errors.append(f'Изменён исходный класс {name}')
    for name, definition in old['attributes'].items():
        if schema['attributes'].get(name) != definition:
            errors.append(f'Изменён исходный атрибут {name}')
    for name, definition in old['properties'].items():
        current = copy.deepcopy(properties.get(name, {}))
        if name == 'хранитсяВ':
            if not set(as_list(definition['domain'])) <= set(as_list(current.get('domain', []))):
                errors.append('Удалена исходная область хранитсяВ')
            current['domain'] = definition['domain']
        if current != definition:
            errors.append(f'Изменено исходное отношение {name}')
    rule_map = {r['id']: r for r in schema['rules']}
    if len(rule_map) != len(schema['rules']):
        errors.append('Повторяются идентификаторы правил')
    for rule in old['rules']:
        if rule_map.get(rule['id']) != rule:
            errors.append(f'Изменено исходное правило {rule["id"]}')
    if set(extension['rules']) != RULE_IDS or len(extension['rules']) != len(RULE_IDS):
        errors.append('Расширение должно декларировать Т1–Т5 ровно по одному разу')
    rules = [rule_map[r] for r in extension['rules']]
    for rule in rules:
        check = rule['check']
        prop = properties.get(check['property'])
        if (check['kind'] != 'some' or check['class'] not in classes
                or check['target'] not in classes or prop is None):
            errors.append(f'{rule["id"]}: некорректное объявление правила some')
            continue
        if (not any(c in ancestors(check['class'], classes) for c in as_list(prop['domain']))
                or not any(c in ancestors(check['target'], classes) for c in as_list(prop['range']))):
            errors.append(f'{rule["id"]}: правило несовместимо с domain/range')
    extension_classes = extension['classes']
    if len(set(extension_classes)) != len(extension_classes) or any(c not in classes for c in extension_classes):
        errors.append('Некорректный список классов расширения')
    if errors:
        return {'status': 'FAIL', 'errors': errors}

    def resolve(target):
        candidates = [root / target, root / (target + '.md')]
        candidates += [root / folder / (target + '.md') for folder in ('notes', 'tasks', 'ontology')]
        found = {p.resolve() for p in candidates if p.exists()}
        if len(found) != 1:
            raise ValueError(f'Не существует или неоднозначна ссылка: {target}')
        path = found.pop()
        if not path.is_relative_to(root):
            raise ValueError(f'Ссылка за пределами проекта: {target}')
        if path.suffix == '.md' and path.parent.name in ('notes', 'tasks', 'ontology'):
            note_type, _ = read_note(path)
            if note_type not in classes:
                raise ValueError(f'Цель не является экземпляром известного класса: {target}')
            return note_type
        return 'Файл'  # Includes the established data/raw/ directory convention.

    selected = {root / 'notes' / (n + '.md') for n in extension['instance_notes'] + extension.get('context_notes', [])}
    for path in (root / 'notes').glob('*.md'):
        note_type, _ = read_note(path)
        if note_type in classes and any(c in ancestors(note_type, classes) for c in extension_classes):
            selected.add(path)
    nodes = {}
    for path in sorted(selected):
        name = str(path.relative_to(root))
        note_type, text = read_note(path)
        if note_type not in classes:
            errors.append(f'{name}: неизвестный type {note_type}')
            continue
        section = re.search(r'^## Связи\s*\n(.*?)(?=^## |\Z)', text, re.S | re.M)
        if not section:
            errors.append(f'{name}: отсутствует раздел Связи')
            continue
        links, types = [], []
        for line in section[1].splitlines():
            if not line.strip():
                continue
            match = re.fullmatch(r'- ([^:]+): \[\[([^\]]+)\]\]', line)
            if not match:
                errors.append(f'{name}: неверная связь {line}')
                continue
            relation, target = match.groups()
            if relation == 'тип':
                types.append(target)
            else:
                try:
                    links.append((relation, resolve(target), target))
                except ValueError as exc:
                    errors.append(f'{name}: {exc}')
        if types != [note_type] or not (root / 'ontology' / (note_type + '.md')).is_file():
            errors.append(f'{name}: связь тип не соответствует type или отсутствует заметка класса')
        nodes[name] = {'type': note_type, 'links': links}
    violations, counts = evaluate(schema, nodes, rules)
    errors.extend(violations)
    instance_counts = {c: sum(c in ancestors(n['type'], classes) for n in nodes.values()) for c in extension_classes}
    return {'status': 'FAIL' if errors else 'PASS', 'errors': errors,
            'scope': 'schema extension, declared/context notes and all notes of extension classes; not full historical graph',
            'instances_by_class': instance_counts, 'instances_checked_per_rule': counts,
            'notes_checked': len(nodes), 'original_classes_preserved': len(old['classes']),
            'original_rules_preserved': len(old['rules'])}


def self_test():
    classes = {c: {'parent': None} for c in ('Файл', 'ASV', 'МодельКлассификатора', 'Образец', 'Результат')}
    classes.update({c: {'parent': 'Результат'} for c in ('ТаксономическоеНазначение', 'ТаксономическийПрофиль')})
    specifications = [('МодельКлассификатора', 'хранитсяВ', 'Файл'), ('ASV', 'хранитсяВ', 'Файл'),
                      ('ТаксономическоеНазначение', 'назначеноASV', 'ASV'),
                      ('ТаксономическоеНазначение', 'полученоМоделью', 'МодельКлассификатора'),
                      ('ТаксономическийПрофиль', 'характеризует', 'Образец')]
    properties = {'хранитсяВ': {'domain': ['МодельКлассификатора', 'ASV'], 'range': 'Файл'},
                  'назначеноASV': {'domain': 'ТаксономическоеНазначение', 'range': 'ASV'},
                  'полученоМоделью': {'domain': 'Результат', 'range': 'МодельКлассификатора'},
                  'характеризует': {'domain': 'Результат', 'range': 'Образец'}}
    rules = [{'id': f'Т{i}', 'check': {'kind': 'some', 'class': c, 'property': p, 'target': t}}
             for i, (c, p, t) in enumerate(specifications, 1)]
    nodes = {c: {'type': c, 'links': []} for c, _, _ in specifications}
    for c, p, t in specifications:
        nodes[c]['links'].append((p, t, 'synthetic target'))
    schema = {'classes': classes, 'properties': properties}
    assert not evaluate(schema, nodes, rules)[0], 'correct graph must pass'
    missing = copy.deepcopy(nodes)
    missing['ТаксономическоеНазначение']['links'] = [x for x in missing['ТаксономическоеНазначение']['links'] if x[0] != 'полученоМоделью']
    assert any(e.startswith('Т4:') for e in evaluate(schema, missing, rules)[0]), 'missing model must fail Т4'
    wrong = copy.deepcopy(nodes)
    wrong['ТаксономическийПрофиль']['links'] = [('характеризует', 'Файл', 'wrong target')]
    faults = evaluate(schema, wrong, rules)[0]
    assert any('range' in e for e in faults) and any(e.startswith('Т5:') for e in faults), 'wrong target must fail range and Т5'
    print('self-test: PASS/0 — clean graph, missing model, wrong profile target; in-memory only')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        result = verify(args.root.resolve())
    except (OSError, UnicodeError, ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
        result = {'status': 'FAIL', 'errors': [str(exc)]}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())

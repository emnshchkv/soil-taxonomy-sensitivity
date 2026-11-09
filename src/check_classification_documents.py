"""Check current report values and links against saved machine evidence."""
from qiime_runtime import selected_python
selected_python()

import argparse
import json
import re
from datetime import datetime
from pathlib import Path
import yaml
from classification_logging import write_json

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / 'results/full_classification/repair_20261109'
NOTES = ['Регрессии проверки классификации', 'Валидация готовых таксономий',
         'Последовательности ASV первого запуска DADA2', 'Проверка договоров классификации',
         'Полная классификация ASV 8 ноября', 'Таксономии трёх моделей',
         'Сверка готовности классификаторов', 'Доказательства готовности классификаторов']

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    errors = []
    def check(condition, message):
        if not condition:
            errors.append(message)
    summary = json.loads((BASE / 'full_verification_retry/summary.json').read_text())
    models = json.loads((BASE / 'model_provenance.json').read_text())
    report_path = ROOT / 'reports/Готовность к классификации.md'
    report = report_path.read_text()
    table = (ROOT / 'results/classifier_readiness/classifier_check_table.md').read_text()
    versions = (ROOT / 'data/refs/VERSIONS.md').read_text()
    check(summary['status'] == 'PASS' and len(summary['results']) == 3, 'full verification')
    elapsed = []
    for row in summary['results']:
        check(row['rows'] == 3374 and row['status'] == 'PASS', 'rows/status')
        check(row['files'][2]['path'] in report, 'result path missing in report')
        runtime = row['execution']['runtime']
        seconds = (datetime.fromisoformat(runtime['end']) -
                   datetime.fromisoformat(runtime['start'])).total_seconds()
        elapsed.append(seconds)
        check(f'{seconds:.6f}' in report, 'duration missing')
        for action in ('validate', 'export'):
            record = json.loads((BASE / 'full_verification_retry' /
                                 f"{row['method']}_{action}.command.json").read_text())
            check(record['exit_code'] == 0, 'CLI failure')
            check((BASE / 'full_verification_retry' / record['stdout']).is_file(), 'stdout missing')
            check((BASE / 'full_verification_retry' / record['stderr']).is_file(), 'stderr missing')
    check(f'{sum(elapsed):.6f}' in report, 'total duration')
    for model in models:
        for text in (report, table):
            check(model['uuid'] in text and model['sha256'] in text, 'model identity missing')
        check(model['sha256'] in versions, 'VERSIONS checksum')
    check('1505' in report and '60078' in report and '2.51%' in report, 'Kr5 limitation')
    check('INCOMPLETE/2' in report and 'не сохранены' in report, 'limitations missing')
    for command in ('project_check', 'regressions', 'preservation'):
        check(json.loads((BASE / f'{command}.command.json').read_text())['exit_code'] == 0,
              command + ' failed')
    schema = yaml.safe_load((ROOT / 'ontology.yaml').read_text())
    properties = schema['properties']
    for name in schema['classes']:
        check((ROOT / 'ontology' / (name + '.md')).is_file(), 'missing class note: ' + name)
    pages = [ROOT / 'notes' / (name + '.md') for name in NOTES] + [report_path]
    all_pages = [p for folder in ('notes', 'tasks', 'ontology', 'reports')
                 for p in (ROOT / folder).glob('*.md')]
    def resolve(target):
        direct = ROOT / target
        if direct.is_file():
            return direct
        if direct.with_suffix('.md').is_file():
            return direct.with_suffix('.md')
        matches = [p for p in all_pages if p.stem == target]
        return matches[0] if len(matches) == 1 else None
    def kind(path):
        match = re.search(r'^type: (.+)$', path.read_text(), re.M)
        return match.group(1) if match else None
    for page in pages:
        source_kind = kind(page)
        for relation, target in re.findall(r'^- ([^:\n]+): \[\[([^\]]+)\]\]', page.read_text(), re.M):
            destination = resolve(target)
            check(destination is not None, f'{page.name}: broken/ambiguous {target}')
            if relation == 'тип':
                check(target == source_kind, f'{page.name}: wrong type link')
                continue
            check(relation in properties, f'{page.name}: unknown relation {relation}')
            if relation not in properties or destination is None:
                continue
            rule = properties[relation]
            domain = rule['domain'] if isinstance(rule['domain'], list) else [rule['domain']]
            check(source_kind in domain, f'{page.name}: wrong domain {relation}')
            if rule['range'] != 'Файл':
                check(kind(destination) == rule['range'], f'{page.name}: wrong range {relation}')
    write_json(args.output, {'passed': not errors, 'errors': errors,
                            'models': len(models), 'graph_pages': len(pages),
                            'total_classification_seconds': sum(elapsed)})
    for error in errors:
        print('FAIL:', error)
    print(f'documents: {"PASS" if not errors else "FAIL"}/{int(bool(errors))}')
    return int(bool(errors))

if __name__ == '__main__':
    raise SystemExit(main())

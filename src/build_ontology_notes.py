"""Create missing class notes from ontology.yaml; preserve existing notes."""
from qiime_runtime import selected_python
selected_python()

from pathlib import Path
import yaml

def main():
    root = Path(__file__).resolve().parent.parent
    schema = yaml.safe_load((root / 'ontology.yaml').read_text())
    directory = root / 'ontology'
    directory.mkdir(exist_ok=True)
    created = []
    for name, item in schema['classes'].items():
        path = directory / (name + '.md')
        if path.exists():
            continue
        relations, attributes = [], []
        for relation, rule in schema['properties'].items():
            domain = rule['domain'] if isinstance(rule['domain'], list) else [rule['domain']]
            if name in domain:
                relations.append(relation)
        for attribute, rule in schema['attributes'].items():
            domain = rule['domain'] if isinstance(rule['domain'], list) else [rule['domain']]
            if name in domain:
                attributes.append(attribute + ' (' + ', '.join(rule['values']) + ')')
        parent = '- подклассОт: [[' + item['parent'] + ']]\n' if item.get('parent') else ''
        text = (f"---\nkind: класс\n---\n# {name}\n\n{item['definition']}\n\n"
                f"Отношения от этого класса: {', '.join(relations) or 'нет'}.\n"
                f"Атрибуты: {'; '.join(attributes) or 'нет'}.\n\n## Связи\n{parent}")
        with path.open('x') as f:
            f.write(text)
        created.append(str(path.relative_to(root)))
    print('Created:', ', '.join(created) or 'none')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

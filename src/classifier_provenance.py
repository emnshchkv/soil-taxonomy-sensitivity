"""Extract small provenance records from existing models, without loading or exporting them."""
from qiime_runtime import selected_python
selected_python()

import argparse
import json
import zipfile
from pathlib import Path
import yaml
from classification_contract import Qza, MODELS, ProvenanceLoader
from classification_logging import sha256, write_json

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path, required=True, help='New JSON file')
    args = parser.parse_args()
    root = args.root.resolve()
    result = []
    for method, filename in MODELS.items():
        artifact = Qza(root / 'data/refs' / filename, 'TaxonomicClassifier')
        provenance = artifact.provenance()
        sources = []
        with zipfile.ZipFile(artifact.path) as z:
            for member in z.namelist():
                if not member.endswith('/action/action.yaml'):
                    continue
                action = yaml.load(z.read(member), Loader=ProvenanceLoader)['action']
                if action.get('type') == 'import' or action.get('action') in (
                        'get_silva_data', 'filter_seqs_length_by_taxon',
                        'reverse_transcribe', 'cull_seqs', 'dereplicate'):
                    sources.append({'member': member, 'action': action})
        result.append({
            'method': method, 'path': str(artifact.path.relative_to(root)),
            'uuid': artifact.uuid, 'sha256': sha256(artifact.path),
            'training_action': provenance['action'],
            'framework_version': provenance['environment']['framework']['version'],
            'feature_classifier_version': provenance['environment']['plugins']['feature-classifier']['version'],
            'sklearn_version_file': json.loads(artifact.text('data/sklearn_version.json')),
            'source_records': sources,
            'interpretation_limit': 'Имена импортов и фильтры не подтверждают подготовку до импорта, регион или праймеры.',
        })
    write_json(args.output, result)
    print(f'model provenance: PASS/0; {args.output}')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

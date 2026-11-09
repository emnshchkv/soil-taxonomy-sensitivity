"""Read-only contracts for classification artifacts; never loads classifier models."""
import csv
import math
import re
import uuid
import zipfile
from pathlib import Path
import yaml

MODELS = {
    'silva138': 'silva-138-99-nb-classifier.qza',
    'gtdb': 'gtdb_classifier_r220.qza',
    'gg2': 'gg_2024.09.full-length.nb.sklearn-1.4.2.qza',
}

class Violation(ValueError):
    pass

def require(condition, code, detail):
    if not condition:
        raise Violation(f'{code}: {detail}')

class ProvenanceLoader(yaml.SafeLoader):
    """QIIME tags are data, never Python objects."""

def tagged(loader, tag, node):
    if isinstance(node, yaml.ScalarNode):
        return loader.construct_scalar(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    return loader.construct_mapping(node)

ProvenanceLoader.add_multi_constructor('!', tagged)

class Qza:
    def __init__(self, path, expected_type):
        self.path = Path(path)
        with zipfile.ZipFile(self.path) as z:
            names = z.namelist()
            require(len(names) == len(set(names)), 'ZIP_DUPLICATE', str(path))
            roots = [n for n in names if re.fullmatch(r'[^/]+/metadata.yaml', n)]
            require(len(roots) == 1, 'QZA_ROOT', str(path))
            self.root = roots[0].split('/')[0]
            self.metadata = yaml.safe_load(z.read(roots[0]))
        self.uuid = self.metadata.get('uuid')
        require(self.uuid == self.root and str(uuid.UUID(self.uuid)) == self.uuid,
                'QZA_UUID', str(path))
        require(self.metadata.get('type') == expected_type, 'QZA_TYPE', str(path))

    def text(self, member):
        with zipfile.ZipFile(self.path) as z:
            return z.read(f'{self.root}/{member}').decode('utf-8')

    def provenance(self):
        return yaml.load(self.text('provenance/action/action.yaml'), Loader=ProvenanceLoader)

def rows(text, required):
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith('#')]
    parsed = list(csv.reader(lines, delimiter='\t'))
    require(bool(parsed), 'TABLE_EMPTY', 'missing header')
    header = parsed[0]
    require(len(header) == len(set(header)), 'TABLE_HEADER', 'duplicate column')
    require(set(required) <= set(header), 'TABLE_HEADER', str(required))
    result = []
    for row in parsed[1:]:
        require(len(row) == len(header), 'TABLE_ROW', str(row[:1]))
        result.append(dict(zip(header, row)))
    return result

def unique(records, key):
    result = {}
    for row in records:
        name = row[key]
        require(bool(name.strip()), 'EMPTY_ID', key)
        require(name not in result, 'DUPLICATE_ID', name)
        result[name] = row
    return result

def sequences(artifact):
    result, current = {}, None
    for line in artifact.text('data/dna-sequences.fasta').splitlines():
        if not line.strip():
            continue
        if line.startswith('>'):
            tokens = line[1:].split()
            require(bool(tokens), 'EMPTY_ID', 'FASTA')
            current = tokens[0]
            require(current not in result, 'DUPLICATE_ID', current)
            result[current] = ''
        else:
            require(current is not None, 'FASTA_FORMAT', 'sequence before ID')
            result[current] += line.strip().upper()
    for key, seq in result.items():
        require(bool(re.fullmatch('[ACGTRYSWKMBDHVN]+', seq)), 'FASTA_SEQUENCE', key)
    require(bool(result), 'FASTA_EMPTY', str(artifact.path))
    return result

def taxonomy(artifact=None, text=None):
    if artifact is not None:
        text = artifact.text('data/taxonomy.tsv')
    table = unique(rows(text, ['Feature ID', 'Taxon', 'Confidence']), 'Feature ID')
    require(bool(table), 'TABLE_EMPTY', 'taxonomy')
    result = {}
    for key, row in table.items():
        value = float(row['Confidence'])
        require(math.isfinite(value) and 0 <= value <= 1, 'CONFIDENCE', key)
        require(bool(row['Taxon'].strip()), 'TAXON_EMPTY', key)
        result[key] = (row['Taxon'], value)
    return result

def flattened(items):
    require(isinstance(items, list), 'PROVENANCE_FORMAT', 'expected list')
    result = {}
    for item in items:
        require(isinstance(item, dict), 'PROVENANCE_FORMAT', 'expected mapping')
        for key, value in item.items():
            require(key not in result, 'PROVENANCE_DUPLICATE', key)
            result[key] = value
    return result

def check_action(artifact, reads, model):
    provenance = artifact.provenance()
    action = provenance['action']
    require(action.get('type') == 'method' and action.get('action') == 'classify_sklearn'
            and action.get('plugin') == 'environment:plugins:feature-classifier',
            'ACTION', str(artifact.path))
    inputs = flattened(action['inputs'])
    parameters = flattened(action['parameters'])
    require(inputs.get('reads') == reads.uuid, 'READS_UUID', str(artifact.path))
    require(inputs.get('classifier') == model.uuid, 'CLASSIFIER_UUID', str(artifact.path))
    require(parameters.get('confidence') == 0.7, 'PARAMETER', 'confidence')
    require(parameters.get('n_jobs') == 1, 'PARAMETER', 'n_jobs')
    return provenance

def check_stats(root):
    metadata = unique(rows((root / 'data/meta/metadata.tsv').read_text(),
                           ['sample-id', 'sample_type']), 'sample-id')
    require(len(metadata) == 11, 'SAMPLE_COUNT', 'metadata')
    expected = set(metadata)
    columns = ['input', 'filtered', 'denoised', 'non-chimeric']
    saved = rows((root / 'results/exact_match_reproduction/provenance/samples.tsv').read_text(),
                 ['run', 'sample-id', *columns])
    saved = unique([r for r in saved if r['run'] == '1'], 'sample-id')
    artifact = Qza(root / 'work/qza/dada2_stats.qza', 'SampleData[DADA2Stats]')
    actual = unique(rows(artifact.text('data/stats.tsv'), ['sample-id', *columns]), 'sample-id')
    require(set(saved) == set(actual) == expected, 'SAMPLE_SET', 'expected all 11 samples')
    for sid in expected:
        for column in columns:
            a, b = saved[sid][column], actual[sid][column]
            require(bool(re.fullmatch(r'\d+', a)) and bool(re.fullmatch(r'\d+', b)),
                    'COUNT_DOMAIN', f'{sid}/{column}')
            require(int(a) == int(b), 'STATS_VALUE', f'{sid}/{column}')
    return len(expected)

def read_inputs(root):
    rep = Qza(root / 'work/qza/rep_seq.qza', 'FeatureData[Sequence]')
    seqs = sequences(rep)
    require(len(seqs) == 3374, 'ASV_COUNT', str(len(seqs)))
    models = {name: Qza(root / 'data/refs' / file, 'TaxonomicClassifier')
              for name, file in MODELS.items()}
    return rep, seqs, models

def readiness(root):
    count = check_stats(root)
    rep, seqs, models = read_inputs(root)
    probe = Qza(root / 'work/classifier_probe/probe5_seqs.qza', 'FeatureData[Sequence]')
    probe_seqs = sequences(probe)
    selected = sorted(seqs)[:5]
    require(set(probe_seqs) == set(selected), 'PROBE_ID_SET', str(probe.path))
    require(all(probe_seqs[key] == seqs[key] for key in selected),
            'PROBE_SEQUENCE', str(probe.path))
    results = []
    for method, model in models.items():
        tax = Qza(root / f'work/classifier_probe/taxonomy_probe_{method}.qza',
                  'FeatureData[Taxonomy]')
        table = taxonomy(tax)
        require(set(table) == set(selected), 'TAXONOMY_ID_SET', method)
        check_action(tax, probe, model)
        results.append({'method': method, 'rows': len(table), 'uuid': tax.uuid})
    return {'samples': count, 'asv': len(seqs), 'probe_ids': selected, 'models': results}

def full_contract(root):
    rep, seqs, models = read_inputs(root)
    result = []
    for method, model in models.items():
        tax = Qza(root / f'work/qza/taxonomy_sklearn_{method}.qza', 'FeatureData[Taxonomy]')
        table = taxonomy(tax)
        require(set(table) == set(seqs), 'TAXONOMY_ID_SET', method)
        provenance = check_action(tax, rep, model)
        result.append((method, tax, table, provenance, rep, model))
    return result

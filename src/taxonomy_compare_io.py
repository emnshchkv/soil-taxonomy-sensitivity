"""Strict input, output schema and source identity for the comparison.

No computation of biological metrics lives in this module.
"""
from pathlib import Path
import hashlib
import io
import json
import math
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parent.parent
MODELS = ('sklearn_silva138', 'sklearn_gtdb', 'sklearn_gg2')
RANKS = ('тип', 'класс', 'порядок', 'семейство', 'род')
PAIRS = ((MODELS[0], MODELS[1]), (MODELS[0], MODELS[2]), (MODELS[1], MODELS[2]))
CATEGORIES = ('both_match', 'both_diff', 'only_a', 'only_b', 'neither')
SAMPLES = ('K_S63', 'Kr1_S80', 'Kr2_S81', 'Kr3_S82', 'Kr4_S83', 'Kr5_S84',
           'LUB2_S92', 'LUB3_S93', 'LUB4_S94', 'LUB5_S95', 'LUB6_S96')
INPUT_FILES = {'ASV_table': 'work/qza/ASV_table.qza', 'rep_seq': 'work/qza/rep_seq.qza',
               'metadata': 'data/meta/metadata.tsv', 'VERSIONS': 'data/refs/VERSIONS.md',
               **{m: f'work/qza/taxonomy_{m}.qza' for m in MODELS}}
CODE_FILES = ('src/complete_taxonomy_comparison.py', 'src/repro_defects.py', 'src/taxonomy_compare_io.py', 'src/taxonomy_compare_core.py',
              'src/taxonomy_comparison.py', 'src/check_taxonomy_comparison.py',
              'src/taxonomy_compare_acceptance.py', 'src/qiime_runtime.py',
              'src/classification_contract.py', 'src/classification_logging.py',
              'tests/test_taxonomy_comparison.py', 'tests/taxonomy_comparison_fixtures.py')
SCHEMAS = {
    'assignments.tsv': ['feature_id', 'soil_asv'] + [c for m in MODELS for c in
        [f'{m}_taxon', *[c for r in RANKS for c in
            [f'{m}_{r}', f'{m}_{r}_key', f'{m}_{r}_assigned']]]],
    'coverage.tsv': ['model', 'rank', 'assigned_asv', 'unassigned_asv', 'total_asv',
                     'assigned_asv_frac', 'unassigned_asv_frac'],
    'coverage_per_sample.tsv': ['model', 'sample', 'rank', 'assigned_reads', 'unassigned_reads',
                               'total_reads', 'assigned_reads_frac', 'unassigned_reads_frac'],
    'pairwise_match.tsv': ['model_a', 'model_b', 'rank', *CATEGORIES, 'total_asv',
                          *[f'{c}_frac' for c in CATEGORIES], 'both_assigned_asv', 'cond_match_asv_frac'],
    'pairwise_match_per_sample.tsv': ['model_a', 'model_b', 'sample', 'rank',
        *[f'{c}_reads' for c in CATEGORIES], 'total_reads',
        *[f'{c}_reads_frac' for c in CATEGORIES], 'both_assigned_reads', 'cond_match_reads_frac'],
    'profiles.tsv': ['model', 'sample', 'rank', 'key', 'count', 'total_reads', 'frac'],
    'distances.tsv': ['model_a', 'model_b', 'sample', 'rank', 'distance',
                      'assigned_reads_frac_a', 'assigned_reads_frac_b'],
    'summary.tsv': ['model_a', 'model_b', 'rank', 'n_samples', 'min', 'median', 'max'],
    'manifest.tsv': ['item', 'path', 'sha256', 'uuid', 'type', 'code_version'],
}
RULES = ('Literal rank-prefix normalization; trim surrounding whitespace; preserve suffixes.\n'
         'Service labels: empty, Unassigned, Unknown_Family.\n'
         'Lineage keys retain known ancestors through the target rank; unknown target = UNCLASSIFIED.\n'
         'No synonyms, filtering, rarefaction or pooling of sample read counts.\n'
         'ASV denominator: positive sum in ten soils; sample denominator: all reads in that soil.\n'
         'Conditional fractions use separate ASV and read denominators; zero denominator = NA.\n'
         'Distance = half L1 over the union of lineage keys including UNCLASSIFIED.\n')

def issue(code, detail, missing=False):
    return {'code': code, 'detail': str(detail), 'missing': missing}

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def _git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, check=False)

def code_identity():
    hashes = {p: sha256(ROOT / p) for p in CODE_FILES}
    commit = _git('rev-parse', 'HEAD').stdout.decode().strip()
    clean = bool(commit) and all(_git('show', f'{commit}:{p}').returncode == 0 and
        hashlib.sha256(_git('show', f'{commit}:{p}').stdout).hexdigest() == hashes[p] for p in CODE_FILES)
    snapshot = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return {'version': commit if clean else 'snapshot-' + snapshot, 'commit': commit,
            'clean': clean, 'files': hashes}

def verify_code_identity(record):
    errors = []
    try:
        actual = code_identity()
        if record['files'] != actual['files']:
            errors.append(issue('CODE_CONTENT', 'recorded source hashes differ from executing source'))
        if record['clean']:
            commit = record['commit']
            if record['version'] != commit or len(commit) != 40:
                errors.append(issue('CODE_VERSION', 'expected full commit SHA'))
            for p in CODE_FILES:
                got = _git('show', f'{commit}:{p}')
                if got.returncode or hashlib.sha256(got.stdout).hexdigest() != record['files'][p]:
                    errors.append(issue('CODE_COMMIT', p))
        else:
            expected = 'snapshot-' + hashlib.sha256(json.dumps(record['files'], sort_keys=True).encode()).hexdigest()
            if record['version'] != expected:
                errors.append(issue('CODE_SNAPSHOT', 'snapshot digest mismatch'))
    except (KeyError, TypeError, OSError, ValueError) as e:
        errors.append(issue('CODE_FORMAT', e))
    return errors

def load_dataset(root):
    """Return dataset and all observed violations, without mutating either source."""
    import biom
    import h5py
    from classification_contract import Qza, rows, unique, sequences, taxonomy, require
    root = Path(root)
    d = {'table': {}, 'samples': [], 'soil_samples': [], 'sequences': {},
         'taxonomies': {}, 'inputs': {}}
    errors = []
    for item, rel in INPUT_FILES.items():
        path = root / rel
        if not path.is_file():
            errors.append(issue('INPUT_MISSING', rel, True))
            continue
        try:
            rec = {'path': rel, 'sha256': sha256(path), 'uuid': '', 'type': 'text'}
            if item == 'metadata':
                meta = unique(rows(path.read_text(), ['sample-id', 'sample_type']), 'sample-id')
                require(set(meta) == set(SAMPLES), 'SAMPLE_SET', 'metadata must contain the fixed 11 IDs')
                for sid, row in meta.items():
                    require(row['sample_type'] == ('negative_control' if sid == 'K_S63' else 'soil'),
                            'SAMPLE_ROLE', sid)
                d['soil_samples'] = sorted(s for s in meta if s != 'K_S63')
            elif item == 'VERSIONS':
                require(bool(path.read_text().strip()), 'VERSIONS_EMPTY', rel)
            else:
                typ = 'FeatureTable[Frequency]' if item == 'ASV_table' else (
                    'FeatureData[Sequence]' if item == 'rep_seq' else 'FeatureData[Taxonomy]')
                q = Qza(path, typ)
                rec.update(uuid=q.uuid, type=typ)
                with zipfile.ZipFile(path) as z:
                    require(all(n.startswith(q.root + '/') for n in z.namelist()), 'QZA_ROOT', rel)
                    if item == 'ASV_table':
                        with h5py.File(io.BytesIO(z.read(q.root + '/data/feature-table.biom')), 'r') as f:
                            table = biom.Table.from_hdf5(f)
                        fids = list(table.ids(axis='observation'))
                        sids = list(table.ids(axis='sample'))
                        require(len(fids) == 3374 and len(set(fids)) == 3374, 'ASV_COUNT', len(fids))
                        require(len(sids) == 11 and set(sids) == set(SAMPLES), 'SAMPLE_SET', sids)
                        out = {}
                        for fid, vals in zip(fids, table.iter_data(axis='observation')):
                            row = {}
                            for sid, value in zip(sids, vals):
                                value = float(value)
                                require(math.isfinite(value) and value >= 0 and value.is_integer(),
                                        'COUNT_DOMAIN', f'{fid}/{sid}/{value}')
                                row[sid] = int(value)
                            require(sum(row.values()) > 0, 'ASV_EMPTY', fid)
                            out[fid] = row
                        for sid in sids:
                            require(sum(r[sid] for r in out.values()) > 0, 'SAMPLE_EMPTY', sid)
                        d['table'], d['samples'] = out, sorted(sids)
                    elif item == 'rep_seq':
                        d['sequences'] = sequences(q)
                        require(len(d['sequences']) == 3374, 'ASV_COUNT', 'FASTA')
                        require(all(len(s) == 180 for s in d['sequences'].values()), 'SEQUENCE_LENGTH', 'expected 180 nt')
                    else:
                        d['taxonomies'][item] = {k: v[0] for k, v in taxonomy(q).items()}
                        require(len(d['taxonomies'][item]) == 3374, 'ASV_COUNT', item)
            d['inputs'][item] = rec
        except Exception as e:
            errors.append(issue('INPUT_INVALID', f'{rel}: {type(e).__name__}: {e}'))
    sets = [set(v) for v in [d['table'], d['sequences'], *d['taxonomies'].values()] if v]
    if sets and any(s != sets[0] for s in sets):
        errors.append(issue('ASV_SET', 'table, FASTA and taxonomies must have the same 3374 IDs'))
    return d, errors

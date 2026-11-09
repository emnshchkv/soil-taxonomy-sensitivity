"""Reproducible private-input regressions. Only temporary copies are mutated.

The two legacy cases execute isolated functions from 8f9832f (not its full CLI).
All current cases exercise the same contract used by the two verification CLIs.
Actual QIIME validate/export are tested separately against unmodified real QZA.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from qiime_runtime import selected_python
selected_python()

import argparse
import ast
import csv
import io
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from classification_contract import MODELS, readiness, full_contract
from classification_logging import record_command, sha256, write_json

PROJECT = Path(__file__).resolve().parents[1]
CASES = [
    ('baseline_ready', 'ready', 0, 'PASS'),
    ('baseline_full', 'full', 0, 'PASS'),
    ('missing_sample', 'ready', 1, 'SAMPLE_SET'),
    ('duplicate_stats', 'ready', 1, 'DUPLICATE_ID'),
    ('duplicate_saved', 'ready', 1, 'DUPLICATE_ID'),
    ('duplicate_metadata', 'ready', 1, 'DUPLICATE_ID'),
    ('duplicate_rep', 'ready', 1, 'DUPLICATE_ID'),
    ('duplicate_probe', 'ready', 1, 'DUPLICATE_ID'),
    ('changed_sequence', 'ready', 1, 'PROBE_SEQUENCE'),
    ('wrong_probe_ids', 'ready', 1, 'TAXONOMY_ID_SET'),
    ('changed_probe_uuid', 'ready', 1, 'READS_UUID'),
    ('wrong_classifier', 'ready', 1, 'CLASSIFIER_UUID'),
    ('full_duplicate', 'full', 1, 'DUPLICATE_ID'),
    ('full_wrong_ids', 'full', 1, 'TAXONOMY_ID_SET'),
    ('full_nan', 'full', 1, 'CONFIDENCE'),
    ('full_negative', 'full', 1, 'CONFIDENCE'),
    ('full_over_one', 'full', 1, 'CONFIDENCE'),
    ('full_empty_taxon', 'full', 1, 'TAXON_EMPTY'),
    ('full_type', 'full', 1, 'QZA_TYPE'),
    ('full_missing', 'full', 1, 'FileNotFoundError'),
    ('full_changed_rep_uuid', 'full', 1, 'READS_UUID'),
    ('full_changed_model_uuid', 'full', 1, 'CLASSIFIER_UUID'),
    ('full_wrong_parameter', 'full', 1, 'PARAMETER'),
    ('full_wrong_action', 'full', 1, 'ACTION'),
    ('reorder_ready', 'ready', 0, 'PASS'),
    ('reorder_full', 'full', 0, 'PASS'),
    ('confidence_bounds', 'full', 0, 'PASS'),
    ('legacy_duplicate_stats', 'legacy', 0, 'LEGACY_FALSE_PASS'),
    ('legacy_duplicate_probe', 'legacy', 0, 'LEGACY_FALSE_PASS'),
]

def paths():
    small = ['data/meta/metadata.tsv',
             'results/exact_match_reproduction/provenance/samples.tsv',
             'work/qza/dada2_stats.qza', 'work/qza/rep_seq.qza',
             'work/classifier_probe/probe5_seqs.qza']
    for method in MODELS:
        small += [f'work/classifier_probe/taxonomy_probe_{method}.qza',
                  f'work/qza/taxonomy_sklearn_{method}.qza']
    models = ['data/refs/' + file for file in MODELS.values()]
    return small, models

def copy_inputs(source, dest):
    small, models = paths()
    for rel in small:
        p = dest / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / rel, p)
    for rel in models:
        p = dest / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.symlink_to(source / rel)

def edit_qza(path, member, change):
    assert not path.is_symlink(), 'never modify linked originals'
    with zipfile.ZipFile(path) as z:
        entries = [(n, z.read(n)) for n in z.namelist()]
    root = next(n.split('/')[0] for n, _ in entries
                if re.fullmatch(r'[^/]+/metadata.yaml', n))
    target = root + '/' + member
    assert sum(n == target for n, _ in entries) == 1
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, data in entries:
            z.writestr(n, change(data.decode()).encode() if n == target else data)

def change_uuid(path):
    assert not path.is_symlink()
    new = '00000000-0000-4000-8000-000000000001'
    with zipfile.ZipFile(path) as z:
        entries = [(n, z.read(n)) for n in z.namelist()]
    root = next(n.split('/')[0] for n, _ in entries
                if re.fullmatch(r'[^/]+/metadata.yaml', n))
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, data in entries:
            if n == root + '/metadata.yaml':
                data = data.decode().replace(root, new).encode()
            z.writestr(new + n[len(root):], data)

def duplicate_row(text):
    lines = text.splitlines()
    row = next(x for x in lines[1:] if x.strip() and not x.startswith('#'))
    return text.rstrip('\n') + '\n' + row + '\n'

def duplicate_fasta(text):
    return text.rstrip('\n') + '\n>' + text.split('>')[1].rstrip('\n') + '\n'

def tax_value(text, col, value):
    lines = list(csv.reader(text.splitlines(), delimiter='\t'))
    lines[1][col] = value
    out = io.StringIO()
    csv.writer(out, delimiter='\t', lineterminator='\n').writerows(lines)
    return out.getvalue()

def reorder_table(text):
    lines = list(csv.reader(text.splitlines(), delimiter='\t'))
    header, rest = lines[0], lines[1:]
    comments = [r for r in rest if r[0].startswith('#')]
    data = [r for r in rest if not r[0].startswith('#')]
    out = io.StringIO()
    order = lambda row: row[:1] + row[1:][::-1]
    csv.writer(out, delimiter='\t', lineterminator='\n').writerows(
        [order(header)] + [order(r) for r in comments + data[::-1]])
    return out.getvalue()

def reorder_fasta(text):
    return ''.join('>' + record for record in text.split('>')[1:][::-1])

def mutate(root, case):
    stats = root / 'work/qza/dada2_stats.qza'
    rep = root / 'work/qza/rep_seq.qza'
    probe = root / 'work/classifier_probe/probe5_seqs.qza'
    small_tax = root / 'work/classifier_probe/taxonomy_probe_silva138.qza'
    full_tax = root / 'work/qza/taxonomy_sklearn_silva138.qza'
    fasta = 'data/dna-sequences.fasta'
    action = 'provenance/action/action.yaml'
    if case == 'missing_sample':
        def remove(text):
            lines = text.splitlines()
            victim = next(x for x in lines[1:] if not x.startswith('#'))
            return '\n'.join(x for x in lines if x != victim) + '\n'
        edit_qza(stats, 'data/stats.tsv', remove)
    elif case in ('duplicate_stats', 'legacy_duplicate_stats'):
        edit_qza(stats, 'data/stats.tsv', duplicate_row)
    elif case == 'duplicate_saved':
        p = root / 'results/exact_match_reproduction/provenance/samples.tsv'
        p.write_text(duplicate_row(p.read_text()))
    elif case == 'duplicate_metadata':
        p = root / 'data/meta/metadata.tsv'
        p.write_text(duplicate_row(p.read_text()))
    elif case == 'duplicate_rep':
        edit_qza(rep, fasta, duplicate_fasta)
    elif case in ('duplicate_probe', 'legacy_duplicate_probe'):
        edit_qza(probe, fasta, duplicate_fasta)
    elif case == 'changed_sequence':
        def replace_base(text):
            lines = text.splitlines()
            lines[1] = ('C' if lines[1][0] != 'C' else 'A') + lines[1][1:]
            return '\n'.join(lines) + '\n'
        edit_qza(probe, fasta, replace_base)
    elif case == 'wrong_probe_ids':
        edit_qza(small_tax, 'data/taxonomy.tsv', lambda t: tax_value(t, 0, 'different-id'))
    elif case == 'changed_probe_uuid':
        change_uuid(probe)
    elif case == 'wrong_classifier':
        edit_qza(small_tax, action, lambda t: re.sub(
            r'(classifier: )[0-9a-f-]+', r'\g<1>00000000-0000-4000-8000-000000000001', t))
    elif case == 'full_duplicate':
        edit_qza(full_tax, 'data/taxonomy.tsv', duplicate_row)
    elif case == 'full_wrong_ids':
        edit_qza(full_tax, 'data/taxonomy.tsv', lambda t: tax_value(t, 0, 'different-id'))
    elif case in ('full_nan', 'full_negative', 'full_over_one', 'full_empty_taxon'):
        col, value = {'full_nan': (2, 'nan'), 'full_negative': (2, '-0.1'),
                      'full_over_one': (2, '1.1'), 'full_empty_taxon': (1, ' ')}[case]
        edit_qza(full_tax, 'data/taxonomy.tsv', lambda t: tax_value(t, col, value))
    elif case == 'full_type':
        edit_qza(full_tax, 'metadata.yaml', lambda t: t.replace('FeatureData[Taxonomy]', 'FeatureData[Sequence]'))
    elif case == 'full_missing':
        full_tax.rename(full_tax.with_suffix('.missing'))
    elif case == 'full_changed_rep_uuid':
        change_uuid(rep)
    elif case == 'full_changed_model_uuid':
        # A tiny stand-in for the metadata-only model reader, never unpickled.
        p = root / 'data/refs' / MODELS['silva138']
        p.unlink()  # removes only this temporary symlink, never its target
        new = '00000000-0000-4000-8000-000000000001'
        with zipfile.ZipFile(p, 'w') as z:
            z.writestr(new + '/metadata.yaml', f'uuid: {new}\ntype: TaxonomicClassifier\n')
    elif case == 'full_wrong_parameter':
        edit_qza(full_tax, action, lambda t: t.replace('confidence: 0.7', 'confidence: 0.8'))
    elif case == 'full_wrong_action':
        edit_qza(full_tax, action, lambda t: t.replace('action: classify_sklearn', 'action: import'))
    elif case == 'reorder_ready':
        edit_qza(stats, 'data/stats.tsv', reorder_table)
        edit_qza(rep, fasta, reorder_fasta)
        edit_qza(probe, fasta, reorder_fasta)
        for method in MODELS:
            edit_qza(root / f'work/classifier_probe/taxonomy_probe_{method}.qza',
                     'data/taxonomy.tsv', reorder_table)
    elif case == 'reorder_full':
        edit_qza(rep, fasta, reorder_fasta)
        for method in MODELS:
            edit_qza(root / f'work/qza/taxonomy_sklearn_{method}.qza',
                     'data/taxonomy.tsv', reorder_table)
    elif case == 'confidence_bounds':
        edit_qza(full_tax, 'data/taxonomy.tsv', lambda t: tax_value(t, 2, '0'))
        other = root / 'work/qza/taxonomy_sklearn_gtdb.qza'
        edit_qza(other, 'data/taxonomy.tsv', lambda t: tax_value(t, 2, '1'))

def run_case(project, case, scope):
    with tempfile.TemporaryDirectory(prefix='classification-case-') as temp:
        root = Path(temp)
        copy_inputs(project, root)
        mutate(root, case)
        if scope == 'legacy':
            source = subprocess.check_output(
                ['git', 'show', '8f9832f:src/verify_classifier_readiness.py'],
                cwd=project, text=True)
            names = {'_read_qza_data', 'check_dada2_stats', 'check_probe5_ids'}
            module = ast.Module(body=[n for n in ast.parse(source).body
                                     if isinstance(n, ast.FunctionDef) and n.name in names],
                                type_ignores=[])
            ns = {'ROOT': root, 'Path': Path, 'zipfile': zipfile}
            exec(compile(module, '8f9832f:isolated_functions', 'exec'), ns)
            function = 'check_dada2_stats' if case.endswith('stats') else 'check_probe5_ids'
            passed = ns[function]()
            print('LEGACY_FALSE_PASS' if passed else 'legacy detected violation')
            return 0 if passed else 1
        try:
            value = readiness(root) if scope == 'ready' else full_contract(root)
            print(f'{case}: PASS; scope={scope}')
            return 0
        except Exception as error:
            print(f'{case}: FAIL: {type(error).__name__}: {error}')
            return 1

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=PROJECT)
    parser.add_argument('--log-dir', type=Path)
    parser.add_argument('--case', choices=[c[0] for c in CASES])
    args = parser.parse_args()
    project = args.project.resolve()
    if args.case:
        case, scope, _, _ = next(c for c in CASES if c[0] == args.case)
        return run_case(project, case, scope)
    if args.log_dir is None:
        parser.error('--log-dir must be a new directory')
    out = args.log_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    small, models = paths()
    before = {p: sha256(project / p) for p in small + models}
    results = []
    for case, scope, expected, marker in CASES:
        result = record_command(out, case, [sys.executable, '-B', str(Path(__file__).resolve()),
                                '--project', str(project), '--case', case], project)
        stdout = (out / result['stdout']).read_text()
        passed = result['exit_code'] == expected and marker in stdout
        results.append({'case': case, 'scope': scope, 'expected_exit': expected,
                        'actual_exit': result['exit_code'], 'expected_marker': marker,
                        'passed': passed})
    after = {p: sha256(project / p) for p in small + models}
    unchanged = before == after
    write_json(out / 'input_hashes.json', {'before': before, 'after': after, 'unchanged': unchanged})
    success = unchanged and all(row['passed'] for row in results)
    write_json(out / 'summary.json', {'passed': success, 'count': len(results), 'cases': results})
    print(f'regressions: {"PASS" if success else "FAIL"}/{0 if success else 1}; {len(results)} cases')
    return 0 if success else 1

if __name__ == '__main__':
    raise SystemExit(main())

#!/usr/bin/env python3
"""Independent, read-only acceptance of the article-number transfer (not biology).

The expected selectors below are fixed independently of the supplied numbers map.
No generator module is imported. Source arithmetic uses Decimal with its own
median implementation; missing/invalid outputs always produce FAIL/1.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CMP = 'results/taxonomy_comparison/finish_20261116/main/'
FAM = 'results/taxonomy_barplots/prokaryote_family_20261118_final/'
FIG = 'figures/taxonomy/family/'
META = 'data/meta/metadata.tsv'
CHOICE = 'memory/denoise_choice.md'
STATS = 'results/exact_match_reproduction/provenance/samples.tsv'
MODELS = ('sklearn_silva138', 'sklearn_gtdb', 'sklearn_gg2')
TAGS = ('Silva', 'Gtdb', 'GgTwo')
RANKS = ('тип', 'класс', 'порядок', 'семейство', 'род')
PAIRS = ((MODELS[0], MODELS[1]), (MODELS[0], MODELS[2]), (MODELS[1], MODELS[2]))
PAIR_TAGS = ('SilvaGtdb', 'SilvaGgTwo', 'GtdbGgTwo')
CS = 'comparison_all_post_dada2'
FS = 'family_prokaryote_selected'
COV_COLS = ['scope', 'model', 'rank', 'assigned_asv', 'total_asv', 'assigned_asv_pct',
            'read_coverage_min_pct', 'read_coverage_median_pct', 'read_coverage_max_pct', 'n_samples']
PAIR_COLS = ['scope', 'model_a', 'model_b', 'rank', 'distance_min', 'distance_median', 'distance_max',
             'n_samples', 'both_match', 'both_diff', 'only_a', 'only_b', 'neither', 'total_asv',
             'both_assigned_asv', 'conditional_asv_pct', 'conditional_read_min_pct',
             'conditional_read_median_pct', 'conditional_read_max_pct', 'conditional_defined_samples']
FAM_COLS = ['scope', 'model', 'sample', 'sample_type', 'post_dada2_reads', 'selected_reads',
            'excluded_reads', 'selected_bacteria_reads', 'selected_archaea_reads', 'named_family_reads',
            'unclassified_family_reads', 'selected_pct_of_post_dada2', 'unclassified_pct_of_selected', 'zero_selected']
NUM_COLS = ['id', 'scope', 'source_file', 'row_key', 'column', 'operation', 'numerator',
            'denominator', 'n', 'unit', 'raw_value', 'digits', 'display_value']
FIG_COLS = ['scope', 'figure_id', 'model', 'path', 'sha256', 'source_table', 'denominator', 'display_threshold_pct']
OUTPUTS = {'comparison_coverage.tsv', 'comparison_pairs.tsv', 'family_sample_counts.tsv',
           'figure_sources.tsv', 'numbers.tsv', 'numbers.tex', 'comparison_main.tex', 'family_main.tex', 'README.md'}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest() if hasattr(hashlib, 'file_digest') else _digest(stream)


def _digest(stream):
    h = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1 << 20), b''):
        h.update(chunk)
    return h.hexdigest()


def rows(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        result = list(reader)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'{path}: invalid header')
        if any(None in row or None in row.values() for row in result):
            raise ValueError(f'{path}: malformed TSV')
        return reader.fieldnames, result


def js(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def percent(n, d):
    return None if not int(d) else Decimal(int(n)) * 100 / Decimal(int(d))


def median(values):
    values = sorted(values)
    if not values:
        return None
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def extent(values):
    values = [v for v in values if v is not None]
    return (min(values), median(values), max(values)) if values else (None, None, None)


def display(value, digits):
    if value is None:
        return 'NA'
    return format(Decimal(value).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_EVEN), f'.{digits}f')


def index(data, keys):
    answer = {}
    for row in data:
        key = tuple(row[k] for k in keys)
        if key in answer:
            raise ValueError(f'duplicate source key {key}')
        answer[key] = row
    return answer


def expected(root):
    """Reconstruct all values from the fixed source files, never from output selectors."""
    take = lambda path: rows(root / path)[1]
    meta = take(META)
    samples = sorted(r['sample-id'] for r in meta)
    soils = sorted(r['sample-id'] for r in meta if r['sample_type'] == 'soil')
    roles = {r['sample-id']: r['sample_type'] for r in meta}
    wanted_soils = [f'Kr{i}_S{79+i}' for i in range(1, 6)] + [f'LUB{i}_S{90+i}' for i in range(2, 7)]
    if soils != wanted_soils or samples != sorted(wanted_soils + ['K_S63']) or roles.get('K_S63') != 'negative_control':
        raise ValueError('source metadata: unexpected sample identities or roles')
    coverage = index(take(CMP + 'coverage.tsv'), ('model', 'rank'))
    per_coverage = index(take(CMP + 'coverage_per_sample.tsv'), ('model', 'rank', 'sample'))
    pair = index(take(CMP + 'pairwise_match.tsv'), ('model_a', 'model_b', 'rank'))
    per_pair = index(take(CMP + 'pairwise_match_per_sample.tsv'), ('model_a', 'model_b', 'rank', 'sample'))
    distances = index(take(CMP + 'distances.tsv'), ('model_a', 'model_b', 'rank', 'sample'))
    filters = index(take(FAM + 'filter_counts.tsv'), ('model', 'sample'))
    profiles = take(FAM + 'family_profiles.tsv')
    assignments = take(CMP + 'assignments.tsv')
    stats = index(take(STATS), ('run', 'sample-id'))
    fm = json.loads((root / FAM / 'prokaryotic_family_manifest.json').read_text())
    choice_text = (root / CHOICE).read_text()
    choice = {p: int(re.search(r'^' + p + r': (\d+)$', choice_text, re.M).group(1)) for p in ('trim-left', 'trunc-len')}
    cov_rows, pair_rows, fam_rows = [], [], []
    for model in MODELS:
        for rank in RANKS:
            row = coverage[model, rank]
            nums = [int(per_coverage[model, rank, sample]['assigned_reads']) for sample in soils]
            dens = [int(per_coverage[model, rank, sample]['total_reads']) for sample in soils]
            low, med, high = extent([percent(a, b) for a, b in zip(nums, dens)])
            cov_rows.append(dict(zip(COV_COLS, (CS, model, rank, int(row['assigned_asv']), int(row['total_asv']),
                percent(row['assigned_asv'], row['total_asv']), low, med, high, len(soils)))))
    for a, b in PAIRS:
        for rank in RANKS:
            row = pair[a, b, rank]
            d = [Decimal(distances[a, b, rank, s]['distance']) for s in soils]
            conditions = [percent(per_pair[a, b, rank, s]['both_match_reads'], per_pair[a, b, rank, s]['both_assigned_reads']) for s in soils]
            low, med, high = extent(conditions)
            pair_rows.append(dict(zip(PAIR_COLS, (CS, a, b, rank, min(d), median(d), max(d), len(soils),
                *(int(row[k]) for k in ('both_match', 'both_diff', 'only_a', 'only_b', 'neither', 'total_asv', 'both_assigned_asv')),
                percent(row['both_match'], row['both_assigned_asv']), low, med, high, sum(v is not None for v in conditions)))))
    for model in MODELS:
        for sample in samples:
            row = filters[model, sample]
            selected = int(row['selected_reads'])
            p = [r for r in profiles if r['model'] == model and r['sample'] == sample]
            if sum(int(r['count']) for r in p) != selected:
                raise ValueError(f'source profile mass {model}/{sample}')
            unknown = sum(int(r['count']) for r in p if r['key'] == 'UNCLASSIFIED')
            total = int(row['total_reads'])
            if total != int(stats['1', sample]['non-chimeric']):
                raise ValueError(f'source total differs from first DADA2 run: {model}/{sample}')
            fam_rows.append(dict(zip(FAM_COLS, (FS, model, sample, roles[sample], total, selected,
                int(row['excluded_reads']), int(row['selected_Bacteria_reads']), int(row['selected_Archaea_reads']),
                selected - unknown, unknown, percent(selected, total), percent(unknown, selected), int(selected == 0)))))
    context = dict(meta=meta, soils=soils, samples=samples, roles=roles, coverage=coverage, per_coverage=per_coverage,
                   pair=pair, per_pair=per_pair, distances=distances, filters=filters, assignments=assignments,
                   stats=stats, fm=fm, choice=choice, cov_rows=cov_rows, pair_rows=pair_rows, fam_rows=fam_rows)
    return context


class Audit:
    def __init__(self):
        self.issues = []
        self.checks = 0

    def require(self, condition, message):
        self.checks += 1
        if not condition:
            self.issues.append(message)

    def number(self, observed, wanted, label):
        if wanted is None:
            self.require(observed == 'NA', label + ': expected NA')
        elif isinstance(wanted, int):
            self.require(observed == str(wanted), label + ': integer mismatch')
        else:
            try:
                found = Decimal(observed)
                good = found.is_finite() and abs(found - wanted) <= Decimal('5e-12') * max(Decimal(1), abs(wanted))
            except Exception:
                good = False
            self.require(good, label + ': numeric mismatch')

    def table(self, directory, name, columns, data, keys):
        header, actual = rows(directory / name)
        self.require(header == columns, name + ': wrong columns/order')
        self.require(len(actual) == len(data), name + ': wrong row count')
        observed = index(actual, keys)
        wanted = index(data, keys)
        self.require(set(observed) == set(wanted), name + ': wrong keys/composition')
        for key in set(observed) & set(wanted):
            for field, value in wanted[key].items():
                label = f'{name} {key} {field}'
                if isinstance(value, (int, Decimal)) or value is None:
                    self.number(observed[key].get(field, ''), value, label)
                else:
                    self.require(observed[key].get(field) == value, label + ': mismatch')


def expected_numbers(c):
    result = []

    def add(name, scope, source, selector, column, operation, value, unit, n, numerator='', denominator=''):
        digits = 2 if unit == 'percent' else 6 if unit == 'distance' else 0
        serialize = lambda x: js(x) if isinstance(x, (list, dict)) else str(x)
        result.append(dict(zip(NUM_COLS, (name, scope, source, js(selector), serialize(column), operation,
            serialize(numerator), serialize(denominator), int(n), unit, value, digits, display(value, digits)))))

    for suffix, role in (('All', None), ('Soil', 'soil'), ('Control', 'negative_control')):
        count = sum(role is None or row['sample_type'] == role for row in c['meta'])
        add('ProjectSamples' + suffix, 'project', META, {} if role is None else {'sample_type': role},
            'sample-id', 'count_rows', count, 'samples', count)
    for suffix, soil in (('All', False), ('Soil', True)):
        count = sum(not soil or row['soil_asv'] == '1' for row in c['assignments'])
        add('ProjectAsv' + suffix, 'project', CMP + 'assignments.tsv', {'soil_asv': '1'} if soil else {},
            'feature_id', 'count_rows', count, 'ASV', count)
    for name, field in (('ProjectTrimLeft', 'trim-left'), ('ProjectTruncLen', 'trunc-len')):
        add(name, 'project', CHOICE, {}, field, 'value', c['choice'][field], 'nt', 1)
    add('ProjectAsvLength', 'project', CHOICE, {}, ['trunc-len', 'trim-left'], 'subtract',
        c['choice']['trunc-len'] - c['choice']['trim-left'], 'nt', 2, c['choice']['trunc-len'], c['choice']['trim-left'])
    kr = c['stats']['1', 'Kr5_S84']; sel = {'run': '1', 'sample-id': 'Kr5_S84'}
    for name, field in (('ProjectKrFiveInput', 'input'), ('ProjectKrFiveRetained', 'non-chimeric')):
        add(name, 'project', STATS, sel, field, 'value', int(kr[field]), 'reads', 1)
    add('ProjectKrFiveRetentionPct', 'project', STATS, sel, ['non-chimeric', 'input'], 'ratio_pct',
        percent(kr['non-chimeric'], kr['input']), 'percent', 1, kr['non-chimeric'], kr['input'])
    add('ProjectModels', 'project', CMP + 'coverage.tsv', {}, 'model', 'count_distinct',
        len({key[0] for key in c['coverage']}), 'models', len(c['coverage']))
    add('SelectedDisplayThresholdPct', FS, FAM + 'prokaryotic_family_manifest.json', {}, 'parameters.min_percent',
        'value', Decimal(c['fm']['parameters']['min_percent']), 'percent', 1)
    for model, tag in zip(MODELS, TAGS):
        sel = {'model': model, 'rank': 'семейство'}
        row = c['coverage'][model, 'семейство']
        add(f'Cmp{tag}FamilyAsvPct', CS, CMP + 'coverage.tsv', sel, ['assigned_asv', 'total_asv'], 'ratio_pct',
            percent(row['assigned_asv'], row['total_asv']), 'percent', 1, row['assigned_asv'], row['total_asv'])
        nums = [int(c['per_coverage'][model, 'семейство', s]['assigned_reads']) for s in c['soils']]
        dens = [int(c['per_coverage'][model, 'семейство', s]['total_reads']) for s in c['soils']]
        add(f'Cmp{tag}FamilyReadMedianPct', CS, CMP + 'coverage_per_sample.tsv', dict(sel, sample=c['soils']),
            ['assigned_reads', 'total_reads'], 'median_ratio_pct', median([percent(n, d) for n, d in zip(nums, dens)]),
            'percent', len(nums), nums, dens)
        count = sum(int(c['filters'][model, s]['selected_reads']) == 0 for s in c['soils'])
        add(f'Selected{tag}ZeroSoilSamples', FS, FAM + 'filter_counts.tsv',
            {'model': model, 'sample_type': 'soil', 'zero_selected': '1'}, 'sample', 'count_rows', count, 'samples', count)
    for (a, b), tag in zip(PAIRS, PAIR_TAGS):
        for rank, name in (('семейство', 'Family'), ('род', 'Genus')):
            values = [Decimal(c['distances'][a, b, rank, s]['distance']) for s in c['soils']]
            add(f'Cmp{tag}{name}DistanceMedian', CS, CMP + 'distances.tsv',
                {'model_a': a, 'model_b': b, 'rank': rank, 'sample': c['soils']}, 'distance', 'median_decimal',
                median(values), 'distance', len(values))
        sel = {'model_a': a, 'model_b': b, 'rank': 'род'}
        row = c['pair'][a, b, 'род']
        add(f'Cmp{tag}GenusConditionalAsvPct', CS, CMP + 'pairwise_match.tsv', sel,
            ['both_match', 'both_assigned_asv'], 'ratio_pct', percent(row['both_match'], row['both_assigned_asv']),
            'percent', 1, row['both_match'], row['both_assigned_asv'])
        selected = [s for s in c['soils'] if int(c['per_pair'][a, b, 'род', s]['both_assigned_reads']) > 0]
        nums = [int(c['per_pair'][a, b, 'род', s]['both_match_reads']) for s in selected]
        dens = [int(c['per_pair'][a, b, 'род', s]['both_assigned_reads']) for s in selected]
        add(f'Cmp{tag}GenusConditionalReadMedianPct', CS, CMP + 'pairwise_match_per_sample.tsv', dict(sel, sample=selected),
            ['both_match_reads', 'both_assigned_reads'], 'median_ratio_pct', median([percent(n, d) for n, d in zip(nums, dens)]),
            'percent', len(nums), nums, dens)
        add(f'Cmp{tag}GenusDefinedSamples', CS, CMP + 'pairwise_match_per_sample.tsv', dict(sel, sample=selected),
            'sample', 'count_rows', len(selected), 'samples', len(selected))
    if len(result) != 37:
        raise ValueError('checker internal macro inventory must have 37 entries')
    return result


def expected_figures(root, c):
    result = []
    for suffix, model in zip(('overview', 'silva', 'gtdb', 'gg2'), ('all', *MODELS)):
        fid = 'prokaryotic_family_' + suffix
        path = FIG + fid + '.pdf'
        if digest(root / path) != digest(root / FAM / 'figures' / (fid + '.pdf')):
            raise ValueError('figure copy differs from validated source: ' + fid)
        result.append(dict(zip(FIG_COLS, (FS, fid, model, path, digest(root / path), FAM + 'family_barplot_data.tsv',
            'selected_reads(model,sample)', Decimal(c['fm']['parameters']['min_percent'])))))
    return result


def check_tex(audit, out, c, numbers):
    text = (out / 'numbers.tex').read_text()
    # Every non-comment line must be exactly one declaration. No hidden extra code.
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('%')]
    declarations = []
    for line in lines:
        match = re.fullmatch(r'\\newcommand\{\\([A-Za-z]+)\}\{([^{}]+)\}', line)
        audit.require(match is not None, 'numbers.tex: invalid declaration')
        if match:
            declarations.append(match.groups())
    audit.require(len(declarations) == 37 and len(dict(declarations)) == 37, 'numbers.tex: macro count/duplicates')
    audit.require(dict(declarations) == {r['id']: r['display_value'] for r in numbers}, 'numbers.tex: actual values/identities differ')
    coverage = index(c['cov_rows'], ('model', 'rank'))
    pairs = index(c['pair_rows'], ('model_a', 'model_b', 'rank'))
    comp = ['% scope: ' + CS, '% Read coverage median (percent); D median; all post-DADA2 reads.',
            r'\begin{tabular}{lrrrrrr}', r'\hline',
            r'Rank & \multicolumn{3}{c}{Read coverage (\%)} & \multicolumn{3}{c}{D} \\',
            r' & SILVA & GTDB & GG2 & SILVA--GTDB & SILVA--GG2 & GTDB--GG2 \\', r'\hline']
    for rank, label in zip(RANKS, ('Phylum', 'Class', 'Order', 'Family', 'Genus')):
        cells = [label] + [display(coverage[m, rank]['read_coverage_median_pct'], 2) for m in MODELS]
        cells += [display(pairs[a, b, rank]['distance_median'], 6) for a, b in PAIRS]
        comp.append(' & '.join(cells) + r' \\')
    comp += [r'\hline', r'\end{tabular}']
    audit.require((out / 'comparison_main.tex').read_text().splitlines() == comp, 'comparison_main.tex: labels/values/scope/layout mismatch')
    family = index(c['fam_rows'], ('model', 'sample'))
    fam = ['% scope: ' + FS, '% Reads: post-DADA2; N: selected reads; U: UNCLASSIFIED / N (percent); N=0 gives NA.',
           r'\begin{tabular}{lrrrrrrr}', r'\hline',
           r'Sample & Reads & \multicolumn{2}{c}{SILVA} & \multicolumn{2}{c}{GTDB} & \multicolumn{2}{c}{GG2} \\',
           r' & & N & U (\%) & N & U (\%) & N & U (\%) \\', r'\hline']
    for sample in c['soils']:
        cells = [sample.replace('_', r'\_'), str(family[MODELS[0], sample]['post_dada2_reads'])]
        for model in MODELS:
            row = family[model, sample]
            cells += [str(row['selected_reads']), display(row['unclassified_pct_of_selected'], 2)]
        fam.append(' & '.join(cells) + r' \\')
    fam += [r'\hline', r'\end{tabular}']
    audit.require((out / 'family_main.tex').read_text().splitlines() == fam, 'family_main.tex: labels/values/scope/layout mismatch')


def source_inventory():
    paths = {META, CHOICE, STATS, 'research/paper_results_contract.md', FIG + 'manifest.json'}
    paths.update(CMP + name + '.tsv' for name in ('assignments', 'coverage', 'coverage_per_sample', 'pairwise_match', 'pairwise_match_per_sample', 'distances'))
    paths.update(FAM + name + '.tsv' for name in ('filter_counts', 'family_profiles', 'asv_mask', 'family_barplot_data'))
    paths.update((FAM + 'prokaryotic_family_manifest.json', FAM + 'verification.json'))
    paths.update(FIG + 'prokaryotic_family_' + suffix + '.pdf' for suffix in ('overview', 'silva', 'gtdb', 'gg2'))
    return paths


def check_manifest(audit, out, root):
    manifest = json.loads((out / 'manifest.json').read_text())
    audit.require(manifest.get('schema_version') == 1, 'manifest: schema_version')
    audit.require(manifest.get('manuscript_status') == 'absent_not_checked', 'manifest: manuscript must be absent_not_checked')
    expected_sets = {'sources': source_inventory(), 'code': {'src/make_paper_results.py', 'src/paper_results_core.py'}, 'outputs': OUTPUTS}
    for name, required in expected_sets.items():
        mapping = manifest.get(name, {})
        audit.require(isinstance(mapping, dict) and set(mapping) == required, 'manifest: ' + name + ' inventory')
        if not isinstance(mapping, dict):
            continue
        base = out if name == 'outputs' else root
        for path in required:
            audit.require(mapping.get(path) == digest(base / path), 'manifest: ' + name + ' SHA256 ' + path)
    audit.require(not re.search(r'/Users/|/private/|/tmp/', (out / 'manifest.json').read_text()), 'manifest: absolute paths')


def audit_bundle(out_dir, source_root=ROOT):
    out, root = Path(out_dir).resolve(), Path(source_root).resolve()
    audit = Audit()
    try:
        with localcontext() as ctx:
            ctx.prec = 70
            before = {p: digest(root / p) for p in source_inventory()}
            c = expected(root)
            audit.table(out, 'comparison_coverage.tsv', COV_COLS, c['cov_rows'], ('model', 'rank'))
            audit.table(out, 'comparison_pairs.tsv', PAIR_COLS, c['pair_rows'], ('model_a', 'model_b', 'rank'))
            audit.table(out, 'family_sample_counts.tsv', FAM_COLS, c['fam_rows'], ('model', 'sample'))
            audit.table(out, 'figure_sources.tsv', FIG_COLS, expected_figures(root, c), ('figure_id',))
            numbers = expected_numbers(c)
            audit.table(out, 'numbers.tsv', NUM_COLS, numbers, ('id',))
            check_tex(audit, out, c, numbers)
            check_manifest(audit, out, root)
            audit.require(before == {p: digest(root / p) for p in before}, 'sources changed during verification')
    except Exception as exc:
        audit.issues.append(f'{type(exc).__name__}: {exc}')
    return {'status': 'FAIL' if audit.issues else 'PASS', 'checks': audit.checks, 'issues': audit.issues,
            'manuscript_status': 'absent_not_checked', 'scope': 'article-number transfer; biological validity not assessed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', required=True, type=Path)
    parser.add_argument('--source-root', default=ROOT, type=Path)
    parser.add_argument('--report', type=Path, help='New JSON report path; never overwrites an existing report.')
    args = parser.parse_args()
    result = audit_bundle(args.out_dir, args.source_root)
    encoded = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if args.report:
        with args.report.open('x', encoding='utf-8') as stream:
            stream.write(encoded)
    print(encoded, end='')
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())


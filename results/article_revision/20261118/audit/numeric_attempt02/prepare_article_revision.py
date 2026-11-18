"""Export requested article diagnostics from existing tables and fixed QZA.

This does not run denoising, classification or the comparison pipeline. Existing
results are read only. Existing export files are never overwritten. The independent
checker lives in check_article_revision.py and does not import this module.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import io
import json
import re
from pathlib import Path
from statistics import median
import zipfile

import biom
import h5py

ROOT = Path(__file__).resolve().parents[1]
MAIN = 'results/taxonomy_comparison/finish_20261116/main'
FAMILY = 'results/taxonomy_barplots/prokaryote_family_20261118_final'
PAIRS = [('sklearn_silva138', 'sklearn_gtdb', 'SilvaGtdb', 'SILVA--GTDB'),
         ('sklearn_silva138', 'sklearn_gg2', 'SilvaGgTwo', 'SILVA--Greengenes2'),
         ('sklearn_gtdb', 'sklearn_gg2', 'GtdbGgTwo', 'GTDB--Greengenes2')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(relative):
    with (ROOT / relative).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path, records):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def fixed(value, decimals=2):
    return format(Decimal(str(value)).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_EVEN), 'f')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    out = args.out_dir.resolve()
    output_names = ['kr_filter_reasons.tsv', 'genus_conditional.tsv', 'genus_conditional_per_sample.tsv',
                    'control_summary.tsv', 'control_overlap.tsv', 'control_shared_asv.tsv',
                    'method_parameters.tsv', 'numbers.tsv', 'numbers.tex', 'kr_filter_table.tex',
                    'genus_denominators_table.tex', 'summary.json', 'manifest.json']
    existing = [name for name in output_names if (out / name).exists()]
    if existing:
        raise SystemExit('Refusing to overwrite existing exports: ' + ', '.join(existing))
    source_paths = [f'{MAIN}/assignments.tsv', f'{MAIN}/pairwise_match.tsv',
                    f'{MAIN}/pairwise_match_per_sample.tsv', f'{FAMILY}/filter_counts.tsv',
                    'results/read_counts_R1.tsv', 'work/qza/ASV_table.qza',
                    'work/qza/dada2_stats.qza', 'data/meta/metadata.tsv',
                    'src/prepare_article_revision.py', 'checks/verify_denoising.py', 'src/check_exact_match.py']
    sources = {p: sha(ROOT / p) for p in source_paths}
    assignments = {r['feature_id']: r for r in rows(f'{MAIN}/assignments.tsv')}
    with zipfile.ZipFile(ROOT / 'work/qza/ASV_table.qza') as archive:
        member = next(n for n in archive.namelist() if n.endswith('/data/feature-table.biom'))
        with h5py.File(io.BytesIO(archive.read(member)), 'r') as stream:
            table = biom.Table.from_hdf5(stream)
    feature_ids = list(table.ids(axis='observation'))
    sample_ids = list(table.ids(axis='sample'))
    values = table.matrix_data.toarray()
    assert (values >= 0).all() and (values.astype(int) == values).all()
    counts = {fid: {sample: int(value) for sample, value in zip(sample_ids, row)}
              for fid, row in zip(feature_ids, values)}
    assert set(counts) == set(assignments)
    metadata = {r['sample-id']: r['sample_type'] for r in rows('data/meta/metadata.tsv')}
    assert set(metadata) == set(sample_ids)
    soil = sorted(s for s in sample_ids if metadata[s] == 'soil')
    assert [s for s in sample_ids if metadata[s] == 'negative_control'] == ['K_S63']
    soil_asv = {f for f, row in counts.items() if sum(row[s] for s in soil)}
    assert soil_asv == {f for f, row in assignments.items() if row['soil_asv'] == '1'}
    filters = rows(f'{FAMILY}/filter_counts.tsv')
    kr = []
    for row in filters:
        if row['model'] == 'sklearn_silva138' and row['sample'].startswith('Kr'):
            record = {'sample': row['sample'], **{k: int(row[k]) for k in
                      ['total_reads', 'Eukaryota_reads', 'unknown_reads', 'organelle_reads', 'selected_reads']}}
            assert record['total_reads'] == sum(record[k] for k in
                   ['Eukaryota_reads', 'unknown_reads', 'organelle_reads', 'selected_reads'])
            assert record['total_reads'] == sum(c[record['sample']] for c in counts.values())
            kr.append(record)

    pooled = rows(f'{MAIN}/pairwise_match.tsv')
    saved_samples = rows(f'{MAIN}/pairwise_match_per_sample.tsv')
    genus = []
    per_sample = []
    for a, b, suffix, label in PAIRS:
        raw = next(r for r in pooled if (r['model_a'], r['model_b'], r['rank']) == (a, b, 'род'))
        both = {f for f in soil_asv if assignments[f][a + '_род_assigned'] == '1'
                and assignments[f][b + '_род_assigned'] == '1'}
        matching = {f for f in both if assignments[f][a + '_род_key']
                    == assignments[f][b + '_род_key']}
        assert len(both) == int(raw['both_assigned_asv'])
        assert len(matching) == int(raw['both_match'])
        pair_samples = []
        for sample in soil:
            saved = next(r for r in saved_samples if
                         (r['model_a'], r['model_b'], r['rank'], r['sample']) == (a, b, 'род', sample))
            both_asv = sum(counts[f][sample] > 0 for f in both)
            both_reads = sum(counts[f][sample] for f in both)
            matching_reads = sum(counts[f][sample] for f in matching)
            assert both_reads == int(saved['both_assigned_reads'])
            assert matching_reads == int(saved['both_match_reads'])
            assert (both_reads == 0) == (saved['cond_match_reads_frac'] == 'NA')
            record = dict(model_a=a, model_b=b, rank='род', sample=sample,
                          observed_asv=sum(c[sample] > 0 for c in counts.values()),
                          both_assigned_asv=both_asv, both_assigned_reads=both_reads,
                          both_match_asv=sum(counts[f][sample] > 0 for f in matching),
                          both_match_reads=matching_reads,
                          conditional_read_fraction=saved['cond_match_reads_frac'])
            pair_samples.append(record)
        per_sample.extend(pair_samples)
        defined = [r for r in pair_samples if r['both_assigned_reads'] > 0]
        genus.append(dict(model_a=a, model_b=b, rank='род',
                          pooled_both_assigned_asv=int(raw['both_assigned_asv']),
                          pooled_both_match_asv=int(raw['both_match']),
                          pooled_conditional_asv_fraction=raw['cond_match_asv_frac'],
                          defined_samples=len(defined),
                          conditional_read_median=format(median(float(r['conditional_read_fraction']) for r in defined), '.17g'),
                          conditional_sample_asv_min=min(r['both_assigned_asv'] for r in defined),
                          conditional_sample_asv_max=max(r['both_assigned_asv'] for r in defined),
                          conditional_sample_asv_median=median(r['both_assigned_asv'] for r in defined)))

    control_ids = {f for f, row in counts.items() if row['K_S63'] > 0}
    shared = control_ids & soil_asv
    input_control = next(r for r in rows('results/read_counts_R1.tsv') if r['sample-id'] == 'K_S63')
    with zipfile.ZipFile(ROOT / 'work/qza/dada2_stats.qza') as archive:
        member = next(n for n in archive.namelist() if n.endswith('/data/stats.tsv'))
        lines = [line for line in archive.read(member).decode().splitlines() if not line.startswith('#q2:')]
        dada = next(r for r in csv.DictReader(lines, delimiter='\t') if r['sample-id'] == 'K_S63')
    control_reads = sum(row['K_S63'] for row in counts.values())
    assert control_reads == int(dada['non-chimeric'])
    assert int(input_control['reads_R1']) == int(dada['input'])
    control = [dict(sample='K_S63', input_reads=int(input_control['reads_R1']),
                    retained_reads=control_reads, observed_asv=len(control_ids),
                    shared_soil_asv=len(shared), control_only_asv=len(control_ids - soil_asv),
                    shared_control_reads=sum(counts[f]['K_S63'] for f in shared),
                    total_soil_reads=sum(row[s] for row in counts.values() for s in soil),
                    shared_soil_reads=sum(counts[f][s] for f in shared for s in soil))]
    control[0]['shared_soil_read_fraction'] = format(control[0]['shared_soil_reads'] / control[0]['total_soil_reads'], '.17g')
    overlap = [dict(sample=s, total_reads=sum(c[s] for c in counts.values()),
                    observed_asv=sum(c[s] > 0 for c in counts.values()),
                    shared_control_asv=sum(counts[f][s] > 0 for f in control_ids),
                    shared_control_reads=sum(counts[f][s] for f in control_ids)) for s in soil]
    shared_ids = [dict(feature_id=f, control_reads=counts[f]['K_S63'],
                       soil_reads=sum(counts[f][s] for s in soil),
                       soil_samples=';'.join(s for s in soil if counts[f][s])) for f in sorted(shared)]

    thresholds = []
    for path in ['checks/verify_denoising.py', 'src/check_exact_match.py']:
        code = (ROOT / path).read_text()
        found = re.findall(r'0\s*<\s*100\s*\*\s*kept\s*<\s*(\d+)\s*\*\s*n', code)
        assert len(found) == 1, (path, found)
        thresholds.append(int(found[0]))
    assert thresholds[0] == thresholds[1]
    parameters = [dict(parameter='retention_threshold_pct', value=thresholds[0],
                       source='checks/verify_denoising.py;src/check_exact_match.py',
                       locator='0 < 100 * kept < threshold * n; soil samples')]
    out.mkdir(parents=True, exist_ok=True)
    tables = {'kr_filter_reasons.tsv': kr, 'genus_conditional.tsv': genus,
              'genus_conditional_per_sample.tsv': per_sample, 'control_summary.tsv': control,
              'control_overlap.tsv': overlap, 'control_shared_asv.tsv': shared_ids,
              'method_parameters.tsv': parameters}
    for name, records in tables.items():
        write_tsv(out / name, records)
    records = []
    def add(name, value, table_name, locator, rounding='integer'):
        records.append(dict(macro='Revision' + name, value=str(value),
                            source=str((out / table_name).relative_to(ROOT)), locator=locator, rounding=rounding))
    for row in kr:
        sample = row['sample']
        for field, suffix in [('total_reads', 'TotalReads'), ('Eukaryota_reads', 'EukaryotaReads'),
                              ('unknown_reads', 'UnknownReads'), ('organelle_reads', 'OrganelleReads'),
                              ('selected_reads', 'SelectedReads')]:
            kr_name = {'Kr1': 'KrOne', 'Kr2': 'KrTwo', 'Kr3': 'KrThree', 'Kr4': 'KrFour', 'Kr5': 'KrFive'}[sample.split('_')[0]]
            add(kr_name + 'Silva' + suffix, row[field], 'kr_filter_reasons.tsv', sample + ';' + field)
    for (a, b, suffix, label), row in zip(PAIRS, genus):
        for field, ending in [('pooled_both_assigned_asv', 'PooledBothAssignedAsv'),
                              ('pooled_both_match_asv', 'PooledBothMatchAsv'), ('defined_samples', 'ConditionalDefinedSamples'),
                              ('conditional_sample_asv_min', 'ConditionalSampleAsvMin'),
                              ('conditional_sample_asv_max', 'ConditionalSampleAsvMax'),
                              ('conditional_sample_asv_median', 'ConditionalSampleAsvMedian')]:
            add(suffix + 'Genus' + ending, row[field], 'genus_conditional.tsv', a + ';' + b + ';' + field)
        add(suffix + 'GenusConditionalReadMedianPct', fixed(Decimal(row['conditional_read_median']) * 100),
            'genus_conditional.tsv', a + ';' + b + ';conditional_read_median*100', 'HALF_EVEN;2 decimals')
    for field, suffix in [('input_reads', 'InputReads'), ('retained_reads', 'RetainedReads'),
                          ('observed_asv', 'ObservedAsv'), ('shared_soil_asv', 'SharedSoilAsv'),
                          ('control_only_asv', 'OnlyAsv'), ('shared_control_reads', 'SharedControlReads'),
                          ('total_soil_reads', 'TotalSoilReads'), ('shared_soil_reads', 'SharedSoilReads')]:
        add('Control' + suffix, control[0][field], 'control_summary.tsv', 'K_S63;' + field)
    add('ControlSharedSoilReadPct', fixed(Decimal(control[0]['shared_soil_read_fraction']) * 100),
        'control_summary.tsv', 'K_S63;shared_soil_read_fraction*100', 'HALF_EVEN;2 decimals')
    add('RetentionThresholdPct', thresholds[0], 'method_parameters.tsv', 'retention_threshold_pct;value')
    write_tsv(out / 'numbers.tsv', records)
    (out / 'numbers.tex').write_text('% Generated from diagnostic TSVs; immutable input tables and QZA.\n' +
        ''.join('\\newcommand{\\' + r['macro'] + '}{' + r['value'] + '}\n' for r in records))
    kr_tex = '\\begin{tabular}{lrrrrr}\n\\toprule\nОбразец & Всего & Eukaryota & Домен не задан & Органеллы & Осталось \\\\\n\\midrule\n'
    for row in kr:
        kr_tex += row['sample'].replace('_', '\\_') + ' & ' + ' & '.join(str(row[k]) for k in
                  ['total_reads', 'Eukaryota_reads', 'unknown_reads', 'organelle_reads', 'selected_reads']) + ' \\\\\n'
    (out / 'kr_filter_table.tex').write_text(kr_tex + '\\bottomrule\n\\end{tabular}\n')
    genus_tex = '\\begin{tabular}{lrrr}\n\\toprule\nПара моделей & Общие ASV & Образцы & ASV в образце \\\\\n\\midrule\n'
    for (_, _, _, label), row in zip(PAIRS, genus):
        genus_tex += f"{label} & {row['pooled_both_assigned_asv']} & {row['defined_samples']} & {row['conditional_sample_asv_min']}--{row['conditional_sample_asv_max']} \\\\\n"
    (out / 'genus_denominators_table.tex').write_text(genus_tex + '\\bottomrule\n\\end{tabular}\n')
    summary = dict(scope='Read-only export and diagnostic aggregation, no new scientific pipeline run.',
                   kr_filter_reasons=kr, genus_conditional=genus, control=control[0], method_parameters=parameters,
                   definitions={'conditional_asv': 'ASVs with nonzero frequency and genus assigned by both models; complete lineage matching is case-sensitive after trimming whitespace.',
                                'conditional_defined_samples': 'Samples with nonzero both_assigned_reads; zero denominators remain NA and are excluded from the median.',
                                'pooled_both_assigned_asv': 'Distinct ASVs across all ten soil samples, not a sum across sample-specific denominators.',
                                'control_overlap': 'ASV identifier set intersection, not proof of contamination or its source; no ASVs removed.',
                                'kr_filter': 'Exclusive categories: Eukaryota, unknown domain, explicit organelle labels within Bacteria/Archaea, selected Bacteria/Archaea.'})
    (out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    assert sources == {p: sha(ROOT / p) for p in source_paths}, 'Sources changed during export'
    manifest = dict(sources=sources, outputs={name: sha(out / name) for name in sorted(output_names) if name != 'manifest.json'},
                    derivation='Saved TSV transfer; additional control intersections and per-sample ASV denominators from fixed post-DADA2 BIOM.',
                    primary_consistency_checks='BIOM totals vs DADA2/control input; sample roles; soil ASV set; pooled/per-sample agreement counts against saved comparison; mutually exclusive Kr partition.')
    (out / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'EXPORTED', 'out_dir': str(out), 'macros': len(records),
                      'control': control[0], 'genus': genus}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()

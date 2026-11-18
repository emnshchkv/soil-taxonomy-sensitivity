#!/usr/bin/env python3
"""Independent raw-QZA-to-family audit, AI4Science-агент.

Only artifact reading is shared with the project IO module. No generator,
renderer or taxonomy arithmetic is imported. Mutations use temporary copies.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
from taxonomy_compare_io import load_dataset, MODELS

ROOT = Path(__file__).resolve().parent.parent
UNKNOWN, OTHER = 'UNCLASSIFIED', 'OTHER_ASSIGNED'
SERVICE = {'', 'Unassigned', 'Unknown_Family'}
DOMAINS = ('Bacteria', 'Archaea', 'Eukaryota', 'unknown')
MANIFEST = 'prokaryotic_family_manifest.json'
TOL = Fraction(1, 10 ** 12)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def indexed(rows, fields, label):
    result = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        need(key not in result, 'Duplicate ' + label + ': ' + str(key))
        result[key] = row
    return result


def exact(row, expected, label):
    for field, value in expected.items():
        need(row.get(field) == str(value), label + ': ' + field + ' mismatch')


def fraction(text, count, total, label):
    need(text == 'NA' if total == 0 else text != 'NA' and
         abs(Fraction(text) - Fraction(count, total)) <= TOL, label + ': fraction/NA mismatch')


def taxon(raw):
    fields = {}
    for token in raw.strip().split(';'):
        token = token.strip()
        if not token or token == 'Unassigned':
            continue
        match = re.fullmatch(r'([dpcofgs])__(.*)', token)
        numeric = re.fullmatch(r'D_([0-6])__(.*)', token)
        need(bool(match or numeric), 'Unrecognized raw taxonomy token: ' + token)
        code, value = (match[1], match[2]) if match else ('dpcofgs'[int(numeric[1])], numeric[2])
        need(code not in fields, 'Duplicate raw taxonomy rank: ' + raw)
        fields[code] = value.strip()
    domain = fields.get('d') or 'unknown'
    need(domain in DOMAINS, 'Unexpected domain requires review: ' + domain)
    markers = sorted({value.lower() for value in fields.values()
                      if value.lower() in {'mitochondria', 'chloroplast'}})
    reason = ('eukaryota' if domain == 'Eukaryota' else 'unknown_domain' if domain == 'unknown'
              else 'organelle' if markers else 'selected')
    family = fields.get('f', '')
    raw_key = UNKNOWN if family in SERVICE else ';'.join(code + '__' + fields[code]
              for code in 'dpcof' if fields.get(code, '') not in SERVICE)
    unknown_reason = ('missing' if 'f' not in fields else 'service' if family in SERVICE
                      else 'uncultured' if family.casefold() == 'uncultured' else '')
    return domain, markers, reason, raw_key, UNKNOWN if unknown_reason else raw_key, family, unknown_reason


def audit(out_dir, source_root=ROOT):
    out, root = Path(out_dir).resolve(), Path(source_root).resolve()
    checks = 0
    try:
        manifest = json.loads((out / MANIFEST).read_text())
        params = manifest['parameters']; threshold = Fraction(params['threshold_fraction'])
        need(manifest['schema'] == 'prokaryotic-family-barplots-v1', 'Unexpected output schema')
        need(params['domain_policy'] == 'per-model explicit Bacteria or Archaea, excluding explicit organelle tokens', 'Wrong domain policy')
        need(0 < threshold <= 1 and Fraction(params['min_percent']) / 100 == threshold, 'Inconsistent display threshold')
        need(params['family_unknown_labels'] == ['', 'Unassigned', 'Unknown_Family', 'uncultured'] and
             params['family_uncultured_casefold'] is True, 'Wrong unknown-family policy')
        need(params['rank'] == 'семейство' and params['models'] == list(MODELS), 'Wrong rank/model order')
        data, issues = load_dataset(root)
        need(not issues, 'Raw source violations: ' + json.dumps(issues, ensure_ascii=False))
        soil, samples, table = data['soil_samples'], data['samples'], data['table']
        roles = {sample: 'soil' if sample in soil else 'negative_control' for sample in samples}
        need(params['soil_samples'] == soil and params['all_samples'] == samples and len(soil) == 10, 'Wrong sample sets/order')
        sources = {row['path']: row['sha256'] for row in data['inputs'].values()}
        need(manifest['sources'] == sources and manifest['inputs_unchanged'] is True, 'Source identity/preservation mismatch')
        need({'src/plot_prokaryotic_families.py', 'src/render_prokaryotic_families.py',
              'src/check_prokaryotic_families.py', 'src/taxonomy_compare_io.py',
              'src/taxonomy_compare_core.py', 'src/classification_contract.py'} <= set(manifest['code_sha256']), 'Incomplete code identity')
        for path, digest in manifest['code_sha256'].items():
            need(sha(ROOT / path) == digest, 'Code hash changed: ' + path)
        for path, digest in manifest['outputs'].items():
            need(sha(out / path) == digest, 'Output hash changed: ' + path)
        mask = indexed(read(out / 'asv_mask.tsv'), ('model', 'feature_id'), 'ASV mask')
        need(set(mask) == {(model, fid) for model in MODELS for fid in table}, 'Incomplete ASV mask')
        details = {}
        for model in MODELS:
            for fid, raw in data['taxonomies'][model].items():
                domain, markers, reason, raw_key, key, label, unknown_reason = taxon(raw)
                counts = table[fid]; soil_reads = sum(counts[sample] for sample in soil)
                details[model, fid] = (domain, markers, reason, key)
                row = mask[model, fid]
                exact(row, dict(raw_taxon=raw, domain=domain, selected=int(reason == 'selected'),
                      exclusion_reason=reason, family_key=key, raw_family_key=raw_key,
                      raw_family_label=label, unknown_family_reason=unknown_reason,
                      soil_asv=int(soil_reads > 0), all_reads=sum(counts.values()), soil_reads=soil_reads,
                      control_reads=sum(counts[s] for s in samples if roles[s] == 'negative_control')), 'ASV mask ' + fid)
                need(json.loads(row['organelle_markers_json']) == markers, 'Wrong organelle union markers')
                checks += 1
        losses = indexed(read(out / 'filter_counts.tsv'), ('model', 'sample'), 'filter counts')
        actual_profiles = indexed(read(out / 'family_profiles.tsv'), ('model', 'sample', 'key'), 'family profile')
        memberships = indexed(read(out / 'family_membership.tsv'), ('model', 'sample', 'key'), 'family membership')
        need(set(losses) == {(m, s) for m in MODELS for s in samples}, 'Missing sample/control filtering audit')
        totals, profiles, displays, colored, expected_keys = {}, {}, {}, {m: set() for m in MODELS}, set()
        threshold_cases = Counter()
        for model in MODELS:
            for sample in samples:
                counts = {fid: row[sample] for fid, row in table.items()}
                total = sum(counts.values())
                domains = {domain: sum(n for fid, n in counts.items() if details[model, fid][0] == domain) for domain in DOMAINS}
                selected = [fid for fid in table if details[model, fid][2] == 'selected']
                organelles = [fid for fid in table if details[model, fid][2] == 'organelle']
                n = sum(counts[fid] for fid in selected); totals[model, sample] = n
                organelle_reads = sum(counts[fid] for fid in organelles)
                need(sum(domains.values()) == total and domains['Bacteria'] + domains['Archaea'] == n + organelle_reads, 'Filtering partition does not conserve mass')
                expected = dict(sample_type=roles[sample], total_reads=total,
                    total_observed_asv=sum(value > 0 for value in counts.values()),
                    **{d + '_reads': value for d, value in domains.items()},
                    organelle_reads=organelle_reads, organelle_observed_asv=sum(counts[fid] > 0 for fid in organelles),
                    selected_reads=n, selected_observed_asv=sum(counts[fid] > 0 for fid in selected), excluded_reads=total - n,
                    selected_Bacteria_reads=sum(counts[fid] for fid in selected if details[model, fid][0] == 'Bacteria'),
                    selected_Archaea_reads=sum(counts[fid] for fid in selected if details[model, fid][0] == 'Archaea'), zero_selected=int(n == 0))
                exact(losses[model, sample], expected, 'Filter counts ' + model + '/' + sample)
                fraction(losses[model, sample]['selected_fraction'], n, total, 'Filter retained fraction')
                values = Counter({UNKNOWN: 0})
                for fid in selected:
                    if counts[fid]: values[details[model, fid][3]] += counts[fid]
                need(sum(values.values()) == n, 'Full family profile mass mismatch')
                profiles[model, sample] = values; displays[model, sample] = Counter()
                for key, count in values.items():
                    pkey = model, sample, key; expected_keys.add(pkey)
                    category = UNKNOWN if key == UNKNOWN else key if n and count and Fraction(count, n) >= threshold else OTHER
                    displays[model, sample][category] += count
                    if key != UNKNOWN:
                        threshold_cases['at_or_above' if category == key else 'below'] += 1
                        if n and Fraction(count, n) == threshold: threshold_cases['exact_boundary'] += 1
                    if sample in soil and category not in (UNKNOWN, OTHER): colored[model].add(category)
                    need(pkey in actual_profiles and pkey in memberships, 'Missing family or threshold-membership row')
                    exact(actual_profiles[pkey], dict(rank='семейство', sample_type=roles[sample], count=count, total_reads=n), 'Family profile')
                    exact(memberships[pkey], dict(sample_type=roles[sample], count=count, total_reads=n, display_category=category), 'Per-sample threshold membership')
                    fraction(actual_profiles[pkey]['fraction'], count, n, 'Family profile')
                    fraction(memberships[pkey]['fraction'], count, n, 'Family membership'); checks += 1
        need(set(actual_profiles) == expected_keys == set(memberships), 'Extra/missing full family keys or samples')
        legend = indexed(read(out / 'family_keys.tsv'), ('model', 'category'), 'legend family')
        segments = indexed(read(out / 'family_barplot_data.tsv'), ('model', 'sample', 'category'), 'bar segment')
        expected_legend, expected_segments, colors, labels = set(), set(), {}, {}
        percent = format(float(threshold * 100), '.12g')
        for model in MODELS:
            categories = colored[model] | {OTHER, UNKNOWN}
            panel_orders = []
            for category in categories:
                lkey = model, category; expected_legend.add(lkey)
                need(lkey in legend, 'Panel legend omits a family passing the threshold')
                row = legend[lkey]; kind = 'other_assigned' if category == OTHER else 'unclassified' if category == UNKNOWN else 'selected'
                full_key = '' if category == OTHER else category
                exact(row, dict(category_kind=kind, full_key=full_key), 'Legend category')
                order = int(row['order']); panel_orders.append(order)
                need(re.fullmatch(r'#[0-9a-fA-F]{6}', row['color']), 'Invalid legend colour')
                if kind == 'selected':
                    fields = dict(token.split('__', 1) for token in category.split(';'))
                    need(set(fields) <= set('dpcof') and fields.get('f') not in SERVICE | {'uncultured'}, 'Legend is not a known family-level key')
                    need(fields['f'] in row['label'], 'Legend omits the actual family label')
                    need(row['color'] not in ('#686868', '#dddddd'), 'Family uses reserved grey')
                    need(category not in colors or colors[category] == row['color'], 'Same full key changes colour')
                    need(category not in labels or labels[category] == row['label'], 'Same full key changes label')
                    colors[category], labels[category] = row['color'], row['label']
                else:
                    exact(row, dict(color='#686868' if category == OTHER else '#dddddd',
                          label=f'Семейства <{percent}%' if category == OTHER else 'Семейство не определено'), 'Neutral legend')
                for sample in soil:
                    skey = model, sample, category; expected_segments.add(skey)
                    need(skey in segments, 'Missing category, including zero-valued legend members')
                    count = displays[model, sample].get(category, 0); n = totals[model, sample]
                    exact(segments[skey], dict(sample_type='soil', rank='семейство', category_kind=kind,
                          count=count, total_reads=n, order=order, label=row['label'], color=row['color'], full_key=full_key), 'Segment ' + str(skey))
                    fraction(segments[skey]['fraction'], count, n, 'Segment fraction'); checks += 1
            need(sorted(panel_orders) == list(range(len(categories))), 'Panel legend orders are not unique/contiguous')
            for sample in soil:
                need(sum(int(segments[model, sample, category]['count']) for category in categories) == totals[model, sample], 'Displayed mass mismatch')
        need(set(legend) == expected_legend and set(segments) == expected_segments, 'Extra displayed category or control in soil plots')
        need(len(set(colors.values())) == len(colors) and len(set(labels.values())) == len(labels), 'Distinct family keys have indistinguishable colours or labels')
        need(manifest['colors'] == colors, 'Manifest full-key colours differ')
        need(manifest['nonempty_soil_samples'] == {m: [s for s in soil if totals[m, s]] for m in MODELS}, 'Wrong nonempty soil set')
        required = {'asv_mask.tsv', 'filter_counts.tsv', 'family_profiles.tsv', 'family_membership.tsv', 'family_keys.tsv', 'family_barplot_data.tsv'}
        required |= {f'figures/prokaryotic_family_{name}.{ext}' for name in ('overview', 'silva', 'gtdb', 'gg2') for ext in ('pdf', 'png', 'svg')}
        need(required <= set(manifest['outputs']), 'Missing required output hashes')
        need(all(sha(root / path) == digest for path, digest in sources.items()), 'Inputs changed during verification')
        return dict(status='PASS', checks=checks, scope='Per-model prokaryotic family profiles and per-sample display threshold',
                    threshold_fraction=str(threshold), threshold_cases=dict(threshold_cases), soil_bars=30,
                    zero_selected_soil_samples={m: [s for s in soil if not totals[m, s]] for m in MODELS}, errors=[])
    except (ValueError, OSError, KeyError, TypeError, ZeroDivisionError) as error:
        return dict(status='FAIL', checks=checks, errors=[str(error)])


def self_test(out, root):
    baseline = audit(out, root)
    need(baseline['status'] == 'PASS', 'Self-test baseline failed: ' + str(baseline['errors']))
    data_name = 'family_barplot_data.tsv'
    originals = {name: (out / name).read_bytes() for name in (MANIFEST, data_name)}
    with tempfile.TemporaryDirectory(prefix='prokaryotic-family-selftest-') as directory:
        copy = Path(directory) / 'outputs'; shutil.copytree(out, copy)
        rows = read(copy / data_name)
        donor = next(row for row in rows if int(row['count']) > 0)
        receiver = next(row for row in rows if row['model'] == donor['model'] and row['sample'] == donor['sample'] and row['category'] != donor['category'])
        for row, change in ((donor, -1), (receiver, 1)):
            row['count'] = str(int(row['count']) + change)
            row['fraction'] = format(int(row['count']) / int(row['total_reads']), '.17g')
        with (copy / data_name).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
        manifest = json.loads((copy / MANIFEST).read_text())
        manifest['outputs'][data_name] = sha(copy / data_name)
        (copy / MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
        mutated = audit(copy, root)
        need(mutated['status'] == 'FAIL' and any('Segment ' in error and 'count mismatch' in error for error in mutated['errors']), 'Semantic corruption was not detected')
        for name, content in originals.items(): (copy / name).write_bytes(content)
        restored = audit(copy, root)
        need(restored['status'] == 'PASS', 'Restored copy did not pass')
    need(all((out / name).read_bytes() == content for name, content in originals.items()), 'Original results changed')
    return dict(status='PASS', baseline='PASS', mutated='FAIL', restored='PASS', detected=mutated['errors'],
                mutation='Moved one read between two categories of one bar, preserving N; recomputed fractions and refreshed the temporary output hash.',
                model=donor['model'], sample=donor['sample'], isolation='Only a temporary output copy changed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--json-out', type=Path)
    args = parser.parse_args()
    try:
        result = self_test(args.out_dir.resolve(), args.source_root.resolve()) if args.self_test else audit(args.out_dir, args.source_root)
    except (ValueError, OSError, KeyError, StopIteration) as error:
        result = dict(status='FAIL', errors=[str(error)])
    if args.json_out:
        with args.json_out.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2); stream.write('\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())

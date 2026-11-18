#!/usr/bin/env python3
"""Independent source-to-barplot verification, AI4Science-агент.

No generator or renderer arithmetic is imported. Reads the source ASV table,
raw taxonomies, and metadata again. The self-test changes only a temporary copy.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile

from taxonomy_compare_io import load_dataset, MODELS, RANKS

ROOT = Path(__file__).resolve().parent.parent
ASSIGNMENTS = 'results/taxonomy_comparison/finish_20261116/main/assignments.tsv'
RANK_CODE = dict(zip(RANKS, ('p', 'c', 'o', 'f', 'g')))
FILE_RANK = dict(zip(RANKS, ('phylum', 'class', 'order', 'family', 'genus')))
DOMAINS = ('Bacteria', 'Archaea', 'Eukaryota', 'unknown', 'unexpected')
TOL = Fraction(1, 10 ** 12)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def read(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def unique(rows, fields, name):
    result = {}
    for row in rows:
        key = tuple(row[f] for f in fields)
        require(key not in result, 'Duplicate ' + name + ': ' + str(key))
        result[key] = row
    return result


def exact(row, expected, context):
    for key, value in expected.items():
        require(row.get(key) == str(value), context + ': ' + key + ' mismatch')


def frac(text, numerator, denominator, context):
    if denominator == 0:
        require(text == 'NA', context + ': zero denominator must be NA')
    else:
        require(text != 'NA' and abs(Fraction(text) - Fraction(numerator, denominator)) <= TOL,
                context + ': fraction mismatch')


def raw_fields(raw):
    text = raw.strip()
    if text in ('', 'Unassigned'):
        return {}, 'unknown', '', 'empty_or_Unassigned'
    fields = {}
    for piece in text.split(';'):
        token = piece.strip()
        if not token:
            continue
        match = re.fullmatch(r'([dpcofgs])__(.*)', token)
        if not match:
            indexed = re.fullmatch(r'D_([0-6])__(.*)', token)
            if indexed:
                rank, value = 'dpcofgs'[int(indexed[1])], indexed[2].strip()
            else:
                return {}, 'unexpected', token, 'unexpected_unprefixed_token'
        else:
            rank, value = match[1], match[2].strip()
        require(rank not in fields, 'Repeated taxonomy rank in ' + raw)
        fields[rank] = value
    domain = fields.get('d', '')
    if not domain:
        return fields, 'unknown', '', 'missing_or_empty_domain'
    return (fields, domain, domain, 'recognized_domain') if domain in DOMAINS[:3] else (
        fields, 'unexpected', domain, 'unexpected_nonempty_domain_token')


def lineage(fields, rank):
    target = RANK_CODE[rank]
    service = {'', 'Unassigned', 'Unknown_Family'}
    if fields.get(target, '') in service:
        return 'UNCLASSIFIED'
    prefix = 'dpcofgs'[:'dpcofgs'.index(target) + 1]
    return ';'.join(code + '__' + fields[code] for code in prefix if fields.get(code, '') not in service)


def audit(out_dir, source_root=ROOT):
    out = Path(out_dir).resolve(); root = Path(source_root).resolve()
    errors = []; checks = 0
    try:
        manifest = json.loads((out / 'taxonomy_barplot_manifest.json').read_text())
        params = manifest['parameters']
        require(manifest['schema'] == 'bacterial-taxonomy-barplots-v1', 'Unexpected manifest schema')
        require(params['domain_policy'] == 'per-model explicit Bacteria', 'Incorrect bacterial policy')
        require(params['organelle_exclusion'] == ['Mitochondria', 'Chloroplast'], 'Incorrect labelled-organelle policy')
        ranks = params['ranks']; top_n = params['top_n']
        require(1 <= top_n <= 20 and isinstance(top_n, int), 'Invalid top_n')
        require(1 <= len(ranks) <= 5 and len(set(ranks)) == len(ranks) and set(ranks) <= set(RANKS), 'Invalid ranks')
        d, issues = load_dataset(root)
        require(not issues, 'Source input violations: ' + json.dumps(issues, ensure_ascii=False))
        assignment_rows = read(root / ASSIGNMENTS)
        assignments = unique(assignment_rows, ['feature_id'], 'source ASV')
        require({k[0] for k in assignments} == set(d['table']), 'Assignment/ASV ID mismatch')
        metadata = unique(read(root / 'data/meta/metadata.tsv'), ['sample-id'], 'source sample')
        roles = {k[0]: v['sample_type'] for k, v in metadata.items()}
        soil = sorted(s for s, role in roles.items() if role == 'soil')
        require(soil == d['soil_samples'] and len(soil) == 10, 'Incorrect soil sample set')
        require(params['soil_samples'] == soil and params['all_samples'] == d['samples'], 'Manifest sample order mismatch')
        require(params['models'] == list(MODELS), 'Manifest model order mismatch')
        expected_sources = {x['path']: x['sha256'] for x in d['inputs'].values()}
        expected_sources[ASSIGNMENTS] = sha(root / ASSIGNMENTS)
        require(manifest['sources'] == expected_sources, 'Scientific source hashes do not match')
        require(manifest['inputs_unchanged'] is True, 'Missing input-preservation declaration')
        for rel, hash_ in manifest['code_sha256'].items():
            require(sha(ROOT / rel) == hash_, 'Code hash changed: ' + rel)
        for rel, hash_ in manifest['outputs'].items():
            require(sha(out / rel) == hash_, 'Artifact hash changed: ' + rel)

        # Reconstruct domain masks directly from source QZA raw Taxon strings.
        actual_mask = unique(read(out / 'bacterial_mask.tsv'), ['feature_id', 'model'], 'mask row')
        require(set(actual_mask) == {(f, m) for f in d['table'] for m in MODELS}, 'Incomplete domain mask')
        domains = {}; parsed = {}; organelles = {}; selected = {m: [] for m in MODELS}
        for fid, counts in d['table'].items():
            source = assignments[fid,]
            soil_reads = sum(counts[s] for s in soil)
            require(source['soil_asv'] == str(int(soil_reads > 0)), 'Source soil_asv mismatch')
            for model in MODELS:
                raw = d['taxonomies'][model][fid]
                require(source[model + '_taxon'] == raw, 'Saved raw Taxon differs from QZA')
                fields, domain, token, reason = raw_fields(raw)
                require(domain != 'unexpected', 'Unexpected nonempty domain token requires review')
                domains[fid, model] = domain; parsed[fid, model] = fields
                found = sorted({value for value in fields.values() if value.casefold() in {'mitochondria', 'chloroplast'}})
                organelles[fid, model] = found
                is_selected = domain == 'Bacteria' and not found
                if is_selected:
                    selected[model].append(fid)
                exact(actual_mask[fid, model], dict(feature_id=fid, model=model, soil_asv=int(soil_reads > 0),
                      domain=domain, domain_token=token, domain_reason=reason,
                      domain_bacteria=int(domain == 'Bacteria'), organelle_flag=int(bool(found)),
                      organelle_tokens_json=json.dumps(found, ensure_ascii=False), selected=int(is_selected),
                      raw_taxon=raw, all_reads=sum(counts.values()), soil_reads=soil_reads,
                      control_reads=sum(counts[s] for s in roles if roles[s] == 'negative_control')), 'Mask ' + fid)
                checks += 1
        losses = unique(read(out / 'domain_by_sample.tsv'), ['model', 'sample'], 'domain sample')
        require(set(losses) == {(m, s) for m in MODELS for s in d['samples']}, 'Incomplete loss audit including control')
        totals = {}
        for model in MODELS:
            for sample in d['samples']:
                counts = {f: c[sample] for f, c in d['table'].items()}
                total = sum(counts.values()); retained = sum(counts[f] for f in selected[model])
                totals[model, sample] = retained
                expected = dict(model=model, sample=sample, sample_type=roles[sample], original_retained_reads=total,
                                original_observed_asv=sum(n > 0 for n in counts.values()), selected_bacterial_reads=retained,
                                excluded_reads=total - retained, zero_bacterial=int(retained == 0))
                for domain in DOMAINS:
                    expected[domain + '_reads'] = sum(n for f, n in counts.items() if domains[f, model] == domain)
                    expected[domain + '_observed_asv'] = sum(n > 0 for f, n in counts.items() if domains[f, model] == domain)
                organelle_reads = sum(n for f, n in counts.items() if domains[f, model] == 'Bacteria' and organelles[f, model])
                organelle_asv = sum(n > 0 for f, n in counts.items() if domains[f, model] == 'Bacteria' and organelles[f, model])
                expected.update(domain_excluded_reads=total - expected['Bacteria_reads'],
                                bacterial_organelle_excluded_reads=organelle_reads,
                                bacterial_organelle_excluded_observed_asv=organelle_asv,
                                selected_bacterial_observed_asv=expected['Bacteria_observed_asv'] - organelle_asv)
                require(retained + organelle_reads == expected['Bacteria_reads'], 'Domain/organelle selection partition mismatch')
                exact(losses[model, sample], expected, 'Domain audit ' + sample + '/' + model)
                frac(losses[model, sample]['retained_fraction'], retained, total, 'Retained fraction')
                frac(losses[model, sample]['excluded_fraction'], total - retained, total, 'Excluded fraction')
                checks += 1
        expected_summary = {}
        for scope, fids in [('all', list(d['table'])), ('soil', [f for f, c in d['table'].items() if sum(c[s] for s in soil) > 0])]:
            expected_summary[scope] = {'asv': len(fids),
                                      'per_model': {m: dict(Counter(domains[f, m] for f in fids)) for m in MODELS},
                                      'bacterial_organelle_excluded_asv': {m: sum(domains[f, m] == 'Bacteria' and bool(organelles[f, m]) for f in fids) for m in MODELS},
                                      'selected_bacterial_asv': {m: sum(domains[f, m] == 'Bacteria' and not organelles[f, m] for f in fids) for m in MODELS},
                                      'common_bacteria_asv': sum(all(domains[f, m] == 'Bacteria' for m in MODELS) for f in fids)}
        require(json.loads((out / 'domain_audit.json').read_text()) == expected_summary, 'ASV domain summary mismatch')
        require(manifest['asv_domain_summary'] == expected_summary, 'Manifest ASV domain summary mismatch')
        nonempty = {m: [s for s in soil if totals[m, s] > 0] for m in MODELS}
        require(manifest['nonempty_bacterial_soil_samples'] == nonempty, 'Nonempty sample set mismatch')
        if params['audit_only']:
            return {'status': 'PASS', 'scope': 'domain audit only', 'checks': checks, 'errors': [], 'warnings': []}

        # Independently parse full lineages from raw Taxon, then aggregate source counts.
        expected_profiles = {}
        profile_counts = {}
        for model in MODELS:
            for rank in ranks:
                for sample in d['samples']:
                    values = Counter({'UNCLASSIFIED': 0})
                    for fid in selected[model]:
                        key = lineage(parsed[fid, model], rank)
                        source = assignments[fid,]
                        require(source[model + '_' + rank + '_key'] == key, 'Source lineage key disagrees with raw Taxon')
                        require(source[model + '_' + rank + '_assigned'] == str(int(key != 'UNCLASSIFIED')), 'Source rank status mismatch')
                        n = d['table'][fid][sample]
                        if n:
                            values[key] += n
                    require(sum(values.values()) == totals[model, sample], 'Bacterial profile mass mismatch')
                    profile_counts[model, sample, rank] = values
                    for key, count in values.items():
                        expected_profiles[model, sample, rank, key] = count
        actual_profiles = unique(read(out / 'bacterial_profiles.tsv'), ['model', 'sample', 'rank', 'key'], 'bacterial profile')
        require(set(actual_profiles) == set(expected_profiles), 'Filtered profile key set mismatch')
        for (model, sample, rank, key), count in expected_profiles.items():
            row = actual_profiles[model, sample, rank, key]
            exact(row, {'count': count, 'total_reads': totals[model, sample], 'sample_type': roles[sample]}, 'Filtered profile')
            frac(row['fraction'], count, totals[model, sample], 'Filtered profile'); checks += 1

        category_rows = unique(read(out / 'taxonomy_barplot_keys.tsv'), ['model', 'rank', 'category'], 'category mapping')
        segment_rows = unique(read(out / 'taxonomy_barplot_data.tsv'), ['model', 'sample', 'rank', 'category'], 'plotted segment')
        expected_categories = set(); expected_segments = set(); rank_colors = {}; rank_labels = {}
        for rank in ranks:
            rank_colors[rank] = {}; rank_labels[rank] = {}
            for model in MODELS:
                valid = nonempty[model]
                keys = {key for s in soil for key in profile_counts[model, s, rank]} - {'UNCLASSIFIED'}
                means = {key: sum((Fraction(profile_counts[model, s, rank].get(key, 0), totals[model, s]) for s in valid), Fraction(0)) / len(valid) for key in keys} if valid else {}
                top = sorted(keys, key=lambda k: (-means[k], k))[:top_n] if valid else []
                other = sorted(keys - set(top))
                definitions = [(key, 'selected', [key]) for key in top]
                definitions += [('OTHER_ASSIGNED', 'other_assigned', other), ('UNCLASSIFIED', 'unclassified', ['UNCLASSIFIED'])]
                for order, (key, kind, members) in enumerate(definitions):
                    ckey = (model, rank, key); expected_categories.add(ckey)
                    require(ckey in category_rows, 'Missing expected selected/aggregate category')
                    row = category_rows[ckey]
                    exact(row, {'category_kind': kind, 'order': order, 'members_count': len(members), 'nonempty_soil_n': len(valid)}, 'Category')
                    require(json.loads(row['full_keys_json']) == members, 'Other/selected lineage membership mismatch')
                    mean = sum((Fraction(sum(profile_counts[model, s, rank].get(k, 0) for k in members), totals[model, s]) for s in valid), Fraction(0)) / len(valid) if valid else None
                    if mean is None:
                        require(row['mean_fraction'] == row['mean_fraction_exact'] == 'NA', 'Empty-model mean must be NA')
                    else:
                        require(Fraction(row['mean_fraction_exact']) == mean and abs(Fraction(row['mean_fraction']) - mean) <= TOL, 'Mean including zero absences mismatch')
                    color = row['color']; label = row['label']
                    require(bool(re.fullmatch(r'#[0-9a-fA-F]{6}', color)), 'Invalid hex color')
                    if kind == 'selected':
                        fields = dict(x.split('__', 1) for x in key.split(';'))
                        require(fields['d'] in label and fields[RANK_CODE[rank]] in label, 'Legend hides domain or terminal taxon name')
                        require(color not in ('#686868', '#dddddd'), 'Assigned key uses an aggregate color')
                        require(key not in rank_colors[rank] or rank_colors[rank][key] == color, 'Same full key has different colors')
                        require(key not in rank_labels[rank] or rank_labels[rank][key] == label, 'Same key has inconsistent labels')
                        rank_colors[rank][key] = color; rank_labels[rank][key] = label
                    else:
                        exact(row, {'color': '#686868' if kind == 'other_assigned' else '#dddddd',
                                    'label': 'Другие назначенные' if kind == 'other_assigned' else 'Без назначения'}, 'Aggregate legend')
                    for sample in soil:
                        skey = (model, sample, rank, key); expected_segments.add(skey)
                        require(skey in segment_rows, 'Missing plotted category, including zero category')
                        count = sum(profile_counts[model, sample, rank].get(k, 0) for k in members)
                        actual = segment_rows[skey]
                        exact(actual, {'count': count, 'total_reads': totals[model, sample], 'category_kind': kind,
                                       'order': order, 'label': label, 'color': color}, 'Segment ' + str(skey))
                        frac(actual['fraction'], count, totals[model, sample], 'Segment fraction'); checks += 1
            require(len(set(rank_colors[rank].values())) == len(rank_colors[rank]), 'Different full keys share colors')
            require(len(set(rank_labels[rank].values())) == len(rank_labels[rank]), 'Different full keys have indistinguishable labels')
        require(set(category_rows) == expected_categories, 'Unexpected category mapping rows')
        require(set(segment_rows) == expected_segments, 'Unexpected plotted samples or categories')
        require(manifest['colors'] == rank_colors, 'Manifest color map mismatch')
        required = {'bacterial_mask.tsv', 'domain_by_sample.tsv', 'domain_audit.json', 'bacterial_profiles.tsv',
                    'taxonomy_barplot_keys.tsv', 'taxonomy_barplot_data.tsv'}
        required |= {'figures/taxonomy_' + FILE_RANK[r] + '.' + ext for r in ranks for ext in ('pdf', 'png', 'svg')}
        require(required <= set(manifest['outputs']), 'Missing generated artifact hashes')
        return {'status': 'PASS', 'scope': 'per-model bacterial barplots', 'checks': checks,
                'ranks': ranks, 'top_n': top_n, 'soil_bars': len(MODELS) * len(soil) * len(ranks),
                'zero_bacterial_soil_samples': {m: [s for s in soil if not totals[m, s]] for m in MODELS},
                'errors': [], 'warnings': ['Model-specific ASV masks and denominators differ; this is not a common-input comparison.']}
    except (ValueError, OSError, KeyError, TypeError, AssertionError, ZeroDivisionError) as exc:
        errors.append(str(exc) or type(exc).__name__)
    return {'status': 'FAIL', 'checks': checks, 'errors': errors, 'warnings': []}


def self_test(out, root):
    before = audit(out, root)
    require(before['status'] == 'PASS', 'Self-test baseline failed: ' + str(before['errors']))
    manifest_name = 'taxonomy_barplot_manifest.json'; data_name = 'taxonomy_barplot_data.tsv'
    originals = {name: (out / name).read_bytes() for name in (manifest_name, data_name)}
    with tempfile.TemporaryDirectory(prefix='bacterial-barplot-selftest-') as tmp:
        copy = Path(tmp) / 'outputs'; shutil.copytree(out, copy)
        values = read(copy / data_name)
        donor = next(r for r in values if int(r['count']) > 0 and int(r['total_reads']) > 0)
        receiver = next(r for r in values if all(r[k] == donor[k] for k in ('model', 'sample', 'rank')) and r['category'] != donor['category'])
        for row, change in ((donor, -1), (receiver, 1)):
            row['count'] = str(int(row['count']) + change)
            row['fraction'] = format(int(row['count']) / int(row['total_reads']), '.17g')
        with (copy / data_name).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(values[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(values)
        manifest = json.loads((copy / manifest_name).read_text())
        manifest['outputs'][data_name] = sha(copy / data_name)
        (copy / manifest_name).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
        bad = audit(copy, root)
        require(bad['status'] == 'FAIL' and any('Segment ' in x and 'count mismatch' in x for x in bad['errors']),
                'Corruption was not caught by semantic source-count comparison')
        for name, content in originals.items():
            (copy / name).write_bytes(content)
        restored = audit(copy, root)
        require(restored['status'] == 'PASS', 'Restoration did not pass')
    require(all((out / name).read_bytes() == content for name, content in originals.items()), 'Original outputs changed')
    return {'status': 'PASS', 'baseline': 'PASS', 'mutated': 'FAIL', 'restored': 'PASS',
            'mutation': 'Moved one read between two categories within one bar, preserving total; recomputed fractions and refreshed the temporary artifact hash.',
            'sample': donor['sample'], 'model': donor['model'], 'rank': donor['rank'],
            'detected': bad['errors'], 'isolation': 'Only a temporary output copy changed; original outputs and scientific inputs were untouched.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--json-out', type=Path)
    args = parser.parse_args()
    try:
        result = self_test(args.out_dir.resolve(), args.source_root.resolve()) if args.self_test else audit(args.out_dir, args.source_root)
    except (ValueError, OSError, KeyError) as exc:
        result = {'status': 'FAIL', 'errors': [str(exc)]}
    if args.json_out:
        args.json_out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result['status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()

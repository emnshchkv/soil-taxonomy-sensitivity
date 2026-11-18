#!/usr/bin/env python3
"""Standalone bacterial barplots from saved assignments and ASV counts.

AI4Science-агент. No denoising, classification, or source-file modification.
Each model selects its own explicit Bacteria assignments, excluding lineages
with exact Mitochondria/Chloroplast tokens at any rank; the selected inputs
and denominators therefore differ between models. This is not a consensus mask.
"""
import argparse
from collections import Counter, defaultdict
import colorsys
import csv
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

from taxonomy_compare_io import load_dataset, MODELS, RANKS

ROOT = Path(__file__).resolve().parent.parent
ASSIGNMENTS = 'results/taxonomy_comparison/finish_20261116/main/assignments.tsv'
RANK_CODES = dict(zip(RANKS, ('p', 'c', 'o', 'f', 'g')))
RANK_EN = dict(zip(('phylum', 'class', 'order', 'family', 'genus'), RANKS))
DOMAIN_VALUES = ('Bacteria', 'Archaea', 'Eukaryota', 'unknown', 'unexpected')
OTHER = 'OTHER_ASSIGNED'
UNKNOWN = 'UNCLASSIFIED'
OTHER_COLOR = '#686868'
UNKNOWN_COLOR = '#dddddd'
PALETTE = ('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
           '#e377c2', '#17becf', '#bcbd22', '#aec7e8', '#ffbb78', '#98df8a',
           '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#9edae5', '#dbdb8d',
           '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173', '#3182bd')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def fraction_string(count, total):
    return format(count / total, '.17g') if total else 'NA'


def get_domain(raw):
    tokens = [x.strip() for x in raw.split(';') if x.strip()]
    if not tokens or tokens == ['Unassigned']:
        return 'unknown', '', 'empty_or_Unassigned'
    values = []
    for token in tokens:
        match = re.fullmatch(r'(?:d__|D_0__)(.*)', token)
        if match:
            values.append(match[1].strip())
        elif not re.fullmatch(r'(?:[pcofgs]__|D_[1-6]__).*', token):
            return 'unexpected', token, 'unexpected_unprefixed_token'
    if len(values) > 1:
        return 'unexpected', ';'.join(values), 'duplicate_domain_tokens'
    if not values or not values[0]:
        return 'unknown', '', 'missing_or_empty_domain'
    value = values[0]
    return (value, value, 'recognized_domain') if value in DOMAIN_VALUES[:3] else (
        'unexpected', value, 'unexpected_nonempty_domain_token')


def get_organelle_tokens(raw):
    found = set()
    for piece in raw.split(';'):
        token = piece.strip()
        match = re.fullmatch(r'(?:[dpcofgs]__|D_[0-6]__)(.*)', token)
        if match:
            value = match[1].strip()
            if value.casefold() in ('mitochondria', 'chloroplast'):
                found.add(value)
    return sorted(found)


def source_data(root):
    dataset, errors = load_dataset(root)
    if errors:
        raise ValueError(json.dumps(errors, ensure_ascii=False))
    rows = read_tsv(root / ASSIGNMENTS)
    assignments = {row['feature_id']: row for row in rows}
    if len(assignments) != len(rows) or set(assignments) != set(dataset['table']):
        raise ValueError('Assignment IDs do not match unique ASV table IDs')
    metadata = read_tsv(root / 'data/meta/metadata.tsv')
    roles = {row['sample-id']: row['sample_type'] for row in metadata}
    if set(roles) != set(dataset['samples']) or len(roles) != len(metadata):
        raise ValueError('Metadata samples are not unique or differ from ASV table')
    soil = dataset['soil_samples']
    if len(soil) != 10:
        raise ValueError('Expected ten soil samples')
    domains = {}; organelles = {}
    mask = []
    for fid, row in assignments.items():
        counts = dataset['table'][fid]
        soil_reads = sum(counts[s] for s in soil)
        if row['soil_asv'] not in ('0', '1') or int(row['soil_asv']) != int(soil_reads > 0):
            raise ValueError('soil_asv disagrees with frequencies: ' + fid)
        for model in MODELS:
            raw = row[model + '_taxon']
            if raw != dataset['taxonomies'][model][fid]:
                raise ValueError('Saved raw Taxon differs from QZA: ' + fid + '/' + model)
            domain, token, reason = get_domain(raw)
            if domain == 'unexpected':
                raise ValueError('Unexpected domain must be reviewed: ' + fid + '/' + raw)
            domains[fid, model] = domain
            organelles[fid, model] = get_organelle_tokens(raw)
            mask.append(dict(feature_id=fid, model=model, soil_asv=row['soil_asv'],
                             domain=domain, domain_token=token, domain_reason=reason,
                             domain_bacteria=int(domain == 'Bacteria'),
                             organelle_flag=int(bool(organelles[fid, model])),
                             organelle_tokens_json=json.dumps(organelles[fid, model], ensure_ascii=False),
                             selected=int(domain == 'Bacteria' and not organelles[fid, model]), raw_taxon=raw,
                             all_reads=sum(counts.values()), soil_reads=soil_reads,
                             control_reads=sum(counts[s] for s in roles if roles[s] == 'negative_control')))
    return dataset, assignments, roles, domains, organelles, mask


def domain_audit(dataset, roles, domains, organelles):
    rows = []
    for model in MODELS:
        for sample in dataset['samples']:
            counts = {fid: values[sample] for fid, values in dataset['table'].items()}
            total = sum(counts.values())
            row = dict(model=model, sample=sample, sample_type=roles[sample],
                       original_retained_reads=total, original_observed_asv=sum(v > 0 for v in counts.values()))
            for domain in DOMAIN_VALUES:
                fids = [fid for fid in counts if domains[fid, model] == domain]
                row[domain + '_reads'] = sum(counts[fid] for fid in fids)
                row[domain + '_observed_asv'] = sum(counts[fid] > 0 for fid in fids)
            organelle_reads = sum(n for fid, n in counts.items()
                                 if domains[fid, model] == 'Bacteria' and organelles[fid, model])
            organelle_asv = sum(n > 0 for fid, n in counts.items()
                               if domains[fid, model] == 'Bacteria' and organelles[fid, model])
            selected = row['Bacteria_reads'] - organelle_reads
            row.update(domain_excluded_reads=total - row['Bacteria_reads'],
                       bacterial_organelle_excluded_reads=organelle_reads,
                       bacterial_organelle_excluded_observed_asv=organelle_asv,
                       selected_bacterial_reads=selected,
                       selected_bacterial_observed_asv=row['Bacteria_observed_asv'] - organelle_asv,
                       excluded_reads=total - selected,
                       retained_fraction=fraction_string(selected, total),
                       excluded_fraction=fraction_string(total - selected, total),
                       zero_bacterial=int(selected == 0))
            assert sum(row[d + '_reads'] for d in DOMAIN_VALUES) == total
            rows.append(row)
    summary = {}
    for label, fids in [('all', list(dataset['table'])),
                        ('soil', [f for f, counts in dataset['table'].items()
                                  if sum(counts[s] for s in dataset['soil_samples']) > 0])]:
        summary[label] = {'asv': len(fids),
                          'per_model': {m: dict(Counter(domains[f, m] for f in fids)) for m in MODELS},
                          'bacterial_organelle_excluded_asv': {m: sum(domains[f, m] == 'Bacteria' and bool(organelles[f, m]) for f in fids) for m in MODELS},
                          'selected_bacterial_asv': {m: sum(domains[f, m] == 'Bacteria' and not organelles[f, m] for f in fids) for m in MODELS},
                          'common_bacteria_asv': sum(all(domains[f, m] == 'Bacteria' for m in MODELS) for f in fids)}
    return rows, summary


def labels_for(keys, rank):
    parsed = {key: dict(x.split('__', 1) for x in key.split(';')) for key in keys}
    bases = {key: parsed[key].get('d', '?') + ': ' + parsed[key][RANK_CODES[rank]] for key in keys}
    labels = dict(bases)
    for base in set(bases.values()):
        collision = [k for k in keys if bases[k] == base]
        if len(collision) < 2:
            continue
        ancestors = list('dpcofg')[:list('dpcofg').index(RANK_CODES[rank])]
        # Add the nearest distinguishing ancestor(s), keeping the taxon name.
        suffix = {key: [] for key in collision}
        for code in reversed(ancestors):
            values = {parsed[key].get(code, '?') for key in collision}
            if len(values) > 1:
                for key in collision:
                    suffix[key].append(code + '__' + parsed[key].get(code, '?'))
                if len({tuple(suffix[k]) for k in collision}) == len(collision):
                    break
        for key in collision:
            labels[key] = base + ' (' + '; '.join(suffix[key]) + ')'
    if len(set(labels.values())) != len(keys):
        raise ValueError('Could not disambiguate displayed labels for ' + rank)
    return labels


def colors_for(keys):
    colors = {}
    for i, key in enumerate(sorted(keys)):
        if i < len(PALETTE):
            color = PALETTE[i]
        else:
            rgb = colorsys.hls_to_rgb((i * 0.61803398875) % 1, .43 + .08 * (i % 2), .7)
            color = '#' + ''.join(f'{round(c * 255):02x}' for c in rgb)
        if color in set(colors.values()) | {OTHER_COLOR, UNKNOWN_COLOR}:
            raise ValueError('Color collision; reduce top-n')
        colors[key] = color
    return colors


def build_profiles(dataset, assignments, domains, organelles, roles, rank_order, top_n):
    soil = dataset['soil_samples']
    totals = {(m, s): sum(counts[s] for fid, counts in dataset['table'].items()
                         if domains[fid, m] == 'Bacteria' and not organelles[fid, m])
              for m in MODELS for s in dataset['samples']}
    grouped = defaultdict(lambda: defaultdict(int))
    for model in MODELS:
        for fid, counts in dataset['table'].items():
            if domains[fid, model] != 'Bacteria' or organelles[fid, model]:
                continue
            for rank in rank_order:
                key = assignments[fid][model + '_' + rank + '_key']
                assigned = assignments[fid][model + '_' + rank + '_assigned']
                if (assigned == '0') != (key == UNKNOWN) or assigned not in ('0', '1'):
                    raise ValueError('Rank key/status disagreement')
                if key != UNKNOWN and not key.startswith('d__Bacteria;'):
                    raise ValueError('Non-bacterial lineage within bacterial mask: ' + key)
                for sample, count in counts.items():
                    if count:
                        grouped[model, sample, rank][key] += count
    profiles = []
    for model in MODELS:
        for sample in dataset['samples']:
            for rank in rank_order:
                count_map = grouped[model, sample, rank]
                count_map.setdefault(UNKNOWN, 0)
                total = totals[model, sample]
                assert sum(count_map.values()) == total
                for key, count in sorted(count_map.items()):
                    profiles.append(dict(model=model, sample=sample, sample_type=roles[sample], rank=rank,
                                         key=key, count=count, total_reads=total,
                                         fraction=fraction_string(count, total)))
    means = {}
    selected = {}
    valid_samples = {}
    for rank in rank_order:
        for model in MODELS:
            valid = [s for s in soil if totals[model, s] > 0]
            valid_samples[model, rank] = valid
            keys = set().union(*(set(grouped[model, s, rank]) for s in soil)) - {UNKNOWN}
            averages = {key: sum((Fraction(grouped[model, s, rank].get(key, 0), totals[model, s])
                                          for s in valid), Fraction(0)) / len(valid) for key in keys} if valid else {}
            means[model, rank] = averages
            selected[model, rank] = sorted(averages, key=lambda k: (-averages[k], k))[:top_n]
    key_rows = []
    segments = []
    color_maps = {}
    for rank in rank_order:
        union = set().union(*(set(selected[m, rank]) for m in MODELS))
        labels = labels_for(union, rank)
        color_maps[rank] = colors_for(union)
        for model in MODELS:
            top = selected[model, rank]
            other_members = sorted(set(means[model, rank]) - set(top))
            category_defs = [(key, 'selected', labels[key], color_maps[rank][key], [key]) for key in top]
            category_defs += [(OTHER, 'other_assigned', 'Другие назначенные', OTHER_COLOR, other_members),
                              (UNKNOWN, 'unclassified', 'Без назначения', UNKNOWN_COLOR, [UNKNOWN])]
            for order, (key, kind, label, color, members) in enumerate(category_defs):
                valid = valid_samples[model, rank]
                mean = sum((sum(grouped[model, s, rank].get(k, 0) for k in members) /
                            Fraction(totals[model, s]) for s in valid), Fraction(0)) / len(valid) if valid else None
                key_rows.append(dict(model=model, rank=rank, category=key, category_kind=kind, order=order,
                                     label=label, color=color, full_keys_json=json.dumps(members, ensure_ascii=False),
                                     members_count=len(members), nonempty_soil_n=len(valid),
                                     mean_fraction='NA' if mean is None else format(float(mean), '.17g'),
                                     mean_fraction_exact='NA' if mean is None else str(mean)))
                for sample in soil:
                    count = sum(grouped[model, sample, rank].get(k, 0) for k in members)
                    segments.append(dict(model=model, sample=sample, rank=rank, category=key,
                                         category_kind=kind, count=count, total_reads=totals[model, sample],
                                         fraction=fraction_string(count, totals[model, sample]), order=order,
                                         label=label, color=color))
    return profiles, key_rows, segments, color_maps, valid_samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--top-n', type=int, default=8)
    parser.add_argument('--ranks', nargs='+', default=list(RANKS), choices=list(RANKS) + list(RANK_EN))
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--force', action='store_true', help='Explicitly permit replacing generated files in an existing output directory')
    args = parser.parse_args()
    if not 1 <= args.top_n <= 20:
        parser.error('--top-n must be between 1 and 20')
    rank_order = [RANK_EN.get(r, r) for r in args.ranks]
    if len(set(rank_order)) != len(rank_order):
        parser.error('Duplicate ranks')
    root = args.source_root.resolve(); out = args.out_dir.resolve()
    if out.exists() and not args.force:
        parser.error('Output already exists; use a new directory or explicitly --force')
    dataset, assignments, roles, domains, organelles, mask = source_data(root)
    source_records = {record['path']: record['sha256'] for record in dataset['inputs'].values()}
    source_records[ASSIGNMENTS] = sha(root / ASSIGNMENTS)
    domain_rows, summary = domain_audit(dataset, roles, domains, organelles)
    result = None if args.audit_only else build_profiles(dataset, assignments, domains, organelles, roles, rank_order, args.top_n)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'bacterial_mask.tsv', mask)
    write_tsv(out / 'domain_by_sample.tsv', domain_rows)
    (out / 'domain_audit.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    generated = [out / name for name in ('bacterial_mask.tsv', 'domain_by_sample.tsv', 'domain_audit.json')]
    color_maps = {}; valid = {}
    if result is not None:
        profiles, key_rows, segments, color_maps, valid = result
        write_tsv(out / 'bacterial_profiles.tsv', profiles)
        write_tsv(out / 'taxonomy_barplot_keys.tsv', key_rows)
        write_tsv(out / 'taxonomy_barplot_data.tsv', segments)
        generated += [out / name for name in ('bacterial_profiles.tsv', 'taxonomy_barplot_keys.tsv', 'taxonomy_barplot_data.tsv')]
        from render_taxonomy_barplots import render_barplots
        renderer_rows = [dict(row, fraction=row['fraction'] if row['fraction'] == 'NA' else float(row['fraction'])) for row in segments]
        generated += render_barplots(renderer_rows, dataset['soil_samples'], list(MODELS), rank_order, out)
    if any(sha(root / rel) != value for rel, value in source_records.items()):
        raise ValueError('Scientific input changed during generation')
    git = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True, check=False)
    code_paths = ['src/plot_taxonomy_profiles.py', 'src/check_taxonomy_barplots.py', 'src/render_taxonomy_barplots.py']
    manifest = {
        'schema': 'bacterial-taxonomy-barplots-v1', 'author': 'AI4Science-агент',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'argv': [sys.executable, *sys.argv], 'working_directory': str(Path.cwd()),
        'code_commit': git.stdout.strip(),
        'code_sha256': {p: sha(ROOT / p) for p in code_paths if (ROOT / p).exists()},
        'sources': source_records, 'inputs_unchanged': True,
        'parameters': {'domain_policy': 'per-model explicit Bacteria',
                       'organelle_exclusion': ['Mitochondria', 'Chloroplast'],
                       'organelle_rule': 'case-insensitive exact complete token after rank prefix and whitespace removal, any rank; each ASV counted once',
                       'top_n': args.top_n,
                       'ranks': rank_order, 'models': list(MODELS), 'soil_samples': dataset['soil_samples'],
                       'all_samples': dataset['samples'], 'audit_only': args.audit_only},
        'selection_rule': 'For each model and rank independently: descending mean fraction over all nonempty bacterial soil samples, absent keys contribute zero; ties use the complete lineage key lexicographically.',
        'normalization': 'Count divided by the per-model selected bacterial reads of that sample after domain and exact labelled-organelle exclusions. Zero selected reads => NA, sample retained. Eukaryota, Archaea and unknown-domain ASVs excluded; rank UNCLASSIFIED inside remaining explicit Bacteria retained. Unlabelled organelles cannot be identified by this rule.',
        'comparison_limit': 'The ASV masks and denominators differ by model; these plots are not a common-input comparison or biological ground truth.',
        'other_assigned': 'Remaining assigned complete keys are grouped for display only; membership differs between model/rank panels and is listed in taxonomy_barplot_keys.tsv.',
        'color_rule': 'Within a rank, identical full keys share a color; different full keys have different colors even when their terminal labels match.',
        'colors': color_maps,
        'nonempty_bacterial_soil_samples': {m: [s for s in dataset['soil_samples'] if next(x for x in domain_rows if x['model'] == m and x['sample'] == s)['selected_bacterial_reads'] > 0] for m in MODELS},
        'asv_domain_summary': summary,
        'outputs': {str(p.relative_to(out)): sha(p) for p in generated},
    }
    (out / 'taxonomy_barplot_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'GENERATED', 'output': str(out), 'audit_only': args.audit_only,
                      'mask_rows': len(mask), 'domain_sample_rows': len(domain_rows),
                      'profiles': 0 if result is None else len(result[0]),
                      'segments': 0 if result is None else len(result[2]),
                      'nonempty_bacterial_soil_samples': manifest['nonempty_bacterial_soil_samples'],
                      'files': len(generated)}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()

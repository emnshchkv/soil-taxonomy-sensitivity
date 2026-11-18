#!/usr/bin/env python3
"""Standalone family profiles of per-model Bacteria and Archaea assignments.

Exclude explicit organelle tokens. Apply the display threshold within each
sample; retain complete family counts and ASV selection evidence separately.
Never overwrite an output directory or change a scientific input.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import colorsys
import csv
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
from taxonomy_compare_io import MODELS, load_dataset, sha256
from taxonomy_compare_core import parse_taxon, lineage_key

ROOT = Path(__file__).resolve().parent.parent
RANK = 'семейство'
UNKNOWN, OTHER = 'UNCLASSIFIED', 'OTHER_ASSIGNED'
UNKNOWN_COLOR, OTHER_COLOR = '#dddddd', '#686868'
DOMAINS = ('Bacteria', 'Archaea', 'Eukaryota', 'unknown')
PREFIXES = dict(enumerate('dpcofgs'))
PALETTE = ('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
           '#e377c2', '#17becf', '#bcbd22', '#aec7e8', '#ffbb78', '#98df8a',
           '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#9edae5', '#dbdb8d',
           '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173', '#3182bd')


def write_tsv(path, rows):
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def fraction_string(count, total):
    return format(count / total, '.17g') if total else 'NA'


def taxon_details(raw):
    # Numeric rank prefixes are purely structural aliases; values stay literal.
    canonical = ';'.join(re.sub(r'^D_([0-6])__', lambda m: PREFIXES[int(m[1])] + '__',
                               part.strip()) for part in raw.split(';'))
    parsed = parse_taxon(canonical, ';')
    raw_family_key = lineage_key(parsed, RANK) or UNKNOWN
    family_values = [m[1].strip() for token in canonical.split(';')
                     if (m := re.fullmatch(r'f__(.*)', token))]
    raw_family_label = family_values[0] if family_values else ''
    if not family_values:
        unknown_family_reason = 'missing'
    elif raw_family_label in ('', 'Unassigned', 'Unknown_Family'):
        unknown_family_reason = 'service'
    elif raw_family_label.casefold() == 'uncultured':
        unknown_family_reason = 'uncultured'
    else:
        unknown_family_reason = ''
    family_key = UNKNOWN if unknown_family_reason else raw_family_key
    domain = parsed.domain or 'unknown'
    if domain not in DOMAINS:
        raise ValueError('Unrecognized domain requires review: ' + raw)
    markers = set()
    for token in canonical.split(';'):
        match = re.fullmatch(r'[dpcofgs]__(.*)', token)
        if match and match[1].strip().lower() in ('mitochondria', 'chloroplast'):
            markers.add(match[1].strip().lower())
    if domain == 'Eukaryota':
        reason = 'eukaryota'
    elif domain == 'unknown':
        reason = 'unknown_domain'
    elif markers:
        reason = 'organelle'
    else:
        reason = 'selected'
    return (domain, sorted(markers), reason, family_key,
            raw_family_label, raw_family_key, unknown_family_reason)


def source_data(root):
    dataset, issues = load_dataset(root)
    if issues:
        raise ValueError(json.dumps(issues, ensure_ascii=False))
    if len(dataset['soil_samples']) != 10:
        raise ValueError('Expected ten soil samples')
    roles = {s: 'soil' if s in dataset['soil_samples'] else 'negative_control'
             for s in dataset['samples']}
    details, mask = {}, []
    for model in MODELS:
        for fid, raw in sorted(dataset['taxonomies'][model].items()):
            record = taxon_details(raw)
            domain, markers, reason, family, raw_label, raw_key, family_reason = record
            counts = dataset['table'][fid]
            soil_reads = sum(counts[s] for s in dataset['soil_samples'])
            details[model, fid] = record
            mask.append(dict(feature_id=fid, model=model, raw_taxon=raw, domain=domain,
                             organelle_markers_json=json.dumps(markers), selected=int(reason == 'selected'),
                             exclusion_reason=reason, family_key=family, raw_family_label=raw_label,
                             raw_family_key=raw_key, unknown_family_reason=family_reason,
                             soil_asv=int(soil_reads > 0),
                             all_reads=sum(counts.values()), soil_reads=soil_reads,
                             control_reads=sum(counts[s] for s in roles if roles[s] == 'negative_control')))
    return dataset, roles, details, mask


def build_profiles(dataset, roles, details, threshold):
    """Full profiles and per-sample threshold decisions, including the control."""
    grouped = defaultdict(lambda: defaultdict(int))
    for (model, fid), record in details.items():
        _, _, reason, family = record[:4]
        if reason == 'selected':
            for sample, count in dataset['table'][fid].items():
                if count:
                    grouped[model, sample][family] += count
    losses, profiles, memberships = [], [], []
    totals, colored = {}, defaultdict(set)
    for model in MODELS:
        for sample in dataset['samples']:
            counts = {f: c[sample] for f, c in dataset['table'].items()}
            total = sum(counts.values())
            domains = {d: sum(c for f, c in counts.items() if details[model, f][0] == d) for d in DOMAINS}
            organelles = {f for f in counts if details[model, f][2] == 'organelle'}
            selected = {f for f in counts if details[model, f][2] == 'selected'}
            selected_reads = sum(counts[f] for f in selected)
            organelle_reads = sum(counts[f] for f in organelles)
            if sum(domains.values()) != total or domains['Bacteria'] + domains['Archaea'] != selected_reads + organelle_reads:
                raise ValueError('Filtering does not conserve counts')
            totals[model, sample] = selected_reads
            losses.append(dict(model=model, sample=sample, sample_type=roles[sample],
                               total_reads=total, total_observed_asv=sum(c > 0 for c in counts.values()),
                               **{d + '_reads': domains[d] for d in DOMAINS},
                               organelle_reads=organelle_reads,
                               organelle_observed_asv=sum(counts[f] > 0 for f in organelles),
                               selected_reads=selected_reads, selected_observed_asv=sum(counts[f] > 0 for f in selected),
                               excluded_reads=total - selected_reads,
                               selected_Bacteria_reads=sum(counts[f] for f in selected if details[model, f][0] == 'Bacteria'),
                               selected_Archaea_reads=sum(counts[f] for f in selected if details[model, f][0] == 'Archaea'),
                               selected_fraction=fraction_string(selected_reads, total), zero_selected=int(selected_reads == 0)))
            count_map = grouped[model, sample]
            count_map.setdefault(UNKNOWN, 0)
            if sum(count_map.values()) != selected_reads:
                raise ValueError('Family profile does not conserve selected reads')
            for key, count in sorted(count_map.items()):
                fraction = fraction_string(count, selected_reads)
                profiles.append(dict(model=model, sample=sample, sample_type=roles[sample], rank=RANK,
                                     key=key, count=count, total_reads=selected_reads, fraction=fraction))
                above = bool(selected_reads and count and Fraction(count, selected_reads) >= threshold)
                category = UNKNOWN if key == UNKNOWN else key if above else OTHER
                memberships.append(dict(model=model, sample=sample, sample_type=roles[sample], key=key,
                                        count=count, total_reads=selected_reads, fraction=fraction,
                                        display_category=category))
                if key != UNKNOWN and above and roles[sample] == 'soil':
                    colored[model].add(key)
    return losses, profiles, memberships, totals, colored


def labels_for(keys):
    parsed = {k: dict(t.split('__', 1) for t in k.split(';')) for k in keys}
    bases = {k: parsed[k]['d'] + ': ' + parsed[k]['f'] for k in keys}
    labels = dict(bases)
    for base in sorted(set(bases.values())):
        collision = [k for k in keys if bases[k] == base]
        if len(collision) < 2:
            continue
        suffix = {k: [] for k in collision}
        for ancestor in 'ocpd':
            if len({parsed[k].get(ancestor, '?') for k in collision}) > 1:
                for key in collision:
                    suffix[key].append(ancestor + '__' + parsed[key].get(ancestor, '?'))
                if len({tuple(suffix[k]) for k in collision}) == len(collision):
                    break
        for key in collision:
            labels[key] += ' (' + '; '.join(suffix[key]) + ')'
    if len(set(labels.values())) != len(keys):
        raise ValueError('Family labels remain ambiguous')
    return labels


def colors_for(keys):
    colors = {}
    for i, key in enumerate(sorted(keys)):
        if i < len(PALETTE):
            color = PALETTE[i]
        else:
            rgb = colorsys.hls_to_rgb((i * .61803398875) % 1, .43 + .08 * (i % 2), .7)
            color = '#' + ''.join(f'{round(c * 255):02x}' for c in rgb)
        if color in set(colors.values()) | {OTHER_COLOR, UNKNOWN_COLOR}:
            raise ValueError('Color collision')
        colors[key] = color
    return colors


def display_tables(dataset, roles, memberships, totals, colored, min_percent):
    all_keys = set().union(*(colored[m] for m in MODELS))
    labels, colors = labels_for(all_keys), colors_for(all_keys)
    membership = defaultdict(list)
    for row in memberships:
        membership[row['model'], row['sample']].append(row)
    key_rows, segments = [], []
    for model in MODELS:
        definitions = [(key, 'selected', labels[key], colors[key], key) for key in sorted(colored[model])]
        definitions += [(OTHER, 'other_assigned', f'Семейства <{min_percent}%', OTHER_COLOR, ''),
                        (UNKNOWN, 'unclassified', 'Семейство не определено', UNKNOWN_COLOR, UNKNOWN)]
        for order, (category, kind, label, color, full_key) in enumerate(definitions):
            key_rows.append(dict(model=model, category=category, category_kind=kind, order=order,
                                 label=label, color=color, full_key=full_key))
            for sample in dataset['soil_samples']:
                count = sum(r['count'] for r in membership[model, sample] if r['display_category'] == category)
                total = totals[model, sample]
                segments.append(dict(model=model, sample=sample, sample_type=roles[sample], rank=RANK,
                                     category=category, category_kind=kind, count=count, total_reads=total,
                                     fraction=fraction_string(count, total), order=order, label=label,
                                     color=color, full_key=full_key))
        for sample in dataset['soil_samples']:
            if sum(r['count'] for r in segments if r['model'] == model and r['sample'] == sample) != totals[model, sample]:
                raise ValueError('Display categories do not conserve selected reads')
    return key_rows, segments, colors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--min-percent', default='2', help='Family threshold within each sample; inclusive, default 2')
    args = parser.parse_args()
    try:
        percent = Fraction(args.min_percent)
        if not 0 < percent <= 100:
            raise ValueError('outside (0, 100]')
    except (ValueError, ZeroDivisionError) as exc:
        parser.error('Invalid --min-percent: ' + str(exc))
    threshold = percent / 100
    label_percent = format(float(percent), '.12g')
    root, out = args.source_root.resolve(), args.out_dir.absolute()
    if out.exists() or out.is_symlink():
        parser.error('Output already exists; choose a new directory')
    dataset, roles, details, mask = source_data(root)
    sources = {r['path']: r['sha256'] for r in dataset['inputs'].values()}
    code_paths = ['src/plot_prokaryotic_families.py', 'src/render_prokaryotic_families.py',
                  'src/check_prokaryotic_families.py', 'src/taxonomy_compare_io.py',
                  'src/taxonomy_compare_core.py', 'src/classification_contract.py']
    code_hashes = {p: sha256(ROOT / p) for p in code_paths}
    losses, profiles, memberships, totals, colored = build_profiles(dataset, roles, details, threshold)
    keys, segments, colors = display_tables(dataset, roles, memberships, totals, colored, label_percent)
    out.mkdir(parents=True, exist_ok=False)
    tables = {'asv_mask.tsv': mask, 'filter_counts.tsv': losses, 'family_profiles.tsv': profiles,
              'family_membership.tsv': memberships, 'family_keys.tsv': keys, 'family_barplot_data.tsv': segments}
    generated = []
    for name, rows in tables.items():
        write_tsv(out / name, rows)
        generated.append(out / name)
    from render_prokaryotic_families import render_family_profiles
    generated += render_family_profiles(segments, dataset['soil_samples'], list(MODELS), out, threshold=threshold)
    if any(sha256(root / p) != digest for p, digest in sources.items()):
        raise ValueError('Scientific input changed during generation')
    if any(sha256(ROOT / p) != digest for p, digest in code_hashes.items()):
        raise ValueError('Generator or checking code changed during generation')
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True, check=False)
    manifest = dict(schema='prokaryotic-family-barplots-v1', author='AI4Science-агент',
                    created_utc=datetime.now(timezone.utc).isoformat(), argv=[sys.executable, *sys.argv],
                    working_directory=str(Path.cwd()), source_root=str(root), code_commit=commit.stdout.strip(),
                    code_sha256=code_hashes, sources=sources, inputs_unchanged=True,
                    parameters=dict(domain_policy='per-model explicit Bacteria or Archaea, excluding explicit organelle tokens',
                                    min_percent=str(percent), threshold_fraction=str(threshold), rank=RANK,
                                    family_unknown_labels=['', 'Unassigned', 'Unknown_Family', 'uncultured'],
                                    family_uncultured_casefold=True,
                                    models=list(MODELS), soil_samples=dataset['soil_samples'], all_samples=dataset['samples']),
                    organelle_rule='Any rank: entire trimmed token value lower() equals mitochondria or chloroplast; union once; excluded within Bacteria or Archaea.',
                    family_rule='Literal lineage through family using taxonomy_compare_core. Empty/Unassigned/Unknown_Family and terminal family uncultured (casefold) give UNCLASSIFIED without excluding ASVs. Other ancestor labels stay literal. Raw family label/key retained in asv_mask. No synonyms or family inference from genus.',
                    display_rule='Each sample separately: count/N >= threshold gives a colored family; smaller assigned families go to OTHER_ASSIGNED in that sample. UNCLASSIFIED remains separate. Zero N => NA.',
                    normalization='N is the selected prokaryotic reads of this model and sample; masks differ across models. These are model assignments, not verified biological identity.',
                    legend_rule='Union of colored family keys in soil samples of each model. Same full key has the same color across models; gray display categories are exceptions.',
                    colors=colors,
                    nonempty_soil_samples={m: [s for s in dataset['soil_samples'] if totals[m, s] > 0] for m in MODELS},
                    outputs={str(p.relative_to(out)): sha256(p) for p in generated})
    (out / 'prokaryotic_family_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(dict(status='GENERATED', output=str(out), mask_rows=len(mask), profile_rows=len(profiles),
                          soil_figures=12, files=len(generated), nonempty_soil_samples=manifest['nonempty_soil_samples']),
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()

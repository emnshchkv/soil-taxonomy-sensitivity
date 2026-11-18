#!/usr/bin/env python3
"""Trace Kr taxonomic disagreements to saved ASVs, raw QZA and frequencies.

AI4Science-агент. Diagnostic summaries only: no sequence upload, classification,
denoising, change to original files, or inference of biological ground truth.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import re
import subprocess
import sys

from classification_contract import Qza, taxonomy
from taxonomy_compare_io import load_dataset, MODELS, RANKS, PAIRS

ROOT = Path(__file__).resolve().parent.parent
OLD = 'results/taxonomy_barplots/bacteria_20261118'
RANK_CODES = dict(zip(RANKS, ('p', 'c', 'o', 'f', 'g')))
DEPTHS = ('unassigned', 'domain', 'phylum', 'class', 'order', 'family')
FULL_DEPTHS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
SERVICE = {'', 'Unassigned', 'Unknown_Family'}
PLACEHOLDER = re.compile(r'^(?:uncultured|unidentified|unclassified|unknown|environmental|metagenom)', re.I)
TOL = Fraction(1, 10 ** 12)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    require(bool(rows), 'No rows to write: ' + str(path))
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def pct(n, d):
    return format(100 * n / d, '.12g') if d else 'NA'


def frac_ok(actual, n, d):
    return actual == 'NA' if not d else actual != 'NA' and abs(Fraction(actual) - Fraction(n, d)) <= TOL


def parse(raw):
    fields = {}
    if raw.strip() not in ('', 'Unassigned'):
        for token in raw.split(';'):
            token = token.strip()
            if not token:
                continue
            match = re.fullmatch(r'([dpcofgs])__(.*)', token)
            if not match:
                match2 = re.fullmatch(r'D_([0-6])__(.*)', token)
                require(match2 is not None, 'Unexpected taxonomy token: ' + token)
                code, value = 'dpcofgs'[int(match2[1])], match2[2].strip()
            else:
                code, value = match[1], match[2].strip()
            require(code not in fields, 'Duplicated taxonomy rank in ' + raw)
            fields[code] = value
    domain = fields.get('d', '')
    domain = 'unknown' if domain in SERVICE else domain
    require(domain in ('Bacteria', 'Archaea', 'Eukaryota', 'unknown'), 'Unexpected domain ' + domain)
    assigned = {code: fields.get(code, '') not in SERVICE for code in 'dpcofgs'}
    present = [i for i, code in enumerate('dpcofgs') if assigned[code]]
    upto = [i for i in present if i <= 4]
    organelles = sorted({value for value in fields.values() if value.casefold() in ('mitochondria', 'chloroplast')})
    family = fields.get('f', '')
    old = domain == 'Bacteria' and not organelles
    new = domain in ('Bacteria', 'Archaea') and not organelles
    return {'fields': fields, 'domain': domain, 'organelle_tokens': organelles,
            'deepest_upto_family': DEPTHS[max(upto) + 1] if upto else 'unassigned',
            'deepest_full': FULL_DEPTHS[max(present)] if present else 'unassigned',
            'gaps_before_deepest': [FULL_DEPTHS[i] for i in range(max(present)) if not assigned['dpcofgs'[i]]] if present else [],
            'family_literal': family, 'family_assigned': assigned['f'],
            'family_placeholder': bool(family and PLACEHOLDER.search(family)),
            'old_selected': bool(old), 'new_selected': bool(new)}


def key_at(info, rank):
    target = RANK_CODES[rank]
    fields = info['fields']
    if fields.get(target, '') in SERVICE:
        return 'UNCLASSIFIED'
    return ';'.join(code + '__' + fields[code] for code in 'dpcofgs'[:'dpcofgs'.index(target) + 1]
                    if fields.get(code, '') not in SERVICE)


def reason(info, new=False):
    if info['domain'] == 'unknown':
        return 'excluded_unknown_domain'
    if info['domain'] == 'Eukaryota':
        return 'excluded_Eukaryota'
    if info['domain'] == 'Archaea' and not new:
        return 'excluded_Archaea'
    if info['organelle_tokens']:
        return 'excluded_labelled_organelle'
    if not info['family_assigned']:
        return 'selected_family_UNCLASSIFIED'
    return 'selected_family_literal_placeholder' if info['family_placeholder'] else 'selected_family_named_literal'


def independent_old_check(root, dataset, info):
    """Recompute old selection and profile/segment counts directly from QZA data."""
    base = root / OLD
    masks = read_tsv(base / 'bacterial_mask.tsv')
    seen = set(); checks = 0
    for row in masks:
        fid, model = row['feature_id'], row['model']
        require((fid, model) not in seen, 'Duplicate old mask row'); seen.add((fid, model))
        expected = info[fid, model]
        require(row['raw_taxon'] == dataset['taxonomies'][model][fid], 'Old raw taxonomy changed')
        require(row['domain'] == expected['domain'], 'Old domain mismatch')
        require(int(row['selected']) == int(expected['old_selected']), 'Old selection mismatch')
        require(int(row['organelle_flag']) == bool(expected['organelle_tokens']), 'Old organelle flag mismatch')
        checks += 4
    require(seen == {(f, m) for f in dataset['table'] for m in MODELS}, 'Old mask incomplete')
    count_maps = defaultdict(Counter); totals = {}
    for model in MODELS:
        chosen = [f for f in dataset['table'] if info[f, model]['old_selected']]
        for sample in dataset['samples']:
            totals[model, sample] = sum(dataset['table'][f][sample] for f in chosen)
            for rank in RANKS:
                counts = count_maps[model, sample, rank]
                counts['UNCLASSIFIED'] = 0
                for fid in chosen:
                    n = dataset['table'][fid][sample]
                    if n:
                        counts[key_at(info[fid, model], rank)] += n
                require(sum(counts.values()) == totals[model, sample], 'Recomputed old profile mass mismatch')
    expected_keys = {(m, s, r, k) for (m, s, r), counts in count_maps.items() for k in counts}
    seen = set()
    for row in read_tsv(base / 'bacterial_profiles.tsv'):
        m, s, r, k = (row[x] for x in ('model', 'sample', 'rank', 'key'))
        require((m, s, r, k) not in seen, 'Duplicate old profile row'); seen.add((m, s, r, k))
        require(int(row['count']) == count_maps[m, s, r][k], 'Old filtered profile count mismatch')
        require(int(row['total_reads']) == totals[m, s], 'Old profile denominator mismatch')
        require(frac_ok(row['fraction'], int(row['count']), totals[m, s]), 'Old profile fraction mismatch')
        checks += 3
    require(seen == expected_keys, 'Old profile keys mismatch')
    categories = {}
    for row in read_tsv(base / 'taxonomy_barplot_keys.tsv'):
        triple = row['model'], row['rank'], row['category']
        require(triple not in categories, 'Duplicate old category row')
        categories[triple] = json.loads(row['full_keys_json'])
    seen = set(); bar_mass = Counter()
    for row in read_tsv(base / 'taxonomy_barplot_data.tsv'):
        m, s, r, c = (row[x] for x in ('model', 'sample', 'rank', 'category'))
        require(s in dataset['soil_samples'], 'Control or unexpected sample in old soil plot')
        require((m, s, r, c) not in seen, 'Duplicate old segment'); seen.add((m, s, r, c))
        expected = sum(count_maps[m, s, r].get(k, 0) for k in categories[m, r, c])
        require(int(row['count']) == expected, 'Old plotted count mismatch')
        require(int(row['total_reads']) == totals[m, s], 'Old plotted denominator mismatch')
        require(frac_ok(row['fraction'], expected, totals[m, s]), 'Old plotted fraction/NA mismatch')
        bar_mass[m, s, r] += expected; checks += 3
    require(seen == {(m, s, r, c) for m, r, c in categories for s in dataset['soil_samples']}, 'Old plotted segment set incomplete')
    require(all(n == totals[m, s] for (m, s, r), n in bar_mass.items()), 'Old plotted mass not conserved')
    return {'status': 'PASS', 'checks': checks, 'mask_rows': len(masks), 'profiles': len(expected_keys),
            'bars': len(bar_mass), 'method': 'Independent raw-QZA domain/organelle/rank parsing and ASV-frequency aggregation; no old generator import.'}


def run(root, out):
    require(not out.exists(), 'Output already exists; use a new directory')
    d, errors = load_dataset(root)
    require(not errors, 'Invalid source dataset: ' + json.dumps(errors, ensure_ascii=False))
    kr = [s for s in d['samples'] if s.startswith('Kr')]
    lub = [s for s in d['samples'] if s.startswith('LUB')]
    require(len(kr) == 5 and len(lub) == 5, 'Expected five Kr and five LUB soils')
    metadata = {r['sample-id']: r['sample_type'] for r in read_tsv(root / 'data/meta/metadata.tsv')}
    raw_conf = {}; tax = {}; info = {}
    source_hashes = {r['path']: r['sha256'] for r in d['inputs'].values()}
    for model in MODELS:
        q = Qza(root / f'work/qza/taxonomy_{model}.qza', 'FeatureData[Taxonomy]')
        tax[model] = taxonomy(q)
        raw_rows = list(csv.DictReader(io.StringIO(q.text('data/taxonomy.tsv')), delimiter='\t'))
        for row in raw_rows:
            fid = row['Feature ID']
            raw_conf[fid, model] = row['Confidence']
            require(float(row['Confidence']) == tax[model][fid][1], 'Confidence parser mismatch')
            require(row['Taxon'] == d['taxonomies'][model][fid], 'Raw Taxon mismatch')
            info[fid, model] = parse(row['Taxon'])
    kr_ids = sorted(f for f, counts in d['table'].items() if sum(counts[s] for s in kr) > 0)
    summaries = []
    reasons = []
    reason_names = ('excluded_unknown_domain', 'excluded_Eukaryota', 'excluded_Archaea', 'excluded_labelled_organelle',
                    'selected_family_UNCLASSIFIED', 'selected_family_literal_placeholder', 'selected_family_named_literal')
    for model in MODELS:
        for sample in d['samples']:
            counts = {f: c[sample] for f, c in d['table'].items()}
            total = sum(counts.values())
            old_n = sum(n for f, n in counts.items() if info[f, model]['old_selected'])
            new_n = sum(n for f, n in counts.items() if info[f, model]['new_selected'])
            row = dict(model=model, sample=sample, sample_type=metadata[sample], total_reads=total,
                       observed_asv=sum(n > 0 for n in counts.values()))
            for domain in ('Bacteria', 'Archaea', 'Eukaryota', 'unknown'):
                row['domain_' + domain + '_reads'] = sum(n for f, n in counts.items() if info[f, model]['domain'] == domain)
            row.update(old_selected_reads=old_n, new_selected_reads=new_n, returned_Archaea_reads=new_n - old_n,
                       old_zero=int(old_n == 0), new_zero=int(new_n == 0))
            for version, selected_n in [('old', old_n), ('new', new_n)]:
                family_unknown = sum(n for f, n in counts.items() if info[f, model][version + '_selected'] and not info[f, model]['family_assigned'])
                placeholder = sum(n for f, n in counts.items() if info[f, model][version + '_selected'] and info[f, model]['family_assigned'] and info[f, model]['family_placeholder'])
                row[version + '_family_UNCLASSIFIED_reads'] = family_unknown
                row[version + '_family_literal_placeholder_reads'] = placeholder
                row[version + '_family_named_literal_reads'] = selected_n - family_unknown - placeholder
                row[version + '_family_UNCLASSIFIED_pct_selected'] = pct(family_unknown, selected_n)
                row[version + '_selected_pct_original'] = pct(selected_n, total)
                counts_by_reason = Counter()
                asv_by_reason = Counter()
                for f, n in counts.items():
                    if n:
                        why = reason(info[f, model], new=version == 'new')
                        counts_by_reason[why] += n; asv_by_reason[why] += 1
                require(sum(counts_by_reason.values()) == total, 'Filter reason partition fails')
                for why in reason_names:
                    reasons.append(dict(model=model, sample=sample, sample_type=metadata[sample], filter=version,
                                        reason=why, reads=counts_by_reason[why], asv=asv_by_reason[why], total_reads=total,
                                        pct_original=pct(counts_by_reason[why], total)))
            summaries.append(row)
    all_asv = []
    for fid in kr_ids:
        counts = d['table'][fid]
        row = dict(feature_id=fid, kr_reads_total=sum(counts[s] for s in kr), **{s: counts[s] for s in kr},
                   K_S63=counts['K_S63'], LUB_reads_total=sum(counts[s] for s in lub),
                   LUB_present_samples=';'.join(s for s in lub if counts[s]), sequence=d['sequences'][fid])
        for model in MODELS:
            parsed = info[fid, model]
            values = {'Taxon': tax[model][fid][0], 'Confidence_raw': raw_conf[fid, model],
                      'domain': parsed['domain'], 'deepest_upto_family': parsed['deepest_upto_family'],
                      'deepest_full': parsed['deepest_full'], 'gaps_before_deepest': json.dumps(parsed['gaps_before_deepest']),
                      'family_literal': parsed['family_literal'], 'family_key': key_at(parsed, 'семейство'),
                      'family_assigned_by_contract': int(parsed['family_assigned']),
                      'family_placeholder_flag': int(parsed['family_placeholder']),
                      'organelle_tokens': json.dumps(parsed['organelle_tokens']),
                      'old_selected': int(parsed['old_selected']), 'new_selected': int(parsed['new_selected']),
                      'old_reason': reason(parsed), 'new_reason': reason(parsed, new=True)}
            row.update({model + '_' + key: value for key, value in values.items()})
        all_asv.append(row)
    by_fid = {row['feature_id']: row for row in all_asv}
    top_rows = []; coverage = {}
    for sample in kr:
        observed = sorted((f for f in kr_ids if d['table'][f][sample]), key=lambda f: (-d['table'][f][sample], f))
        total = sum(d['table'][f][sample] for f in observed); cumulative = 0; selected_count = 0
        for order, fid in enumerate(observed, 1):
            n = d['table'][fid][sample]; cumulative += n; selected_count += 1
            top_rows.append(dict(sample=sample, top_order=order, reads=n, total_reads=total, pct_original=pct(n, total),
                                 cumulative_reads=cumulative, cumulative_pct=pct(cumulative, total), **by_fid[fid]))
            if cumulative * 10 >= total * 9 and order >= min(5, len(observed)):
                break
        coverage[sample] = dict(top_asv=selected_count, reads=cumulative, total=total, pct=pct(cumulative, total))
        require(cumulative * 10 >= total * 9 and selected_count >= min(5, len(observed)), 'Top-ASV coverage rule failed')
    depth_rows = []; transitions = []; family_tokens = []
    for sample in kr:
        total = sum(c[sample] for c in d['table'].values())
        for model in MODELS:
            depth_counter = Counter(); depth_asv = Counter(); old_depth = Counter(); new_depth = Counter()
            family_counts = defaultdict(lambda: [0, 0, 0, 0])
            for fid in kr_ids:
                n = d['table'][fid][sample]
                if not n:
                    continue
                p = info[fid, model]; depth = p['deepest_upto_family']
                depth_counter[depth] += n; depth_asv[depth] += 1
                old_depth[depth] += n * p['old_selected']; new_depth[depth] += n * p['new_selected']
                k = (p['family_literal'] or '<absent>', int(p['family_assigned']), int(p['family_placeholder']))
                values = family_counts[k]; values[0] += n; values[1] += 1
                values[2] += n * p['old_selected']; values[3] += n * p['new_selected']
            require(sum(depth_counter.values()) == total, 'Classification-depth mass mismatch')
            for depth in DEPTHS:
                depth_rows.append(dict(model=model, sample=sample, deepest_upto_family=depth, reads=depth_counter[depth],
                                       asv=depth_asv[depth], total_reads=total, pct_original=pct(depth_counter[depth], total),
                                       old_selected_reads=old_depth[depth], new_selected_reads=new_depth[depth]))
            for (token, assigned, placeholder), values in sorted(family_counts.items()):
                family_tokens.append(dict(model=model, sample=sample, family_literal=token,
                                          assigned_by_current_contract=assigned, placeholder_flag=placeholder,
                                          reads=values[0], asv=values[1], old_selected_reads=values[2],
                                          new_selected_reads=values[3], pct_original=pct(values[0], total)))
        for a, b in PAIRS:
            grouped = defaultdict(lambda: [0, 0, []])
            for fid in kr_ids:
                n = d['table'][fid][sample]
                if not n:
                    continue
                p, q = info[fid, a], info[fid, b]
                key = p['domain'], p['deepest_upto_family'], q['domain'], q['deepest_upto_family']
                grouped[key][0] += n; grouped[key][1] += 1; grouped[key][2].append(fid)
            require(sum(x[0] for x in grouped.values()) == total, 'Cross-model transition mass mismatch')
            for key, (n, asv, ids) in sorted(grouped.items(), key=lambda x: (-x[1][0], x[0])):
                transitions.append(dict(sample=sample, model_a=a, model_b=b, domain_a=key[0], deepest_a=key[1],
                                        domain_b=key[2], deepest_b=key[3], reads=n, asv=asv, total_reads=total,
                                        pct_original=pct(n, total), feature_ids=';'.join(sorted(ids))))
    old_paths = sorted(p for p in (root / OLD).rglob('*') if p.is_file())
    old_hashes = {str(p.relative_to(root)): sha(p) for p in old_paths}
    independent = independent_old_check(root, d, info)
    out.mkdir(parents=True); logs = out / 'logs'; logs.mkdir()
    checker_argv = [sys.executable, str(root / 'src/check_taxonomy_barplots.py'), '--out-dir', str(root / OLD)]
    started = datetime.now(timezone.utc).isoformat()
    with (logs / 'old_checker.stdout.txt').open('w') as stdout, (logs / 'old_checker.stderr.txt').open('w') as stderr:
        check = subprocess.run(checker_argv, cwd=root, stdout=stdout, stderr=stderr, check=False)
    (logs / 'old_checker.command.json').write_text(json.dumps({'argv': checker_argv, 'start_utc': started,
        'end_utc': datetime.now(timezone.utc).isoformat(), 'exit_code': check.returncode,
        'stdout': 'old_checker.stdout.txt', 'stderr': 'old_checker.stderr.txt'}, ensure_ascii=False, indent=2) + '\n')
    require(check.returncode == 0, 'Existing old barplot checker did not pass; see logs')
    outputs = {'asv_diagnostics.tsv': all_asv, 'top_asv_by_kr.tsv': top_rows,
               'sample_filter_summary.tsv': summaries, 'filter_reasons.tsv': reasons,
               'classification_depth.tsv': depth_rows, 'cross_model_transitions.tsv': transitions,
               'family_literal_tokens.tsv': family_tokens}
    for name, rows in outputs.items():
        write_tsv(out / name, rows)
    implementation = Path(sys.prefix) / f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages/q2_feature_classifier/_skl.py'
    code_text = implementation.read_text()
    require('accepted_cum_prob' in code_text and '1.0 - cum_prob' in code_text, 'Confidence implementation context changed')
    confidence_note = {'package_version': importlib.metadata.version('q2-feature-classifier'),
                       'local_file': str(implementation), 'sha256': sha(implementation), 'lines': [107, 138],
                       'semantics': 'The implementation sums class probabilities for each successive node and stops when node probability is below confidence. For an accepted prefix, saved Confidence is the cumulative probability of the last accepted rank; for Unassigned it is 1-cum_prob. A d__Bacteria Confidence is not a family-level confidence.'}
    require(all(sha(root / rel) == value for rel, value in source_hashes.items()), 'Scientific input changed')
    require(all(sha(root / rel) == value for rel, value in old_hashes.items()), 'Old barplot artifact changed')
    verification = {'status': 'PASS', 'kr_asv': len(kr_ids), 'original_samples': len(d['samples']),
                    'top_asv_coverage': coverage, 'old_checker_exit': check.returncode,
                    'independent_old_check': independent, 'scientific_inputs_unchanged': True,
                    'old_outputs_unchanged': True, 'mass_checks': 'ASV/sample totals, old profiles and segments, reason partitions, depth partitions and every pairwise transition conserve original frequencies.',
                    'limitations': ['Assignments are model outputs, not biological ground truth.',
                                    'Family depth follows the existing literal contract; placeholder labels are flagged separately and do not alter profiles.',
                                    'Known exact organelle labels can be excluded; unlabelled organelles cannot be identified by this rule.',
                                    'Sequences were read locally and not uploaded or reclassified.']}
    (out / 'verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2) + '\n')
    manifest = {'author': 'AI4Science-агент', 'argv': [sys.executable, *sys.argv],
                'created_utc': datetime.now(timezone.utc).isoformat(), 'script_sha256': sha(Path(__file__)),
                'sources': source_hashes, 'old_outputs_before': old_hashes,
                'confidence_semantics': confidence_note,
                'depth_rule': 'Deepest explicit non-service rank up to family, independently of missing ancestors; full deepest rank and gaps also recorded.',
                'placeholder_rule': 'Report-only case-insensitive prefix flag: uncultured, unidentified, unclassified, unknown, environmental, metagenom. Does not change the original assigned-rank contract.',
                'top_rule': 'For each Kr sample: descending read count, ties by ASV ID; stop after >=90% original retained reads and at least min(5,observed_ASV) rows.',
                'output_sha256': {name: sha(out / name) for name in [*outputs, 'verification.json']}}
    (out / 'audit_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    kr_summaries = [r for r in summaries if r['sample'] in kr]
    compact = [{k: r[k] for k in ('model', 'sample', 'total_reads', 'old_selected_reads', 'new_selected_reads',
               'returned_Archaea_reads', 'old_family_UNCLASSIFIED_pct_selected', 'new_family_UNCLASSIFIED_pct_selected',
               'new_family_literal_placeholder_reads')} for r in kr_summaries]
    return {'status': 'PASS', 'output': str(out), 'kr_asv': len(kr_ids), 'top_coverage': coverage,
            'samples': compact, 'old_filter_verification': independent,
            'placeholder_tokens': sorted({r['family_literal'] for r in family_tokens if r['placeholder_flag']})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    args = parser.parse_args(); started = datetime.now(timezone.utc).isoformat()
    out = args.out_dir.resolve()
    try:
        result = run(args.source_root.resolve(), out); exit_code = 0
    except (ValueError, OSError, KeyError) as exc:
        result = {'status': 'FAIL', 'errors': [str(exc)]}; exit_code = 1
    text = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if (out / 'logs').is_dir():
        (out / 'logs/audit.stdout.txt').write_text(text)
        (out / 'logs/audit.command.json').write_text(json.dumps({'argv': [sys.executable, *sys.argv],
            'start_utc': started, 'end_utc': datetime.now(timezone.utc).isoformat(), 'exit_code': exit_code}, ensure_ascii=False, indent=2) + '\n')
    print(text, end=''); raise SystemExit(exit_code)


if __name__ == '__main__':
    main()

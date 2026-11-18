#!/usr/bin/env python3
"""Isolated CLI mutation tests for check_paper_results; never edit source results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def snapshot(paths):
    return {str(p): {'sha256': sha(p), 'mode': p.stat().st_mode, 'size': p.stat().st_size} for p in paths}


def edit_tsv(folder, name, change):
    path = folder / name
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        columns, rows = reader.fieldnames, list(reader)
    change(rows)
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(stream.getvalue(), encoding='utf-8')


def first(rows, **selector):
    return next(r for r in rows if all(r[k] == v for k, v in selector.items()))


def patch_text(folder, name, before, after):
    path = folder / name
    text = path.read_text()
    if before not in text:
        raise ValueError(f'test precondition: {before!r} absent from {name}')
    path.write_text(text.replace(before, after, 1))


def refresh(folder):
    path = folder / 'manifest.json'
    manifest = json.loads(path.read_text())
    manifest['outputs'] = {p: sha(folder / p) for p in manifest['outputs']}
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')


def mutation_cases():
    def wrong_macro(folder):
        patch_text(folder, 'numbers.tex', r'\newcommand{\ProjectSamplesAll}{11}', r'\newcommand{\ProjectSamplesAll}{12}')

    def wrong_model_selector(folder):
        def change(rows):
            row = first(rows, id='CmpSilvaFamilyAsvPct')
            selector = json.loads(row['row_key'])
            selector['model'] = 'sklearn_gtdb'
            row['row_key'] = json.dumps(selector, ensure_ascii=False, sort_keys=True, separators=(',', ':'))
        edit_tsv(folder, 'numbers.tsv', change)

    def wrong_source_selector(folder):
        def change(rows):
            row = first(rows, id='CmpSilvaFamilyAsvPct')
            row['source_file'] = 'results/taxonomy_barplots/prokaryote_family_20261118_final/filter_counts.tsv'
        edit_tsv(folder, 'numbers.tsv', change)

    def wrong_selected_denominator(folder):
        def change(rows):
            row = first(rows, model='sklearn_silva138', sample='Kr1_S80')
            row['unclassified_pct_of_selected'] = str(Decimal(row['unclassified_family_reads']) * 100 / Decimal(row['post_dada2_reads']))
        edit_tsv(folder, 'family_sample_counts.tsv', change)

    def wrong_conditional_denominator(folder):
        def change(rows):
            row = first(rows, id='CmpSilvaGtdbGenusConditionalAsvPct')
            row['column'] = '["both_match","total_asv"]'
            row['denominator'] = '3309'
            row['raw_value'] = str(Decimal(row['numerator']) * 100 / Decimal(row['denominator']))
        edit_tsv(folder, 'numbers.tsv', change)

    def na_to_zero(folder):
        def change(rows):
            row = next(r for r in rows if r['zero_selected'] == '1')
            if row['unclassified_pct_of_selected'] != 'NA':
                raise ValueError('test precondition: N=0 must give NA')
            row['unclassified_pct_of_selected'] = '0'
        edit_tsv(folder, 'family_sample_counts.tsv', change)

    def one_read_transfer(folder):
        def change(rows):
            row = next(r for r in rows if int(r['named_family_reads']) > 1 and int(r['unclassified_family_reads']) > 0)
            row['named_family_reads'] = str(int(row['named_family_reads']) - 1)
            row['unclassified_family_reads'] = str(int(row['unclassified_family_reads']) + 1)
            row['unclassified_pct_of_selected'] = str(Decimal(row['unclassified_family_reads']) * 100 / Decimal(row['selected_reads']))
        edit_tsv(folder, 'family_sample_counts.tsv', change)

    def wrong_count(folder):
        def change(rows):
            row = first(rows, id='ProjectSamplesSoil')
            row['raw_value'] = row['display_value'] = '11'
            row['n'] = '11'
        edit_tsv(folder, 'numbers.tsv', change)
        patch_text(folder, 'numbers.tex', r'\newcommand{\ProjectSamplesSoil}{10}', r'\newcommand{\ProjectSamplesSoil}{11}')

    def wrong_soil_role(folder):
        def change(rows):
            first(rows, model='sklearn_silva138', sample='K_S63')['sample_type'] = 'soil'
        edit_tsv(folder, 'family_sample_counts.tsv', change)

    def wrong_soil_composition(folder):
        def change(rows):
            row = first(rows, id='CmpSilvaFamilyReadMedianPct')
            selector = json.loads(row['row_key'])
            selector['sample'][0] = 'K_S63'
            row['row_key'] = json.dumps(selector, ensure_ascii=False, sort_keys=True, separators=(',', ':'))
        edit_tsv(folder, 'numbers.tsv', change)

    def missing_soil_row(folder):
        def change(rows):
            rows.remove(first(rows, model='sklearn_silva138', sample='Kr1_S80'))
        edit_tsv(folder, 'family_sample_counts.tsv', change)

    def wrong_comparison_labels(folder):
        patch_text(folder, 'comparison_main.tex', ' & SILVA & GTDB & GG2 &', ' & GTDB & SILVA & GG2 &')

    def wrong_family_labels(folder):
        patch_text(folder, 'family_main.tex', r'\multicolumn{2}{c}{SILVA} & \multicolumn{2}{c}{GTDB}', r'\multicolumn{2}{c}{GTDB} & \multicolumn{2}{c}{SILVA}')

    def wrong_manuscript_status(folder):
        path = folder / 'manifest.json'
        data = json.loads(path.read_text())
        data['manuscript_status'] = 'checked'
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')

    return {fn.__name__: fn for fn in (wrong_macro, wrong_model_selector, wrong_source_selector,
        wrong_selected_denominator, wrong_conditional_denominator, na_to_zero, one_read_transfer,
        wrong_count, wrong_soil_role, wrong_soil_composition, wrong_comparison_labels,
        wrong_family_labels, wrong_manuscript_status, missing_soil_row)}


def run_suite(out_dir, source_root=ROOT):
    out, root = Path(out_dir).resolve(), Path(source_root).resolve()
    manifest = json.loads((out / 'manifest.json').read_text())
    protected = list(out.rglob('*'))
    protected = [p for p in protected if p.is_file()]
    protected += [root / p for group in ('sources', 'code') for p in manifest[group]]
    before = snapshot(sorted(set(protected)))
    evidence = []

    def check(name, folder, wanted):
        command = [sys.executable, str(ROOT / 'src/check_paper_results.py'), '--out-dir', str(folder), '--source-root', str(root)]
        run = subprocess.run(command, capture_output=True, text=True, cwd=root)
        try:
            report = json.loads(run.stdout)
        except ValueError:
            report = {'status': 'INVALID_CHECKER_OUTPUT', 'stdout': run.stdout}
        ok = run.returncode == wanted and report.get('status') == ('PASS' if wanted == 0 else 'FAIL')
        if wanted == 1:
            ok = ok and bool(report.get('issues')) and any('SHA256' not in s for s in report.get('issues', []))
        evidence.append({'case': name, 'command': command, 'expected_exit_code': wanted,
                         'actual_exit_code': run.returncode, 'accepted': ok, 'result': report, 'stderr': run.stderr})
        return ok

    with tempfile.TemporaryDirectory(prefix='paper-results-mutations-') as temp:
        copy = Path(temp) / 'bundle'
        shutil.copytree(out, copy)
        baseline_ok = check('baseline', copy, 0)
        if baseline_ok:
            pristine = {p.relative_to(copy): p.read_bytes() for p in copy.rglob('*') if p.is_file()}
            with localcontext() as ctx:
                ctx.prec = 70
                for name, mutate in mutation_cases().items():
                    for path, content in pristine.items():
                        (copy / path).write_bytes(content)
                    mutate(copy)
                    refresh(copy)  # Every mutated asset has a valid temporary output SHA.
                    check(name, copy, 1)
            for path, content in pristine.items():
                (copy / path).write_bytes(content)
            check('restored', copy, 0)
    unchanged = before == snapshot(sorted(set(protected)))
    return {'status': 'PASS' if evidence and all(case['accepted'] for case in evidence) and unchanged else 'FAIL',
            'case_count': len(evidence), 'cases': evidence, 'original_files_unchanged': unchanged,
            'protected_files': len(before), 'isolation': 'Mutations only in TemporaryDirectory; source files remain read-only.',
            'hash_policy': 'Temporary output hashes refreshed after every mutation; semantic rejection required.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', required=True, type=Path)
    parser.add_argument('--source-root', default=ROOT, type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    try:
        result = run_suite(args.out_dir, args.source_root)
    except Exception as exc:
        result = {'status': 'FAIL', 'error': f'{type(exc).__name__}: {exc}'}
    text = json.dumps(result, ensure_ascii=False, indent=2) + '\n'
    if args.report:
        with args.report.open('x', encoding='utf-8') as stream:
            stream.write(text)
    print(text, end='')
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())

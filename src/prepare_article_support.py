"""Copy-edit table headings and export additional saved result scalars.

No scientific processing. New output directory only; deterministic outputs.
Run in qiime2-amplicon-2025.7 from the project root.
"""
import argparse
import csv
import hashlib
import json
from decimal import Decimal, ROUND_HALF_EVEN
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)
    bundle = ROOT / 'results/paper_step4/20261118'
    substitutions = {
        'comparison_main.tex': {'Rank &': 'Ранг &', 'Read coverage': 'Покрытие ридов', 'Phylum &': 'Тип &', 'Class &': 'Класс &', 'Order &': 'Порядок &', 'Family &': 'Семейство &', 'Genus &': 'Род &'},
        'family_main.tex': {'Sample & Reads &': 'Образец & Риды &'},
    }
    sources = {}
    for name, mapping in substitutions.items():
        p = bundle / name
        value = p.read_text()
        for old, new in mapping.items():
            assert old in value
            value = value.replace(old, new)
        (args.out_dir / name).write_text(value)
        sources[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    records = []
    def add(macro, value, path, locator, rounding):
        records.append(dict(macro=macro,value=str(value),source=path,locator=locator,rounding=rounding))
        sources[path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    path = 'results/paper_step4/validation_20261118/numeric_verification.json'
    add('ArticleTransferChecks', json.loads((ROOT/path).read_text())['checks'], path, 'checks', 'integer')
    path = 'results/paper_step4/validation_20261118/numeric_tests.json'
    tests = json.loads((ROOT/path).read_text())
    add('ArticleTransferCases', tests['case_count'], path, 'case_count', 'integer')
    add('ArticleDetectedMutations', sum(c['expected_exit_code']==1 and c['actual_exit_code']==1 for c in tests['cases']), path, 'cases: expected_exit_code=actual_exit_code=1; count', 'integer')
    path = 'results/paper_step4/20261118/comparison_pairs.tsv'
    rows = list(csv.DictReader((ROOT/path).open(), delimiter='\t'))
    row = next(r for r in rows if (r['model_a'],r['model_b'],r['rank'])==('sklearn_silva138','sklearn_gtdb','семейство'))
    for field,suffix in [('distance_min','Min'),('distance_max','Max')]:
        value = format(Decimal(row[field]).quantize(Decimal('.000001'),rounding=ROUND_HALF_EVEN),'f')
        add('ArticleSilvaGtdbFamilyDistance'+suffix,value,path,'sklearn_silva138;sklearn_gtdb;семейство;'+field,'HALF_EVEN;6 decimals')
    (args.out_dir/'extra_numbers.tex').write_text('% Additional saved scalars; no scientific recomputation.\n'+''.join('\\newcommand{\\'+r['macro']+'}{'+r['value']+'}\n' for r in records))
    with (args.out_dir/'extra_numbers.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['macro','value','source','locator','rounding'],delimiter='\t');writer.writeheader();writer.writerows(records)
    (args.out_dir/'manifest.json').write_text(json.dumps({'sources':sources,'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.out_dir.iterdir())}},ensure_ascii=False,indent=2)+'\n')
    print('Exported translated headings and five scalars from saved results.')

if __name__ == '__main__':
    main()

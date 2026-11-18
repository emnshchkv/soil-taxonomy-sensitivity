"""CLI rejection tests in isolated manuscript copies; scientific sources read-only.

Run after the complete draft exists. Each case gets a fresh copy. The test never
edits a source reached through an external TeX input or a symlink.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT/'manuscript/ai4science_ysda2026'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory():
    paths = list(ARTICLE.rglob('*'))
    for relative in ('results/paper_step4/20261118', 'results/article_step5/20261118'):
        paths.extend((ROOT/relative).glob('*'))
    return {str(p.relative_to(ROOT)):digest(p) for p in paths if p.is_file()}


def replace_one(path, old, new):
    text = path.read_text()
    if old not in text:
        raise ValueError(f'mutation anchor missing in {path.name}: {old!r}')
    path.write_text(text.replace(old,new,1))


def source_with(article, needle):
    found = [p for p in article.rglob('*.tex') if needle in p.read_text()]
    if len(found) != 1:
        raise ValueError(f'expected one local source with {needle!r}, found {found}')
    return found[0]


def mutate_table(article):
    # Only the relocated temporary copy is edited; relative source inputs in the
    # rest of this case still point to the protected original exports.
    source = ROOT/'results/article_step5/20261118/family_main.tex'
    temporary = article/'mutated_family_main.tex'
    shutil.copyfile(source,temporary)
    replace_one(temporary, '12351', '12352')
    needle = '../../results/article_step5/20261118/family_main.tex'
    replace_one(source_with(article,needle),needle,'mutated_family_main.tex')


def mutate_claim(article):
    path = article/'claims.tsv'
    with path.open(newline='') as f:
        reader = csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
    rows[0]['source']='results/nonexistent-claim-evidence.tsv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,help='new JSON report path')
    args=parser.parse_args()
    before=inventory()
    cases=[]
    specs=[
        ('baseline',0,None,True),
        ('final_baseline',0,None,False),
        ('numeric_macro_override',1,lambda a:replace_one(a/'paper.tex','\\begin{document}','\\renewcommand{\\ProjectSamplesSoil}{999}\n\\begin{document}'),True),
        ('literal_result_number',1,lambda a:replace_one(a/'sections/problem.tex','\\ProjectSamplesSoil{}','999'),True),
        ('wrong_existing_macro',1,lambda a:replace_one(a/'sections/problem.tex','\\ProjectSamplesSoil{}','\\ProjectSamplesAll{}'),True),
        ('unknown_citation',1,lambda a:replace_one(source_with(a,'\\cite{bokulich2018}'),'\\cite{bokulich2018}','\\cite{invented2099}'),True),
        ('missing_mandatory_section',1,lambda a:replace_one(source_with(a,'\\section{Ограничения}'),'\\section{Ограничения}',''),True),
        ('changed_actual_table_cell',1,mutate_table,True),
        ('overview_band_gap_or_overlap',1,lambda a:replace_one(source_with(a,'766.03bp'),'766.03bp','765.03bp'),True),
        ('unmapped_short_prose',1,lambda a:replace_one(a/'sections/problem.tex','\\label{sec:problem}','\\label{sec:problem}\n\nПокрытие составляет 999 процентов.\n'),True),
        ('missing_claim_source',1,mutate_claim,True),
        ('running_title_warning',1,lambda a:(a/'paper.log').write_bytes(b'Title exceeds size limitations for running head.\n'+(a/'paper.log').read_bytes()),False),
        ('missing_pdf_final',2,lambda a:(a/'paper.pdf').unlink(missing_ok=True),False),
        ('restored',0,None,True),
    ]
    intended = {
        'numeric_macro_override':'numeric macro overridden',
        'literal_result_number':'literal result numbers outside verified macros',
        'wrong_existing_macro':'macro list does not match actual claim',
        'unknown_citation':'unknown citations',
        'missing_mandatory_section':'eight required sections are missing or out of order',
        'changed_actual_table_cell':'actual numeric table differs',
        'overview_band_gap_or_overlap':'overview bands leave a gap, overlap, or omit source content',
        'unmapped_short_prose':'unmapped prose paragraph',
        'missing_claim_source':'unsupported/missing source',
        'running_title_warning':'running title exceeds template size limit',
        'missing_pdf_final':'missing build output: paper.pdf',
    }
    with tempfile.TemporaryDirectory(prefix='article-rejection-') as tmp:
        for name,expected,mutation,draft in specs:
            case_root=Path(tmp)/name
            article=case_root/'manuscript/ai4science_ysda2026'
            shutil.copytree(ARTICLE,article)
            # These are read-only logical inputs: mutations are strictly limited
            # to `article`, never a path through these symlinks.
            for directory in ('results','figures'):
                (case_root/directory).symlink_to(ROOT/directory,target_is_directory=True)
            if mutation:
                mutation(article)
            command=[sys.executable,str(ROOT/'src/check_article.py'),'--root',str(ROOT),'--article-dir',str(article)]
            if draft:
                command.append('--draft')
            process=subprocess.run(command,capture_output=True,text=True)
            try:
                result=json.loads(process.stdout)
            except json.JSONDecodeError:
                result={'invalid_output':process.stdout[-3000:]}
            messages=result.get('errors',[])+result.get('incomplete',[])
            intended_issue=name not in intended or any(intended[name] in m for m in messages)
            cases.append({'case':name,'expected_exit_code':expected,'actual_exit_code':process.returncode,
                          'detected_intended_issue':intended_issue,
                          'accepted':process.returncode==expected and intended_issue,'result':result,'stderr':process.stderr})
    after=inventory()
    preserved=before==after
    report={'status':'PASS' if preserved and all(c['accepted'] for c in cases) else 'FAIL',
            'case_count':len(cases),'cases':cases,'source_files_preserved':preserved,
            'preserved_file_count':len(before),'changed_paths':sorted(k for k in set(before)|set(after) if before.get(k)!=after.get(k))}
    encoded=json.dumps(report,ensure_ascii=False,indent=2)+'\n'
    if args.report:
        with args.report.open('x') as f:f.write(encoded)
    print(encoded,end='')
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':
    sys.exit(main())

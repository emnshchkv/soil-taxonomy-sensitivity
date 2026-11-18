"""CLI rejection tests in isolated manuscript copies; scientific sources read-only.

Run after the complete draft exists. Each case gets a fresh copy. The test never
edits a source reached through an external TeX input or a symlink.
--draft-only runs source checks while explicitly reporting skipped PDF cases.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT/'paper/ai4science_ysda2026'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory():
    paths = list(ARTICLE.rglob('*'))
    for relative in ('results/paper_step4/20261118', 'results/article_step5/20261118', 'results/article_revision/20261118'):
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


def crop_model_figure(article):
    path=source_with(article,'prokaryotic_family_silva.pdf')
    replace_one(path,'\\includegraphics[','\\includegraphics[trim={1bp 0 0 0},clip,')


def duplicate_model_figure(article):
    replace_one(article/'paper.tex','\\end{document}',
                '\\includegraphics[width=\\textwidth]{../../figures/taxonomy/family/prokaryotic_family_silva.pdf}\n\\end{document}')


def mutate_revision_table(article):
    source=ROOT/'results/article_revision/20261118/kr_filter_table.tex'
    target=article/source.name
    shutil.copyfile(source,target)
    # Even a wrong sample identifier is an error in the diagnostic table.
    replace_one(target,'Kr1','Kr9')
    needle='../../results/article_revision/20261118/kr_filter_table.tex'
    replace_one(source_with(article,needle),needle,target.name)


def mutate_claim(article, claim_id=None):
    path = article/'claims.tsv'
    with path.open(newline='') as f:
        reader = csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
    row = rows[0] if claim_id is None else next(r for r in rows if r['claim_id'] == claim_id)
    row['source']='results/nonexistent-claim-evidence.tsv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)


def mark_written(article, relative='sections/problem.tex', optional=False):
    """Change declared authorship only in a synthetic temporary manuscript."""
    path=article/relative
    text=path.read_text()
    if '\\porozhdeno{' not in text:
        raise ValueError(f'authorship fixture anchor missing in {relative}')
    mark='\\napisano[автор]{' if optional else '\\napisano{'
    path.write_text(text.replace('\\porozhdeno{',mark))


def mark_all_written(article):
    # This tests validation of an entirely human-declared manuscript. It does
    # not reclassify the real article's prose or claim to establish authorship.
    for path in article.rglob('*.tex'):
        text=path.read_text()
        if '\\porozhdeno{' in text:
            path.write_text(text.replace('\\porozhdeno{','\\napisano{'))


def conflicting_authorship(article):
    path=article/'sections/problem.tex'
    path.write_text('\\napisano{'+path.read_text()+'}\n')


def submission_notice_cases():
    """Exercise extracted PDF-text notice rules without altering a real PDF."""
    spec=importlib.util.spec_from_file_location('article_notice_checker',ROOT/'src/check_article.py')
    checker=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    legacy='Локальный обезличенный черновик AI4Science-YSDA 2026. Статья не подана.'
    ready='Обезличенная версия для подачи на AI4Science-YSDA 2026.'
    examples=[
        ('notice_local_draft',legacy,0,None),
        ('notice_ready_to_submit',ready,0,None),
        ('notice_ready_wrapped',ready.replace(' ','\n'),0,None),
        ('notice_missing','Анонимные авторы',1,'PDF must identify'),
        ('notice_draft_without_status','Локальный обезличенный черновик AI4Science-YSDA 2026.',1,'PDF must identify'),
        ('notice_false_submission_ru',ready+' Работа подана.',1,'incorrectly claims'),
        ('notice_false_submission_article',ready+' Статья была подана.',1,'incorrectly claims'),
        ('notice_false_submission_en',ready+' Submitted to AI4Science-YSDA 2026.',1,'incorrectly claims'),
        ('notice_conflicting_status',legacy+' Работа подана.',1,'incorrectly claims'),
    ]
    cases=[]
    for name,text,expected,intended in examples:
        audit=checker.Audit()
        checker.check_submission_notice(audit,text)
        actual=1 if audit.errors else 0
        detected=intended is None or any(intended in error for error in audit.errors)
        cases.append({'case':name,'kind':'extracted-PDF-text notice unit case',
                      'expected_exit_code':expected,'actual_exit_code':actual,
                      'detected_intended_issue':detected,'accepted':actual==expected and detected,
                      'result':{'errors':audit.errors},'stderr':''})
    return cases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,help='new JSON report path')
    parser.add_argument('--draft-only',action='store_true',help='run source-only cases and explicitly skip final PDF checks')
    args=parser.parse_args()
    before=inventory()
    cases=[]
    specs=[
        ('baseline',0,None,True),
        ('final_baseline',0,None,False),
        ('written_mark',0,mark_written,True),
        ('written_optional_argument',0,lambda a:mark_written(a,optional=True),True),
        ('entirely_written',0,mark_all_written,True),
        ('written_literal_result_number',1,lambda a:(mark_written(a),replace_one(a/'sections/problem.tex','\\ProjectSamplesSoil{}','999')),True),
        ('written_unknown_citation',1,lambda a:(mark_written(a,'sections/related.tex'),replace_one(a/'sections/related.tex','\\cite{bokulich2018}','\\cite{invented2099}')),True),
        ('written_missing_claim_source',1,lambda a:(mark_written(a),mutate_claim(a,'P01')),True),
        ('missing_authorship_mark',1,lambda a:replace_one(a/'sections/problem.tex','\\porozhdeno{\\Claim{P01}','\\textit{\\Claim{P01}'),True),
        ('empty_written_mark',1,lambda a:replace_one(a/'sections/problem.tex','\\porozhdeno{\\Claim{P01}','\\napisano{}\\textit{\\Claim{P01}'),True),
        ('source_mark_without_authorship',1,lambda a:replace_one(a/'sections/related.tex','\\porozhdeno{\\Claim{L01}','\\istochnik{\\Claim{L01}'),True),
        ('conflicting_authorship_marks',1,conflicting_authorship,True),
        ('numeric_macro_override',1,lambda a:replace_one(a/'paper.tex','\\begin{document}','\\renewcommand{\\ProjectSamplesSoil}{999}\n\\begin{document}'),True),
        ('literal_result_number',1,lambda a:replace_one(a/'sections/problem.tex','\\ProjectSamplesSoil{}','999'),True),
        ('wrong_existing_macro',1,lambda a:replace_one(a/'sections/problem.tex','\\ProjectSamplesSoil{}','\\ProjectSamplesAll{}'),True),
        ('unknown_citation',1,lambda a:replace_one(source_with(a,'\\cite{bokulich2018}'),'\\cite{bokulich2018}','\\cite{invented2099}'),True),
        ('missing_mandatory_section',1,lambda a:replace_one(source_with(a,'\\section{Ограничения}'),'\\section{Ограничения}',''),True),
        ('changed_actual_table_cell',1,mutate_table,True),
        ('cropped_model_figure',1,crop_model_figure,True),
        ('duplicated_model_figure',1,duplicate_model_figure,True),
        ('changed_revision_table',1,mutate_revision_table,True),
        ('wrong_model_identity',1,lambda a:replace_one(a/'models.tex','gtdb_classifier_r220.qza','gtdb_classifier_r999.qza'),True),
        ('unmapped_short_prose',1,lambda a:replace_one(a/'sections/problem.tex','\\label{sec:problem}','\\label{sec:problem}\n\nПокрытие составляет 999 процентов.\n'),True),
        ('missing_claim_source',1,mutate_claim,True),
        ('running_title_warning',1,lambda a:(a/'paper.log').write_bytes(b'Title exceeds size limitations for running head.\n'+(a/'paper.log').read_bytes()),False),
        ('missing_pdf_final',2,lambda a:(a/'paper.pdf').unlink(missing_ok=True),False),
        ('restored',0,None,True),
    ]
    intended = {
        'written_literal_result_number':'literal result numbers outside verified macros',
        'written_unknown_citation':'unknown citations',
        'written_missing_claim_source':'P01: unsupported/missing source',
        'missing_authorship_mark':'P01: not enclosed in authorship mark',
        'empty_written_mark':'P01: not enclosed in authorship mark',
        'source_mark_without_authorship':'L01: not enclosed in authorship mark',
        'conflicting_authorship_marks':'conflicting authorship marks',
        'numeric_macro_override':'numeric macro overridden',
        'literal_result_number':'literal result numbers outside verified macros',
        'wrong_existing_macro':'macro list does not match actual claim',
        'unknown_citation':'unknown citations',
        'missing_mandatory_section':'eight required sections are missing or out of order',
        'changed_actual_table_cell':'actual numeric table differs',
        'cropped_model_figure':'single-page figure must preserve its full source',
        'duplicated_model_figure':'model PDF must appear exactly once',
        'changed_revision_table':'actual revision table differs',
        'wrong_model_identity':'model filename missing/different',
        'unmapped_short_prose':'unmapped prose paragraph',
        'missing_claim_source':'unsupported/missing source',
        'running_title_warning':'running title exceeds template size limit',
        'missing_pdf_final':'missing build output: paper.pdf',
    }
    skipped_cases=[name for name,_,_,draft in specs if args.draft_only and not draft]
    if args.draft_only:
        specs=[spec for spec in specs if spec[3]]
    with tempfile.TemporaryDirectory(prefix='article-rejection-') as tmp:
        for name,expected,mutation,draft in specs:
            case_root=Path(tmp)/name
            article=case_root/'paper/ai4science_ysda2026'
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
    if not args.draft_only:
        cases.extend(submission_notice_cases())
    after=inventory()
    preserved=before==after
    report={'status':'PASS' if preserved and all(c['accepted'] for c in cases) else 'FAIL',
            'scope':'source-only CLI regressions; final PDF and notice unit cases not checked' if args.draft_only else 'source/final-PDF CLI regressions and extracted-PDF-text notice unit cases',
            'draft_only':args.draft_only,'skipped_cases':skipped_cases,
            'case_count':len(cases),'cases':cases,'source_files_preserved':preserved,
            'preserved_file_count':len(before),'changed_paths':sorted(k for k in set(before)|set(after) if before.get(k)!=after.get(k))}
    encoded=json.dumps(report,ensure_ascii=False,indent=2)+'\n'
    if args.report:
        with args.report.open('x') as f:f.write(encoded)
    print(encoded,end='')
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':
    sys.exit(main())

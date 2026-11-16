"""Check the saved report, comparison graph and preserved protocol/history."""
from pathlib import Path
import re
import subprocess
import sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'src'))
from qiime_runtime import selected_python
selected_python()
import yaml
from report_taxonomy_comparison import render, DEFAULT_OUTPUT
from complete_taxonomy_comparison import protection


def main():
    issues=[]
    if DEFAULT_OUTPUT.read_text()!=render():issues.append('Report differs from saved evidence')
    original=subprocess.check_output(['git','show','e6dec7b:METHOD_COMPARISON.md'],cwd=ROOT)
    if not (ROOT/'METHOD_COMPARISON.md').read_bytes().startswith(original):issues.append('Original method changed')
    if subprocess.check_output(['git','rev-parse','hw2^{commit}'],cwd=ROOT,text=True).strip()!='26841271b01926c56204eb76043beb3c04c8c6fa':issues.append('hw2 changed')
    protected=protection(ROOT/'results/taxonomy_comparison/finish_20261116/protected_before.json')
    issues+=protected['violations']
    schema=yaml.safe_load((ROOT/'ontology.yaml').read_text())
    pages={}
    by_name={}
    edges=[]
    for folder in ('notes','tasks','reports','ontology'):
        for p in (ROOT/folder).glob('*.md'):
            text=p.read_text()
            meta=yaml.safe_load(text.split('---',2)[1]) if text.startswith('---') else {}
            pages[p]=(meta or {},text)
            by_name.setdefault(p.stem,[]).append(p)
            for relation,target in re.findall(r'^- ([^:\n]+): \[\[([^\]]+)\]\]',text,re.M):
                edges.append((p,relation,target))
    titles=['Метод буквального сравнения таксономий','Набор готовых таксономий для сравнения',
      'Проверка расчёта сравнения на синтетических данных','Результаты сравнения таксономических профилей',
      'Доказательства воспроизведения сравнения','Сравнение готовых таксономий по фиксированным рангам',
      'Сверка сравнения с чистой копией','Сравнение таксономических профилей',
      'Завершить этап сравнения в отдельной ветке','Зафиксировать договор сравнения и исходное состояние',
      'Согласовать расчёт сравнения и независимую проверку','Выполнить и воспроизвести сравнение таксономий',
      'Оформить проверенные результаты сравнения','Провести сравнение таксономических профилей',
      'Зафиксировать методику сравнения таксономий','Реализовать расчёт и проверку сравнения таксономий',
      'Провести сравнение и воспроизведение таксономий','Оформить отчёт и завершить этап сравнения',
      'Завершить исправления и этап 4 сравнения таксономий']
    def resolve(target):
        p=ROOT/target
        if p.exists():return p
        if not p.suffix and p.with_suffix('.md').is_file():return p.with_suffix('.md')
        if len(by_name.get(target,[]))==1:return by_name[target][0]
        return None
    for title in titles:
        matches=by_name.get(title,[])
        if len(matches)!=1:issues.append('Missing/ambiguous page '+title);continue
        p=matches[0];meta,text=pages[p];typ=meta.get('type')
        if typ not in schema['classes']:issues.append('Unknown type '+title)
        if f'- тип: [[{typ}]]' not in text:issues.append('Missing type edge '+title)
        for target in re.findall(r'\[\[([^\]]+)\]\]',text):
            if resolve(target) is None:issues.append('Missing link '+title+' -> '+target)
        for origin,relation,target in edges:
            if origin!=p or relation=='тип':continue
            definition=schema['properties'].get(relation)
            if not definition:issues.append('Unknown relation '+relation);continue
            domain=definition['domain'];domain=[domain] if isinstance(domain,str) else domain
            if typ not in domain:issues.append('Wrong domain '+title+'/'+relation)
            dest=resolve(target)
            dest_type=pages.get(dest,({},''))[0].get('type','Файл')
            if dest is not None and dest_type!=definition['range']:issues.append('Wrong range '+title+'/'+relation+'/'+target)
        if typ=='Задача':
            if meta.get('status')!='выполнена' or '## Итог' not in text:issues.append('Unfinished task '+title)
            incoming=[(a,r) for a,r,t in edges if t==title and r in ('решает','выполняет','подзадачаОт')]
            if not incoming:issues.append('Task lacks completion edge '+title)
            for child,r in incoming:
                if r=='подзадачаОт' and pages[child][0].get('status') not in ('выполнена','принята'):
                    issues.append('Unfinished child '+child.stem)
    for issue in issues:print(issue)
    print(f'{"FAIL/1" if issues else "PASS/0"} — report, {len(titles)} graph pages, protocol, hw2, {protected["checked"]} protected files')
    return 1 if issues else 0

if __name__=='__main__':raise SystemExit(main())

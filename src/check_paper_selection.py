import csv,json,re
from pathlib import Path
root=Path.cwd()
read=lambda p:list(csv.DictReader((root/p).open(),delimiter='\t'))
folder='results/taxonomy_barplots/prokaryote_family_20261118_final/'
fc={(r['model'],r['sample']):r for r in read(folder+'filter_counts.tsv')}
models=['sklearn_silva138','sklearn_gtdb','sklearn_gg2']
kr=['Kr1_S80','Kr2_S81','Kr3_S82','Kr4_S83','Kr5_S84']
expected=[[38,0,2,249,0],[10142,8013,9147,11518,1413],[10114,8098,8974,12013,1465]]
for m,counts in zip(models,expected):
    assert [int(fc[m,s]['selected_reads']) for s in kr]==counts
assert [int(fc[models[0],s]['total_reads']) for s in kr]==[12351,8904,9985,12293,1505]
bar={(r['model'],r['sample'],r['category']):r for r in read(folder+'family_barplot_data.tsv')}
for m,u,n in zip(models,[246,11515,12010],[249,11518,12013]):
    r=bar[m,'Kr4_S83','UNCLASSIFIED']; assert int(r['count'])==u and int(r['total_reads'])==n and n-u==3
for m in models[1:]:
    for s in [kr[i] for i in [0,1,2,4]]:
        r=bar[m,s,'UNCLASSIFIED']; assert r['count']==r['total_reads']
for m,u,n in zip(models,[2496,5273,5468],[19084,19090,19086]):
    r=bar[m,'LUB5_S95','UNCLASSIFIED']; assert int(r['count'])==u and int(r['total_reads'])==n
r=fc[models[0],'Kr2_S81']; assert [int(r[c]) for c in ['Eukaryota_reads','unknown_reads','organelle_reads']]==[7780,1113,11]
r=next(r for r in read('results/taxonomy_barplots/kr_diagnostics_20261118/asv_diagnostics.tsv') if r['feature_id']=='121fc77a2e0cf578c1369ab32511cda3')
assert (r['Kr2_S81'],r['Kr3_S82'])==('1002','1138')
assert 'Didymellaceae' in r[models[0]+'_Taxon'] and 'Ascochyta' in r[models[0]+'_Taxon']
assert all(r[m+'_Taxon']=='d__Bacteria' for m in models[1:])
links=[]
for doc in ['research/results_selection.md','research/paper_results_contract.md']:
    for target in re.findall(r'\[\[([^\]]+)\]\]',(root/doc).read_text()):
        assert (root/target).exists(),target
        links.append(target)
result={'status':'PASS','source_claim_groups':5,'links_checked':len(links),'source_scope':'selection documents and existing source rows; no manuscript','documents':['research/results_selection.md','research/paper_results_contract.md']}
p=root/'results/paper_step4/preflight_20261118/selection_validation.json';p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False))

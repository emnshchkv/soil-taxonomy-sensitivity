"""Isolated CLI corruption cases for independent article-diagnostic verification."""
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

ROOT=Path(__file__).resolve().parents[1]
BUNDLE=ROOT/"results/article_revision/20261118"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table_change(directory, name, field, value):
    path=directory/name
    with path.open() as f:
        reader=csv.DictReader(f,delimiter="\t")
        fields=reader.fieldnames
        rows=list(reader)
    rows[0][field]=value
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def text_change(directory, name, old, new):
    path=directory/name
    value=path.read_text()
    if old not in value:
        raise ValueError("mutation anchor missing: "+old)
    path.write_text(value.replace(old,new,1))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report",type=Path)
    args=parser.parse_args()
    before={p.name:sha(p) for p in BUNDLE.iterdir() if p.is_file()}
    specs=[
        ("baseline",None,None),
        ("wrong_control_overlap",lambda p:table_change(p,"control_summary.tsv","shared_soil_asv","13"),"control_summary.tsv: value differs"),
        ("wrong_sample_denominator",lambda p:table_change(p,"genus_conditional_per_sample.tsv","both_assigned_asv","999"),"genus_conditional_per_sample.tsv: value differs"),
        ("wrong_pooled_denominator",lambda p:table_change(p,"genus_conditional.tsv","pooled_both_assigned_asv","999"),"genus_conditional.tsv: value differs"),
        ("wrong_kr_partition",lambda p:table_change(p,"kr_filter_reasons.tsv","Eukaryota_reads","999"),"kr_filter_reasons.tsv: value differs"),
        ("wrong_numeric_macro",lambda p:table_change(p,"numbers.tsv","value","999"),"revision scalar value/provenance differs"),
        ("wrong_tex_table",lambda p:text_change(p,"kr_filter_table.tex","12351","12352"),"Kr TeX table differs"),
        ("wrong_quality_threshold",lambda p:table_change(p,"method_parameters.tsv","value","29"),"method_parameters.tsv: value differs"),
        ("restored",None,None),
    ]
    cases=[]
    with tempfile.TemporaryDirectory(prefix="article-revision-tests-") as tmp:
        for name,mutate,intended in specs:
            copy=Path(tmp)/name
            copy.mkdir()
            for path in BUNDLE.iterdir():
                if path.is_file():
                    shutil.copy2(path,copy/path.name)
            if mutate:
                mutate(copy)
                # Let forged manifest hashes agree with changed files: the case
                # must fail the independent value checks, not merely integrity.
                manifest=json.loads((copy/"manifest.json").read_text())
                manifest["outputs"]={n:sha(copy/n) for n in manifest["outputs"]}
                (copy/"manifest.json").write_text(json.dumps(manifest,ensure_ascii=False))
            proc=subprocess.run([sys.executable,str(ROOT/"src/check_article_revision.py"),"--out-dir",str(copy)],
                                capture_output=True,text=True)
            try:
                result=json.loads(proc.stdout)
            except json.JSONDecodeError:
                result={"errors":[proc.stdout[-2000:],proc.stderr[-2000:]]}
            expected=0 if mutate is None else 1
            detected=intended is None or any(intended in err for err in result.get("errors",[]))
            cases.append(dict(case=name,expected_exit_code=expected,actual_exit_code=proc.returncode,
                              detected_intended_issue=detected,accepted=proc.returncode==expected and detected,
                              result=result,stderr=proc.stderr))
    after={p.name:sha(p) for p in BUNDLE.iterdir() if p.is_file()}
    preserved=before==after
    report=dict(status="PASS" if preserved and all(c["accepted"] for c in cases) else "FAIL",
                scope="Temporary diagnostic copies; frozen scientific sources and manuscript untouched.",
                cases=cases,case_count=len(cases),source_files_preserved=preserved)
    text=json.dumps(report,ensure_ascii=False,indent=2)+"\n"
    if args.report:
        with args.report.open("x") as f:
            f.write(text)
    print(text,end="")
    return 0 if report["status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())

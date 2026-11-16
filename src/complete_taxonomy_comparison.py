"""One main comparison plus clean local reproduction; never invokes DADA2/classification.

This command refuses existing main/reproduction outputs. All subprocess streams
and true exit codes are preserved. Real QZA copies remain temporary and private.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'src'))
from qiime_runtime import selected_python
selected_python()
from classification_logging import record_command, write_json
from taxonomy_compare_io import code_identity, load_dataset, sha256, INPUT_FILES
from taxonomy_compare_acceptance import compare_tables


def protection(baseline):
    original = json.loads(baseline.read_text())
    # These are authorized stage4 documents or unrelated instruction text.
    excluded = {'METHOD_COMPARISON.md', 'diary.md', 'AGENT_TAXONOMY_COMPARISON_PROMPT.txt'}
    records, violations = {}, []
    for rel, expected in original.items():
        if rel in excluded:
            continue
        p = ROOT/rel
        if not p.is_file():
            violations.append(rel+': missing')
            continue
        actual = {'sha256': sha256(p), 'mode': stat.S_IMODE(p.stat().st_mode), 'size': p.stat().st_size}
        records[rel] = actual
        if actual != expected:
            violations.append(rel+': hash/mode/size differs')
    return {'checked': len(records), 'violations': violations, 'files': records}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-dir', type=Path, required=True)
    ap.add_argument('--protocol-commit', required=True)
    args=ap.parse_args()
    base=args.run_dir.resolve()
    base.mkdir(parents=True, exist_ok=True)
    logs=base/'execution_logs'
    logs.mkdir(exist_ok=False)
    commands={}
    def run(name, argv, cwd=ROOT, expected=0):
        rec=record_command(logs,name,argv,cwd)
        if rec['exit_code'] != expected:
            raise RuntimeError(f'{name}: exit {rec["exit_code"]}, expected {expected}; see {logs}')
        return logs/f'{name}.command.json'
    identity=code_identity()
    if not identity['clean']:
        raise RuntimeError('Commit all source and tests before actual comparison')
    commit=identity['commit']
    for name in ('main','reproduction','tests_final','historical_replay'):
        if (base/name).exists():
            raise RuntimeError(f'refuse existing output {base/name}')
    before=protection(base/'protected_before.json')
    write_json(base/'protected_start.json',before)
    if before['violations']:
        raise RuntimeError('protected inputs changed before execution')
    data, issues=load_dataset(ROOT)
    if issues:
        raise RuntimeError(str(issues))
    write_json(base/'inputs_before.json',data['inputs'])
    commands['tests']=run('tests',[sys.executable,'-B','tests/test_taxonomy_comparison.py',
                                  '--log-dir',base/'tests_final'])
    run('historical_replay',[sys.executable,'-B','src/repro_defects.py',
                             '--log-dir',base/'historical_replay'])
    commands['main']=run('main',[sys.executable,'-B','src/taxonomy_comparison.py',
                                 '--input-dir',ROOT,'--out-dir',base/'main'])
    commands['tables_main']=run('tables_main',[sys.executable,'-B','src/check_taxonomy_comparison.py',
                                '--input-dir',ROOT,'--out-dir',base/'main','--tables-only'])
    with tempfile.TemporaryDirectory(prefix='taxonomy-clean-reproduction-') as temp:
        clone=Path(temp)/'source'
        run('clone',['git','clone','--local','--no-hardlinks','--no-checkout',ROOT,clone])
        run('checkout',['git','checkout','--detach',commit],clone)
        for rel in INPUT_FILES.values():
            dest=clone/rel
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,dest)
            dest.chmod(0o444)
        # Only private QZA subtrees are made non-writable; no source permission changes.
        (clone/'work/qza').chmod(0o555)
        copied, issues=load_dataset(clone)
        if issues or copied['inputs'] != data['inputs']:
            raise RuntimeError('copied input mismatch')
        def status():
            return subprocess.check_output(['git','status','--porcelain'],cwd=clone,text=True)
        clean_before=status()
        if clean_before:
            raise RuntimeError('clone is not clean')
        commands['reproduction']=run('reproduction',[sys.executable,'-B','src/taxonomy_comparison.py',
                                      '--input-dir',clone,'--out-dir',base/'reproduction'],clone)
        commands['tables_reproduction']=run('tables_reproduction',[sys.executable,'-B',
                        'src/check_taxonomy_comparison.py','--input-dir',clone,
                        '--out-dir',base/'reproduction','--tables-only'],clone)
        copied_after, issues=load_dataset(clone)
        if issues:
            raise RuntimeError(str(issues))
        checkout={'commit':commit,'status_before':clean_before,'status_after':status(),
                  'clone_command_record':str(logs/'clone.command.json'),
                  'checkout_command_record':str(logs/'checkout.command.json'),
                  'inputs_before':copied['inputs'],'inputs_after':copied_after['inputs'],
                  'readonly_inputs':all(not stat.S_IMODE((clone/p).stat().st_mode)&0o222 for p in INPUT_FILES.values())}
        write_json(base/'clean_checkout.json',checkout)
        # TemporaryDirectory must remove the private copy, including its read-only directory.
        (clone/'work/qza').chmod(0o755)
    diff,counts=compare_tables(base/'main',base/'reproduction')
    write_json(base/'reproduction_comparison.json',{'tolerance':1e-12,'violations':diff,'rows':counts})
    if diff:
        raise RuntimeError(str(diff))
    commands['project_check']=run('project_check',[sys.executable,'-B','checks/check.py'])
    after,issues=load_dataset(ROOT)
    if issues or after['inputs']!=data['inputs']:
        raise RuntimeError('project inputs changed')
    write_json(base/'inputs_after.json',after['inputs'])
    protected=protection(base/'protected_before.json')
    write_json(base/'protected_after.json',protected)
    if protected['violations']:
        raise RuntimeError('protected files changed')
    out=base/'main'
    def relative(p):return os.path.relpath(p,out)
    artifacts={}
    for dirname in ('execution_logs','tests_final','historical_replay','reproduction'):
        for p in (base/dirname).rglob('*'):
            if p.is_file():artifacts[relative(p)]=sha256(p)
    for name in ('inputs_before.json','inputs_after.json','clean_checkout.json',
                 'reproduction_comparison.json','protected_before.json','protected_start.json','protected_after.json'):
        p=base/name;artifacts[relative(p)]=sha256(p)
    evidence={'schema_version':1,'code_commit':commit,'protocol_commit':args.protocol_commit,
              'inputs_before':data['inputs'],'inputs_after':after['inputs'],
              'tests_summary':relative(base/'tests_final/summary.json'),
              'commands':{k:relative(v) for k,v in commands.items()},
              'reproduction_dir':'../reproduction','checkout_record':relative(base/'clean_checkout.json'),
              'protected_record':relative(base/'protected_after.json'),'artifacts':artifacts}
    write_json(out/'acceptance.json',evidence)
    run('full_acceptance',[sys.executable,'-B','src/check_taxonomy_comparison.py',
                              '--input-dir',ROOT,'--out-dir',out])
    print(f'PASS/0 — comparison, tests and local reproduction: {base}')

if __name__=='__main__':
    main()

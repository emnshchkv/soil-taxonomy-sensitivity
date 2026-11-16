"""Evidence required for the full comparison verdict (no writes or reruns)."""
import csv
from datetime import datetime
import os
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from taxonomy_compare_io import ROOT, SCHEMAS, CODE_FILES, sha256, issue, verify_code_identity


def compare_tables(first, second):
    """Join on exact keys; tolerate only floating metric columns, never IDs/counts."""
    from check_taxonomy_comparison import ROW_KEYS, read_rows
    errors, counts = [], {}
    for name, header in SCHEMAS.items():
        try:
            tables = []
            for folder in (first, second):
                local_errors=[]
                rows=read_rows(Path(folder)/name,header,local_errors)
                if local_errors:
                    errors.extend(local_errors)
                    raise ValueError('invalid reproduction table')
                mapped = {tuple(r[k] for k in ROW_KEYS[name]): r for r in rows}
                if len(mapped) != len(rows):
                    raise ValueError('duplicate key')
                tables.append(mapped)
            a,b=tables
            if set(a) != set(b):
                raise ValueError('row keys differ')
            for key in a:
                for col in header:
                    va,vb=a[key][col],b[key][col]
                    if va == vb:
                        continue
                    is_float = col.endswith('_frac') or col in {'frac','distance','min','median','max',
                                                               'assigned_reads_frac_a','assigned_reads_frac_b'}
                    if not is_float:
                        raise ValueError(f'{key}/{col}: exact value differs')
                    xa,xb=float(va),float(vb)
                    if not math.isfinite(xa) or not math.isfinite(xb) or abs(xa-xb)>1e-12:
                        raise ValueError(f'{key}/{col}: numeric value differs')
            counts[name]=len(a)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            errors.append(issue('REPRODUCTION', f'{name}: {exc}', isinstance(exc,FileNotFoundError)))
    return errors, counts


def check_evidence(out_dir, dataset, run):
    """Validate test records, committed source, frozen inputs and clean reproduction."""
    out_dir = Path(out_dir)
    evidence = out_dir / 'acceptance.json'
    if not evidence.is_file():
        return [issue('EVIDENCE_MISSING', evidence, True)]
    errors = []
    try:
        e = json.loads(evidence.read_text())
        def get(rel):
            p = (out_dir / rel).resolve()
            if not p.is_file():
                errors.append(issue('EVIDENCE_MISSING', rel, True))
                return None
            expected = e['artifacts'].get(rel)
            if expected != sha256(p):
                errors.append(issue('EVIDENCE_HASH', rel))
            return p
        for rel in e['artifacts']:
            get(rel)
        if not run['code']['clean'] or e['code_commit'] != run['code']['commit']:
            errors.append(issue('EVIDENCE_CODE', 'main calculation must use committed source'))
        code = e['code_commit']
        protocol = e['protocol_commit']
        if not re.fullmatch('[0-9a-f]{40}',code) or not re.fullmatch('[0-9a-f]{40}',protocol) or code==protocol:
            errors.append(issue('PROTOCOL_VERSION','distinct full source and prior protocol commits required'))
        ancestry = subprocess.run(['git', 'merge-base', '--is-ancestor', protocol, code], cwd=ROOT, capture_output=True)
        proto = subprocess.run(['git', 'show', f'{protocol}:METHOD_COMPARISON.md'], cwd=ROOT, capture_output=True)
        if ancestry.returncode or proto.returncode:
            errors.append(issue('PROTOCOL_COMMIT', 'protocol must precede source commit'))
        else:
            text = proto.stdout.decode()
            for item, rec in dataset['inputs'].items():
                if f'`{rec["path"]}` | `{rec["sha256"]}`' not in text:
                    errors.append(issue('PROTOCOL_INPUT', item))
        if e['inputs_before'] != dataset['inputs'] or e['inputs_after'] != dataset['inputs']:
            errors.append(issue('INPUT_PRESERVATION', 'inputs differ from recorded before/after state'))
        test_path = get(e['tests_summary'])
        if test_path:
            tests = json.loads(test_path.read_text())
            if not tests['success'] or tests['tests_run'] < 1 or tests['failures'] or tests['errors'] or tests.get('skipped',0):
                errors.append(issue('SELF_TEST', 'synthetic tests did not all pass'))
            errors.extend(verify_code_identity(tests['code']))
            if tests['code']['files'] != run['code']['files']:
                errors.append(issue('SELF_TEST_CODE', 'tests and calculation use different source'))
            if len(tests['commands']) < 40 or tests['tests_run'] < 70:
                errors.append(issue('SELF_TEST_LOG', 'missing CLI cases'))
            names=[c['name'] for c in tests['commands']]
            metas=[c['metadata_file'] for c in tests['commands']]
            if len(set(names))!=len(names) or len(set(metas))!=len(metas):
                errors.append(issue('SELF_TEST_DUPLICATE','CLI records must be unique'))
            required_cases={'baseline_calculation':0,'baseline_tables':0,
                'test_missing_input':2,'test_missing_summary_dependency':2,
                'test_bad_input_outranks_missing_output':1,'test_frequency_negative':1,
                'test_frequency_fraction':1,'test_frequency_nan':1,'test_frequency_infinity':1,
                'test_duplicate_metadata':1,'test_mutated_distance':1,'test_mutated_summary':1}
            for case,expected in required_cases.items():
                matches=[c for c in tests['commands'] if c['name'].split('_',1)[-1]==case]
                if len(matches)!=1 or matches[0]['expected_exit']!=expected:
                    errors.append(issue('SELF_TEST_COVERAGE',case))
            for command in tests['commands']:
                meta_rel=os.path.relpath(test_path.parent/command['metadata_file'],out_dir)
                meta_path=get(meta_rel)
                if meta_path and json.loads(meta_path.read_text()) != command:
                    errors.append(issue('SELF_TEST_METADATA', meta_rel))
                for stream in ('stdout','stderr'):
                    rel=os.path.relpath(test_path.parent/command[stream],out_dir)
                    stream_path=get(rel)
                    if stream_path and sha256(stream_path) != command[stream+'_sha256']:
                        errors.append(issue('SELF_TEST_STREAM',rel))
                if command['exit_code'] != command['expected_exit']:
                    errors.append(issue('SELF_TEST_CASE', command.get('command','unknown')))
        required = {'tests', 'main', 'tables_main', 'reproduction', 'tables_reproduction', 'project_check'}
        if set(e['commands']) != required:
            errors.append(issue('COMMAND_SET', 'required command records differ'))
        records={}
        scripts={'tests':'test_taxonomy_comparison.py','main':'taxonomy_comparison.py',
                 'tables_main':'check_taxonomy_comparison.py','reproduction':'taxonomy_comparison.py',
                 'tables_reproduction':'check_taxonomy_comparison.py','project_check':'check.py'}
        for label, rel in e['commands'].items():
            p = get(rel)
            if not p:
                continue
            c = json.loads(p.read_text())
            records[label]=c
            argv=c['argv']
            if len(argv)<3 or not Path(argv[0]).name.startswith('python') or argv[1]!='-B' or Path(argv[2]).name!=scripts[label]:
                errors.append(issue('COMMAND_SCRIPT',label))
            start,end=datetime.fromisoformat(c['start_utc']),datetime.fromisoformat(c['end_utc'])
            if start.tzinfo is None or end.tzinfo is None or start>end:
                errors.append(issue('COMMAND_TIME',label))
            if label=='tests':
                try:
                    if (Path(c['cwd'])/argv[argv.index('--log-dir')+1]).resolve() != (out_dir/e['tests_summary']).parent.resolve():
                        errors.append(issue('COMMAND_TESTS_PATH',label))
                except (ValueError,IndexError):
                    errors.append(issue('COMMAND_ARGUMENTS',label))
            if label in ('main','tables_main','reproduction','tables_reproduction'):
                expected_out=out_dir if label in ('main','tables_main') else out_dir/e['reproduction_dir']
                try:
                    if Path(argv[argv.index('--out-dir')+1]).resolve() != expected_out.resolve():
                        errors.append(issue('COMMAND_OUTPUT',label))
                    actual_input=Path(argv[argv.index('--input-dir')+1]).resolve()
                    if label in ('main','tables_main') and actual_input!=Path(run['input_dir']).resolve():
                        errors.append(issue('COMMAND_INPUT',label))
                    if Path(c['cwd']).resolve()!=actual_input:
                        errors.append(issue('COMMAND_CWD',label))
                    if label.startswith('tables_') and '--tables-only' not in argv:
                        errors.append(issue('COMMAND_SCOPE',label))
                except (ValueError,IndexError):
                    errors.append(issue('COMMAND_ARGUMENTS',label))
            if c['exit_code'] != 0 or not c['argv'] or not c['start_utc'] or not c['end_utc']:
                errors.append(issue('COMMAND_RESULT', label))
            for stream in ('stdout', 'stderr'):
                s = (p.parent / c[stream]).resolve()
                relative = str(Path(os.path.relpath(s, out_dir)))
                get(relative)
        if set(records)==required:
            sequence=('tests','main','tables_main','reproduction','tables_reproduction','project_check')
            for a,b in zip(sequence,sequence[1:]):
                if datetime.fromisoformat(records[a]['end_utc']) > datetime.fromisoformat(records[b]['start_utc']):
                    errors.append(issue('COMMAND_ORDER',a+' -> '+b))
        reproduction = (out_dir / e['reproduction_dir']).resolve()
        rpath = get(str(Path(os.path.relpath(reproduction / 'run.json', out_dir))))
        if rpath:
            rr = json.loads(rpath.read_text())
            errors.extend(verify_code_identity(rr['code']))
            for label in ('reproduction','tables_reproduction'):
                if label in records:
                    argv=records[label]['argv']
                    if Path(argv[argv.index('--input-dir')+1]).resolve()!=Path(rr['input_dir']).resolve():
                        errors.append(issue('REPRODUCTION_COMMAND_INPUT',label))
            if set(rr['outputs']) != set(run['outputs']):
                errors.append(issue('REPRODUCTION_OUTPUT_SET','output set differs'))
            if rr['code'] != run['code'] or rr['inputs'] != dataset['inputs']:
                errors.append(issue('REPRODUCTION_IDENTITY', 'reproduction source or inputs differ'))
            for name, expected in rr['outputs'].items():
                p = get(str(Path(os.path.relpath(reproduction / name, out_dir))))
                if p and sha256(p) != expected:
                    errors.append(issue('REPRODUCTION_HASH', name))
            diff, _ = compare_tables(out_dir, reproduction)
            errors.extend(diff)
        p = get(e['checkout_record'])
        if p:
            checkout = json.loads(p.read_text())
            if checkout['commit'] != code or checkout['status_before'] or checkout['status_after']:
                errors.append(issue('CLEAN_CHECKOUT', 'source checkout was not clean'))
            if checkout['inputs_before'] != dataset['inputs'] or checkout['inputs_after'] != dataset['inputs']:
                errors.append(issue('REPRODUCTION_INPUTS', 'copy inputs changed'))
            if not checkout['readonly_inputs']:
                errors.append(issue('REPRODUCTION_PERMISSIONS', 'private copies not read-only'))
        p = get(e['protected_record'])
        if p:
            protection = json.loads(p.read_text())
            if protection['violations'] or protection['checked'] < 83:
                errors.append(issue('PROTECTED_FILES', 'protected input/results changed'))
    except (KeyError, TypeError, ValueError, OSError) as exc:
        errors.append(issue('EVIDENCE_FORMAT', f'{type(exc).__name__}: {exc}', isinstance(exc,FileNotFoundError)))
    return errors

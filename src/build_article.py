"""Build the existing article with pdfLaTeX/BibTeX and save actual commands.

Does not install TeX or modify the system environment. Pass --tex-bin only
when a working TeX distribution is not already on PATH.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tex-bin', type=Path)
    parser.add_argument('--log-dir', type=Path, required=True)
    args = parser.parse_args()
    args.log_dir.mkdir(parents=True, exist_ok=False)
    directory = ROOT/'paper/ai4science_ysda2026'
    env = os.environ.copy()
    if args.tex_bin:
        env['PATH'] = str(args.tex_bin) + os.pathsep + env.get('PATH','')
    tools = {name:shutil.which(name,path=env.get('PATH')) for name in ('pdflatex','bibtex')}
    if not all(tools.values()):
        raise SystemExit('TeX missing; no packages installed: '+repr(tools))
    report = {'started_utc':datetime.now(timezone.utc).isoformat(),'status':'INCOMPLETE','commands':[]}
    commands = [[tools['pdflatex'],'-interaction=nonstopmode','-halt-on-error','paper.tex'],
                [tools['bibtex'],'paper'],
                [tools['pdflatex'],'-interaction=nonstopmode','-halt-on-error','paper.tex'],
                [tools['pdflatex'],'-interaction=nonstopmode','-halt-on-error','paper.tex']]
    for i,command in enumerate(commands,1):
        # TeX may wrap a log line in the middle of a UTF-8 byte sequence.
        # Preserve original bytes instead of making decoding part of the build.
        result = subprocess.run(command,cwd=directory,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        name=f'{i:02d}_{Path(command[0]).name}.log'
        (args.log_dir/name).write_bytes(result.stdout)
        report['commands'].append({'command':command,'cwd':str(directory),'log':name,'exit_code':result.returncode})
        (args.log_dir/'build.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
        if result.returncode:
            report['status']='FAIL'
            (args.log_dir/'build.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
            print('Build failed: '+str(args.log_dir/name));return result.returncode
    report['status']='PASS'
    report['outputs']={f'paper.{ext}':hashlib.sha256((directory/f'paper.{ext}').read_bytes()).hexdigest() for ext in ('pdf','aux','log','bbl')}
    (args.log_dir/'build.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('PASS/0 — pdfLaTeX → BibTeX → pdfLaTeX → pdfLaTeX; '+str(args.log_dir/'build.json'))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

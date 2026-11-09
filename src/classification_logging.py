"""Complete command evidence: unabridged streams, real exit codes and UTC times."""
import hashlib
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def write_json(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, default=str)
        f.write('\n')

def record_command(directory, name, command, cwd):
    directory = Path(directory)
    command = list(map(str, command))
    start = datetime.now(timezone.utc).isoformat()
    stdout = directory / f'{name}.stdout.txt'
    stderr = directory / f'{name}.stderr.txt'
    with stdout.open('xb') as out, stderr.open('xb') as err:
        try:
            run = subprocess.run(command, cwd=cwd, stdout=out, stderr=err, check=False)
            code = run.returncode
        except OSError as error:
            err.write(f'{type(error).__name__}: {error}\n'.encode())
            code = 127
    result = {'argv': command, 'command': shlex.join(command), 'cwd': str(cwd),
              'start_utc': start, 'end_utc': datetime.now(timezone.utc).isoformat(),
              'exit_code': code, 'stdout': stdout.name, 'stderr': stderr.name}
    write_json(directory / f'{name}.command.json', result)
    print(f'{name}: exit={code}; logs={directory}', flush=True)
    return result

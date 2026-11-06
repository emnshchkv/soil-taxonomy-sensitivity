#!/usr/bin/env python3
"""Exercise the real CLI on isolated copies of existing outputs; never run DADA2."""
import argparse
import contextlib
import csv
import io
import json
import os
from pathlib import Path
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from qiime_runtime import selected_python
selected_python()
import check_exact_match as exact
import verify_denoising as api
import biom
import h5py

def snapshot(root):
    paths = set()
    for name in ("data", "work", "memory", "results", "checks", "env"):
        paths.update(p for p in (root / name).rglob("*") if p.is_file())
    return {str(p.relative_to(root)): (api.sha256(p), stat.S_IMODE(p.stat().st_mode))
            for p in sorted(paths)}

def prepare(source, target):
    # Only raw FASTQ and immutable imported archives are shared, and are never mutated.
    (target / "data").mkdir(parents=True)
    (target / "data/raw").symlink_to(source / "data/raw", target_is_directory=True)
    shutil.copytree(source / "data/meta", target / "data/meta")
    shutil.copytree(source / "env", target / "env")
    shutil.copy2(source / "data/raw.sha256", target / "data/raw.sha256")
    for name in ("read_counts_R1.tsv", "quality_by_position.tsv", "quality_by_sample.tsv"):
        (target / "results").mkdir(exist_ok=True)
        shutil.copy2(source / "results" / name, target / "results" / name)
    (target / "memory").mkdir()
    shutil.copy2(source / "memory/denoise_choice.md", target / "memory/denoise_choice.md")
    for directory in api.RUNS:
        (target / directory).mkdir(parents=True, exist_ok=True)
        for product in api.PRODUCTS:
            shutil.copy2(source / directory / product[0], target / directory / product[0])
    for name in ("work/qza/raw_reads_single.qza", "work/qzv/reads_summary.qzv"):
        p = target / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.symlink_to(source / name)
    (target / ".exact-regression-copy").write_text("isolated mutable outputs\n")

@contextlib.contextmanager
def absent(path):
    assert not path.is_symlink()
    backup = path.with_name(path.name + ".test-hidden")
    path.rename(backup)
    try:
        yield
    finally:
        backup.rename(path)

@contextlib.contextmanager
def transfer(root):
    assert (root / ".exact-regression-copy").is_file()
    path = root / api.RUNS[2] / "ASV_table.qza"
    assert not path.is_symlink() and path.stat().st_nlink == 1
    original = path.read_bytes()
    with zipfile.ZipFile(io.BytesIO(original)) as z:
        members = [(info, z.read(info.filename)) for info in z.infolist()]
    entry = next(i for i, (info, _) in enumerate(members)
                 if info.filename.endswith("/data/feature-table.biom"))
    with h5py.File(io.BytesIO(members[entry][1]), "r") as h:
        table = biom.Table.from_hdf5(h)
    values = table.matrix_data.toarray()
    j = list(table.ids("sample")).index("Kr1_S80")
    donor = int(values[:, j].argmax())
    receiver = (donor + 1) % len(values)
    assert values[donor, j] >= 1
    values[donor, j] -= 1
    values[receiver, j] += 1
    changed = biom.Table(values, table.ids("observation"), table.ids("sample"))
    buf = io.BytesIO()
    with h5py.File(buf, "w") as h:
        changed.to_hdf5(h, "isolated intentional read transfer")
    members[entry] = (members[entry][0], buf.getvalue())
    try:
        with zipfile.ZipFile(path, "w") as z:
            for info, payload in members:
                z.writestr(info, payload)
        yield
    finally:
        path.write_bytes(original)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    args = parser.parse_args()
    source, logs = args.project.resolve(), args.log_dir.resolve()
    if source == logs or source in logs.parents:
        parser.error("--log-dir must be outside the protected project during the test")
    logs.mkdir(parents=True, exist_ok=False)
    before = snapshot(source)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "diff", "--", "src/check_exact_match.py",
                                         "tests/test_exact_match_cli.py"], cwd=ROOT))
    rows = []
    with tempfile.TemporaryDirectory(prefix="exact-match-cases-") as tmp:
        root = Path(tmp)
        prepare(source, root)

        def case(label, expected):
            # This test already runs in the selected conda interpreter. Nesting
            # conda activation can fail in compiler activation hooks before Python.
            command = [sys.executable, "-B", str(ROOT / "src/check_exact_match.py"), "--root", str(root)]
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"), timeout=180)
            output = result.stdout + result.stderr
            (logs / (label + ".txt")).write_text(
                f"Commit: {revision}\nTracked code differs from commit: {dirty}\n"
                f"Wrapper SHA256: {api.sha256(ROOT / 'src/check_exact_match.py')}\n"
                f"Command: {shlex.join(command)}\nExpected exit: {expected}\n"
                + output + f"\nexit_code={result.returncode}\n")
            ok = result.returncode == expected
            if label in ("complete", "restored"):
                ok = ok and sum(line.startswith("runs=") for line in output.splitlines()) == 3
                ok = ok and "exact-match: PASS/0" in output
            if label.startswith("transfer"):
                ok = ok and "REPEAT_DIFFERENCE" in output and "(0, 2, 1, 0)" in output
            rows.append((label, expected, result.returncode, ok, label + ".txt"))
            print(f"{label}: expected={expected} actual={result.returncode} matched={ok}", flush=True)

        case("complete", 0)
        with absent(root / "results/read_counts_R1.tsv"):
            case("missing_counts", 2)
        with absent(root / "results/quality_by_sample.tsv"):
            case("missing_quality", 2)
        with absent(root / api.RUNS[2]):
            case("missing_run", 2)
        with transfer(root):
            case("transfer", 1)
            with absent(root / "results/read_counts_R1.tsv"):
                case("transfer_and_missing", 1)
        case("restored", 0)

    after = snapshot(source)
    (logs / "protected_before.json").write_text(json.dumps(before, indent=2) + "\n")
    (logs / "protected_after.json").write_text(json.dumps(after, indent=2) + "\n")
    with (logs / "summary.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(("case", "expected_exit", "actual_exit", "matched", "log"))
        writer.writerows(rows)
    print(f"Protected source files unchanged: {before == after}; files={len(before)}")
    return 0 if before == after and all(row[3] for row in rows) else 1

if __name__ == "__main__":
    sys.exit(main())

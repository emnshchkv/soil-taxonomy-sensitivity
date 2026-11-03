#!/usr/bin/env python3
"""Read-only integration checks; private artifacts are optional explicit inputs."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import stat
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from qiime_runtime import selected_python
selected_python()

spec = importlib.util.spec_from_file_location("verifier", ROOT/"checks/verify_denoising.py")
api = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = api
spec.loader.exec_module(api)


def snapshot(root):
    files = set()
    for directory in ("data", "work", "memory", "results", "tests/fixtures"):
        files.update(p for p in (root/directory).rglob("*") if p.is_file())
    files.update(root/p for p in ("METHOD.md", "checks/check.py", "ontology.yaml"))
    return {str(p.relative_to(root)): (api.sha256(p),stat.S_IMODE(p.stat().st_mode))
            for p in sorted(files)}


def safety():
    before = snapshot(ROOT)
    run = subprocess.run([sys.executable,"-B",str(ROOT/"checks/verify_denoising.py"),"--self-test"],
                         cwd=ROOT, text=True, capture_output=True)
    after = snapshot(ROOT)
    if run.returncode:
        print(run.stdout, run.stderr)
    assert run.returncode == 0, run.returncode
    changed = sorted(name for name in before.keys() | after.keys() if before.get(name) != after.get(name))
    assert not changed, f"project files changed during self-test: {changed}"
    print(f"self-test isolation: PASS; {len(before)} files, SHA256 and modes unchanged")
    # Regression: the obsolete self-test overwrote these paths even when populated.
    with tempfile.TemporaryDirectory(prefix="verifier-sentinels-") as tmp:
        root = Path(tmp)
        sentinels = ["memory/denoise_choice.md"] + [
            f"{d}/{p[0]}" for d in api.RUNS for p in api.PRODUCTS]
        for name in sentinels:
            path = root/name
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(b"SENTINEL: must never be overwritten\n")
        before = {p:api.sha256(root/p) for p in sentinels}
        run = subprocess.run([sys.executable,"-B",str(ROOT/"checks/verify_denoising.py"),
                              "--self-test","--root",str(root)], cwd=root, capture_output=True)
        assert run.returncode == 0, run.stderr.decode()
        assert before == {p:api.sha256(root/p) for p in sentinels}
        print("populated foreign work/memory sentinels: PASS; 10 files unchanged")


def real_archives(directory):
    import biom
    import pandas as pd
    from Bio import SeqIO
    from qiime2 import Artifact
    before = {p[0]:api.sha256(directory/p[0]) for p in api.PRODUCTS}
    observed = api.load_run(directory)
    artifacts = [Artifact.load(str(directory/p[0])) for p in api.PRODUCTS]
    for artifact,product in zip(artifacts,api.PRODUCTS):
        artifact.validate(level="max")
        assert str(artifact.type) == product[1]
    canonical = artifacts[0].view(biom.Table)
    assert set(canonical.ids("sample")) == observed.samples
    assert set(canonical.ids("observation")) == observed.features
    expected = {(f,s):int(canonical.get_value_by_ids(f,s))
                for f in canonical.ids("observation") for s in canonical.ids("sample")
                if canonical.get_value_by_ids(f,s)}
    assert expected == observed.matrix
    with tempfile.TemporaryDirectory(prefix="verifier-official-export-") as tmp:
        for artifact,product in zip(artifacts[1:],api.PRODUCTS[1:]):
            destination = Path(tmp)/product[2]
            artifact.export_data(str(destination))
        dna = {r.id:str(r.seq) for r in SeqIO.parse(Path(tmp)/"representative_sequences/dna-sequences.fasta","fasta")}
        assert dna == observed.dna
        frame = pd.read_csv(Path(tmp)/"denoising_stats/stats.tsv",sep="\t",index_col=0)
        frame = frame.loc[~frame.index.str.startswith("#q2:")]
        expected_stats = {str(s):tuple(int(row[c]) for c in api.STATS) for s,row in frame.iterrows()}
        assert expected_stats == observed.stats
    counts = {s:row[0] for s,row in expected_stats.items()}
    errors = api.validate(api.Dataset({1:observed},counts,"K_S63",(20,200),api.INPUT_ID))
    assert not errors, errors
    assert before == {p[0]:api.sha256(directory/p[0]) for p in api.PRODUCTS}
    print(f"real QZA vs QIIME Artifact.validate/view/export: PASS; {len(dna)} ASV, {len(counts)} samples")
    print("reader check only; archived outputs, no new DADA2 execution")
    print("real archive SHA256:",json.dumps(before,sort_keys=True))


def real_quality(archive):
    with tempfile.TemporaryDirectory(prefix="verifier-real-quality-") as tmp:
        root = Path(tmp)
        for name in ("data","env"):
            (root/name).symlink_to(ROOT/name, target_is_directory=True)
        (root/"work").mkdir()
        for name in ("qza","qzv"):
            (root/"work"/name).symlink_to(ROOT/"work"/name,target_is_directory=True)
        (root/"results").mkdir()
        (root/"results/read_counts_R1.tsv").symlink_to(ROOT/"results/read_counts_R1.tsv")
        with tarfile.open(archive,"r") as tar:
            for name in ("quality_by_position.tsv","quality_by_sample.tsv"):
                payload = tar.extractfile("results/"+name).read()
                (root/"results"/name).write_bytes(payload)
                print(f"reference {name} SHA256={hashlib.sha256(payload).hexdigest()}")
        out = api.Audit()
        api.check_quality(root,out)
        assert out.exit_code == 0, (out.errors,out.incomplete)
        path = root/"results/quality_by_position.tsv"
        text = path.read_text()
        rows = list(csv.DictReader(io.StringIO(text),delimiter="\t"))
        print(f"real saved profiles vs independent FASTQ histogram: PASS; {len(rows)} positions, 11 samples")
        rows[0]["q50"] = str(float(rows[0]["q50"])+1)
        with path.open("w") as f:
            writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
        out = api.Audit()
        api.check_quality(root,out)
        assert [v.code for v in out.errors] == ["QUALITY_VALUE"], out.errors
        assert out.exit_code == 1
        print("real saved median +1 in temporary copy: FAIL/1 QUALITY_VALUE (expected)")
        path.write_text(text)
        out = api.Audit(); api.check_quality(root,out)
        assert out.exit_code == 0, (out.errors,out.incomplete)
        print("restored real profile in temporary copy: PASS/0")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-run",type=Path)
    parser.add_argument("--quality-archive",type=Path)
    args = parser.parse_args()
    safety()
    if args.real_run:
        real_archives(args.real_run)
    if args.quality_archive:
        real_quality(args.quality_archive)
    print("integration checks: PASS/0")

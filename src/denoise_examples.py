"""Versioned, invented QIIME-layout examples for METHOD.md. Mutations are confined to temporary directories."""
import argparse
import copy
import sys
import gzip
import io
import tempfile
import uuid
import zipfile
from pathlib import Path

if __name__ == "__main__":
    sys.dont_write_bytecode = True
    from qiime_runtime import selected_python
    selected_python()

import biom
import h5py
import numpy as np
import yaml


CONTROL, SOIL = "SYN_water", "SYN_soil"
COUNTS = {CONTROL: 10, SOIL: 12}
DNA = {"a": "ACGT", "b": "TGCA"}
INPUT = "ff500000-0000-4000-8000-000000000001"
ORIGINAL = {CONTROL: {"a": 2, "b": 3}, SOIL: {"a": 5, "b": 4}}
STATS = {CONTROL: (10, 8, 7, 5), SOIL: (12, 11, 10, 9)}
SETTINGS = dict(trim_left=1, trunc_len=5, max_ee=3.0, min_fold_parent_over_abundance=2.0,
                n_threads=6, trunc_q=2, pooling_method="independent", chimera_method="consensus",
                allow_one_off=False, n_reads_learn=1000000, hashed_feature_ids=True, retain_all_samples=True)


def stable_id(key):
    return str(uuid.uuid5(uuid.NAMESPACE_URL, "synthetic.example/denoise/" + key))


def hdf5_table(matrix, features, samples):
    values = np.array([[matrix[s][f] for s in samples] for f in features])
    t = biom.Table(values, features, samples)
    buffer = io.BytesIO()
    with h5py.File(buffer, "w") as h:
        t.to_hdf5(h, "synthetic verifier fixture")
    return buffer.getvalue()


def archive(path, semantic, output, member, payload, run_number, settings, mutation="clean"):
    """QIIME 2 archive member names, BIOM 2.1 and marked synthetic execution metadata."""
    aid = stable_id(f"{run_number}/{output}")
    prov = dict(execution=dict(uuid=stable_id(str(run_number))),
                environment=dict(framework=dict(version="2025.7.0"), plugins=dict(dada2=dict(version="2025.7.0"))),
                action=dict(type="method", plugin="environment:plugins:dada2", action="denoise_single",
                            inputs=[{"demultiplexed_seqs": INPUT}], parameters=[{k:v} for k,v in settings.items()],
                            **{"output-name":output}))
    if run_number == 3:
        if mutation == "version":
            prov["environment"]["plugins"]["dada2"]["version"] = "2024.10.0"
        elif mutation == "input":
            prov["action"]["inputs"] = [{"demultiplexed_seqs": stable_id("wrong-input")}]
        elif mutation == "mixed_execution" and output == "representative_sequences":
            prov["execution"]["uuid"] = stable_id("wrong-execution")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(aid+"/metadata.yaml", yaml.safe_dump(dict(uuid=aid, type=semantic)))
        z.writestr(aid+"/provenance/action/action.yaml", yaml.safe_dump(prov))
        z.writestr(aid+"/data/"+member, payload)


def fixture(root, mutation="clean"):
    for number in (1,2,3):
        directory = root / str(number)
        directory.mkdir()
        features = list(DNA)
        samples = list(COUNTS)
        dna = dict(DNA)
        matrix = copy.deepcopy(ORIGINAL)
        stats = dict(STATS)
        settings = dict(SETTINGS)
        execution = number
        if number == 3:
            if mutation == "transfer":
                matrix[SOIL]["a"] += 1
                matrix[SOIL]["b"] -= 1
            elif mutation == "order":
                features.reverse()
                samples.reverse()
            elif mutation == "reuse":
                execution = 1
            elif mutation == "length":
                dna["a"] += "A"
            elif mutation == "balance":
                matrix[SOIL]["a"] += 1
            elif mutation == "param":
                settings["trim_left"] = 0
            elif mutation == "increasing":
                stats[SOIL] = (12, 13, 10, 9)
            elif mutation == "low_retention":
                stats[SOIL] = (12, 11, 10, 3)
                matrix[SOIL] = {"a":2, "b":1}
        raw = hdf5_table(matrix, features, samples)
        fasta = "".join(f">{f}\n{dna[f]}\n" for f in features)
        header = "sample-id\tinput\tfiltered\tdenoised\tnon-chimeric\n#q2:types\tnumeric\tnumeric\tnumeric\tnumeric\n"
        lines = ["\t".join([s,*map(str,stats[s])]) for s in samples]
        if number == 3 and mutation == "duplicate":
            lines.append(lines[0])
        rows = header + "\n".join(lines) + "\n"
        archive(directory/"ASV_table.qza","FeatureTable[Frequency]","table","feature-table.biom",
                raw, execution, settings, mutation)
        archive(directory/"rep_seq.qza","FeatureData[Sequence]","representative_sequences","dna-sequences.fasta",
                fasta, execution, settings, mutation)
        archive(directory/"dada2_stats.qza","SampleData[DADA2Stats]","denoising_stats","stats.tsv",
                rows, execution, settings, mutation)


def example_audit(api, mutation, verbose=False):
    with tempfile.TemporaryDirectory(prefix="verifier-synthetic-") as tmp:
        root = Path(tmp)
        fixture(root, "balance" if mutation == "fail_and_missing" else mutation)
        out = api.Audit()
        paths = ("1","2","3")
        if mutation == "missing":
            paths = ("1","2","absent")
        elif mutation == "fail_and_missing":
            paths = ("absent","2","3")
        api.read_repeats(root, out, COUNTS, (1,5), INPUT, CONTROL, paths, verbose)
        return out


def exercise_properties(api, check):
    """Six contract families on genuinely different matrix dimensions."""
    with tempfile.TemporaryDirectory(prefix="verifier-properties-") as tmp:
        root = Path(tmp)
        fixture(root)
        initial = {i: api.load_run(root/str(i)) for i in (1, 2, 3)}
    families = ("axis order", "feature renaming", "balanced transfer", "per-sample balance",
                "stage monotonicity", "count domain", "purity")
    results = {name: [] for name in families}
    for size in (1, 2, 4, 7):
        samples = [CONTROL] + [f"SYN_soil_{i}" for i in range(size)]
        features = [f"f{i}" for i in range(size+1)]
        dna = {f: "".join("ACGT"[(i >> shift) & 3] for shift in (6,4,2,0))
               for i,f in enumerate(features)}
        matrix = {(f,s): i+j+1 for i,f in enumerate(features) for j,s in enumerate(samples)}
        sums = {s:sum(n for (f,sid),n in matrix.items() if sid == s) for s in samples}
        counts = {s:n+3 for s,n in sums.items()}
        runs = copy.deepcopy(initial)
        for run in runs.values():
            run.samples, run.features = set(samples), set(features)
            run.dna, run.matrix = dict(dna), dict(matrix)
            run.stats = {s:(n+3,n+2,n+1,n) for s,n in sums.items()}
        def codes():
            return {v.code for v in api.validate(api.Dataset(runs, counts, CONTROL, (1,5), INPUT))}
        check(f"baseline {len(features)}x{len(samples)}", set(), codes())
        before = copy.deepcopy(runs)
        codes()
        results["purity"].append(runs == before and codes() == codes())
        original = copy.deepcopy(runs[3])
        runs[3].matrix = dict(reversed(list(original.matrix.items())))
        runs[3].dna = dict(reversed(list(original.dna.items())))
        runs[3].stats = dict(reversed(list(original.stats.items())))
        results["axis order"].append(codes() == set())
        runs[3] = copy.deepcopy(original)
        rename = {f:"renamed_"+f for f in features}
        runs[3].features = set(rename.values())
        runs[3].dna = {rename[f]:v for f,v in dna.items()}
        runs[3].matrix = {(rename[f],s):n for (f,s),n in matrix.items()}
        results["feature renaming"].append(codes() == set())
        runs[3] = copy.deepcopy(original)
        sid = samples[-1]
        runs[3].matrix[(features[0],sid)] += 1
        runs[3].matrix[(features[1],sid)] -= 1
        results["balanced transfer"].append(codes() == {"REPEAT_DIFFERENCE"})
        runs[3] = copy.deepcopy(original)
        runs[3].matrix[(features[0],sid)] += 1
        results["per-sample balance"].append("READ_BALANCE" in codes())
        runs[3] = copy.deepcopy(original)
        n, _, den, nc = runs[3].stats[sid]
        runs[3].stats[sid] = (n,n+1,den,nc)
        results["stage monotonicity"].append("STATS_COUNT" in codes())
        domain = []
        for bad in (-1, 0.5, float("nan"), float("inf"), True):
            runs[3] = copy.deepcopy(original)
            runs[3].matrix[(features[0],sid)] = bad
            domain.append("COUNT_VALUE" in codes())
        results["count domain"].append(all(domain))
    for name,outcomes in results.items():
        check("property " + name + " sizes=2x2,3x3,5x5,8x8", 4, sum(outcomes))


def self_test(api):
    failures = []
    def check(name, expected, seen):
        print(f"{'ok' if seen == expected else 'FAIL'} {name}: {seen}")
        if seen != expected:
            failures.append((name, seen))
    snapshot = Path(__file__).resolve().parents[1]/"tests/fixtures/denoise-clean"
    with (snapshot/"SHA256SUMS").open() as f:
        recorded = {name: digest for digest, name in (line.split() for line in f)}
    expected_files = {f"{i}/{p[0]}" for i in (1,2,3) for p in api.PRODUCTS}
    check("versioned fixture files and SHA256", True, set(recorded) == expected_files and
          all(api.sha256(snapshot/name) == recorded[name] for name in expected_files))
    snapshot_check = api.Audit()
    api.read_repeats(snapshot, snapshot_check, COUNTS, (1,5), INPUT, CONTROL, ("1","2","3"), False)
    check("versioned fixture verdict", True, not snapshot_check.errors and not snapshot_check.incomplete)
    cases = {
        "clean": set(), "order": set(), "transfer": {"REPEAT_DIFFERENCE"},
        "duplicate": {"SAMPLE_KEY"}, "reuse": {"REUSED_RUN"},
        "length": {"ASV_LENGTH","REPEAT_DIFFERENCE"}, "balance": {"READ_BALANCE","REPEAT_DIFFERENCE"},
        "version": {"PROVENANCE"}, "input": {"PROVENANCE"}, "mixed_execution": {"PROVENANCE"},
        "fail_and_missing": {"READ_BALANCE","REPEAT_DIFFERENCE","INCOMPLETE"},
        "param": {"PROVENANCE"}, "increasing":{"STATS_COUNT","REPEAT_DIFFERENCE"},
        "missing":{"INCOMPLETE"}, "low_retention":{"INCOMPLETE","REPEAT_DIFFERENCE"},
    }
    for case, expected in cases.items():
        out = example_audit(api, case)
        seen = {v.code for v in out.errors}
        if out.incomplete:
            seen.add("INCOMPLETE")
        check(case, sorted(expected), sorted(seen))
        expected_exit = 1 if expected - {"INCOMPLETE"} else 2 if expected else 0
        check(case + " exit", expected_exit, out.exit_code)
    with tempfile.TemporaryDirectory(prefix="verifier-quality-") as tmp:
        path = Path(tmp)/"two.fastq.gz"
        with gzip.open(path, "wt", encoding="ascii") as f:
            f.write("@x\nA\n+\n+\n@y\nT\n+\n?\n")
        positions, summary = api.statistics(api.quality_histogram(path))
        check("two reads Q10,Q30", (2,15,20,25), tuple(positions["1"].values()))
        check("sample summary", (2,2,20,50), tuple(summary.values()))
        wrong = {"n_bases":"2","q25":"15","q50":"21","q75":"25"}
        out = api.Audit()
        api.compare_values(out, (SOIL,"1"), positions["1"], wrong)
        check("median +1", ["QUALITY_VALUE"], [v.code for v in out.errors])
        with gzip.open(path, "wt", encoding="ascii") as f:
            f.write("@x\nACGT\n+\n???\n")
        out = api.Audit()
        out.read("short quality", api.quality_histogram, path)
        check("malformed FASTQ", ["FASTQ_FORMAT"], [v.code for v in out.errors])
    from verifier_quality_cases import exercise_quality
    exercise_quality(api, check)
    exercise_properties(api, check)
    return failures


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create or repair ONLY a synthetic red/green case.")
    parser.add_argument("action", choices=["create", "repair"])
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    case = args.directory
    marker = case / ".verifier-synthetic"
    if args.action == "create":
        case.mkdir()  # Refuse any existing directory.
        fixture(case, "transfer")
        marker.write_text("synthetic verifier transfer v1\\n")
        print("created: third run soil [6,3], expected [5,4]")
    else:
        if (case.is_symlink() or (case/"3").is_symlink() or marker.is_symlink()
                or marker.read_text() != "synthetic verifier transfer v1\\n"):
            parser.error("expected a case created by this command")
        # Replace just the affected table atomically. Keep the other eight archives intact.
        target = case/"3/ASV_table.qza"
        with tempfile.TemporaryDirectory(dir=case/"3") as tmp:
            repaired = Path(tmp)/target.name
            archive(repaired, "FeatureTable[Frequency]", "table", "feature-table.biom",
                    hdf5_table(ORIGINAL, list(DNA), list(COUNTS)), 3, SETTINGS)
            repaired.replace(target)
        print("repaired: third run soil [5,4]")

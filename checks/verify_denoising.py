#!/usr/bin/env python3
"""Read-only verification of R1 quality and three DADA2 runs (METHOD.md §§3–7)."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import math
import re
import sys
import time
import uuid
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from itertools import combinations
from pathlib import Path

sys.dont_write_bytecode = True
if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from qiime_runtime import selected_python
    selected_python()
import biom
import h5py
import numpy as np
import yaml
from Bio.SeqIO.QualityIO import FastqGeneralIterator

PROJECT = Path(__file__).resolve().parents[1]
INPUT_ID = "88f0322c-49af-48a6-b7d7-03db0b914866"
SAMPLES = {"K_S63", *(f"Kr{i}_S{79+i}" for i in range(1, 6)),
           *(f"LUB{i}_S{90+i}" for i in range(2, 7))}
HASHES = {
    "data/raw.sha256": "aeaaba15a9b6479938af3ddeb86bf86c23ff1ca9270c5ac3c685f1ac5f0300ab",
    "data/meta/metadata.tsv": "41ceed920479044d4fd8aaa75e1dc099593e3545fb649260925c6c7645e1b123",
    "data/meta/fq-manifest-single.tsv": "f5c8a4d06f352e4b8ebf519f9d362043de3a5779ba270d71adc9b84b6f636f17",
    "results/read_counts_R1.tsv": "242ae373be49497594a2a6fb799bdae14425b8c0a09fa271cbe58252f4fa9075",
    "work/qza/raw_reads_single.qza": "994baa25162444d445610afce07217312ed0a8c05b020e4b356d8377bf53e5f0",
    "work/qzv/reads_summary.qzv": "727859cc3cf164b239ac83defad338b30c4624ce18ab84b41330d9ad8d5bb1fb",
    "env/explicit.txt": "a063218673eca51c9266813365781d372a421396bc64238c34d025b43f30ed6b",
}
PARAMS = dict(max_ee=3, min_fold_parent_over_abundance=2, n_threads=6, trunc_q=2,
              pooling_method="independent", chimera_method="consensus", allow_one_off=False,
              n_reads_learn=1000000, hashed_feature_ids=True, retain_all_samples=True)
STATS = ("input", "filtered", "denoised", "non-chimeric")
RUNS = ("work/qza", "work/denoise_repeats/02", "work/denoise_repeats/03")
PRODUCTS = (
    ("ASV_table.qza", "FeatureTable[Frequency]", "table", "feature-table.biom"),
    ("rep_seq.qza", "FeatureData[Sequence]", "representative_sequences", "dna-sequences.fasta"),
    ("dada2_stats.qza", "SampleData[DADA2Stats]", "denoising_stats", "stats.tsv"),
)
TOLERANCE = Decimal("0.0000005")


@dataclass(frozen=True)
class Violation:
    code: str
    run: str
    key: str
    expected: object
    observed: object


class Invalid(ValueError):
    def __init__(self, code, detail):
        super().__init__(detail)
        self.code = code


@dataclass
class Audit:
    errors: list[Violation] = field(default_factory=list)
    incomplete: list[str] = field(default_factory=list)

    @property
    def exit_code(self):
        return 1 if self.errors else 2 if self.incomplete else 0

    def same(self, code, run, key, want, got):
        if want != got:
            self.errors.append(Violation(code, run, str(key), want, got))

    def read(self, target, function, *args):
        try:
            return function(*args)
        except FileNotFoundError as e:
            self.incomplete.append(f"{target}: missing {e.filename or e}")
        except Invalid as e:
            self.errors.append(Violation(e.code, "-", str(target), "valid", str(e)))
        except (zipfile.BadZipFile, UnicodeError, yaml.YAMLError, csv.Error,
                KeyError, TypeError, ValueError, AttributeError, IndexError, EOFError) as e:
            self.errors.append(Violation("SCHEMA", "-", str(target), "readable data", str(e)))
        except OSError as e:
            self.incomplete.append(f"{target}: {e}")
        return None

    def finish(self, stage):
        for v in self.errors:
            print(f"{v.code} run={v.run} key={v.key} expected={v.expected} observed={v.observed}")
        for message in self.incomplete:
            print("INCOMPLETE", message)
        code = self.exit_code
        print(f"{stage}: {('PASS/0', 'FAIL/1', 'INCOMPLETE/2')[code]}")
        return code


def table(text, columns, key, duplicate_code="SAMPLE_KEY"):
    rows = csv.DictReader(io.StringIO(text), delimiter="\t", strict=True)
    fields = rows.fieldnames or []
    if len(set(fields)) != len(fields) or not set(columns) <= set(fields):
        raise Invalid("SCHEMA", f"required {columns}, columns {fields}")
    indexed = {}
    key = (key,) if isinstance(key, str) else key
    for row in rows:
        if any((row.get(k) or "").startswith("#q2:") for k in key):
            continue
        if None in row or any(row.get(k) in (None, "") for k in columns):
            raise Invalid("SCHEMA", "missing or extra field")
        values = tuple(row[k] for k in key)
        idx = values[0] if len(key) == 1 else values
        if idx in indexed:
            raise Invalid(duplicate_code, f"duplicate {idx}")
        indexed[idx] = row
    return indexed


def file_table(path, *args):
    return table(path.read_text(), *args)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs(root, out):
    for name, expected in HASHES.items():
        digest = out.read(name, sha256, root / name)
        if digest is not None:
            out.same("INPUT_HASH", "-", name, expected, digest)
    lines = out.read("raw.sha256", (root / "data/raw.sha256").read_text)
    if lines is not None:
        for line in lines.splitlines():
            if not re.fullmatch(r"[a-f0-9]{64}  data/raw/[\w.-]+\.fastq\.gz", line):
                out.errors.append(Violation("SCHEMA", "-", "raw.sha256", "hash  relative path", line))
                continue
            expected, name = line.split()
            digest = out.read(name, sha256, root / name)
            if digest is not None:
                out.same("INPUT_HASH", "-", name, expected, digest)


def read_qzv(path):
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if re.fullmatch(r"[^/]+/data/per-sample-fastq-counts.tsv", n)]
        if len(names) != 1:
            raise Invalid("SCHEMA", "expected one per-sample count table")
        return table(z.read(names[0]).decode(), ["sample ID", "forward sequence count"], "sample ID")


def input_counts(root, out):
    cols = ("sample-id", "reads_R1", "reads_R2")
    rows = out.read("read_counts", file_table, root / "results/read_counts_R1.tsv", cols, "sample-id")
    if rows is None:
        return None
    out.same("SAMPLE_KEY", "-", "read_counts", SAMPLES, set(rows))
    counts = {}
    for sid, row in rows.items():
        values = [row[c] for c in cols[1:]]
        if not all(re.fullmatch(r"[1-9]\d*", v) for v in values):
            out.errors.append(Violation("READ_COUNT", "-", sid, "positive integers", values))
            continue
        counts[sid] = int(values[0])
        out.same("READ_COUNT", "-", sid, counts[sid], int(values[1]))
    return counts


def quality_histogram(path):
    """Stream validated Phred+33 records. Pool equal lengths within bounded batches."""
    histogram = np.zeros((0, 94), dtype=np.int64)

    def flush(batch):
        nonlocal histogram
        for length in {len(q) for q in batch}:
            if not length:
                raise Invalid("FASTQ_FORMAT", "empty read")
            chunk = b"".join(q for q in batch if len(q) == length)
            scores = np.frombuffer(chunk, dtype=np.uint8).reshape(-1, length)
            if scores.min() < 33 or scores.max() > 126:
                raise Invalid("FASTQ_FORMAT", "invalid Phred+33 or empty read")
            if length > len(histogram):
                histogram = np.pad(histogram, ((0, length - len(histogram)), (0, 0)))
            bins = np.arange(length) * 94 + (scores - 33)
            histogram[:length] += np.bincount(bins.ravel(), minlength=length*94).reshape(length, 94)

    try:
        with gzip.open(path, "rt", encoding="ascii") as f:
            batch = []
            for _, _, qual in FastqGeneralIterator(f):
                batch.append(qual.encode("ascii"))
                if len(batch) >= 512:
                    flush(batch)
                    batch = []
            if batch:
                flush(batch)
    except (ValueError, UnicodeError, gzip.BadGzipFile, EOFError) as e:
        raise Invalid("FASTQ_FORMAT", str(e)) from e
    if not histogram.size:
        raise Invalid("FASTQ_FORMAT", "empty FASTQ")
    return histogram


def statistics(histogram):
    positions = {}
    for pos, row in enumerate(histogram, 1):
        cdf = row.cumsum()
        n = int(cdf[-1])
        values = {}
        for label, numerator in (("q25", 1), ("q50", 2), ("q75", 3)):
            lo, remainder = divmod((n - 1) * numerator, 4)
            a = int(cdf.searchsorted(lo + 1))
            b = int(cdf.searchsorted(min(lo + 2, n)))
            values[label] = Decimal(4*a + remainder*(b-a))/4
        positions[str(pos)] = dict(n_bases=n, **values)
    totals = histogram.sum(axis=0)
    bases = int(totals.sum())
    summary = dict(n_reads=int(histogram[0].sum()), n_bases=bases,
                   mean_q=Decimal(int(totals @ np.arange(94)))/bases,
                   pct_q30=Decimal(int(totals[30:].sum()) * 100)/bases)
    return positions, summary


def compare_values(out, key, expected, observed):
    for name, want in expected.items():
        raw = observed.get(name, "")
        try:
            number = Decimal(raw)
            if isinstance(want, int):
                ok = re.fullmatch(r"\d+", raw) and number == want
            else:
                ok = number.is_finite() and abs(want-number) <= TOLERANCE
        except (InvalidOperation, ValueError):
            ok = False
        if not ok:
            out.errors.append(Violation("QUALITY_VALUE", "-", f"{key}/{name}", want, raw))


def check_quality(root, out, hashes=True):
    if hashes:
        verify_inputs(root, out)
    counts = input_counts(root, out)
    meta = out.read("metadata", file_table, root / "data/meta/metadata.tsv",
                    ("sample-id", "sample_type", "file_r1", "site_code", "site"), "sample-id")
    manifest = out.read("manifest", file_table, root / "data/meta/fq-manifest-single.tsv",
                        ("sample-id", "absolute-filepath"), "sample-id")
    qzv = out.read("reads_summary", read_qzv, root / "work/qzv/reads_summary.qzv")
    pos = out.read("quality_by_position", file_table, root / "results/quality_by_position.tsv",
                   ("sample-id", "position", "n_bases", "q25", "q50", "q75"),
                   ("sample-id", "position"), "QUALITY_KEY")
    overall = out.read("quality_by_sample", file_table, root / "results/quality_by_sample.tsv",
                       ("sample-id", "n_reads", "n_bases", "mean_q", "pct_q30"), "sample-id", "QUALITY_KEY")
    for name, rows in (("metadata", meta), ("manifest", manifest), ("QZV", qzv), ("quality_by_sample", overall)):
        if rows is not None:
            out.same("SAMPLE_KEY", "-", name, SAMPLES, set(rows))
    expected_keys = set()
    complete = True
    for sid in sorted(SAMPLES):
        if meta is None or sid not in meta:
            complete = False
            continue
        row = meta[sid]
        fname = row["file_r1"]
        expected_name = sid + ("_first_launch" if sid.startswith("LUB") else "") + "_R1_001.fastq.gz"
        out.same("SAMPLE_KEY", "-", sid, expected_name, fname)
        out.same("SAMPLE_KEY", "-", sid, "negative_control" if sid == "K_S63" else "soil", row["sample_type"])
        if manifest is not None and sid in manifest:
            out.same("SAMPLE_KEY", "-", sid, fname, Path(manifest[sid]["absolute-filepath"]).name)
        if fname != expected_name:
            complete = False
            continue
        hist = out.read(sid, quality_histogram, root / "data/raw" / fname)
        if hist is None:
            complete = False
            continue
        per_pos, summary = statistics(hist)
        if counts is not None:
            out.same("READ_COUNT", "-", sid, counts.get(sid), summary["n_reads"])
        if qzv is not None and sid in qzv:
            out.same("READ_COUNT", "-", f"{sid}/QZV", str(summary["n_reads"]), qzv[sid]["forward sequence count"])
        if overall is not None and sid in overall:
            compare_values(out, sid, summary, overall[sid])
        for p, expected in per_pos.items():
            key = (sid, p)
            expected_keys.add(key)
            if pos is not None and key in pos:
                compare_values(out, key, expected, pos[key])
            elif pos is not None:
                out.errors.append(Violation("QUALITY_KEY", "-", sid, p, "missing"))
    if complete and pos is not None:
        for extra in sorted(pos.keys() - expected_keys):
            out.errors.append(Violation("QUALITY_KEY", "-", extra, "no extra row", "present"))
    return counts


class ProvenanceLoader(yaml.SafeLoader):
    pass


ProvenanceLoader.add_constructor("!ref", lambda loader, node: loader.construct_scalar(node))
ProvenanceLoader.add_multi_constructor("!", lambda loader, tag, node: loader.construct_scalar(node)
                                      if isinstance(node, yaml.ScalarNode) else loader.construct_sequence(node))


def folded(items):
    result = {}
    for item in items:
        if not isinstance(item, dict) or len(item) != 1 or result.keys() & item.keys():
            raise Invalid("PROVENANCE", "invalid or duplicate action field")
        result.update(item)
    return result


def read_product(path, product):
    _, semantic_type, output, filename = product
    with zipfile.ZipFile(path) as z:
        roots = [n.split("/")[0] for n in z.namelist() if re.fullmatch(r"[^/]+/metadata.yaml", n)]
        if len(roots) != 1 or len(z.namelist()) != len(set(z.namelist())):
            raise Invalid("SCHEMA", "archive must have a single root and unique paths")
        aid = roots[0]
        uuid.UUID(aid)
        def get(member):
            try:
                return z.read(f"{aid}/{member}")
            except KeyError as e:
                raise Invalid("SCHEMA", f"missing archive member {member}") from e
        meta = yaml.safe_load(get("metadata.yaml"))
        if meta.get("uuid") != aid or meta.get("type") != semantic_type:
            raise Invalid("SCHEMA", f"wrong type or root UUID, expected {semantic_type}")
        prov = yaml.load(get("provenance/action/action.yaml"), Loader=ProvenanceLoader)
        payload = get("data/" + filename)
    action = prov["action"]
    execution = prov["execution"]["uuid"]
    uuid.UUID(execution)
    env = prov["environment"]
    versions = (str(env["framework"]["version"]), str(env["plugins"]["dada2"]["version"]))
    return (aid, execution, action["type"], str(action["plugin"]).split(":")[-1],
            action["action"], action["output-name"], folded(action["inputs"]),
            folded(action["parameters"]), versions), payload


def parse_sequences(raw):
    records = {}
    current = None
    for line in raw.decode("ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].strip()
            if not current or current in records:
                raise Invalid("ASV_KEY", f"duplicate/empty FASTA ID: {current}")
            records[current] = ""
        elif line:
            if current is None or not re.fullmatch("[ACGT]+", line):
                raise Invalid("SCHEMA", "invalid representative DNA sequence")
            records[current] += line
    if not records or len(set(records.values())) != len(records):
        raise Invalid("ASV_KEY", "empty FASTA or multiple IDs for the same DNA")
    return records


@dataclass
class Run:
    prov: list
    samples: set
    features: set
    matrix: dict
    dna: dict
    stats: dict


def load_run(directory):
    products = [read_product(directory / p[0], p) for p in PRODUCTS]
    try:
        with h5py.File(io.BytesIO(products[0][1]), "r") as f:
            for axis in ("observation", "sample"):
                values = list(f[axis + "/ids"][:])
                if len(values) != len(set(values)):
                    raise Invalid("ASV_KEY", f"duplicate BIOM {axis} ID")
            t = biom.Table.from_hdf5(f)
    except Invalid:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Invalid("SCHEMA", f"invalid BIOM/HDF5: {e}") from e
    samples = list(t.ids(axis="sample"))
    features = list(t.ids(axis="observation"))
    if len(set(samples)) != len(samples) or len(set(features)) != len(features):
        raise Invalid("ASV_KEY", "duplicate matrix axis")
    matrix = {}
    sparse = t.matrix_data.tocoo()
    for i, j, x in zip(sparse.row, sparse.col, sparse.data):
        if not math.isfinite(x) or x < 0 or x != int(x):
            raise Invalid("COUNT_VALUE", f"{features[i]}/{samples[j]}: {x}")
        if x:
            matrix[(features[i], samples[j])] = int(x)
    text = products[2][1].decode()
    sample_key = text.split("\t", 1)[0]
    if sample_key not in ("sample-id", "#SampleID", "sample ID", "id"):
        raise Invalid("SCHEMA", "invalid DADA2 sample ID header: " + sample_key)
    rows = table(text, (sample_key, *STATS), sample_key)
    stats = {}
    for sid, row in rows.items():
        if not all(re.fullmatch(r"\d+", row[k]) for k in STATS):
            raise Invalid("STATS_COUNT", f"invalid count for {sid}")
        stats[sid] = tuple(int(row[k]) for k in STATS)
    return Run([p[0] for p in products], set(samples), set(features), matrix,
               parse_sequences(products[1][1]), stats)


@dataclass
class Dataset:
    runs: dict
    counts: dict
    control: str
    trim: tuple
    input_uuid: str


def run_metrics(data):
    for a, b in combinations(data.runs, 2):
        r, s = data.runs[a], data.runs[b]
        def by_dna(run):
            return {(run.dna[f], sid): n for (f,sid),n in run.matrix.items() if f in run.dna}
        left, right = by_dna(r), by_dna(s)
        deltas = [abs(left.get(k,0)-right.get(k,0)) for k in left.keys() | right.keys()]
        statdiff = sum(sum(x != y for x,y in zip(r.stats[sid],s.stats[sid]))
                       for sid in r.stats.keys() & s.stats.keys())
        yield (a,b), (len(set(r.dna.values()) ^ set(s.dna.values())),
                       sum(d != 0 for d in deltas), max(deltas, default=0), statdiff)


def count_value(value):
    """BIOM stores floats; accept integral values, never bool, NaN or infinity."""
    return (not isinstance(value, (bool, np.bool_))
            and isinstance(value, (int, float, np.integer, np.floating))
            and math.isfinite(value) and value >= 0 and value == int(value))


def validate(data):
    """Pure contract for the objects that could be read. Availability is handled by I/O."""
    out = Audit()
    if (len(data.trim) != 2 or not all(count_value(x) for x in data.trim)
            or not 0 <= data.trim[0] < data.trim[1] <= 251):
        return [Violation("SCHEMA", "-", "trim", "0 <= L < T <= 251", data.trim)]
    if (not data.counts or data.control not in data.counts
            or not all(count_value(x) and x > 0 for x in data.counts.values())):
        return [Violation("READ_COUNT", "-", "counts", "positive integers and control", data.counts)]
    execs = {}
    comparable = {}
    left, right = data.trim
    params = dict(PARAMS, trim_left=left, trunc_len=right)
    for r, run in data.runs.items():
        safe_to_compare = True
        executions = set()
        ids = set()
        out.same("PROVENANCE", r, "output count", 3, len(run.prov))
        for product, provenance in zip(PRODUCTS, run.prov):
            if len(provenance) != 9:
                out.errors.append(Violation("PROVENANCE", r, "record", "9 fields", provenance))
                continue
            aid, execution, kind, plugin, method, output, inputs, observed, versions = provenance
            try:
                uuid.UUID(aid)
                uuid.UUID(execution)
            except (ValueError, TypeError, AttributeError):
                out.errors.append(Violation("PROVENANCE", r, "UUID", "valid UUIDs", (aid, execution)))
            ids.add(aid)
            executions.add(execution)
            out.same("PROVENANCE", r, "action", ("method","dada2","denoise_single",product[2]),
                     (kind,plugin,method,output))
            out.same("PROVENANCE", r, "input", {"demultiplexed_seqs":data.input_uuid}, inputs)
            out.same("PROVENANCE", r, "versions", ("2025.7.0","2025.7.0"), versions)
            out.same("PROVENANCE", r, "parameters", params, observed)
        out.same("PROVENANCE", r, "artifact UUID count", 3, len(ids))
        out.same("PROVENANCE", r, "execution UUID count", 1, len(executions))
        for execution in executions:
            if execution in execs:
                out.errors.append(Violation("REUSED_RUN", r, "execution UUID", "distinct", execs[execution]))
            execs[execution] = r
        out.same("ASV_KEY", r, "FASTA vs matrix", run.features, set(run.dna))
        out.same("SAMPLE_KEY", r, "table samples", set(data.counts), run.samples)
        out.same("SAMPLE_KEY", r, "stats samples", set(data.counts), set(run.stats))
        if not run.features or len(set(run.dna.values())) != len(run.dna):
            out.errors.append(Violation("ASV_KEY", r, "sequences", "nonempty unique DNA", run.dna))
        for feature, seq in run.dna.items():
            out.same("ASV_LENGTH", r, feature, right-left, len(seq))
            if not re.fullmatch("[ACGT]+", seq):
                out.errors.append(Violation("ASV_KEY", r, feature, "DNA alphabet ACGT", seq))
        sums = Counter()
        for (feature, sid), n in run.matrix.items():
            if feature not in run.features or sid not in run.samples:
                out.errors.append(Violation("ASV_KEY", r, (feature, sid), "matrix axes", "unknown key"))
            if not count_value(n):
                out.errors.append(Violation("COUNT_VALUE", r, (feature, sid), "nonnegative integer", n))
                safe_to_compare = False
                continue
            sums[sid] += n
        for sid, counters in run.stats.items():
            if len(counters) != 4 or not all(count_value(n) for n in counters):
                out.errors.append(Violation("STATS_COUNT", r, sid, "four nonnegative integers", counters))
                safe_to_compare = False
                continue
            inp, filt, den, nc = counters
            out.same("READ_COUNT", r, sid, data.counts.get(sid), inp)
            if not inp >= filt >= den >= nc >= 0 or (sid != data.control and nc == 0):
                out.errors.append(Violation("STATS_COUNT", r, sid, "nonincreasing; soil > 0", run.stats[sid]))
            out.same("READ_BALANCE", r, sid, nc, sums[sid])
        if safe_to_compare:
            comparable[r] = run
    for pair, metric in run_metrics(Dataset(comparable, data.counts, data.control, data.trim, data.input_uuid)):
        if any(metric):
            out.errors.append(Violation("REPEAT_DIFFERENCE", str(pair), "ASV/cells/max_reads/stats", (0,0,0,0), metric))
    return out.errors


def read_choice(path):
    content = path.read_text()
    values = []
    for name in ("trim-left","trunc-len"):
        fields = re.findall(rf"^{name}:\s*(.*?)\s*$", content, re.M)
        if not fields or fields == ["<число>"]:
            raise FileNotFoundError(f"решение человека: {name}")
        if len(fields) != 1 or not fields[0].isdigit():
            raise Invalid("SCHEMA", f"invalid {name}")
        values.append(int(fields[0]))
    if not 0 <= values[0] < values[1] <= 251:
        raise Invalid("SCHEMA", f"invalid trim range: {values}")
    return tuple(values)


def read_repeats(root, out, counts, choice, input_uuid=INPUT_ID, control="K_S63", paths=RUNS, verbose=True):
    runs = {}
    for i, name in enumerate(paths, 1):
        item = out.read(name, load_run, root / name)
        if item is not None:
            runs[i] = item
    if counts is None or choice is None:
        return
    data = Dataset(runs, counts, control, choice, input_uuid)
    out.errors.extend(validate(data))
    for i, r in runs.items():
        for sid, (n, _, _, kept) in sorted(r.stats.items()):
            if n:
                if verbose:
                    nonzero = sum(s == sid for _,s in r.matrix)
                    print(f"run={i} {sid}: input={n} non-chimeric={kept} retention={100*kept/n:.2f}% ASV={nonzero}")
                if sid != control and 0 < 100*kept < 30*n:
                    out.incomplete.append(f"run {i}: {sid} retention below 30%; human review")
    if verbose:
        for pair, metric in run_metrics(data):
            print(f"runs={pair[0]}-{pair[1]} different_sequences/cells/max_delta_reads/stats={metric}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--self-test", action="store_true")
    modes.add_argument("--quality-only", action="store_true")
    modes.add_argument("--example-dir", type=Path, help="read an explicitly synthetic three-run case")
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="project root, input paths are read-only")
    args = parser.parse_args()
    from denoise_examples import self_test, COUNTS, INPUT, CONTROL
    root = args.root.resolve()
    start = time.monotonic()
    out = Audit()
    if args.example_dir:
        read_repeats(args.example_dir, out, COUNTS, (1,5), INPUT, CONTROL, ("1","2","3"))
        return out.finish("synthetic example")
    if not args.quality_only:
        failures = self_test(sys.modules[__name__])
        out.errors.extend(Violation("SELF_TEST", "-", name, "expected verdict", seen) for name,seen in failures)
        if args.self_test:
            return out.finish("self-test")
    counts = check_quality(root, out)
    if not args.quality_only:
        choice = out.read("denoise_choice", read_choice, root/"memory/denoise_choice.md")
        read_repeats(root, out, counts, choice, paths=RUNS)
    print(f"elapsed_seconds={time.monotonic()-start:.2f}")
    return out.finish("quality-only" if args.quality_only else "full")


if __name__ == "__main__":
    sys.path.insert(0, str(PROJECT / "src"))
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"VERIFIER_ERROR FAIL/1: {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)

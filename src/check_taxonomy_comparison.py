"""Read-only independent numerical verification and final acceptance of comparison.

This module does not import the calculation core or runner. Shared code supplies
file readers, the output schema and source identity, never expected metrics.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import Counter
from pathlib import Path

from qiime_runtime import selected_python

selected_python()

from taxonomy_compare_io import (  # noqa: E402
    CATEGORIES, MODELS, PAIRS, RANKS, RULES, SCHEMAS,
    load_dataset, sha256, verify_code_identity,
)

ROOT = Path(__file__).resolve().parent.parent
TOL = 1e-12
PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
UNKNOWN = {"", "Unassigned", "Unknown_Family"}
ROW_KEYS = {
    "assignments.tsv": ("feature_id",),
    "coverage.tsv": ("model", "rank"),
    "coverage_per_sample.tsv": ("model", "sample", "rank"),
    "pairwise_match.tsv": ("model_a", "model_b", "rank"),
    "pairwise_match_per_sample.tsv": ("model_a", "model_b", "sample", "rank"),
    "profiles.tsv": ("model", "sample", "rank", "key"),
    "distances.tsv": ("model_a", "model_b", "sample", "rank"),
    "summary.tsv": ("model_a", "model_b", "rank"),
    "manifest.tsv": ("item",),
}


def issue(code, detail, missing=False):
    return {"code": code, "detail": str(detail), "missing": missing}


def parse_independently(raw):
    """Place labels at explicit domain-through-species positions; never shift ranks."""
    if not isinstance(raw, str):
        raise ValueError("Taxon must be a string")
    labels = [None] * len(PREFIXES)
    if raw.strip() in {"", "Unassigned"}:
        return labels
    seen = set()
    for segment in raw.split(";"):
        segment = segment.strip()
        if not segment:
            continue
        prefix = segment[:3]
        if prefix not in PREFIXES:
            raise ValueError(f"Unexpected rank fragment: {segment!r}")
        if prefix in seen:
            raise ValueError(f"Duplicate rank: {prefix}")
        seen.add(prefix)
        label = segment[3:].strip()
        labels[PREFIXES.index(prefix)] = None if label in UNKNOWN else label
    return labels


def independent_key(labels, rank_index):
    position = rank_index + 1
    if labels[position] is None:
        return None
    return ";".join(PREFIXES[i] + labels[i] for i in range(position + 1)
                    if labels[i] is not None)


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def independent_expected(dataset):
    """Compute reference rows from raw strings and counts with independent loops.

    Pure function: no file access or mutation. The input contract is checked by
    the reader; this also supports small hand-computed datasets in unit tests.
    """
    table = dataset["table"]
    features = sorted(table)
    soils = sorted(dataset["soil_samples"])
    taxa = dataset["taxonomies"]
    soil_features = [fid for fid in features if sum(table[fid][s] for s in soils) > 0]
    if not soil_features:
        raise ValueError("No ASV with positive soil abundance")
    totals = {s: sum(table[fid][s] for fid in features) for s in soils}
    if any(t <= 0 for t in totals.values()):
        raise ValueError("Soil read total must be positive")
    parsed = {m: {fid: parse_independently(taxa[m][fid]) for fid in features} for m in MODELS}
    keys = {m: {fid: [independent_key(parsed[m][fid], i) for i in range(len(RANKS))]
                for fid in features} for m in MODELS}
    outputs = {name: [] for name in SCHEMAS if name != "manifest.tsv"}
    soil_set = set(soil_features)
    for fid in features:
        row = {"feature_id": fid, "soil_asv": int(fid in soil_set)}
        for model in MODELS:
            row[f"{model}_taxon"] = taxa[model][fid]
            for i, rank in enumerate(RANKS):
                row[f"{model}_{rank}"] = parsed[model][fid][i + 1] or ""
                row[f"{model}_{rank}_key"] = keys[model][fid][i] or "UNCLASSIFIED"
                row[f"{model}_{rank}_assigned"] = int(keys[model][fid][i] is not None)
        outputs["assignments.tsv"].append(row)

    profiles = {}
    covered = {}
    for model in MODELS:
        for i, rank in enumerate(RANKS):
            known = sum(keys[model][fid][i] is not None for fid in soil_features)
            n = len(soil_features)
            outputs["coverage.tsv"].append({
                "model": model, "rank": rank, "assigned_asv": known,
                "unassigned_asv": n - known, "total_asv": n,
                "assigned_asv_frac": known / n, "unassigned_asv_frac": (n - known) / n,
            })
            for sample in soils:
                counts = Counter()
                for fid in features:
                    amount = table[fid][sample]
                    if amount:
                        counts[keys[model][fid][i] or "UNCLASSIFIED"] += amount
                total = totals[sample]
                profiles[(model, sample, rank)] = counts
                unassigned = counts.get("UNCLASSIFIED", 0)
                covered[(model, sample, rank)] = (total - unassigned) / total
                outputs["coverage_per_sample.tsv"].append({
                    "model": model, "sample": sample, "rank": rank,
                    "assigned_reads": total - unassigned, "unassigned_reads": unassigned,
                    "total_reads": total, "assigned_reads_frac": (total - unassigned) / total,
                    "unassigned_reads_frac": unassigned / total,
                })
                for key in sorted(counts):
                    outputs["profiles.tsv"].append({
                        "model": model, "sample": sample, "rank": rank, "key": key,
                        "count": counts[key], "total_reads": total, "frac": counts[key] / total,
                    })

    for first, second in PAIRS:
        for i, rank in enumerate(RANKS):
            membership = {}
            for fid in soil_features:
                a, b = keys[first][fid][i], keys[second][fid][i]
                if a is None:
                    category = "neither" if b is None else "only_b"
                elif b is None:
                    category = "only_a"
                else:
                    category = "both_match" if a == b else "both_diff"
                membership[fid] = category
            counts = Counter(membership.values())
            n = len(soil_features)
            both = counts["both_match"] + counts["both_diff"]
            outputs["pairwise_match.tsv"].append({
                "model_a": first, "model_b": second, "rank": rank,
                **{c: counts[c] for c in CATEGORIES}, "total_asv": n,
                **{f"{c}_frac": counts[c] / n for c in CATEGORIES},
                "both_assigned_asv": both, "cond_match_asv_frac": ratio(counts["both_match"], both),
            })
            distances = []
            for sample in soils:
                weighted = Counter()
                for fid, category in membership.items():
                    weighted[category] += table[fid][sample]
                total = totals[sample]
                both_reads = weighted["both_match"] + weighted["both_diff"]
                outputs["pairwise_match_per_sample.tsv"].append({
                    "model_a": first, "model_b": second, "sample": sample, "rank": rank,
                    **{f"{c}_reads": weighted[c] for c in CATEGORIES}, "total_reads": total,
                    **{f"{c}_reads_frac": weighted[c] / total for c in CATEGORIES},
                    "both_assigned_reads": both_reads,
                    "cond_match_reads_frac": ratio(weighted["both_match"], both_reads),
                })
                a = profiles[(first, sample, rank)]
                b = profiles[(second, sample, rank)]
                distance = math.fsum(abs(a.get(k, 0) - b.get(k, 0)) / total
                                     for k in sorted(set(a) | set(b))) / 2
                distances.append(distance)
                outputs["distances.tsv"].append({
                    "model_a": first, "model_b": second, "sample": sample, "rank": rank,
                    "distance": distance,
                    "assigned_reads_frac_a": covered[(first, sample, rank)],
                    "assigned_reads_frac_b": covered[(second, sample, rank)],
                })
            outputs["summary.tsv"].append({
                "model_a": first, "model_b": second, "rank": rank, "n_samples": len(soils),
                "min": min(distances), "median": statistics.median(distances), "max": max(distances),
            })
    return outputs


def read_rows(path, header, problems):
    if not path.is_file():
        problems.append(issue("MISSING_OUTPUT", path.name, True))
        return None
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t", strict=True)
            actual = next(reader, None)
            if actual != list(header):
                problems.append(issue("OUTPUT_HEADER", f"{path.name}: {actual!r}"))
                return None
            rows = []
            for line_no, values in enumerate(reader, 2):
                if len(values) != len(header):
                    problems.append(issue("OUTPUT_ROW_WIDTH", f"{path.name}:{line_no}"))
                    continue
                rows.append(dict(zip(header, values)))
            validate_domains(path.name, rows, problems)
            return rows
    except (OSError, UnicodeError, csv.Error) as exc:
        problems.append(issue("OUTPUT_READ", f"{path.name}: {exc}"))
        return None


def validate_domains(name, rows, problems):
    """Validate intrinsic constraints even when another input is unavailable."""
    seen = set()
    for line_no, row in enumerate(rows, 2):
        key = tuple(row[c] for c in ROW_KEYS[name])
        if key in seen:
            problems.append(issue("OUTPUT_DUPLICATE_KEY", f"{name}:{line_no}: {key}"))
        seen.add(key)
        if any(not value for value in key):
            problems.append(issue("OUTPUT_EMPTY_KEY", f"{name}:{line_no}"))
        for column, value in row.items():
            if column in {"model", "model_a", "model_b"} and value not in MODELS:
                problems.append(issue("OUTPUT_MODEL", f"{name}:{line_no}: {value!r}"))
            if column == "rank" and value not in RANKS:
                problems.append(issue("OUTPUT_RANK", f"{name}:{line_no}: {value!r}"))
            integer = (column in {*CATEGORIES, "count", "n_samples", "soil_asv"}
                       or column.endswith(("_reads", "_asv", "_assigned")))
            numeric = column.endswith("_frac") or column in {"frac", "distance", "min", "median", "max"}
            if integer:
                valid = bool(re.fullmatch(r"0|[1-9][0-9]*", value))
                if valid and (column == "soil_asv" or column.endswith("_assigned")):
                    valid = int(value) in (0, 1)
                if not valid:
                    problems.append(issue("OUTPUT_COUNT_DOMAIN", f"{name}:{line_no}: {column}={value!r}"))
            elif numeric:
                if value == "NA" and column.startswith("cond_match_"):
                    continue
                try:
                    parsed = float(value)
                    valid = math.isfinite(parsed) and 0 <= parsed <= 1
                except ValueError:
                    valid = False
                if not valid:
                    problems.append(issue("OUTPUT_FRACTION_DOMAIN", f"{name}:{line_no}: {column}={value!r}"))


def same_value(actual, expected):
    if expected is None:
        return actual == "NA"
    if isinstance(expected, bool):
        return actual == str(int(expected))
    if isinstance(expected, int):
        return bool(re.fullmatch(r"0|[1-9][0-9]*", actual)) and int(actual) == expected
    if isinstance(expected, float):
        try:
            value = float(actual)
        except (ValueError, TypeError):
            return False
        return math.isfinite(value) and abs(value - expected) <= TOL
    return actual == expected


def compare_rows(name, actual, expected, key_columns, problems, skip_columns=()):
    if actual is None:
        return
    indexed = {}
    for row in actual:
        key = tuple(row[k] for k in key_columns)
        if key in indexed:
            problems.append(issue("OUTPUT_DUPLICATE_KEY", f"{name}: {key}"))
        indexed[key] = row
    wanted = {tuple(str(r[k]) for k in key_columns): r for r in expected}
    if set(indexed) != set(wanted):
        missing = sorted(set(wanted) - set(indexed))
        extra = sorted(set(indexed) - set(wanted))
        problems.append(issue("OUTPUT_KEYS", f"{name}: missing={len(missing)}, extra={len(extra)}; "
                              f"examples={missing[:2] + extra[:2]}"))
    for key in sorted(set(indexed) & set(wanted)):
        row, reference = indexed[key], wanted[key]
        for column in SCHEMAS[name]:
            if column not in skip_columns and not same_value(row[column], reference[column]):
                problems.append(issue("OUTPUT_VALUE", f"{name}: {key}, {column}: "
                                      f"{row[column]!r} != {reference[column]!r}"))


def read_json(path, problems, *, missing_code="MISSING_OUTPUT"):
    if not path.is_file():
        problems.append(issue(missing_code, path.name, True))
        return None
    try:
        def pairs(items):
            obj = {}
            for key, value in items:
                if key in obj:
                    raise ValueError(f"Duplicate JSON key: {key}")
                obj[key] = value
            return obj
        def bad_constant(value):
            raise ValueError(f"Non-finite JSON number: {value}")
        value = json.loads(path.read_text(encoding="utf-8"),
                           object_pairs_hook=pairs, parse_constant=bad_constant)
        if not isinstance(value, dict):
            raise ValueError("Expected a JSON object")
        return value
    except (OSError, UnicodeError, ValueError) as exc:
        problems.append(issue("JSON_FORMAT", f"{path.name}: {exc}"))
        return None


def check_tables(input_dir, out_dir):
    """Validate readable outputs even if another required file is missing."""
    input_dir, out_dir = Path(input_dir), Path(out_dir)
    problems = []
    try:
        dataset, input_problems = load_dataset(input_dir)
    except Exception as exc:
        dataset, input_problems = {}, [issue("INPUT_READ", f"{type(exc).__name__}: {exc}")]
    problems.extend(input_problems)
    actual = {name: read_rows(out_dir / name, header, problems)
              for name, header in SCHEMAS.items()}
    for name in ("README.md", "normalization_rules.txt"):
        path = out_dir / name
        if not path.is_file():
            problems.append(issue("MISSING_OUTPUT", name, True))
            continue
        try:
            contents = path.read_text(encoding="utf-8")
            if not contents.strip():
                problems.append(issue("OUTPUT_EMPTY", name))
            if name == "normalization_rules.txt" and contents != RULES:
                problems.append(issue("NORMALIZATION_RULES", name))
        except (OSError, UnicodeError) as exc:
            problems.append(issue("OUTPUT_READ", f"{name}: {exc}"))
    run = read_json(out_dir / "run.json", problems)
    code_version = None
    if run is not None:
        try:
            problems.extend(verify_code_identity(run["code"]))
            code_version = run["code"]["version"]
            if not isinstance(code_version, str):
                raise ValueError("code.version must be a string")
            if run["inputs"] != dataset.get("inputs") and not input_problems:
                problems.append(issue("RUN_INPUTS", "run inputs differ from actual inputs"))
            names = set(SCHEMAS) | {"README.md", "normalization_rules.txt"}
            recorded = run["outputs"]
            if not isinstance(recorded, dict) or set(recorded) != names:
                problems.append(issue("RUN_OUTPUT_SET", "wrong output manifest keys"))
            if isinstance(recorded, dict):
                for name in names & set(recorded):
                    path = out_dir / name
                    digest = recorded[name]
                    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
                        problems.append(issue("RUN_OUTPUT_HASH_FORMAT", name))
                    elif path.is_file() and digest != sha256(path):
                        problems.append(issue("RUN_OUTPUT_HASH", name))
        except Exception as exc:
            problems.append(issue("RUN_FORMAT", f"{type(exc).__name__}: {exc}"))
    # Detect malformed taxonomy even if another input is missing.
    for model, taxa in dataset.get("taxonomies", {}).items():
        for fid, raw in taxa.items():
            try:
                parse_independently(raw)
            except ValueError as exc:
                problems.append(issue("TAXON_FORMAT", f"{model}/{fid}: {exc}"))
    if not input_problems:
        try:
            expected = independent_expected(dataset)
            for name, rows in expected.items():
                compare_rows(name, actual[name], rows, ROW_KEYS[name], problems)
            manifest = [dict(item=item, **record, code_version=code_version or "")
                        for item, record in sorted(dataset["inputs"].items())]
            compare_rows("manifest.tsv", actual["manifest.tsv"], manifest,
                         ROW_KEYS["manifest.tsv"], problems,
                         skip_columns=("code_version",) if code_version is None else ())
        except Exception as exc:
            problems.append(issue("ORACLE_INPUT", f"{type(exc).__name__}: {exc}"))
    return problems, dataset, run


def verdict(problems):
    if any(not p.get("missing", False) for p in problems):
        return 1
    return 2 if problems else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=ROOT)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "results/taxonomy_comparison/main")
    parser.add_argument("--tables-only", action="store_true",
                        help="Verify numerical tables only; not final stage acceptance")
    args = parser.parse_args(argv)
    try:
        problems, dataset, run = check_tables(args.input_dir, args.out_dir)
        if not problems and not args.tables_only:
            from taxonomy_compare_acceptance import check_evidence
            problems.extend(check_evidence(args.out_dir, dataset, run))
    except Exception as exc:
        problems = [issue("CHECK_ERROR", f"{type(exc).__name__}: {exc}")]
    code = verdict(problems)
    for problem in problems:
        print(f"{problem['code']}: {problem['detail']}")
    label = ("PASS", "FAIL", "INCOMPLETE")[code]
    scope = "tables only; stage acceptance not checked" if args.tables_only else "full stage acceptance"
    print(f"{label}/{code} — {scope}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())

"""Independently check the saved domain audit against raw QZA taxonomies.

Uses the existing read-only QZA loader, but no domain-audit calculation code.
Writes only the requested verification JSON; never changes scientific inputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
from taxonomy_compare_io import ROOT, MODELS, load_dataset, sha256

DOMAINS = ("Bacteria", "Archaea", "Eukaryota", "unknown", "unexpected")


def records(path, key):
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    found = {row[key]: row for row in rows}
    if len(found) != len(rows):
        raise ValueError(f"Duplicate {key}: {path}")
    return found


def domain(raw):
    """Strictly classify the actual observed syntax; fail on unfamiliar syntax."""
    if raw.strip() in ("", "Unassigned"):
        return "unknown", "", "empty_or_Unassigned"
    tokens = [part.strip() for part in raw.split(";") if part.strip()]
    values = [m.group(1).strip() for token in tokens
              if (m := re.fullmatch(r"(?:d__|D_0__)(.*)", token))]
    if len(values) != 1 or values[0] not in DOMAINS[:3]:
        raise ValueError(f"Unfamiliar domain syntax, requires review: {raw!r}")
    if any(not re.fullmatch(r"(?:[dpcofgs]__|D_[0-6]__).*$", t) for t in tokens):
        raise ValueError(f"Unprefixed taxonomy token: {raw!r}")
    return values[0], values[0], "recognized_domain"


def flags(values):
    kinds = set(values)
    common = kinds == {"Bacteria"}
    excluded = bool(kinds & {"Eukaryota", "Archaea"})
    bacterial = "Bacteria" in kinds
    out = {
        "common_bacteria": common,
        "any_eukaryota_or_archaea": excluded,
        "any_bacteria": bacterial,
        "any_bacteria_but_disagreement_or_unknown": bacterial and not common,
        "any_bacteria_with_eukaryota_or_archaea": bacterial and excluded,
        "any_unknown_domain": "unknown" in kinds,
        "any_unexpected_domain": "unexpected" in kinds,
    }
    partition = ("any_eukaryota_or_archaea" if excluded else
                 "common_bacteria" if common else
                 "bacteria_and_unknown_only" if kinds <= {"Bacteria", "unknown"}
                 and bacterial else "other")
    return {k: int(v) for k, v in out.items()}, partition


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path,
                        default=ROOT / "results/taxonomy_barplots/preflight_20261118")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    folder = args.audit_dir.resolve()
    output = args.output or folder / "domain_independent_check.json"
    if output.exists():
        parser.error("Verification output already exists; choose a new --output")
    audit = json.loads((folder / "domain_audit.json").read_text())
    before = {p: sha256(ROOT / p) for p in audit["sources_before"]}
    errors = []
    checked = 0

    def check(ok, detail):
        nonlocal checked
        checked += 1
        if not ok:
            errors.append(detail)

    def equals(actual, expected, detail):
        check(actual == expected, f"{detail}: {actual!r} != {expected!r}")

    try:
        equals(before, audit["sources_before"], "Inputs vs audit start")
        equals(before, audit["sources_after"], "Inputs vs audit end")
        equals(audit["inputs_unchanged"], True, "Recorded unchanged flag")
        for path, digest in audit["outputs"].items():
            equals(sha256(ROOT / path), digest, f"Audit output hash {path}")
        ds, issues = load_dataset(ROOT)
        equals(issues, [], "Read-only raw QZA input validation")
        if issues:
            raise ValueError("Cannot continue with invalid QZA inputs")
        table = ds["table"]
        soils = ds["soil_samples"]
        samples = ds["samples"]
        equals(samples, audit["samples"], "Samples")
        equals(soils, audit["soil_samples"], "Soils")
        saved = records(folder / "domain_assignments.tsv", "feature_id")
        prior = records(ROOT / "results/taxonomy_comparison/finish_20261116/main/assignments.tsv", "feature_id")
        per_sample = records(folder / "domain_by_sample.tsv", "sample_id")
        equals(set(saved), set(table), "Domain audit ASV set")
        equals(set(prior), set(table), "Original comparison ASV set")
        equals(set(per_sample), set(samples), "Domain sample set")
        domains, booleans, partitions, soil_ids = {}, {}, {}, set()
        for fid, counts in table.items():
            row = saved[fid]
            soil_reads = sum(counts[s] for s in soils)
            if soil_reads:
                soil_ids.add(fid)
            values = []
            for model in MODELS:
                raw = ds["taxonomies"][model][fid]
                equals(row[model + "_taxon"], raw, f"Raw QZA taxonomy {model}/{fid}")
                equals(prior[fid][model + "_taxon"], raw, f"Comparison raw taxonomy {model}/{fid}")
                parsed = domain(raw)
                for suffix, value in zip(("domain", "domain_token", "domain_reason"), parsed):
                    equals(row[f"{model}_{suffix}"], value, f"Domain {model}/{fid}/{suffix}")
                domains[model, fid] = parsed[0]
                values.append(parsed[0])
            for field, value in {"soil_asv": int(soil_reads > 0), "all_reads": sum(counts.values()),
                                 "soil_reads": soil_reads, "control_reads": counts["K_S63"]}.items():
                equals(int(row[field]), value, f"ASV counts {fid}/{field}")
            booleans[fid], partitions[fid] = flags(values)
            for flag, value in booleans[fid].items():
                equals(int(row[flag]), value, f"ASV flag {fid}/{flag}")
            equals(row["exclusive_partition"], partitions[fid], f"ASV partition {fid}")
        for key, ids in (("asv_all", set(table)), ("asv_soil", soil_ids)):
            summary = {"asv": len(ids), "per_model": {
                model: {d: sum(domains[model, f] == d for f in ids) for d in DOMAINS}
                for model in MODELS}, "overlapping_flags": {
                flag: sum(booleans[f][flag] for f in ids) for flag in booleans[next(iter(table))]},
                "exclusive_partition": dict(Counter(partitions[f] for f in ids))}
            equals(summary, audit[key], key)
        with (ROOT / "results/exact_match_reproduction/provenance/samples.tsv").open() as stream:
            first = {r["sample-id"]: r for r in csv.DictReader(stream, delimiter="\t") if r["run"] == "1"}
        counts_summary = {}
        for sid in samples:
            row = per_sample[sid]
            total = sum(v[sid] for v in table.values())
            observed = {f for f, v in table.items() if v[sid] > 0}
            expected = {"total_retained_reads": total, "observed_asv": len(observed)}
            equals(row["sample_type"], "soil" if sid in soils else "negative_control", f"Role {sid}")
            equals(total, int(first[sid]["non-chimeric"]), f"DADA2 total {sid}")
            equals(len(observed), int(first[sid]["nonzero_ASV"]), f"DADA2 ASVs {sid}")
            counts_summary[sid] = {"total_reads": total}
            for model in MODELS:
                for d in DOMAINS:
                    ids = {f for f in observed if domains[model, f] == d}
                    expected[f"{model}_{d}_reads"] = sum(table[f][sid] for f in ids)
                    expected[f"{model}_{d}_asv"] = len(ids)
                equals(sum(expected[f"{model}_{d}_reads"] for d in DOMAINS), total, f"Read partition {sid}/{model}")
                equals(sum(expected[f"{model}_{d}_asv"] for d in DOMAINS), len(observed), f"ASV partition {sid}/{model}")
                counts_summary[sid][model] = {d: expected[f"{model}_{d}_reads"] for d in DOMAINS}
                pct = 100 * (1 - expected[f"{model}_Bacteria_reads"] / total)
                check(abs(float(row[f"{model}_bacteria_excluded_pct"]) - pct) < 5e-10, f"Excluded percent {sid}/{model}")
            for flag in booleans[next(iter(table))]:
                ids = {f for f in observed if booleans[f][flag]}
                expected[flag + "_reads"] = sum(table[f][sid] for f in ids)
                expected[flag + "_asv"] = len(ids)
            expected["common_excluded_reads"] = total - expected["common_bacteria_reads"]
            expected["zero_common_bacteria"] = int(expected["common_bacteria_reads"] == 0)
            for field, value in expected.items():
                equals(int(row[field]), value, f"Sample {sid}/{field}")
            for field, value in {"common_excluded_pct": 100 * expected["common_excluded_reads"] / total,
                                 "common_retained_pct": 100 * expected["common_bacteria_reads"] / total}.items():
                check(abs(float(row[field]) - value) < 5e-10, f"Sample percent {sid}/{field}")
        for model in MODELS:
            zeros = [s for s in samples if counts_summary[s][model]["Bacteria"] == 0]
            equals(zeros, audit["zero_per_model_bacteria_samples"][model], f"Zero Bacteria samples {model}")
            excluded = [100 * (1 - counts_summary[s][model]["Bacteria"] / counts_summary[s]["total_reads"]) for s in soils]
            for got, expected in zip(audit["soil_per_model_bacteria_excluded_pct_range"][model], [min(excluded), max(excluded)]):
                check(abs(got - expected) < 5e-10, f"Excluded percent range {model}")
        common_counts = {s: sum(table[f][s] for f in table if booleans[f]["common_bacteria"]) for s in samples}
        equals([s for s in samples if common_counts[s] == 0], audit["zero_common_bacteria_samples"], "Zero common Bacteria samples")
        excluded = [100 * (1 - common_counts[s] / counts_summary[s]["total_reads"]) for s in soils]
        for got, expected in zip(audit["soil_common_excluded_pct_range"], [min(excluded), max(excluded)]):
            check(abs(got - expected) < 5e-10, "Common excluded percent range")
        equals(audit["unexpected_tokens"], {m: {} for m in MODELS}, "No unfamiliar domain syntax")
        equals(audit["status"], "PASS", "Recorded audit status")
        equals(audit["violations"], [], "Recorded audit violations")
        equals([counts_summary[f"Kr{i}_S{79+i}"]["sklearn_silva138"]["Bacteria"] for i in range(1, 6)],
               [40, 11, 11, 249, 0], "SILVA Kr1–Kr5 bacterial read counts")
    except Exception as exc:
        errors.append(f"{type(exc).__name__}: {exc}")
    after = {p: sha256(ROOT / p) for p in before}
    equals(after, before, "Inputs unchanged during independent verification")
    result = {"author": "AI4Science-агент", "created_utc": datetime.now(timezone.utc).isoformat(),
              "status": "FAIL" if errors else "PASS", "checks": checked, "violations": errors,
              "method": "Read raw QZA Taxon and integer ASV frequencies via read-only load_dataset; independently parse domain prefixes with regex and re-sum ASVs/reads. No domain audit computation code imported.",
              "scope": "Verify existing audit; no common mask selected or used for downstream profiles.",
              "command": "conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/check_taxonomy_domain_audit.py",
              "checker_sha256": sha256(Path(__file__)), "sources_before": before, "sources_after": after,
              "inputs_unchanged": before == after, "sample_domain_read_counts": locals().get("counts_summary", {})}
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(f"{result['status']}: {checked} checks; {len(errors)} violations; {output}")
    for error in errors[:10]:
        print(error)
    return int(bool(errors))


if __name__ == "__main__":
    raise SystemExit(main())

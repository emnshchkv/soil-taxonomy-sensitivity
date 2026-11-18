"""Independent read-only verification of the article's diagnostic supplements.

Reads the fixed BIOM CSR matrix directly, without biom.Table or importing any
exporter. The saved scientific comparison/filter tables remain read-only.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import io
import json
from pathlib import Path
import re
from statistics import median
import zipfile

import h5py

ROOT = Path(__file__).resolve().parents[1]
OUT = Path("results/article_revision/20261118")
MAIN = Path("results/taxonomy_comparison/finish_20261116/main")
FAMILY = Path("results/taxonomy_barplots/prokaryote_family_20261118_final")
PAIRS = [("sklearn_silva138", "sklearn_gtdb", "SilvaGtdb", "SILVA--GTDB"),
         ("sklearn_silva138", "sklearn_gg2", "SilvaGgTwo", "SILVA--Greengenes2"),
         ("sklearn_gtdb", "sklearn_gg2", "GtdbGgTwo", "GTDB--Greengenes2")]


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix(root):
    """Decode the BIOM observation-major sparse matrix independently."""
    with zipfile.ZipFile(root / "work/qza/ASV_table.qza") as archive:
        members = [n for n in archive.namelist() if n.endswith("/data/feature-table.biom")]
        if len(members) != 1:
            raise ValueError("expected exactly one BIOM matrix")
        with h5py.File(io.BytesIO(archive.read(members[0])), "r") as handle:
            decode = lambda seq: [v.decode() if isinstance(v, bytes) else str(v) for v in seq]
            samples = decode(handle["sample/ids"][:])
            features = decode(handle["observation/ids"][:])
            data = handle["observation/matrix/data"][:]
            indices = handle["observation/matrix/indices"][:]
            indptr = handle["observation/matrix/indptr"][:]
            result = {s: {} for s in samples}
            for i, feature in enumerate(features):
                for k in range(int(indptr[i]), int(indptr[i+1])):
                    value = data[k]
                    if int(value) != value or value < 0:
                        raise ValueError("noninteger or negative BIOM frequency")
                    if value:
                        result[samples[int(indices[k])]][feature] = int(value)
    return result


def tex_rows(path):
    result = []
    for line in path.read_text().splitlines():
        if "&" in line and re.match(r"(?:Kr[1-5]|LUB[2-6])", line.strip()):
            cells = [s.strip().replace(r"\_", "_") for s in line.split("&")]
            cells[-1] = re.sub(r"\\\\.*", "", cells[-1]).strip()
            result.append(cells)
    return result


def inspect(root=ROOT, out=None):
    root = Path(root).resolve()
    out = Path(out).resolve() if out else root / OUT
    errors, count = [], 0
    def check(condition, message):
        nonlocal count
        count += 1
        if not condition:
            errors.append(message)
    def equal_records(name, expected, key):
        actual = rows(out / name)
        key_fn = lambda row: tuple(row[k] for k in key)
        actual_map = {key_fn(r): r for r in actual}
        expected_map = {key_fn(r): r for r in expected}
        check(len(actual_map) == len(actual) and actual_map.keys() == expected_map.keys(), name + ": row keys differ")
        for row_key in expected_map.keys() & actual_map.keys():
            a, e = actual_map[row_key], expected_map[row_key]
            check(a.keys() == e.keys(), name + ": column schema differs")
            for field in e.keys() & a.keys():
                value = e[field]
                if isinstance(value, (int, float, Decimal)):
                    try:
                        same = abs(Decimal(a[field]) - Decimal(str(value))) <= Decimal("1e-12")
                    except Exception:
                        same = False
                else:
                    same = a[field] == value
                check(same, name + ": value differs " + repr(row_key) + " " + field)
    try:
        values = matrix(root)
        metadata = {r["sample-id"]: r["sample_type"] for r in rows(root / "data/meta/metadata.tsv")}
        check(set(values) == set(metadata), "BIOM samples differ from metadata")
        soil = sorted(s for s in values if metadata[s] == "soil")
        check(len(soil) == 10 and metadata.get("K_S63") == "negative_control", "sample roles/count differ")
        soil_ids = set().union(*(set(values[s]) for s in soil))
        assignments = {r["feature_id"]: r for r in rows(root / MAIN / "assignments.tsv")}
        check(soil_ids == {f for f, r in assignments.items() if r["soil_asv"] == "1"}, "soil ASV flags differ from BIOM")
        all_features = set().union(*(set(v) for v in values.values()))
        check(all_features == set(assignments), "assignment features differ from BIOM")

        kr_fields = ["total_reads", "Eukaryota_reads", "unknown_reads", "organelle_reads", "selected_reads"]
        kr = [{"sample": r["sample"], **{k: int(r[k]) for k in kr_fields}}
              for r in rows(root / FAMILY / "filter_counts.tsv")
              if r["model"] == "sklearn_silva138" and r["sample"].startswith("Kr")]
        check(len(kr) == 5, "expected five SILVA Kr samples")
        for r in kr:
            check(r["total_reads"] == sum(values[r["sample"]].values()), "Kr total differs from BIOM: " + r["sample"])
            check(r["total_reads"] == sum(r[k] for k in kr_fields[1:]), "Kr exclusions do not partition reads: " + r["sample"])
        equal_records("kr_filter_reasons.tsv", kr, ["sample"])

        saved_pooled = rows(root / MAIN / "pairwise_match.tsv")
        saved_sample = rows(root / MAIN / "pairwise_match_per_sample.tsv")
        genus, per_sample = [], []
        for a, b, _, _ in PAIRS:
            assigned = {f for f in soil_ids if assignments[f][a+"_род_assigned"] == assignments[f][b+"_род_assigned"] == "1"}
            matching = {f for f in assigned if assignments[f][a+"_род_key"] == assignments[f][b+"_род_key"]}
            raw = [r for r in saved_pooled if (r["model_a"], r["model_b"], r["rank"]) == (a,b,"род")]
            check(len(raw) == 1, "pooled genus selector is not unique")
            raw = raw[0]
            check(int(raw["both_assigned_asv"]) == len(assigned) and int(raw["both_match"]) == len(matching), "pooled genus counts differ from assignment flags")
            group = []
            for sample in soil:
                present = values[sample]
                both_reads = sum(present.get(f,0) for f in assigned)
                match_reads = sum(present.get(f,0) for f in matching)
                saved = [r for r in saved_sample if (r["model_a"],r["model_b"],r["rank"],r["sample"]) == (a,b,"род",sample)]
                check(len(saved) == 1, "sample genus selector is not unique")
                saved = saved[0]
                check(int(saved["both_assigned_reads"]) == both_reads and int(saved["both_match_reads"]) == match_reads, "sample genus reads differ from BIOM")
                fraction = saved["cond_match_reads_frac"]
                check(fraction == "NA" if not both_reads else abs(float(fraction)-match_reads/both_reads)<1e-12, "conditional fraction/NA differs")
                group.append(dict(model_a=a,model_b=b,rank="род",sample=sample,
                                  observed_asv=len(present),both_assigned_asv=len(assigned & set(present)),
                                  both_assigned_reads=both_reads,both_match_asv=len(matching & set(present)),
                                  both_match_reads=match_reads,conditional_read_fraction=fraction))
            per_sample.extend(group)
            defined = [r for r in group if r["both_assigned_reads"]]
            genus.append(dict(model_a=a,model_b=b,rank="род",pooled_both_assigned_asv=len(assigned),
                              pooled_both_match_asv=len(matching),pooled_conditional_asv_fraction=raw["cond_match_asv_frac"],
                              defined_samples=len(defined),conditional_read_median=median(float(r["conditional_read_fraction"]) for r in defined),
                              conditional_sample_asv_min=min(r["both_assigned_asv"] for r in defined),
                              conditional_sample_asv_max=max(r["both_assigned_asv"] for r in defined),
                              conditional_sample_asv_median=median(r["both_assigned_asv"] for r in defined)))
        equal_records("genus_conditional.tsv", genus, ["model_a","model_b"])
        equal_records("genus_conditional_per_sample.tsv", per_sample, ["model_a","model_b","sample"])

        negative = values["K_S63"]
        shared = set(negative) & soil_ids
        raw_count = [r for r in rows(root / "results/read_counts_R1.tsv") if r["sample-id"] == "K_S63"][0]
        with zipfile.ZipFile(root / "work/qza/dada2_stats.qza") as archive:
            data = archive.read(next(n for n in archive.namelist() if n.endswith("/data/stats.tsv"))).decode()
        stats = next(r for r in csv.DictReader([line for line in data.splitlines() if not line.startswith("#q2:")], delimiter="\t")
                     if r["sample-id"] == "K_S63")
        check(int(raw_count["reads_R1"]) == int(raw_count["reads_R2"]) == int(stats["input"]), "control raw/input counts differ")
        check(sum(negative.values()) == int(stats["non-chimeric"]), "control retained reads differ")
        overlap = [dict(sample=s,total_reads=sum(values[s].values()),observed_asv=len(values[s]),
                        shared_control_asv=len(set(values[s]) & set(negative)),
                        shared_control_reads=sum(values[s].get(f,0) for f in negative)) for s in soil]
        control = dict(sample="K_S63",input_reads=int(raw_count["reads_R1"]),retained_reads=sum(negative.values()),
                       observed_asv=len(negative),shared_soil_asv=len(shared),control_only_asv=len(set(negative)-soil_ids),
                       shared_control_reads=sum(negative[f] for f in shared),
                       total_soil_reads=sum(r["total_reads"] for r in overlap),
                       shared_soil_reads=sum(r["shared_control_reads"] for r in overlap))
        control["shared_soil_read_fraction"] = control["shared_soil_reads"]/control["total_soil_reads"]
        equal_records("control_summary.tsv", [control], ["sample"])
        equal_records("control_overlap.tsv", overlap, ["sample"])
        shared_rows = [dict(feature_id=f,control_reads=negative[f],soil_reads=sum(values[s].get(f,0) for s in soil),
                            soil_samples=";".join(s for s in soil if f in values[s])) for f in sorted(shared)]
        equal_records("control_shared_asv.tsv", shared_rows, ["feature_id"])

        # Validate every scalar selector and its rounding against independently
        # reconstructed records. Explicit macro names prevent silent relabeling.
        expected = {}
        def scalar(name, value, filename, locator, rounding="integer"):
            expected["Revision"+name] = dict(value=str(value),source=str(OUT/filename),locator=locator,rounding=rounding)
        for r in kr:
            name = dict(Kr1="KrOne",Kr2="KrTwo",Kr3="KrThree",Kr4="KrFour",Kr5="KrFive")[r["sample"].split("_")[0]]
            for field,suffix in zip(kr_fields,["TotalReads","EukaryotaReads","UnknownReads","OrganelleReads","SelectedReads"]):
                scalar(name+"Silva"+suffix,r[field],"kr_filter_reasons.tsv",r["sample"]+";"+field)
        for (a,b,name,_),r in zip(PAIRS,genus):
            for field,suffix in [("pooled_both_assigned_asv","PooledBothAssignedAsv"),("pooled_both_match_asv","PooledBothMatchAsv"),
                                 ("defined_samples","ConditionalDefinedSamples"),("conditional_sample_asv_min","ConditionalSampleAsvMin"),
                                 ("conditional_sample_asv_max","ConditionalSampleAsvMax"),("conditional_sample_asv_median","ConditionalSampleAsvMedian")]:
                rounding = "exact median; no rounding" if field=="conditional_sample_asv_median" else "integer"
                scalar(name+"Genus"+suffix,r[field],"genus_conditional.tsv",a+";"+b+";"+field,rounding)
            displayed = format((Decimal(format(r["conditional_read_median"],".17g"))*100).quantize(Decimal(".01"),rounding=ROUND_HALF_EVEN),"f")
            scalar(name+"GenusConditionalReadMedianPct",displayed,"genus_conditional.tsv",a+";"+b+";conditional_read_median*100","HALF_EVEN;2 decimals")
        for field,suffix in [("input_reads","InputReads"),("retained_reads","RetainedReads"),("observed_asv","ObservedAsv"),
                             ("shared_soil_asv","SharedSoilAsv"),("control_only_asv","OnlyAsv"),("shared_control_reads","SharedControlReads"),
                             ("total_soil_reads","TotalSoilReads"),("shared_soil_reads","SharedSoilReads")]:
            scalar("Control"+suffix,control[field],"control_summary.tsv","K_S63;"+field)
        pct = format((Decimal(format(control["shared_soil_read_fraction"],".17g"))*100).quantize(Decimal(".01"),rounding=ROUND_HALF_EVEN),"f")
        scalar("ControlSharedSoilReadPct",pct,"control_summary.tsv","K_S63;shared_soil_read_fraction*100","HALF_EVEN;2 decimals")
        threshold_sources = ["checks/verify_denoising.py", "src/check_exact_match.py"]
        thresholds = []
        for source in threshold_sources:
            matches = re.findall(r"0\s*<\s*100\s*\*\s*kept\s*<\s*(\d+)\s*\*\s*n", (root/source).read_text())
            check(len(matches) == 1, "retention threshold selector differs: "+source)
            thresholds.append(int(matches[0]))
        check(thresholds[0] == thresholds[1], "retention thresholds disagree")
        equal_records("method_parameters.tsv",
                      [dict(parameter="retention_threshold_pct",value=thresholds[0],
                            source=";".join(threshold_sources),locator="0 < 100 * kept < threshold * n; soil samples")],
                      ["parameter"])
        scalar("RetentionThresholdPct",thresholds[0],"method_parameters.tsv","retention_threshold_pct;value")
        actual = rows(out/"numbers.tsv")
        check(len(actual)==len(expected) and {r["macro"] for r in actual}==expected.keys(), "revision macro names/count differ")
        for r in actual:
            check({k:v for k,v in r.items() if k!="macro"} == expected.get(r["macro"]), "revision scalar value/provenance differs: "+r["macro"])
        definitions = re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{([^{}]+)\}",(out/"numbers.tex").read_text())
        check(len(definitions)==len(expected) and dict(definitions)=={k:v["value"] for k,v in expected.items()}, "revision TeX macro values differ")
        check(tex_rows(out/"kr_filter_table.tex") == [[r["sample"]]+[str(r[k]) for k in kr_fields] for r in kr], "Kr TeX table differs")
        denominator_rows = []
        for sample in soil:
            row = [sample]
            for a,b,_,_ in PAIRS:
                result = next(r for r in per_sample if (r["model_a"],r["model_b"],r["sample"])==(a,b,sample))
                row.extend([str(result["both_assigned_asv"]),str(result["both_assigned_reads"])])
            denominator_rows.append(row)
        check(tex_rows(out/"genus_denominators_table.tex") == denominator_rows, "genus TeX table differs")
        summary=json.loads((out/"summary.json").read_text())
        check(summary["kr_filter_reasons"]==kr, "summary Kr values differ")
        # Summary float serialization is compared numerically through TSV above;
        # require that the JSON exposes exactly the same saved records.
        for name,filename in [("genus_conditional","genus_conditional.tsv")]:
            check([{k:str(v) for k,v in r.items()} for r in summary[name]] == rows(out/filename), "summary differs from "+filename)
        check({k:str(v) for k,v in summary["control"].items()} == rows(out/"control_summary.tsv")[0], "summary control differs")
        manifest=json.loads((out/"manifest.json").read_text())
        required_sources = {str(MAIN/name) for name in ("assignments.tsv","pairwise_match.tsv","pairwise_match_per_sample.tsv")}
        required_sources.update({str(FAMILY/"filter_counts.tsv"),"results/read_counts_R1.tsv","work/qza/ASV_table.qza",
                                 "work/qza/dada2_stats.qza","data/meta/metadata.tsv","src/prepare_article_revision.py",
                                 *threshold_sources})
        check(set(manifest["sources"]) == required_sources, "revision input manifest is incomplete or contains unexpected paths")
        for name,value in manifest["sources"].items():
            if name in required_sources:
                check((root/name).is_file() and digest(root/name)==value, "revision source hash differs: "+name)
        for name,value in manifest["outputs"].items():
            check(Path(name).name == name and (out/name).is_file() and digest(out/name)==value, "revision output hash differs: "+name)
    except (OSError,KeyError,ValueError,TypeError,IndexError,csv.Error) as exc:
        errors.append(type(exc).__name__+": "+str(exc))
    return dict(status="FAIL" if errors else "PASS",checks=count,errors=errors,
                scope="Independent diagnostic transfer verification; no biological truth or contamination verdict.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--out-dir",type=Path,default=ROOT/OUT)
    parser.add_argument("--report",type=Path)
    args=parser.parse_args()
    result=inspect(args.root,args.out_dir)
    text=json.dumps(result,ensure_ascii=False,indent=2)+"\n"
    if args.report:
        with args.report.open("x") as f:
            f.write(text)
    print(text,end="")
    return 0 if result["status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())

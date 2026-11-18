"""Deterministic presentation of saved results; no scientific pipeline rerun."""

import csv
import hashlib
import json
import re
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CMP = "results/taxonomy_comparison/finish_20261116/main/"
FAMILY = "results/taxonomy_barplots/prokaryote_family_20261118_final/"
META = "data/meta/metadata.tsv"
CHOICE = "memory/denoise_choice.md"
PROVENANCE = "results/exact_match_reproduction/provenance/samples.tsv"
FAMILY_MANIFEST = FAMILY + "prokaryotic_family_manifest.json"
COMPARISON_SCOPE = "comparison_all_post_dada2"
FAMILY_SCOPE = "family_prokaryote_selected"
MODELS = ("sklearn_silva138", "sklearn_gtdb", "sklearn_gg2")
MODEL_NAMES = ("Silva", "Gtdb", "GgTwo")
RANKS = ("тип", "класс", "порядок", "семейство", "род")
RANK_LABELS = ("Phylum", "Class", "Order", "Family", "Genus")
PAIRS = ((0, 1), (0, 2), (1, 2))
COUNT_CATEGORIES = ("both_match", "both_diff", "only_a", "only_b", "neither")

COVERAGE_COLUMNS = (
    "scope", "model", "rank", "assigned_asv", "total_asv", "assigned_asv_pct",
    "read_coverage_min_pct", "read_coverage_median_pct", "read_coverage_max_pct", "n_samples",
)
PAIR_COLUMNS = (
    "scope", "model_a", "model_b", "rank", "distance_min", "distance_median", "distance_max",
    "n_samples", *COUNT_CATEGORIES, "total_asv", "both_assigned_asv", "conditional_asv_pct",
    "conditional_read_min_pct", "conditional_read_median_pct", "conditional_read_max_pct",
    "conditional_defined_samples",
)
FAMILY_COLUMNS = (
    "scope", "model", "sample", "sample_type", "post_dada2_reads", "selected_reads", "excluded_reads",
    "selected_bacteria_reads", "selected_archaea_reads", "named_family_reads", "unclassified_family_reads",
    "selected_pct_of_post_dada2", "unclassified_pct_of_selected", "zero_selected",
)
FIGURE_COLUMNS = (
    "scope", "figure_id", "model", "path", "sha256", "source_table", "denominator", "display_threshold_pct",
)
NUMBER_COLUMNS = (
    "id", "scope", "source_file", "row_key", "column", "operation", "numerator", "denominator",
    "n", "unit", "raw_value", "digits", "display_value",
)


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(relative):
    with (ROOT / relative).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def indexed(rows, keys):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        if key in result:
            raise ValueError(f"Duplicate source key: {key}")
        result[key] = row
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def ratio(numerator, denominator):
    return Fraction(numerator * 100, denominator) if denominator else None


def median(values):
    values = sorted(values)
    if not values:
        return None
    half = len(values) // 2
    return values[half] if len(values) % 2 else (values[half - 1] + values[half]) / 2


def raw(value):
    if value is None:
        return "NA"
    value = Fraction(value)
    with localcontext() as context:
        context.prec = 60
        return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def display(value, digits):
    if value is None:
        return "NA"
    value = Fraction(value)
    with localcontext() as context:
        context.prec = 60
        number = Decimal(value.numerator) / Decimal(value.denominator)
        rounded = number.quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_EVEN)
        return f"{rounded:.{digits}f}"


def write_tsv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_bundle():
    """Return plain output rows and traceable scalar records from fixed inputs."""
    metadata = read_rows(META)
    metadata_index = indexed(metadata, ("sample-id",))
    all_samples = sorted(key[0] for key in metadata_index)
    soil = sorted(row["sample-id"] for row in metadata if row["sample_type"] == "soil")
    require(len(metadata) == 11 and len(soil) == 10, "Expected eleven samples and ten soils")
    require([row["sample-id"] for row in metadata if row["sample_type"] == "negative_control"] == ["K_S63"],
            "Expected only K_S63 as negative control")
    assignments = read_rows(CMP + "assignments.tsv")
    indexed(assignments, ("feature_id",))
    coverage_source = read_rows(CMP + "coverage.tsv")
    coverage = indexed(coverage_source, ("model", "rank"))
    coverage_per_sample = indexed(read_rows(CMP + "coverage_per_sample.tsv"), ("model", "sample", "rank"))
    pairs = indexed(read_rows(CMP + "pairwise_match.tsv"), ("model_a", "model_b", "rank"))
    pairs_per_sample = indexed(read_rows(CMP + "pairwise_match_per_sample.tsv"), ("model_a", "model_b", "sample", "rank"))
    distances = indexed(read_rows(CMP + "distances.tsv"), ("model_a", "model_b", "sample", "rank"))
    filters = indexed(read_rows(FAMILY + "filter_counts.tsv"), ("model", "sample"))
    profiles = read_rows(FAMILY + "family_profiles.tsv")
    indexed(profiles, ("model", "sample", "rank", "key"))
    family_manifest = json.loads((ROOT / FAMILY_MANIFEST).read_text(encoding="utf-8"))
    threshold = Fraction(family_manifest["parameters"]["min_percent"])
    require(set(coverage) == {(model, rank) for model in MODELS for rank in RANKS}, "Unexpected model/rank set")
    require(set(coverage_per_sample) == {(model, sample, rank) for model in MODELS for sample in soil for rank in RANKS},
            "Unexpected coverage sample set")
    expected_pair_keys = {(MODELS[a], MODELS[b], rank) for a, b in PAIRS for rank in RANKS}
    expected_sample_pair_keys = {(MODELS[a], MODELS[b], sample, rank) for a, b in PAIRS for sample in soil for rank in RANKS}
    require(set(pairs) == expected_pair_keys, "Unexpected pairwise ASV rows")
    require(set(pairs_per_sample) == set(distances) == expected_sample_pair_keys, "Unexpected pairwise sample rows")
    require(set(filters) == {(model, sample) for model in MODELS for sample in all_samples}, "Unexpected family sample set")

    coverage_rows, pair_rows, family_rows, numbers = [], [], [], []

    def number(identifier, scope, source_file, selector, column, operation, value,
               digits=0, unit="samples", numerator="", denominator="", n=1):
        numbers.append(dict(zip(NUMBER_COLUMNS, (
            identifier, scope, source_file, compact(selector), compact(column) if isinstance(column, list) else column,
            operation, compact(numerator) if isinstance(numerator, list) else str(numerator),
            compact(denominator) if isinstance(denominator, list) else str(denominator), str(n), unit,
            raw(value), str(digits), display(value, digits),
        ))))

    for identifier, selector, values in (
        ("ProjectSamplesAll", {}, metadata),
        ("ProjectSamplesSoil", {"sample_type": "soil"}, [row for row in metadata if row["sample_type"] == "soil"]),
        ("ProjectSamplesControl", {"sample_type": "negative_control"}, [row for row in metadata if row["sample_type"] == "negative_control"]),
    ):
        number(identifier, "project", META, selector, "sample-id", "count_rows", len(values), n=len(values))
    for identifier, selector, values in (
        ("ProjectAsvAll", {}, assignments),
        ("ProjectAsvSoil", {"soil_asv": "1"}, [row for row in assignments if row["soil_asv"] == "1"]),
    ):
        number(identifier, "project", CMP + "assignments.tsv", selector, "feature_id", "count_rows", len(values), unit="ASV", n=len(values))
    choice_text = (ROOT / CHOICE).read_text(encoding="utf-8")
    parameters = {}
    for field, identifier in (("trim-left", "ProjectTrimLeft"), ("trunc-len", "ProjectTruncLen")):
        matches = re.findall(rf"^{field}:\s*(\d+)\s*$", choice_text, re.MULTILINE)
        require(len(matches) == 1, f"Expected exactly one {field}")
        parameters[field] = int(matches[0])
        number(identifier, "project", CHOICE, {}, field, "value", parameters[field], unit="nt")
    number("ProjectAsvLength", "project", CHOICE, {}, ["trunc-len", "trim-left"], "subtract",
           parameters["trunc-len"] - parameters["trim-left"], unit="nt",
           numerator=parameters["trunc-len"], denominator=parameters["trim-left"], n=2)
    kr5_key = {"run": "1", "sample-id": "Kr5_S84"}
    kr5 = indexed(read_rows(PROVENANCE), ("run", "sample-id"))[("1", "Kr5_S84")]
    for identifier, column in (("ProjectKrFiveInput", "input"), ("ProjectKrFiveRetained", "non-chimeric")):
        number(identifier, "project", PROVENANCE, kr5_key, column, "value", int(kr5[column]), unit="reads")
    number("ProjectKrFiveRetentionPct", "project", PROVENANCE, kr5_key, ["non-chimeric", "input"],
           "ratio_pct", ratio(int(kr5["non-chimeric"]), int(kr5["input"])), digits=2, unit="percent",
           numerator=int(kr5["non-chimeric"]), denominator=int(kr5["input"]))
    number("ProjectModels", "project", CMP + "coverage.tsv", {}, "model", "count_distinct",
           len({row["model"] for row in coverage_source}), unit="models", n=len(coverage_source))
    number("SelectedDisplayThresholdPct", FAMILY_SCOPE, FAMILY_MANIFEST, {}, "parameters.min_percent", "value", threshold,
           digits=2, unit="percent")

    for model, model_name in zip(MODELS, MODEL_NAMES):
        for rank in RANKS:
            source = coverage[(model, rank)]
            numerator, denominator = int(source["assigned_asv"]), int(source["total_asv"])
            sample_rows = [coverage_per_sample[(model, sample, rank)] for sample in soil]
            numerators = [int(row["assigned_reads"]) for row in sample_rows]
            denominators = [int(row["total_reads"]) for row in sample_rows]
            require(all(denominators), "Zero post-DADA2 soil denominator")
            values = [ratio(num, den) for num, den in zip(numerators, denominators)]
            coverage_rows.append(dict(zip(COVERAGE_COLUMNS, (
                COMPARISON_SCOPE, model, rank, numerator, denominator, raw(ratio(numerator, denominator)),
                raw(min(values)), raw(median(values)), raw(max(values)), len(values),
            ))))
            if rank == "семейство":
                number(f"Cmp{model_name}FamilyAsvPct", COMPARISON_SCOPE, CMP + "coverage.tsv", {"model": model, "rank": rank},
                       ["assigned_asv", "total_asv"], "ratio_pct", ratio(numerator, denominator), digits=2,
                       unit="percent", numerator=numerator, denominator=denominator)
                number(f"Cmp{model_name}FamilyReadMedianPct", COMPARISON_SCOPE, CMP + "coverage_per_sample.tsv",
                       {"model": model, "rank": rank, "sample": soil}, ["assigned_reads", "total_reads"], "median_ratio_pct",
                       median(values), digits=2, unit="percent", numerator=numerators, denominator=denominators, n=len(values))
        zero_samples = [sample for sample in soil if filters[(model, sample)]["zero_selected"] == "1"]
        number(f"Selected{model_name}ZeroSoilSamples", FAMILY_SCOPE, FAMILY + "filter_counts.tsv",
               {"model": model, "sample_type": "soil", "zero_selected": "1"}, "sample", "count_rows", len(zero_samples), n=len(zero_samples))

    for a, b in PAIRS:
        model_a, model_b = MODELS[a], MODELS[b]
        pair_name = MODEL_NAMES[a] + MODEL_NAMES[b]
        for rank in RANKS:
            selector = {"model_a": model_a, "model_b": model_b, "rank": rank}
            source = pairs[(model_a, model_b, rank)]
            values = [Fraction(distances[(model_a, model_b, sample, rank)]["distance"]) for sample in soil]
            sample_rows = [(sample, pairs_per_sample[(model_a, model_b, sample, rank)]) for sample in soil]
            defined = [(sample, row) for sample, row in sample_rows if int(row["both_assigned_reads"]) != 0]
            numerators = [int(row["both_match_reads"]) for _, row in defined]
            denominators = [int(row["both_assigned_reads"]) for _, row in defined]
            conditional = [ratio(num, den) for num, den in zip(numerators, denominators)]
            numerator, denominator = int(source["both_match"]), int(source["both_assigned_asv"])
            pair_rows.append(dict(zip(PAIR_COLUMNS, (
                COMPARISON_SCOPE, model_a, model_b, rank, raw(min(values)), raw(median(values)), raw(max(values)), len(values),
                *(int(source[key]) for key in COUNT_CATEGORIES), int(source["total_asv"]), denominator,
                raw(ratio(numerator, denominator)), raw(min(conditional) if conditional else None),
                raw(median(conditional)), raw(max(conditional) if conditional else None), len(conditional),
            ))))
            if rank in ("семейство", "род"):
                rank_name = "Family" if rank == "семейство" else "Genus"
                number(f"Cmp{pair_name}{rank_name}DistanceMedian", COMPARISON_SCOPE, CMP + "distances.tsv",
                       {**selector, "sample": soil}, "distance", "median_decimal", median(values), digits=6, unit="distance", n=len(values))
            if rank == "род":
                number(f"Cmp{pair_name}GenusConditionalAsvPct", COMPARISON_SCOPE, CMP + "pairwise_match.tsv",
                       selector, ["both_match", "both_assigned_asv"], "ratio_pct", ratio(numerator, denominator), digits=2,
                       unit="percent", numerator=numerator, denominator=denominator)
                defined_selector = {**selector, "sample": [sample for sample, _ in defined]}
                number(f"Cmp{pair_name}GenusConditionalReadMedianPct", COMPARISON_SCOPE, CMP + "pairwise_match_per_sample.tsv",
                       defined_selector, ["both_match_reads", "both_assigned_reads"], "median_ratio_pct", median(conditional), digits=2,
                       unit="percent", numerator=numerators, denominator=denominators, n=len(conditional))
                number(f"Cmp{pair_name}GenusDefinedSamples", COMPARISON_SCOPE, CMP + "pairwise_match_per_sample.tsv",
                       defined_selector, "sample", "count_rows", len(conditional), n=len(conditional))

    profile_groups = {}
    for row in profiles:
        require(row["rank"] == "семейство", "Non-family row in family profile")
        profile_groups.setdefault((row["model"], row["sample"]), []).append(row)
    for model in MODELS:
        for sample in all_samples:
            source = filters[(model, sample)]
            sample_type = metadata_index[(sample,)]["sample_type"]
            require(source["sample_type"] == sample_type, "Family sample role mismatch")
            selected, total = int(source["selected_reads"]), int(source["total_reads"])
            profile = profile_groups.get((model, sample), [])
            require(sum(int(row["count"]) for row in profile) == selected, "Family profile sum does not equal selected reads")
            require(all(int(row["total_reads"]) == selected and row["sample_type"] == sample_type for row in profile),
                    "Family profile denominator or role mismatch")
            unknown = sum(int(row["count"]) for row in profile if row["key"] == "UNCLASSIFIED")
            require(int(source["zero_selected"]) == int(selected == 0), "Family zero_selected mismatch")
            family_rows.append(dict(zip(FAMILY_COLUMNS, (
                FAMILY_SCOPE, model, sample, sample_type, total, selected, int(source["excluded_reads"]),
                int(source["selected_Bacteria_reads"]), int(source["selected_Archaea_reads"]), selected - unknown, unknown,
                raw(ratio(selected, total)), raw(ratio(unknown, selected)), int(selected == 0),
            ))))
    require(len(numbers) == 37 and len({row["id"] for row in numbers}) == 37, "Expected exactly 37 unique macros")
    require(all(re.fullmatch(r"[A-Za-z]+", row["id"]) for row in numbers), "Macro names must contain ASCII letters only")

    figure_rows = []
    for suffix, model in (("overview", "all"), ("silva", MODELS[0]), ("gtdb", MODELS[1]), ("gg2", MODELS[2])):
        identifier = "prokaryotic_family_" + suffix
        path = "figures/taxonomy/family/" + identifier + ".pdf"
        figure_rows.append(dict(zip(FIGURE_COLUMNS, (
            FAMILY_SCOPE, identifier, model, path, sha256(ROOT / path), FAMILY + "family_barplot_data.tsv",
            "selected_reads(model,sample)", raw(threshold),
        ))))
    return coverage_rows, pair_rows, family_rows, figure_rows, numbers, soil


def latex_tables(coverage_rows, pair_rows, family_rows, soil):
    coverage = indexed(coverage_rows, ("model", "rank"))
    pairs = indexed(pair_rows, ("model_a", "model_b", "rank"))
    families = indexed(family_rows, ("model", "sample"))
    comparison = [
        "% scope: comparison_all_post_dada2",
        "% Read coverage median (percent); D median; all post-DADA2 reads.",
        r"\begin{tabular}{lrrrrrr}", r"\hline",
        r"Rank & \multicolumn{3}{c}{Read coverage (\%)} & \multicolumn{3}{c}{D} \\",
        r" & SILVA & GTDB & GG2 & SILVA--GTDB & SILVA--GG2 & GTDB--GG2 \\",
        r"\hline",
    ]
    for rank, label in zip(RANKS, RANK_LABELS):
        values = [display(Fraction(coverage[(model, rank)]["read_coverage_median_pct"]), 2) for model in MODELS]
        values += [display(Fraction(pairs[(MODELS[a], MODELS[b], rank)]["distance_median"]), 6) for a, b in PAIRS]
        comparison.append(" & ".join([label, *values]) + r" \\")
    comparison += [r"\hline", r"\end{tabular}"]
    family = [
        "% scope: family_prokaryote_selected",
        "% Reads: post-DADA2; N: selected reads; U: UNCLASSIFIED / N (percent); N=0 gives NA.",
        r"\begin{tabular}{lrrrrrrr}", r"\hline",
        r"Sample & Reads & \multicolumn{2}{c}{SILVA} & \multicolumn{2}{c}{GTDB} & \multicolumn{2}{c}{GG2} \\",
        r" & & N & U (\%) & N & U (\%) & N & U (\%) \\",
        r"\hline",
    ]
    for sample in soil:
        totals = {families[(model, sample)]["post_dada2_reads"] for model in MODELS}
        require(len(totals) == 1, "Post-DADA2 reads differ across models")
        values = [sample.replace("_", r"\_"), str(next(iter(totals)))]
        for model in MODELS:
            row = families[(model, sample)]
            value = row["unclassified_pct_of_selected"]
            values += [str(row["selected_reads"]), "NA" if value == "NA" else display(Fraction(value), 2)]
        family.append(" & ".join(values) + r" \\")
    family += [r"\hline", r"\end{tabular}"]
    return "\n".join(comparison) + "\n", "\n".join(family) + "\n"


README = """# Числовой комплект шага 4

Дата договора: 2026-11-18. Контракт: research/paper_results_contract.md.
Рукопись: absent_not_checked. Здесь перенесены сохранённые результаты; DADA2,
классификация, научные сравнения и семейственные рисунки не выполнялись заново.

comparison_coverage.tsv: 15 строк, ключ model/rank. Целые assigned_asv/total_asv
из comparison main/coverage.tsv; read_coverage — проценты assigned_reads/total_reads
из coverage_per_sample.tsv отдельно по десяти почвам, затем min/median/max.
comparison_pairs.tsv: 15 строк, ключ model_a/model_b/rank. Пять категорий и общий
ASV-знаменатель из pairwise_match.tsv. Условный процент — both_match/both_assigned.
По ридам источник pairwise_match_per_sample.tsv; нулевой both_assigned_reads
исключается только из условной сводки; число определённых почв показано отдельно.
Расстояния — сохранённый distance из distances.tsv, сводка по десяти почвам.
Область comparison_all_post_dada2 использует все риды после DADA2 и прежнюю
буквальную обработку таксонов; uncultured не переопределяется.

family_sample_counts.tsv: 33 строки, ключ model/sample; sample_type берётся из
metadata.tsv и сверяется с исходниками. Источник filter_counts.tsv, а UNCLASSIFIED
выбирается по точному key из family_profiles.tsv. post_dada2_reads = total_reads
filter_counts; selected_reads = N. named_family_reads = N − UNCLASSIFIED.
selected_pct_of_post_dada2 = 100*N/post_dada2_reads; unclassified_pct_of_selected
= 100*UNCLASSIFIED/N. При N=0 последняя доля NA, не ноль. Контроль K_S63 остаётся
в подробной таблице, но не входит в почвенные сводки или основной TeX-фрагмент.
Область family_prokaryote_selected — модельный отбор Bacteria + Archaea с
исключением органелл/Eukaryota/неизвестного домена; точное f__uncultured здесь
считается неопределённым семейством. Это назначения моделей, не биологическая истина.

figure_sources.tsv регистрирует четыре существующих PDF figures/taxonomy/family;
файлы не копируются. Source_table — family_barplot_data.tsv. Знаменатель каждой
панели — selected_reads(model,sample), порог отображения — parameters.min_percent
prokaryotic_family_manifest.json. В manifest.json также зарегистрированы маска,
полные профили, сохранённая проверка и манифест переноса рисунков.

numbers.tsv содержит 37 скаляров и их полное происхождение. row_key — JSON
селектор: строка означает точное равенство, массив sample — заданный набор строк
в указанном порядке. column — поле либо JSON-массив полей числителя/знаменателя.
operation: count_rows, count_distinct, value, subtract, ratio_pct,
median_ratio_pct или median_decimal. numerator/denominator у ratio_pct — целые
значения; у median_ratio_pct — массивы целых в порядке sample; у subtract — два
операнда. n — число исходных выбранных строк (для subtract — два параметра).
raw_value: точная рациональная арифметика с десятичной сериализацией precision=60;
для печати Decimal ROUND_HALF_EVEN: проценты 2 знака, D 6, целые 0. NA сохраняется.
Медиана чётного числа значений — среднее двух центральных.

numbers.tex определяет 37 ASCII-команд newcommand. comparison_main.tex:
пять рангов, медианы покрытия и D; family_main.tex: десять почв, объём до отбора,
N и процент UNCLASSIFIED для каждой модели. Это простые tabular-фрагменты,
не рукопись и не проверка окончательной вёрстки. Исходные русские ранги TSV
соответствуют Phylum/Class/Order/Family/Genus в TeX.

manifest.json фиксирует SHA256 всех использованных источников, двух новых
файлов генератора и девяти выходов. Сам manifest не включает собственный хеш.
Пути относительны корню проекта или каталогу комплекта; время и абсолютные
рабочие пути не сохраняются, поэтому повторный запуск побайтно воспроизводим.

Запуск из корня проекта, выходной каталог должен отсутствовать:
conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/make_paper_results.py --out-dir results/paper_step4/NEW

Генерация сама не является независимой проверкой. Приёмка выполняется отдельно
src/check_paper_results.py; отсутствие рукописи остаётся absent_not_checked.
"""


def generate(out_dir):
    """Create a new bundle, refusing an existing destination before any write."""
    out_dir = Path(out_dir)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    if out_dir.exists():
        raise FileExistsError(f"Output directory already exists: {out_dir}")
    coverage, pairs, family, figures, numbers, soil = build_bundle()
    comparison_tex, family_tex = latex_tables(coverage, pairs, family, soil)
    sources = [META, CHOICE, PROVENANCE, "research/paper_results_contract.md"]
    sources += [CMP + name + ".tsv" for name in (
        "assignments", "coverage", "coverage_per_sample", "pairwise_match", "pairwise_match_per_sample", "distances",
    )]
    sources += [FAMILY + name + ".tsv" for name in ("filter_counts", "family_profiles", "asv_mask", "family_barplot_data")]
    sources += [FAMILY_MANIFEST, FAMILY + "verification.json", "figures/taxonomy/family/manifest.json"]
    sources += [row["path"] for row in figures]
    code = ["src/make_paper_results.py", "src/paper_results_core.py"]
    source_hashes = {path: sha256(ROOT / path) for path in sorted(sources)}
    code_hashes = {path: sha256(ROOT / path) for path in sorted(code)}
    out_dir.mkdir(parents=True, exist_ok=False)
    for name, columns, rows in (
        ("comparison_coverage.tsv", COVERAGE_COLUMNS, coverage),
        ("comparison_pairs.tsv", PAIR_COLUMNS, pairs),
        ("family_sample_counts.tsv", FAMILY_COLUMNS, family),
        ("figure_sources.tsv", FIGURE_COLUMNS, figures),
        ("numbers.tsv", NUMBER_COLUMNS, numbers),
    ):
        write_tsv(out_dir / name, columns, rows)
    macros = "% Generated from numbers.tsv; manuscript_status: absent_not_checked\n"
    macros += "".join("\\newcommand{\\" + row["id"] + "}{" + row["display_value"] + "}\n" for row in numbers)
    for name, text in (("numbers.tex", macros), ("comparison_main.tex", comparison_tex),
                       ("family_main.tex", family_tex), ("README.md", README)):
        (out_dir / name).write_text(text, encoding="utf-8", newline="\n")
    manifest = {
        "schema_version": 1,
        "manuscript_status": "absent_not_checked",
        "sources": source_hashes,
        "code": code_hashes,
        "outputs": {path.name: sha256(path) for path in sorted(out_dir.iterdir())},
    }
    require(source_hashes == {path: sha256(ROOT / path) for path in sorted(sources)}, "Sources changed during generation")
    require(code_hashes == {path: sha256(ROOT / path) for path in sorted(code)}, "Generator changed during generation")
    (out_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"status": "generated", "rows": {"comparison_coverage": len(coverage), "comparison_pairs": len(pairs),
            "family_sample_counts": len(family), "figure_sources": len(figures), "numbers": len(numbers)},
            "manuscript_status": "absent_not_checked"}

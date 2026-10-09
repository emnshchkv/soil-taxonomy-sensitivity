#!/usr/bin/env python3
"""Критерий готовности этапа 1: чтение ридов, импорт, денойзинг, классификация (без сравнения методов).

Запуск из корня проекта, в окружении qiime2-amplicon-2025.7 (нужны biom и h5py):
    python3 checks/check.py
Код возврата: 0 — все проверки пройдены, 1 — есть провалы (печатаются списком).
Скрипт ничего не изменяет. Папка checks/ закрыта от записи: chmod -R a-w checks
(поправить проверку самому: chmod -R u+w checks).

Сверки с независимыми источниками (агент не может их подделать, не подделав сами данные):
  - число ридов: read_counts_R1.tsv  <->  demux summarize  <->  вход DADA2 (stats);
  - параметры денойзинга: memory/denoise_choice.md  <->  provenance таблицы ASV  <->  длина ASV;
  - прочитанные ДНК: сумма по образцу в таблице ASV  ==  non-chimeric из stats;
  - какая база использована: uuid классификатора в provenance таксономии  <->  uuid файла в data/refs/.
"""
import csv, hashlib, io, re, sys, zipfile
from pathlib import Path

ROOT = Path(".")
EXPECTED_SAMPLES, NEG = 11, "K_S63"
METHODS = {  # метод -> файл классификатора в data/refs/
    "sklearn_silva138": "silva-138-99-nb-classifier.qza",
    "sklearn_gtdb": "gtdb_classifier_r220.qza",
    "sklearn_gg2": "gg_2024.09.full-length.nb.sklearn-1.4.2.qza",
}
fails, info = [], []
UUID = r"([0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12})"

def fail(msg): fails.append(msg)

def tsv(text, skip_hash=True):
    lines = [l for l in text.splitlines() if l.strip() and not (skip_hash and l.startswith("#"))]
    return list(csv.DictReader(lines, delimiter="\t"))

def read_table(path, cols):
    p = ROOT / path
    if not p.is_file() or p.stat().st_size == 0:
        fail(f"{path}: файла нет или он пустой"); return []
    rows = tsv(p.read_text())
    miss = [c for c in cols if not rows or c not in rows[0]]
    if miss:
        fail(f"{path}: нет столбцов {miss}"); return []
    return rows

def to_int(path, row, col):
    try: return int(row[col])
    except (ValueError, TypeError, KeyError):
        fail(f"{path}: не целое число в столбце {col}: {row.get(col)!r}"); return None

class Qza:
    def __init__(self, path, expect_type):
        self.path, self.ok = path, False
        p = ROOT / path
        if not p.is_file():
            fail(f"{path}: файла нет"); return
        try: self.z = zipfile.ZipFile(p)
        except zipfile.BadZipFile:
            fail(f"{path}: не архив QIIME 2"); return
        self.root = self.z.namelist()[0].split("/")[0]
        self.uuid = self.root
        meta = self.text("metadata.yaml")
        m = re.search(r"^type:\s*(.+?)\s*$", meta or "", re.M)
        t = m.group(1).strip("'\"") if m else None
        if t != expect_type:
            fail(f"{path}: тип {t!r}, ожидался {expect_type!r}"); return
        self.ok = True
    def raw(self, name):
        try: return self.z.read(f"{self.root}/{name}")
        except KeyError:
            fail(f"{self.path}: внутри нет {name}"); return None
    def text(self, name):
        b = self.raw(name); return b.decode() if b is not None else None
    def provenance(self): return self.text("provenance/action/action.yaml") or ""
    def find(self, suffix):
        return next((n for n in self.z.namelist() if n.endswith(suffix)), None)

def refs_uuid(name):
    p = ROOT / "data/refs" / name
    if not p.is_file(): fail(f"data/refs/{name}: файла нет"); return None
    return zipfile.ZipFile(p).namelist()[0].split("/")[0]

# ---------- 1. метаданные, манифест, неизменность и защита сырых данных ----------
meta = read_table("data/meta/metadata.tsv", ["sample-id", "sample_type"])
ids = {r["sample-id"] for r in meta}
types = [r["sample_type"] for r in meta]
if meta and (len(ids) != EXPECTED_SAMPLES or types.count("negative_control") != 1
             or types.count("soil") != EXPECTED_SAMPLES - 1 or NEG not in ids):
    fail(f"metadata.tsv: ожидается {EXPECTED_SAMPLES} образцов: {NEG} (negative_control) и остальные soil")
soil = {r["sample-id"] for r in meta if r["sample_type"] == "soil"}

man = read_table("data/meta/fq-manifest-single.tsv", ["sample-id", "absolute-filepath"])
if man:
    if {r["sample-id"] for r in man} != ids: fail("манифест и metadata.tsv: разные наборы образцов")
    for r in man:
        f = r["absolute-filepath"]
        if "_R1_" not in f or "_R2_" in f: fail(f"манифест: не R1-файл: {f}")
        if not Path(f).is_file(): fail(f"манифест: файла нет: {f}")

sums = ROOT / "data/raw.sha256"
if not sums.is_file():
    fail("data/raw.sha256 нет (shasum -a 256 data/raw/*.fastq.gz > data/raw.sha256)")
else:
    for line in sums.read_text().splitlines():
        if not line.strip(): continue
        h, name = line.split(None, 1)
        p = ROOT / name.strip().lstrip("*")
        if not p.is_file(): fail(f"raw.sha256: файла нет: {p}"); continue
        d = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""): d.update(chunk)
        if d.hexdigest() != h: fail(f"контрольная сумма не совпала (файл изменён): {p}")

for p in [ROOT / "checks", ROOT / "checks/check.py", ROOT / "data/raw", *sorted((ROOT / "data/raw").glob("*.fastq.gz"))]:
    if p.exists() and p.stat().st_mode & 0o222:
        fail(f"{p} не закрыт от записи (chmod -R a-w {p.parts[0]}/{p.parts[1] if p.parts[0]=='data' and len(p.parts)>1 else ''})".replace("/)", ")"))

# ---------- 2. число ридов ----------
rc = read_table("results/read_counts_R1.tsv",
                ["sample-id", "reads_R1", "reads_R2", "mean_len_R1", "min_len_R1", "max_len_R1"])
reads = {}
if rc:
    if {r["sample-id"] for r in rc} != ids: fail("read_counts_R1.tsv: набор образцов не совпадает с metadata.tsv")
    for r in rc:
        a, b = to_int("read_counts_R1.tsv", r, "reads_R1"), to_int("read_counts_R1.tsv", r, "reads_R2")
        if a is None or b is None: continue
        reads[r["sample-id"]] = a
        if a <= 0: fail(f"read_counts_R1.tsv: 0 ридов у {r['sample-id']}")
        if a != b: fail(f"read_counts_R1.tsv: R1 != R2 у {r['sample-id']} ({a} vs {b})")
        try:
            lo, mean, hi = float(r["min_len_R1"]), float(r["mean_len_R1"]), float(r["max_len_R1"])
            if not lo <= mean <= hi: fail(f"read_counts_R1.tsv: нарушено min<=mean<=max у {r['sample-id']}")
        except ValueError: fail(f"read_counts_R1.tsv: длины не числа у {r['sample-id']}")

# ---------- 3. импорт и demux summarize ----------
imp = Qza("work/qza/raw_reads_single.qza", "SampleData[SequencesWithQuality]")
if imp.ok:
    mf = imp.text("data/MANIFEST")
    if mf is not None:
        rows = [l.split(",") for l in mf.splitlines() if l.strip() and not l.startswith("#")][1:]
        if {r[0] for r in rows} != ids: fail("raw_reads_single.qza: образцы не совпадают с metadata.tsv")
        if any(r[2].strip() != "forward" for r in rows): fail("raw_reads_single.qza: есть не-forward риды (нужен только R1)")
sm = Qza("work/qzv/reads_summary.qzv", "Visualization")
if sm.ok:
    n = sm.find("per-sample-fastq-counts.tsv")
    if not n: fail("reads_summary.qzv: нет per-sample-fastq-counts.tsv")
    else:
        got = {}
        for r in tsv(sm.z.read(n).decode()):
            vals = list(r.values())
            try: got[vals[0]] = int(float(vals[1]))
            except (ValueError, IndexError): fail(f"reads_summary.qzv: не разобрана строка {r}")
        for s, c in reads.items():
            if got.get(s) != c: fail(f"reads_summary.qzv: ридов у {s} = {got.get(s)}, в read_counts_R1.tsv = {c}")

# ---------- 4. денойзинг ----------
dc = (ROOT / "memory/denoise_choice.md")
trim = trunc = None
if not dc.is_file(): fail("memory/denoise_choice.md нет (решение человека: trim-left и trunc-len)")
else:
    t = dc.read_text()
    m1, m2 = re.search(r"^trim-left:\s*(\d+)", t, re.M), re.search(r"^trunc-len:\s*(\d+)", t, re.M)
    if not (m1 and m2): fail("denoise_choice.md: нужны строки «trim-left: <число>» и «trunc-len: <число>»")
    else: trim, trunc = int(m1.group(1)), int(m2.group(1))

tab = Qza("work/qza/ASV_table.qza", "FeatureTable[Frequency]")
rep = Qza("work/qza/rep_seq.qza", "FeatureData[Sequence]")
st = Qza("work/qza/dada2_stats.qza", "SampleData[DADA2Stats]")
Qza("work/qzv/dada2_stats.qzv", "Visualization")

if tab.ok and trim is not None:
    pv = tab.provenance()
    if "denoise_single" not in pv: fail("ASV_table.qza: получена не через dada2 denoise-single")
    for name, want in (("trim_left", trim), ("trunc_len", trunc)):
        m = re.search(rf"^\s*-?\s*{name}:\s*(\d+)", pv, re.M)
        if not m: fail(f"ASV_table.qza: в provenance нет параметра {name}")
        elif int(m.group(1)) != want: fail(f"ASV_table.qza: {name} = {m.group(1)}, в denoise_choice.md = {want}")

asv_ids = None
if rep.ok:
    fa = rep.text("data/dna-sequences.fasta") or ""
    seqs, cur = {}, None
    for l in fa.splitlines():
        if l.startswith(">"): cur = l[1:].strip(); seqs[cur] = ""
        elif cur: seqs[cur] += l.strip()
    asv_ids = set(seqs)
    if not seqs: fail("rep_seq.qza: нет последовательностей")
    if trunc:
        bad = {len(s) for s in seqs.values()} - {trunc - trim}
        if bad: fail(f"rep_seq.qza: длины ASV {sorted(bad)}, ожидалась {trunc - trim} (trunc-len − trim-left)")

nonchim = {}
if st.ok:
    for r in tsv(st.text("data/stats.tsv") or ""):
        sid = r.get("sample-id") or r.get("#SampleID") or list(r.values())[0]
        try: v = {k: int(r[k]) for k in ("input", "filtered", "denoised", "non-chimeric")}
        except (KeyError, ValueError): fail(f"dada2_stats.qza: не разобрана строка {r}"); continue
        nonchim[sid] = v["non-chimeric"]
        if not v["input"] >= v["filtered"] >= v["denoised"] >= v["non-chimeric"]:
            fail(f"dada2_stats.qza: ридов не убывает по этапам у {sid}: {v}")
        if sid in reads and v["input"] != reads[sid]:
            fail(f"dada2_stats.qza: input = {v['input']} у {sid}, в read_counts_R1.tsv = {reads[sid]}")
        if sid in soil and v["non-chimeric"] == 0: fail(f"dada2_stats.qza: у почвенного образца {sid} 0 ридов после денойзинга")
        if v["input"]: info.append(f"{sid}: сохранено {100*v['non-chimeric']/v['input']:.1f}% ридов")
    if set(nonchim) != ids: fail("dada2_stats.qza: образцы не совпадают с metadata.tsv")

if tab.ok:
    try:
        import biom, h5py
        raw = tab.raw("data/feature-table.biom")
        if raw is not None:
            t = biom.Table.from_hdf5(h5py.File(io.BytesIO(raw), "r"))
            tsum = dict(zip(t.ids("sample"), t.sum("sample")))
            for sid, nc in nonchim.items():
                if int(round(tsum.get(sid, 0))) != nc:
                    fail(f"ASV_table.qza: сумма по {sid} = {int(round(tsum.get(sid, 0)))}, non-chimeric в stats = {nc}")
            if asv_ids is not None and set(t.ids("observation")) != asv_ids:
                fail("ASV_table.qza и rep_seq.qza: разные наборы ASV")
    except ImportError:
        fail("нет модулей biom/h5py: запустите проверку в окружении qiime2-amplicon-2025.7")

# ---------- 5. классификация ----------
for method, ref in METHODS.items():
    tx = Qza(f"work/qza/taxonomy_{method}.qza", "FeatureData[Taxonomy]")
    if not tx.ok: continue
    rows = tsv(tx.text("data/taxonomy.tsv") or "")
    got = [r.get("Feature ID") for r in rows]
    if len(got) != len(set(got)): fail(f"taxonomy_{method}.qza: повторяющиеся Feature ID")
    if asv_ids is not None and set(got) != asv_ids: fail(f"taxonomy_{method}.qza: ASV не совпадают с rep_seq.qza")
    for r in rows:
        try:
            if not 0 <= float(r["Confidence"]) <= 1: fail(f"taxonomy_{method}.qza: Confidence вне [0,1]: {r}"); break
        except (KeyError, ValueError): fail(f"taxonomy_{method}.qza: нет числового Confidence: {r}"); break
        if not r.get("Taxon", "").strip(): fail(f"taxonomy_{method}.qza: пустой Taxon: {r}"); break
    pv = tx.provenance()
    if "classify_sklearn" not in pv: fail(f"taxonomy_{method}.qza: получена не через classify-sklearn")
    m = re.search(rf"-\s*classifier:\s*{UUID}", pv)
    want = refs_uuid(ref)
    if not m: fail(f"taxonomy_{method}.qza: в provenance нет classifier")
    elif want and m.group(1) != want: fail(f"taxonomy_{method}.qza: использован не {ref} (uuid не совпадает)")
    m = re.search(rf"-\s*reads:\s*{UUID}", pv)
    if rep.ok and m and m.group(1) != rep.uuid: fail(f"taxonomy_{method}.qza: классифицированы не ASV из rep_seq.qza")
    m = re.search(r"^\s*-?\s*confidence:\s*([0-9.]+)", pv, re.M)
    if m and abs(float(m.group(1)) - 0.7) > 1e-9: fail(f"taxonomy_{method}.qza: confidence = {m.group(1)}, нужен порог по умолчанию 0.7")

# ---------- 6. базы и журнал ----------
vf = ROOT / "data/refs/VERSIONS.md"
if not vf.is_file(): fail("data/refs/VERSIONS.md нет")
else:
    vt = vf.read_text()
    for ref in METHODS.values():
        p = ROOT / "data/refs" / ref
        if not p.is_file(): continue
        d = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""): d.update(chunk)
        if ref not in vt or d.hexdigest() not in vt:
            fail(f"VERSIONS.md: нет имени файла или верной sha256 для {ref}")
dr = ROOT / "diary.md"
if not dr.is_file() or "Итог проверки:" not in dr.read_text():
    fail("diary.md: нет ни одной записи с полем «Итог проверки:»")

if fails:
    print("ПРОВЕРКА НЕ ПРОЙДЕНА:"); [print(" -", f) for f in fails]; sys.exit(1)
print("ПРОВЕРКА ПРОЙДЕНА"); [print("  ", i) for i in info]; sys.exit(0)

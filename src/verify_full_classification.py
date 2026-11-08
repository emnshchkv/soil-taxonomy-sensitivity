"""Полная проверка трёх готовых таксономий (AGENT_CLASSIFICATION_FINISH_PROMPT).

Для каждой модели (silva138, gtdb, gg2):
1. qiime tools validate (level=max) — через CLI.
2. Тип FeatureData[Taxonomy], ровно 3374 уникальных ID, точное совпадение
   набора с rep_seq.qza (нет пропущенных или лишних ASV).
3. Confidence — конечные числа в [0,1]; Taxon непустой.
4. Корневое действие classify_sklearn, UUID входа (reads) равен UUID rep_seq.qza,
   UUID классификатора соответствует модели, confidence=0.7, n_jobs=1.
5. Выборочная сверка трёх детерминированно выбранных ASV (первые три по
   лексикографической сортировке ID) с qiime tools export.

Запуск: conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/verify_full_classification.py
Код 0 — все проверки пройдены; 1 — провал.
"""
import math
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import qiime_runtime  # noqa: F401  (выбирает интерпретатор окружения)

ROOT = Path(__file__).resolve().parent.parent

MODEL_UUIDS = {
    "silva138": "70b4b5f4-8fce-40bd-b508-afacbc12a5ed",
    "gtdb": "5d5461cc-6a51-434b-90ab-040f388e4221",
    "gg2": "49ccfb0a-155d-404b-80b4-818d2aeb53b2",
}
MODEL_PATHS = {
    "silva138": "data/refs/silva-138-99-nb-classifier.qza",
    "gtdb": "data/refs/gtdb_classifier_r220.qza",
    "gg2": "data/refs/gg_2024.09.full-length.nb.sklearn-1.4.2.qza",
}
TAX_PATHS = {
    "silva138": "work/qza/taxonomy_sklearn_silva138.qza",
    "gtdb": "work/qza/taxonomy_sklearn_gtdb.qza",
    "gg2": "work/qza/taxonomy_sklearn_gg2.qza",
}
REP_SEQ_UUID = "978af6f8-91bb-4df2-9e19-b395076e5f3c"


def _read_qza_data(path: Path, suffix: str) -> str:
    with zipfile.ZipFile(path) as z:
        data = [n for n in z.namelist() if n.endswith(suffix)]
        if not data:
            raise FileNotFoundError(f"в {path} нет файла *{suffix}")
        return z.read(data[0]).decode()


def _read_root_action(path: Path) -> str:
    with zipfile.ZipFile(path) as z:
        action = [n for n in z.namelist()
                  if n.endswith("provenance/action/action.yaml")][0]
        return z.read(action).decode()


def _parse_action(action_yaml: str) -> dict:
    m = re.search(r"action: (\S+)", action_yaml)
    action = m.group(1) if m else None
    reads = re.search(r"reads: ([0-9a-f-]+)", action_yaml)
    clf = re.search(r"classifier: ([0-9a-f-]+)", action_yaml)
    conf = re.search(r"confidence: ([0-9.]+)", action_yaml)
    njobs = re.search(r"n_jobs: ([0-9]+)", action_yaml)
    return {
        "action": action,
        "reads": reads.group(1) if reads else None,
        "classifier": clf.group(1) if clf else None,
        "confidence": conf.group(1) if conf else None,
        "n_jobs": njobs.group(1) if njobs else None,
    }


def _read_rep_seq_ids() -> set:
    fasta = _read_qza_data(ROOT / "work/qza/rep_seq.qza", "dna-sequences.fasta")
    ids = [l[1:].strip() for l in fasta.splitlines() if l.startswith(">")]
    return set(ids)


def _read_taxonomy(path: Path) -> dict:
    """Читает taxonomy.tsv (из QZA или экспорта): {id: (taxon, confidence)}."""
    if path.suffix == ".qza":
        txt = _read_qza_data(path, "taxonomy.tsv")
    else:
        txt = path.read_text()
    result = {}
    for line in txt.splitlines():
        if line.startswith("Feature ID") or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        result[parts[0]] = (parts[1], parts[2])
    return result


def check_full_taxonomy(method: str) -> bool:
    """Проверяет одну полную таксономию."""
    ok = True
    tax_path = ROOT / TAX_PATHS[method]
    model_path = ROOT / MODEL_PATHS[method]

    # 1. qiime tools validate
    r = subprocess.run(
        ["qiime", "tools", "validate", str(tax_path)],
        capture_output=True, text=True,
    )
    valid = r.returncode == 0 and "appears to be valid" in r.stdout
    if not valid:
        print(f"FAIL: {method}: qiime tools validate не прошёл (код {r.returncode})")
        ok = False

    # 2. Число уникальных ID и совпадение с rep_seq.qza
    tax = _read_taxonomy(tax_path)
    rep_ids = _read_rep_seq_ids()
    tax_ids = set(tax.keys())
    n_unique = len(tax_ids)
    if n_unique != 3374:
        print(f"FAIL: {method}: {n_unique} уникальных ID, ожидалось 3374")
        ok = False
    if tax_ids != rep_ids:
        missing = rep_ids - tax_ids
        extra = tax_ids - rep_ids
        print(f"FAIL: {method}: набор ID не совпадает с rep_seq.qza")
        if missing:
            print(f"  пропущено: {sorted(missing)[:3]}...")
        if extra:
            print(f"  лишние: {sorted(extra)[:3]}...")
        ok = False

    # 3. Confidence конечные в [0,1], Taxon непустой
    conf_bad = []
    taxon_empty = []
    for sid, (taxon, conf_str) in tax.items():
        try:
            conf = float(conf_str)
        except ValueError:
            conf_bad.append(sid)
            continue
        if not math.isfinite(conf) or not (0.0 <= conf <= 1.0):
            conf_bad.append(sid)
        if not taxon or not taxon.strip():
            taxon_empty.append(sid)
    if conf_bad:
        print(f"FAIL: {method}: {len(conf_bad)} ASV с Confidence вне [0,1] или неконечным")
        ok = False
    if taxon_empty:
        print(f"FAIL: {method}: {len(taxon_empty)} ASV с пустым Taxon")
        ok = False

    # 4. Корневое действие и UUID
    action_yaml = _read_root_action(tax_path)
    parsed = _parse_action(action_yaml)
    if parsed["action"] != "classify_sklearn":
        print(f"FAIL: {method}: корневое действие {parsed['action']}, "
              f"ожидалось classify_sklearn")
        ok = False
    if parsed["reads"] != REP_SEQ_UUID:
        print(f"FAIL: {method}: UUID входа {parsed['reads']}, ожидалось {REP_SEQ_UUID}")
        ok = False
    if parsed["classifier"] != MODEL_UUIDS[method]:
        print(f"FAIL: {method}: UUID классификатора {parsed['classifier']}, "
              f"ожидалось {MODEL_UUIDS[method]}")
        ok = False
    if parsed["confidence"] != "0.7":
        print(f"FAIL: {method}: confidence={parsed['confidence']}, ожидалось 0.7")
        ok = False
    if parsed["n_jobs"] != "1":
        print(f"FAIL: {method}: n_jobs={parsed['n_jobs']}, ожидалось 1")
        ok = False

    # 5. Выборочная сверка трёх детерминированно выбранных ASV
    #    (первые три по лексикографической сортировке ID)
    probe3 = sorted(tax_ids)[:3]
    export_dir = ROOT / "work/qzv_export/full_classification" / method
    export_dir.mkdir(parents=True, exist_ok=True)
    r2 = subprocess.run(
        ["qiime", "tools", "export", "--input-path", str(tax_path),
         "--output-path", str(export_dir)],
        capture_output=True, text=True,
    )
    if r2.returncode != 0:
        print(f"FAIL: {method}: qiime tools export не прошёл (код {r2.returncode})")
        ok = False
    else:
        export_tsv = export_dir / "taxonomy.tsv"
        if not export_tsv.exists():
            print(f"FAIL: {method}: нет экспорта {export_tsv}")
            ok = False
        else:
            export_tax = _read_taxonomy(export_tsv)
            for sid in probe3:
                if sid not in export_tax:
                    print(f"FAIL: {method}: {sid} нет в экспорте")
                    ok = False
                    continue
                if export_tax[sid] != tax[sid]:
                    print(f"FAIL: {method}: {sid} в экспорте != в QZA")
                    ok = False

    print(f"Проверка {method}: validate={'PASS' if valid else 'FAIL'}, "
          f"ID={n_unique}, Confidence∈[0,1]={'PASS' if not conf_bad else 'FAIL'}, "
          f"action={parsed['action']}, reads={parsed['reads']}, "
          f"classifier={parsed['classifier']}, confidence={parsed['confidence']}, "
          f"n_jobs={parsed['n_jobs']}, probe3={probe3}")
    return ok


if __name__ == "__main__":
    results = [check_full_taxonomy(m) for m in ["silva138", "gtdb", "gg2"]]
    print("ИТОГ:", "PASS" if all(results) else "FAIL")
    sys.exit(0 if all(results) else 1)
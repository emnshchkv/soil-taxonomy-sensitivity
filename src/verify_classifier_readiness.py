"""Воспроизводимые сверки готовности классификаторов (AGENT_FULL_CLASSIFICATION_PROMPT).

Выполняет три проверки существующих файлов:
1. Статистика DADA2 по 11 образцам: сверка samples.tsv (run=1) с актуальным QZA dada2_stats.qza.
   Требует ровно 11 образцов и полного совпадения набора.
2. Выбор пяти ID и совпадение их последовательностей с актуальным rep_seq.qza.
3. Проверка трёх проб классификации: тип, число ASV, Confidence, UUID входа и
   классификатора в корневом действии, параметры confidence/n_jobs.

Запуск: conda run --no-capture-output -n qiime2-amplicon-2025.7 python3 src/verify_classifier_readiness.py
Код 0 — все проверки пройдены; 1 — провал.
"""
import hashlib
import os
import re
import sys
import zipfile
from pathlib import Path

import biom
from qiime2 import Artifact

import qiime_runtime  # noqa: F401  (выбирает интерпретатор окружения)

ROOT = Path(__file__).resolve().parent.parent

EXPECTED_SAMPLES = [
    "K_S63", "Kr1_S80", "Kr2_S81", "Kr3_S82", "Kr4_S83", "Kr5_S84",
    "LUB2_S92", "LUB3_S93", "LUB4_S94", "LUB5_S95", "LUB6_S96",
]

# UUID классификаторов (из metadata.yaml каждого QZA)
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
PROBE_PATHS = {
    "silva138": "work/classifier_probe/taxonomy_probe_silva138.qza",
    "gtdb": "work/classifier_probe/taxonomy_probe_gtdb.qza",
    "gg2": "work/classifier_probe/taxonomy_probe_gg2.qza",
}
REP_SEQ_UUID = "978af6f8-91bb-4df2-9e19-b395076e5f3c"
PROBE5_SEQS_UUID = "f191f736-c066-42be-aa44-80e29282c3a8"


def _qza_uuid(path: Path) -> str:
    """UUID артефакта из корневого metadata.yaml."""
    with zipfile.ZipFile(path) as z:
        meta = [n for n in z.namelist()
                if n.endswith("metadata.yaml") and "/provenance/" not in n][0]
        return meta.split("/")[0]


def _read_qza_data(path: Path, suffix: str) -> str:
    """Читает первый файл внутри QZA, оканчивающийся на suffix."""
    with zipfile.ZipFile(path) as z:
        data = [n for n in z.namelist() if n.endswith(suffix)]
        if not data:
            raise FileNotFoundError(f"в {path} нет файла *{suffix}")
        return z.read(data[0]).decode()


def _read_root_action(path: Path) -> str:
    """Читает корневое provenance/action/action.yaml артефакта."""
    with zipfile.ZipFile(path) as z:
        action = [n for n in z.namelist()
                  if n.endswith("provenance/action/action.yaml")][0]
        return z.read(action).decode()


def _parse_action(action_yaml: str) -> dict:
    """Извлекает из action.yaml: действие, UUID входов, параметры."""
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


# --- Проверка 1: статистика DADA2 по 11 образцам из актуального QZA ---
def check_dada2_stats() -> bool:
    """Сверка samples.tsv (run=1) с актуальным QZA dada2_stats.qza.

    Читает stats.tsv напрямую из QZA (не из старого экспорта), требует
    ровно 11 образцов и полного совпадения набора с samples.tsv.
    """
    samples_tsv = ROOT / "results/exact_match_reproduction/provenance/samples.tsv"
    qza = ROOT / "work/qza/dada2_stats.qza"
    if not samples_tsv.exists():
        print(f"FAIL: нет {samples_tsv}")
        return False
    if not qza.exists():
        print(f"FAIL: нет {qza}")
        return False

    # samples.tsv (run=1)
    run1 = {}
    with open(samples_tsv) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            if row.get("run") == "1":
                run1[row["sample-id"]] = row

    # актуальный QZA dada2_stats.qza — читаем stats.tsv изнутри
    stats_text = _read_qza_data(qza, "stats.tsv")
    df = {}
    for line in stats_text.splitlines():
        if line.startswith("#") or line.startswith("sample-id"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 6:
            continue
        df[parts[0]] = {
            "input": parts[1], "filtered": parts[2],
            "denoised": parts[4], "non-chimeric": parts[5],
        }

    # ровно 11 образцов и полное совпадение набора
    qza_ids = set(df.keys())
    run1_ids = set(run1.keys())
    if len(run1_ids) != 11:
        print(f"FAIL: samples.tsv run=1 содержит {len(run1_ids)} образцов, ожидалось 11")
        return False
    if qza_ids != run1_ids:
        print(f"FAIL: набор образцов QZA != samples.tsv run=1")
        print(f"  только в QZA: {sorted(qza_ids - run1_ids)}")
        print(f"  только в samples.tsv: {sorted(run1_ids - qza_ids)}")
        return False

    ok = True
    for sid in sorted(run1_ids):
        for col in ["input", "filtered", "denoised", "non-chimeric"]:
            a = int(run1[sid][col])
            b = int(df[sid][col])
            if a != b:
                print(f"FAIL: {sid} {col}: samples.tsv={a} QZA={b}")
                ok = False
    print(f"Проверка 1 (статистика DADA2): {'PASS' if ok else 'FAIL'}, "
          f"образцов={len(run1_ids)} (ровно 11, набор совпадает)")
    return ok


# --- Проверка 2: выбор 5 ID и совпадение последовательностей с rep_seq.qza ---
def check_probe5_ids() -> bool:
    """Выбор первых 5 ID по лексикографической сортировке из актуального rep_seq.qza.

    Проверяет и ID, и их последовательности.
    """
    qza = ROOT / "work/qza/rep_seq.qza"
    if not qza.exists():
        print(f"FAIL: нет {qza}")
        return False

    fasta = _read_qza_data(qza, "dna-sequences.fasta")
    seqs = {}
    cur_id = None
    for line in fasta.splitlines():
        if line.startswith(">"):
            cur_id = line[1:].strip()
            seqs[cur_id] = []
        elif cur_id is not None:
            seqs[cur_id].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}

    # Число записей в fasta (по заголовкам) должно равняться числу уникальных ID
    n_records = sum(1 for l in fasta.splitlines() if l.startswith(">"))
    if n_records != len(seqs):
        print(f"FAIL: rep_seq.qza содержит {n_records} записей, но {len(seqs)} "
              f"уникальных ID — есть дубликаты")
        return False
    if len(seqs) != 3374:
        print(f"FAIL: rep_seq.qza содержит {len(seqs)} ASV, ожидалось 3374")
        return False

    ids_sorted = sorted(seqs.keys())
    probe5 = ids_sorted[:5]
    expected = [
        "0021cfa2357aba7fe9d1b92465676a8d",
        "0027dadcdc187da02f8f422ae4690999",
        "002d2685b94200ea0724b24286fec6ba",
        "006661799e8f204cef635d45144494e9",
        "0078d1c10db41a22b23bd33997e74aaf",
    ]
    ok = (probe5 == expected)
    if not ok:
        print(f"FAIL: первые 5 ID не совпадают")
        print(f"  получено: {probe5}")
        print(f"  ожидалось: {expected}")
        return False

    # Проверка последовательностей: probe5_seqs.qza должен содержать те же
    # ID и те же последовательности, что и rep_seq.qza для этих 5 ASV.
    probe5_qza = ROOT / "work/classifier_probe/probe5_seqs.qza"
    if not probe5_qza.exists():
        print(f"FAIL: нет {probe5_qza}")
        return False
    probe_fasta = _read_qza_data(probe5_qza, "dna-sequences.fasta")
    probe_seqs = {}
    cur_id = None
    for line in probe_fasta.splitlines():
        if line.startswith(">"):
            cur_id = line[1:].strip()
            probe_seqs[cur_id] = []
        elif cur_id is not None:
            probe_seqs[cur_id].append(line.strip())
    probe_seqs = {k: "".join(v) for k, v in probe_seqs.items()}

    if set(probe_seqs.keys()) != set(probe5):
        print(f"FAIL: probe5_seqs.qza ID != первые 5 ID rep_seq.qza")
        return False
    for sid in probe5:
        if probe_seqs[sid] != seqs[sid]:
            print(f"FAIL: последовательность {sid} в probe5_seqs.qza != rep_seq.qza")
            return False

    print(f"Проверка 2 (выбор 5 ID и последовательности): PASS, "
          f"всего ASV={len(seqs)}, 5 ID и их последовательности совпадают")
    return True


# --- Проверка 3: пробы и provenance ---
def check_probes() -> bool:
    """Проверка трёх проб классификации и provenance классификаторов.

    Для каждой пробы: тип FeatureData[Taxonomy], ровно 5 ASV без дубликатов,
    Confidence∈[0,1], корневое действие classify_sklearn, UUID входа (reads)
    равен probe5_seqs.qza, UUID классификатора соответствует модели,
    confidence=0.7, n_jobs=1.
    """
    ok = True
    for method in ["silva138", "gtdb", "gg2"]:
        probe = ROOT / PROBE_PATHS[method]
        model = ROOT / MODEL_PATHS[method]
        if not probe.exists():
            print(f"FAIL: нет пробы {probe}")
            ok = False
            continue
        if not model.exists():
            print(f"FAIL: нет модели {model}")
            ok = False
            continue

        # Тип и число ASV
        tax = Artifact.load(str(probe))
        df = tax.view(pd.DataFrame)
        n = len(df)
        n_unique = df.index.nunique()
        conf_ok = df["Confidence"].astype(float).between(0, 1).all()
        taxon_nonempty = df["Taxon"].notna().all() and (df["Taxon"].astype(str).str.len() > 0).all()
        if n != 5 or n_unique != 5:
            print(f"FAIL: проба {method}: ASV={n}, уникальных={n_unique}, ожидалось 5")
            ok = False
        if not conf_ok:
            print(f"FAIL: проба {method}: Confidence не все в [0,1]")
            ok = False
        if not taxon_nonempty:
            print(f"FAIL: проба {method}: есть пустой Taxon")
            ok = False

        # Корневое действие и UUID входа/классификатора
        action_yaml = _read_root_action(probe)
        parsed = _parse_action(action_yaml)
        if parsed["action"] != "classify_sklearn":
            print(f"FAIL: проба {method}: корневое действие {parsed['action']}, "
                  f"ожидалось classify_sklearn")
            ok = False
        if parsed["reads"] != PROBE5_SEQS_UUID:
            print(f"FAIL: проба {method}: UUID входа {parsed['reads']}, "
                  f"ожидалось {PROBE5_SEQS_UUID}")
            ok = False
        if parsed["classifier"] != MODEL_UUIDS[method]:
            print(f"FAIL: проба {method}: UUID классификатора {parsed['classifier']}, "
                  f"ожидалось {MODEL_UUIDS[method]}")
            ok = False
        if parsed["confidence"] != "0.7":
            print(f"FAIL: проба {method}: confidence={parsed['confidence']}, ожидалось 0.7")
            ok = False
        if parsed["n_jobs"] != "1":
            print(f"FAIL: проба {method}: n_jobs={parsed['n_jobs']}, ожидалось 1")
            ok = False

        # UUID модели в файле совпадает с ожидаемым
        model_uuid = _qza_uuid(model)
        if model_uuid != MODEL_UUIDS[method]:
            print(f"FAIL: модель {method}: UUID {model_uuid}, ожидалось {MODEL_UUIDS[method]}")
            ok = False

        print(f"Проверка 3 ({method}): проба ASV={n}, Confidence∈[0,1]={conf_ok}, "
              f"action={parsed['action']}, reads={parsed['reads']}, "
              f"classifier={parsed['classifier']}, confidence={parsed['confidence']}, "
              f"n_jobs={parsed['n_jobs']}")
    return ok


if __name__ == "__main__":
    import pandas as pd
    results = [
        check_dada2_stats(),
        check_probe5_ids(),
        check_probes(),
    ]
    print("ИТОГ:", "PASS" if all(results) else "FAIL")
    sys.exit(0 if all(results) else 1)
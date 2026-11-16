"""Version 1 synthetic input generator; no project inputs are read.

Logical contents, UUIDs and ZIP timestamps are deterministic. Minimal QZA
containers exercise the BIOM/HDF5, FASTA and taxonomy readers, not QIIME SDK
provenance validation (that has separate real-artifact checks).
"""
from __future__ import annotations

import io
import uuid
import zipfile
from datetime import datetime
from pathlib import Path

import biom
import h5py
import numpy as np

VERSION = "taxonomy-comparison-synthetic-v1"
SAMPLES = ["K_S63", "Kr1_S80", "Kr2_S81", "Kr3_S82", "Kr4_S83", "Kr5_S84",
           "LUB2_S92", "LUB3_S93", "LUB4_S94", "LUB5_S95", "LUB6_S96"]
MODELS = ["sklearn_silva138", "sklearn_gtdb", "sklearn_gg2"]
FEATURES = [f"synthetic_{i:04d}" for i in range(3374)]


def qza(path, kind, members):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    identity = str(uuid.uuid5(uuid.NAMESPACE_URL, VERSION + "/" + path.name))
    formats = {"FeatureTable[Frequency]": "BIOMV210DirFmt",
               "FeatureData[Sequence]": "DNASequencesDirectoryFormat",
               "FeatureData[Taxonomy]": "TSVTaxonomyDirectoryFormat"}
    contents = {"metadata.yaml": f"uuid: {identity}\ntype: {kind}\nformat: {formats.get(kind, 'TestFormat')}\n",
                "VERSION": "QIIME 2\narchive: 6\nframework: 2025.7.0\n"}
    contents.update({"data/" + name: value for name, value in members.items()})
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, value in sorted(contents.items()):
            info = zipfile.ZipInfo(identity + "/" + name, (2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, value)


def matrix():
    values = np.zeros((len(FEATURES), len(SAMPLES)), dtype=float)
    values[:-1, 1:] = 1
    values[:5, 1:] = np.asarray([9, 1, 2, 3, 4])[:, None]
    values[:2, 2] = 0  # Kr2: jointly assigned ASVs exist, but have zero reads.
    values[-1, 0] = 7  # One control-only ASV, excluded from soil denominators.
    return values


def write_table(root, values=None, features=None, samples=None, kind="FeatureTable[Frequency]"):
    values = matrix() if values is None else values
    features = FEATURES if features is None else features
    samples = SAMPLES if samples is None else samples
    table = biom.Table(values, features, samples, input_is_dense=True)
    payload = io.BytesIO()
    with h5py.File(payload, "w") as handle:
        table.to_hdf5(handle, VERSION, creation_date=datetime(2026, 1, 1))
    qza(Path(root) / "work/qza/ASV_table.qza", kind,
        {"feature-table.biom": payload.getvalue()})


def taxon(index, model):
    # At genus rank A/B categories are match, diff, only_a, only_b, neither.
    if index >= 5:
        return "Unassigned"
    genus = {
        MODELS[0]: ["G1", "G2", "GA", None, None],
        MODELS[1]: ["G1", "G9", None, "GB", None],
        MODELS[2]: ["G1", "G2", None, None, None],
    }[model][index]
    family = "F" if index < 4 else "Unknown_Family"
    parts = ["d__Bacteria", "p__P", "c__C", "o__O", "f__" + family,
             "g__" + (genus or ""), "s__"]
    return ("; " if model == MODELS[2] else ";").join(parts)


def sequence(index):
    digits = "".join("ACGT"[(index >> (2 * shift)) & 3] for shift in range(7, -1, -1))
    return digits + "A" * 172


def create_fixture(root):
    root = Path(root)
    write_table(root)
    qza(root / "work/qza/rep_seq.qza", "FeatureData[Sequence]", {
        "dna-sequences.fasta": "".join(f">{fid}\n{sequence(i)}\n" for i, fid in enumerate(FEATURES))})
    for model in MODELS:
        content = "Feature ID\tTaxon\tConfidence\n" + "".join(
            f"{fid}\t{taxon(i, model)}\t0.9\n" for i, fid in enumerate(FEATURES))
        qza(root / f"work/qza/taxonomy_{model}.qza", "FeatureData[Taxonomy]", {"taxonomy.tsv": content})
    (root / "data/meta").mkdir(parents=True, exist_ok=True)
    (root / "data/meta/metadata.tsv").write_text("sample-id\tsample_type\n" + "".join(
        f"{sid}\t{'negative_control' if sid == 'K_S63' else 'soil'}\n" for sid in SAMPLES))
    (root / "data/refs").mkdir(parents=True, exist_ok=True)
    (root / "data/refs/VERSIONS.md").write_text(
        f"# Synthetic fixture {VERSION}\nNo trained models or private inputs are included.\n")
    return root


def alter_member(path, suffix, transform):
    """Rewrite one member inside a disposable synthetic archive."""
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        members = [(info, archive.read(info.filename)) for info in archive.infolist()]
    found = 0
    with zipfile.ZipFile(path, "w") as archive:
        for info, value in members:
            if info.filename.endswith(suffix):
                value = transform(value)
                found += 1
            archive.writestr(info, value)
    if found != 1:
        raise ValueError(f"Expected one {suffix}; found {found}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.destination.exists():
        parser.error("destination must not exist")
    create_fixture(args.destination)
    print(f"Created {VERSION}: {args.destination}")

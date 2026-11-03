"""Small saved-TSV oracle; deliberately independent from histogram calculations."""
import gzip
import tempfile
import zipfile
from pathlib import Path


def exercise_quality(api, check):
    with tempfile.TemporaryDirectory(prefix="verifier-quality-contract-") as tmp:
        root = Path(tmp)
        for name in ("data/meta", "data/raw", "results", "work/qzv"):
            (root/name).mkdir(parents=True)
        meta = "sample-id\tsample_type\tfile_r1\tsite_code\tsite\n"
        manifest = "sample-id\tabsolute-filepath\n"
        counts = "sample-id\treads_R1\treads_R2\n"
        qzv = "sample ID\tforward sequence count\n"
        pos = "sample-id\tposition\tn_bases\tq25\tq50\tq75\n"
        summary = "sample-id\tn_reads\tn_bases\tmean_q\tpct_q30\n"
        # Two bases per read: columns [10,30] and [20,40]. Manual oracle.
        for sid in sorted(api.SAMPLES):
            fname = sid + ("_first_launch" if sid.startswith("LUB") else "") + "_R1_001.fastq.gz"
            with gzip.open(root/"data/raw"/fname, "wt", encoding="ascii") as f:
                f.write("@r1\nAC\n+\n+5\n@r2\nTG\n+\n?I\n")
            role = "negative_control" if sid == "K_S63" else "soil"
            meta += f"{sid}\t{role}\t{fname}\tTEST\tSYNTHETIC\n"
            manifest += f"{sid}\t/irrelevant-machine/{fname}\n"
            counts += f"{sid}\t2\t2\n"
            qzv += f"{sid}\t2\n"
            pos += f"{sid}\t1\t2\t15\t20\t25\n{sid}\t2\t2\t25\t30\t35\n"
            summary += f"{sid}\t2\t4\t25\t50\n"
        for path,text in (("data/meta/metadata.tsv",meta),
                          ("data/meta/fq-manifest-single.tsv",manifest),
                          ("results/read_counts_R1.tsv",counts),
                          ("results/quality_by_position.tsv",pos),
                          ("results/quality_by_sample.tsv",summary)):
            (root/path).write_text(text)
        with zipfile.ZipFile(root/"work/qzv/reads_summary.qzv", "w") as z:
            z.writestr("synthetic/data/per-sample-fastq-counts.tsv",qzv)
        def audit():
            out = api.Audit()
            api.check_quality(root, out, hashes=False)
            return out
        out = audit()
        check("saved profiles manual oracle", ([],0), ([v.code for v in out.errors],out.exit_code))
        positional = root/"results/quality_by_position.tsv"
        positional.write_text(pos.replace("\t15\t20\t25", "\t15\t21\t25",1))
        check("saved median mutation", ["QUALITY_VALUE"], [v.code for v in audit().errors])
        positional.write_text(pos)
        overall = root/"results/quality_by_sample.tsv"
        overall.write_text(summary.replace("\t25\t50", "\t26\t50",1))
        check("saved summary mutation", ["QUALITY_VALUE"], [v.code for v in audit().errors])
        overall.write_text(summary)
        positional.write_text(pos + pos.splitlines()[1] + "\n")
        check("duplicate position", ["QUALITY_KEY"], [v.code for v in audit().errors])
        positional.unlink()
        check("missing profile is incomplete", 2, audit().exit_code)
        overall.write_text(summary.replace("\t25\t50", "\t26\t50",1))
        check("profile error outranks missing file", 1, audit().exit_code)

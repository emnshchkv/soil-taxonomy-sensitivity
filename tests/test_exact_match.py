#!/usr/bin/env python3
"""Public regressions for availability, exact comparisons and the narrow exception."""
import contextlib
import copy
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from qiime_runtime import selected_python
if __name__ == "__main__":
    selected_python()
import check_exact_match as exact
import verify_denoising as api


def dataset():
    counts = {sid: 100 for sid in api.SAMPLES}
    runs = {}
    for n in (1, 2, 3):
        run = api.load_run(ROOT / "tests/fixtures/denoise-clean" / str(n))
        run.samples = set(counts)
        run.dna = {"a": "A" * 180, "b": "C" * 180}
        run.features = set(run.dna)
        run.matrix = {(f, sid): 35 for f in run.features for sid in counts}
        run.stats = {sid: (100, 90, 80, 70) for sid in counts}
        run.prov = [(*p[:6], {"demultiplexed_seqs": api.INPUT_ID},
                     dict(api.PARAMS, trim_left=20, trunc_len=200), p[8]) for p in run.prov]
        runs[n] = run
    return runs, counts


class ExactMatchTests(unittest.TestCase):
    def setUp(self):
        self.runs, self.counts = dataset()

    def assess(self, *, counts="default", choice=(20, 200), incomplete=()):
        out = api.Audit(incomplete=list(incomplete))
        with contextlib.redirect_stdout(io.StringIO()):
            metrics = exact.assess_runs(self.runs, self.counts if counts == "default" else counts,
                                        choice, out)
        return out, metrics

    def test_complete(self):
        out, metrics = self.assess()
        self.assertEqual(out.exit_code, 0)
        self.assertEqual(metrics, {p: (0, 0, 0, 0) for p in exact.PAIRS})

    def test_missing_counts(self):
        self.assertEqual(self.assess(counts=None)[0].exit_code, 2)

    def test_missing_profile(self):
        self.assertEqual(self.assess(incomplete=["quality_by_sample missing"])[0].exit_code, 2)

    def test_missing_run(self):
        del self.runs[3]
        out, metrics = self.assess()
        self.assertEqual(out.exit_code, 2)
        self.assertEqual(set(metrics), {(1, 2)})

    def test_missing_choice(self):
        self.assertEqual(self.assess(choice=None)[0].exit_code, 2)

    def transfer(self):
        self.runs[3].matrix[("a", "Kr1_S80")] += 1
        self.runs[3].matrix[("b", "Kr1_S80")] -= 1

    def test_transfer(self):
        self.transfer()
        out, metrics = self.assess()
        self.assertEqual(out.exit_code, 1)
        self.assertEqual(metrics[(1, 3)], (0, 2, 1, 0))

    def test_transfer_outranks_missing_counts(self):
        self.transfer()
        out, _ = self.assess(counts=None)
        self.assertEqual(out.exit_code, 1)
        self.assertTrue(out.incomplete)

    def test_error_outranks_missing_run(self):
        del self.runs[1]
        self.transfer()
        self.assertEqual(self.assess()[0].exit_code, 1)

    def test_reused_execution(self):
        self.runs[3] = copy.deepcopy(self.runs[1])
        out, _ = self.assess()
        self.assertEqual(out.exit_code, 1)
        self.assertIn("REUSED_RUN", [v.code for v in out.errors])

    def test_all_eleven_required(self):
        sid = "LUB6_S96"
        del self.counts[sid]
        for run in self.runs.values():
            run.samples.remove(sid)
            del run.stats[sid]
            run.matrix = {k: v for k, v in run.matrix.items() if k[1] != sid}
        self.assertEqual(self.assess()[0].exit_code, 1)

    def retention(self, sid, n, kept):
        self.counts[sid] = n
        for run in self.runs.values():
            run.stats[sid] = (n, n, n, kept)
            run.matrix[("a", sid)] = kept // 2
            run.matrix[("b", sid)] = kept - kept // 2

    def test_known_retention_exception(self):
        self.retention("Kr5_S84", 60078, 1505)
        self.assertEqual(self.assess()[0].exit_code, 0)

    def test_new_kr5_loss_requires_review(self):
        self.retention("Kr5_S84", 60078, 1504)
        self.assertEqual(self.assess()[0].exit_code, 2)

    def test_other_soil_loss_requires_review(self):
        self.retention("Kr4_S83", 100, 20)
        self.assertEqual(self.assess()[0].exit_code, 2)

    def test_zero_soil_is_failure(self):
        self.retention("Kr5_S84", 100, 0)
        self.assertEqual(self.assess()[0].exit_code, 1)

    def test_unrelated_incomplete_survives_exception(self):
        self.retention("Kr5_S84", 60078, 1505)
        self.assertEqual(self.assess(incomplete=["missing file"])[0].exit_code, 2)

    def test_choice_read_failure_is_audited(self):
        with patch.object(exact, "check_quality", return_value=self.counts), \
             patch.object(exact, "read_choice", side_effect=FileNotFoundError("choice")), \
             patch.object(exact, "load_run", side_effect=list(self.runs.values())):
            out, runs, metrics = exact.audit(Path("unused"))
        self.assertEqual(out.exit_code, 2)
        self.assertEqual(len(runs), 3)
        self.assertEqual(set(metrics), exact.PAIRS)


if __name__ == "__main__":
    unittest.main(verbosity=2)

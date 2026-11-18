"""Isolated CLI checks; all papers, Git history, and metadata are synthetic."""
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER = Path(__file__).resolve().parents[1] / "src/check_literature_registry.py"


class LiteratureRegistryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="literature-registry-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "project"
        self.root.mkdir()
        self.base = "manuscript/ai4science_ysda2026/"
        self.evidence = "lit/evidence/20261117"
        self.git("init", "-q")
        records, bib, sections = [], [], []
        for i in range(5):
            key, wid, doi = f"work{i}", f"W{i + 1}", f"10.9999/example{i}"
            online, issue = (2023, 2024) if i == 0 else (2020, 2020)
            primary = f"paper/{key}.txt"
            self.write(primary, f"Synthetic primary source {i}\n")
            record = {"key": key, "wid": wid, "doi": doi, "title": f"Example title {i}", "journal": "Example Journal",
                      "issue_year": issue, "openalex_year": online, "discrepancy": f"Online {online}; issue {issue}",
                      "doi_human_confirmation": "confirmed", "primary_check": "verified",
                      "human_confirmation_quote": "Synthetic approval", "human_confirmation_date": "2026-11-17",
                      "local_sha256": self.sha(primary)}
            records.append(record)
            self.json(f"{self.base}evidence/openalex_{wid}_work.json", {"response": {
                "id": f"https://openalex.org/{wid}", "doi": f"https://doi.org/{doi}",
                "display_name": record["title"], "publication_year": online,
                "primary_location": {"source": {"display_name": "Example Journal"}}}})
            bib.append(f"@article{{{key},\n title = {{{{Example}} title {i}}},\n journal = {{Example Journal}},\n year = {{{issue}}},\n doi = {{{doi}}}\n}}")
            sections.append(f"## {key} — Author ({issue})\n{doi}; {wid}.\n"
                            "**Утверждение [источник].** Synthetic claim.\n"
                            "**Происхождение и локатор.** Abstract.\n"
                            "**Отличие проекта [порождено].** Different dataset.\n"
                            f"**Сверка метаданных.** Online {online}; issue {issue}.\n")
        self.json(f"{self.base}evidence/verified_sources.json", records)
        self.write(f"{self.base}paper.bib", "\n\n".join(bib))
        self.write(f"{self.base}paper.tex", "% \\cite{comment_only}\n\\section{Related}\n\\cite{work0,work1}\n\\citep[p. 4]{work2}\n\\textcite{work3,work4}\n")
        self.git("add", ".")
        self.git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "synthetic sources")
        self.commit = self.git("rev-parse", "HEAD").strip()
        transfers = []
        paths = [f"{self.base}paper.bib", f"{self.base}evidence/verified_sources.json"]
        paths.extend(f'{self.base}evidence/openalex_{r["wid"]}_work.json' for r in records)
        for source in paths:
            destination = "lit/references.bib" if source.endswith("paper.bib") else f"{self.evidence}/{Path(source).name}"
            content = (self.root/source).read_bytes()
            self.write(destination, content)
            transfers.append({"source": source, "destination": destination, "sha256": self.sha(destination), "bytes": len(content)})
        self.manifest = {"source_commit": self.commit, "source_base": self.base, "files": transfers,
                         "human_confirmation": {"quote": "Synthetic approval", "date": "2026-11-17"},
                         "primary_files": [{"key": r["key"], "current_path": f'paper/{r["key"]}.txt', "sha256": r["local_sha256"]} for r in records]}
        self.json("lit/transfer_manifest.json", self.manifest)
        self.write("lit/lit.md", "---\nstate: draft\n---\n# Registry\n\n" + "\n".join(sections))
        # Only the synthetic fixture is modified: the article remains in Git and
        # is removed from this fixture's worktree to model the user's paused text.
        import shutil
        shutil.rmtree(self.root/"manuscript")

    def git(self, *args):
        return subprocess.check_output(["git", "-C", str(self.root), *args], text=True, stderr=subprocess.PIPE)

    def write(self, path, content):
        target = self.root/path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content if isinstance(content, bytes) else content.encode())

    def json(self, path, content):
        self.write(path, json.dumps(content, ensure_ascii=False))

    def sha(self, path):
        return hashlib.sha256((self.root/path).read_bytes()).hexdigest()

    def snapshot(self):
        return {str(p.relative_to(self.root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in self.root.rglob("*") if p.is_file()}

    def run_checker(self, *args):
        before = self.snapshot()
        run = subprocess.run([sys.executable, str(CHECKER), "--root", str(self.root), *args],
                             cwd=self.temp.name, capture_output=True, text=True, timeout=30)
        self.assertEqual(before, self.snapshot(), "checker must be read-only")
        self.assertNotIn("Traceback", run.stderr)
        return run.returncode, json.loads(run.stdout)

    def test_saved_citations_and_metadata_pass_with_no_current_manuscript(self):
        code, report = self.run_checker()
        self.assertEqual(code, 0, report)
        self.assertEqual(report["current_manuscript_status"], "absent_not_checked")
        self.assertEqual(report["historical_citation_keys"], [f"work{i}" for i in range(5)])
        self.assertEqual([r["line"] for r in report["historical_citations"]], [3, 4, 5])

    def test_unknown_historical_citation_fails(self):
        self.write(f"{self.base}paper.tex", "\\cite{foreign_work}\n")
        self.git("add", f"{self.base}paper.tex")
        self.git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "foreign citation")
        self.manifest["source_commit"] = self.git("rev-parse", "HEAD").strip()
        self.json("lit/transfer_manifest.json", self.manifest)
        (self.root/f"{self.base}paper.tex").unlink()
        (self.root/self.base).rmdir()
        code, report = self.run_checker()
        self.assertEqual(code, 1, report)
        self.assertTrue(any("historical unknown citation 'foreign_work'" in e for e in report["errors"]))

    def test_doi_substitution_fails_metadata_check(self):
        path = self.root/"lit/references.bib"
        path.write_text(path.read_text().replace("10.9999/example0", "10.9999/wrong"))
        code, report = self.run_checker()
        self.assertEqual(code, 1, report)
        self.assertIn("work0: DOI mismatch", report["errors"])

    def test_explicit_missing_manuscript_is_incomplete(self):
        code, report = self.run_checker("--article-dir", "absent-article")
        self.assertEqual(code, 2, report)
        self.assertEqual(report["current_manuscript_status"], "requested_path_missing")

    def test_present_current_manuscript_is_also_checked(self):
        self.write(f"{self.base}paper.tex", "\\cite{current_unknown}\n")
        code, report = self.run_checker()
        self.assertEqual(code, 1, report)
        self.assertEqual(report["current_manuscript_status"], "checked")
        self.assertTrue(any("current unknown citation 'current_unknown'" in e for e in report["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)

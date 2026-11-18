"""Read-only verification of the transferred literature registry and citation keys.

Uses saved OpenAlex replies and Git objects; no network, article restoration,
TeX execution, or scientific calculation. PASS concerns this registry and the
scanned sources, never the content or readiness of an absent future manuscript.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


EVIDENCE = "lit/evidence/20261117"
CURRENT_ARTICLE = Path("paper/ai4science_ysda2026")
FIELDS = (
    "**Утверждение [источник].**", "**Происхождение и локатор.**",
    "**Отличие проекта [порождено].**", "**Сверка метаданных.**",
)
CITE = re.compile(
    r"\\(?P<command>[Cc]ite[a-zA-Z]*|[Tt]extcite|[Pp]arencite|[Aa]utocite|footcite|supercite|nocite)"
    r"\*?(?:\s*\[[^\]]*\])*\s*\{(?P<keys>[^}]*)\}"
)


def normalize(value: str) -> str:
    return " ".join(value.replace("{", "").replace("}", "").replace("’", "'")
                    .replace("---", "—").casefold().split())


def doi(value: str) -> str:
    return value.removeprefix("https://doi.org/").casefold()


def bibliography(text: str) -> dict[str, dict[str, str]]:
    """Parse the project's literal brace-valued BibTeX fields, including nested braces."""
    entries = {}
    starts = list(re.finditer(r"@article\s*\{\s*([^,\s]+)\s*,", text, re.I))
    for i, match in enumerate(starts):
        key = match.group(1)
        if key in entries:
            raise ValueError(f"duplicate BibTeX key: {key}")
        block = text[match.end():starts[i + 1].start() if i + 1 < len(starts) else len(text)]
        entries[key] = dict(re.findall(r"^\s*(\w+)\s*=\s*\{(.*?)\}\s*,?\s*$", block, re.M | re.S))
    return entries


def citation_locations(text: str, path: str) -> list[dict]:
    # Retain newlines so locations refer to the original source. Literal common
    # cite commands are supported; this intentionally does not expand TeX macros.
    text = re.sub(r"(?<!\\)%[^\n]*", "", text)
    return [{"path": path, "line": text.count("\n", 0, m.start()) + 1,
             "command": m.group("command"),
             "keys": [key.strip() for key in m.group("keys").split(",")]}
            for m in CITE.finditer(text)]


def verify(root: Path, article_dir: Path | None = None) -> dict:
    result = {"status": "INCOMPLETE", "scope": "registry, saved metadata and literal citation keys; not scientific or article acceptance",
              "errors": [], "incomplete": [], "historical_citations": [], "current_citations": [],
              "current_manuscript_status": "absent_not_checked"}
    errors, missing = result["errors"], result["incomplete"]

    def require(ok, message):
        if not ok:
            errors.append(message)

    def read(relative):
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"path outside project: {relative}")
        return path.read_bytes()

    def git(*args):
        call = subprocess.run(["git", "-C", str(root), *args], capture_output=True)
        if call.returncode:
            raise FileNotFoundError(f"Git object unavailable: {' '.join(args)}")
        return call.stdout

    try:
        manifest = json.loads(read("lit/transfer_manifest.json"))
        records = json.loads(read(f"{EVIDENCE}/verified_sources.json"))
        text = read("lit/lit.md").decode("utf-8")
        bib = bibliography(read("lit/references.bib").decode("utf-8"))
        require(bool(re.match(r"\A---\s*\nstate:\s*draft\s*\n", text)), "lit/lit.md: expected state: draft")
        selected = {record["key"]: record for record in records}
        require(len(selected) == len(records) and 5 <= len(selected) <= 10, "selected works must have 5–10 distinct keys")
        sections = {}
        headings = list(re.finditer(r"^## ([a-z0-9]+) — .+$", text, re.M))
        for i, heading in enumerate(headings):
            key = heading.group(1)
            require(key not in sections, f"duplicate registry section: {key}")
            sections[key] = text[heading.end():headings[i + 1].start() if i + 1 < len(headings) else len(text)]
        require(set(sections) == set(selected), "registry section keys differ from selected works")
        require(set(bib) == set(selected), "BibTeX keys differ from selected works")
        commit, base = manifest["source_commit"], manifest["source_base"]
        require(bool(re.fullmatch(r"[0-9a-f]{40}", commit)), "source_commit must be a full commit ID")
        if errors and not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("invalid source commit")
        result["source_commit"] = commit
        destinations = set()
        for item in manifest["files"]:
            current = read(item["destination"])
            original = git("show", f'{commit}:{item["source"]}')
            destinations.add(item["destination"])
            require(hashlib.sha256(current).hexdigest() == item["sha256"], f'hash mismatch: {item["destination"]}')
            require(current == original and len(current) == item["bytes"], f'transfer differs from Git source: {item["destination"]}')
        required_evidence = {"lit/references.bib", f"{EVIDENCE}/verified_sources.json"}
        primary = {item["key"]: item for item in manifest["primary_files"]}
        require(set(primary) == set(selected), "primary file keys differ from selected works")
        confirmation = manifest["human_confirmation"]
        for key, record in selected.items():
            section, fields = sections.get(key, ""), bib.get(key, {})
            for label in FIELDS:
                require(bool(re.search(re.escape(label) + r"\s*\S", section)), f"{key}: missing content after {label}")
            require(record["doi"] in section and record["wid"] in section, f"{key}: missing DOI or OpenAlex ID in registry")
            reply_path = f'{EVIDENCE}/openalex_{record["wid"]}_work.json'
            required_evidence.add(reply_path)
            response = json.loads(read(reply_path))["response"]
            require(response["id"].rstrip("/").endswith("/" + record["wid"]), f"{key}: OpenAlex ID mismatch")
            require(doi(fields.get("doi", "")) == doi(record["doi"]) == doi(response["doi"]), f"{key}: DOI mismatch")
            require(normalize(fields.get("title", "")) == normalize(record["title"]) == normalize(response["display_name"]), f"{key}: title mismatch")
            require(normalize(fields.get("journal", "")) == normalize(record["journal"]) == normalize(response["primary_location"]["source"]["display_name"]), f"{key}: journal mismatch")
            require(fields.get("year") == str(record["issue_year"]), f"{key}: issue year mismatch")
            require(record["openalex_year"] == response["publication_year"], f"{key}: OpenAlex year mismatch")
            if record["issue_year"] != record["openalex_year"]:
                require(bool(record.get("discrepancy")) and all(str(record[field]) in section for field in ("issue_year", "openalex_year")), f"{key}: issue/online year difference not documented")
            require(record["doi_human_confirmation"] == "confirmed" and record["primary_check"] == "verified", f"{key}: missing saved verification")
            require(record["human_confirmation_quote"] == confirmation["quote"] and record["human_confirmation_date"] == confirmation["date"], f"{key}: confirmation provenance mismatch")
            if key in primary:
                item = primary[key]
                require(hashlib.sha256(read(item["current_path"])).hexdigest() == item["sha256"] == record["local_sha256"], f"{key}: primary source hash mismatch")
        require(required_evidence <= destinations, "required bibliography evidence missing from transfer manifest")
        files = git("ls-tree", "-r", "--name-only", commit, "--", base).decode().splitlines()
        historical_tex = [path for path in files if path.endswith(".tex")]
        if not historical_tex:
            missing.append("historical manuscript has no TeX sources")
        for path in historical_tex:
            result["historical_citations"].extend(citation_locations(git("show", f"{commit}:{path}").decode(), path))
        result["historical_tex_files"] = historical_tex
        current_dir = (root / article_dir).resolve() if article_dir is not None else root / CURRENT_ARTICLE
        if current_dir.is_dir():
            current_files = sorted(current_dir.rglob("*.tex"))
            result["current_manuscript_status"] = "checked" if current_files else "present_without_tex"
            if not current_files:
                missing.append("current manuscript directory has no TeX sources")
            for path in current_files:
                result["current_citations"].extend(citation_locations(path.read_text(), str(path)))
        elif article_dir is not None:
            result["current_manuscript_status"] = "requested_path_missing"
            missing.append(f"explicit article directory missing: {current_dir}")
        for scope in ("historical", "current"):
            used = set()
            for location in result[f"{scope}_citations"]:
                for key in location["keys"]:
                    if key == "*" and location["command"] == "nocite":
                        continue
                    used.add(key)
                    require(key in selected, f'{scope} unknown citation {key!r} at {location["path"]}:{location["line"]}')
            result[f"{scope}_citation_keys"] = sorted(used)
        result["selected_keys"] = sorted(selected)
        result["transferred_files_checked"] = len(destinations)
        result["primary_files_checked"] = len(primary)
    except FileNotFoundError as exc:
        missing.append(str(exc))
    except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
        errors.append(f"invalid literature evidence: {exc}")
    result["status"] = "FAIL" if errors else "INCOMPLETE" if missing else "PASS"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--article-dir", type=Path, help="require this current article directory (default inspection: paper/ai4science_ysda2026 when present)")
    args = parser.parse_args()
    result = verify(args.root.resolve(), args.article_dir)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return {"PASS": 0, "FAIL": 1, "INCOMPLETE": 2}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())

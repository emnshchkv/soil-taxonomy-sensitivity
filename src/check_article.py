"""Read-only validation of the actual step-5 manuscript (not biological truth).

Default: require a complete manuscript and built PDF. --draft permits inspection
before PDF production, but missing required manuscript parts remain INCOMPLETE.
The frozen number checker is an independent verifier, never a generator.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_EVEN
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTICLE = Path('manuscript/ai4science_ysda2026')
BUNDLE = Path('results/paper_step4/20261118')
SECTIONS = ('Постановка', 'Место в области', 'Метод', 'Проверка',
            'Ограничения', 'Эксперимент', 'Воспроизводимость', 'Заключение')
FIELDS = ('claim_id', 'section', 'summary', 'basis', 'source', 'locator',
          'macros', 'citations', 'verification')
KINDS = {'design', 'result', 'literature', 'mathematical', 'verification', 'limitation'}
STYLES = ('ai4science_ysda2026.sty', 'ai4science_ysda2026.bst',
          'algorithm.sty', 'algorithmic.sty', 'fancyhdr.sty')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comments(text):
    return re.sub(r'(?<!\\)%[^\n]*', '', text)


def group(text, pos, opening='{', closing='}'):
    while pos < len(text) and text[pos].isspace():
        pos += 1
    if pos >= len(text) or text[pos] != opening:
        raise ValueError(f'expected {opening!r} near {text[pos:pos+50]!r}')
    start = pos + 1
    depth = 1
    pos += 1
    while pos < len(text):
        if text[pos] == '\\':
            pos += 2
            continue
        if text[pos] == opening:
            depth += 1
        elif text[pos] == closing:
            depth -= 1
            if depth == 0:
                return text[start:pos], start, pos + 1
        pos += 1
    raise ValueError('unbalanced TeX group')


def commands(text, name, nargs=1, optional=False):
    """Find balanced arguments, retaining the exact content and spans."""
    for match in re.finditer(r'\\' + re.escape(name) + r'(?![A-Za-z])\*?', text):
        pos = match.end()
        while pos < len(text) and text[pos].isspace():
            pos += 1
        if optional and pos < len(text) and text[pos] == '[':
            _, _, pos = group(text, pos, '[', ']')
        args, starts = [], []
        try:
            for _ in range(nargs):
                value, start, pos = group(text, pos)
                args.append(value)
                starts.append(start)
        except ValueError:
            continue
        yield match.start(), pos, args, starts


def read_rows(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def separated(text):
    return {s.strip() for s in text.split(';') if s.strip()}


def bib_keys(text):
    return re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', text)


def citation_keys(text):
    keys = set()
    for match in re.finditer(r'\\(?:cite|citep|citet|citealp|citeauthor|citeyear)\*?(?:\[[^\]]*\])*\{([^}]+)\}', text):
        keys.update(k.strip() for k in match.group(1).split(','))
    return keys


def table_cells(text):
    body = re.search(r'\\begin\{tabular\}(.*?)\\end\{tabular\}', text, re.S)
    if not body:
        return []
    result = []
    for line in body.group(1).splitlines():
        if '&' not in line:
            continue
        cells = [s.strip() for s in re.sub(r'\\\\.*$', '', line).split('&')]
        # Actual data rows have plain numeric/NA cells after their row label.
        if all(re.fullmatch(r'(?:\d+(?:\.\d+)?|NA)', c) for c in cells[1:]):
            result.append(cells)
    return result


class Audit:
    def __init__(self):
        self.errors, self.incomplete = [], []
        self.count = 0

    def require(self, condition, message):
        self.count += 1
        if not condition:
            self.errors.append(message)

    def need(self, condition, message):
        self.count += 1
        if not condition:
            self.incomplete.append(message)


def inspect(root=ROOT, article_dir=None, draft=False):
    root = Path(root).resolve()
    directory = Path(article_dir).resolve() if article_dir else root / ARTICLE
    audit = Audit()
    details = {'mode': 'draft' if draft else 'final', 'article_directory': str(directory),
               'scope': 'actual TeX, numeric transfer, claims structure, citations, assets and PDF; semantic review is separately recorded; no biological accuracy verdict'}
    all_files = {}

    def resolve(name, suffix=''):
        p = directory / name
        if suffix and not p.suffix:
            p = p.with_suffix(suffix)
        return p.resolve()

    def expand(path, stack=()):
        path = path.resolve()
        if path in stack:
            audit.errors.append('cyclic TeX input: ' + str(path)); return ''
        if not path.is_file():
            audit.incomplete.append('missing TeX input: ' + str(path)); return ''
        text = comments(path.read_text(encoding='utf-8'))
        all_files[path] = text
        pieces, cursor = [], 0
        inputs = sorted(list(commands(text, 'input')) + list(commands(text, 'include')))
        for start, end, args, _ in inputs:
            audit.require('\\' not in args[0], 'dynamic input path is unsupported: ' + args[0])
            pieces.append(text[cursor:start])
            pieces.append(expand(resolve(args[0], '.tex'), stack + (path,)))
            cursor = end
        pieces.append(text[cursor:])
        return ''.join(pieces)

    try:
        paper = directory / 'paper.tex'
        audit.need(paper.is_file(), 'paper.tex is missing')
        if not paper.is_file():
            return finish(audit, details)
        text = expand(paper)
        begin = text.find('\\begin{document}')
        end = text.find('\\end{document}')
        audit.need(begin >= 0 and end > begin, 'complete document environment is missing')
        preamble = text[:begin] if begin >= 0 else ''
        body = text[begin + len('\\begin{document}'):end] if end > begin >= 0 else text
        details['tex_files'] = [str(p) for p in all_files]
        audit.require(bool(re.search(r'\\usepackage\s*\{ai4science_ysda2026\}', preamble)), 'anonymous template package required without options')
        audit.require(not re.search(r'\\usepackage\s*\[[^]]*(?:final|preprint)[^]]*\]\s*\{ai4science_ysda2026\}', preamble), 'non-anonymous package option')
        audit.require(not re.search(r'\\(?:geometry|setstretch|fontsize)\b|\\(?:setlength|addtolength)\s*\{?\\(?:textwidth|textheight|oddsidemargin|evensidemargin|baselineskip)', preamble), 'template layout override')
        for filename in STYLES:
            copy, original = directory / filename, root / 'ai4science_ysda2026_template' / filename
            audit.need(copy.is_file(), 'missing template copy: ' + filename)
            audit.require(original.is_file() and copy.is_file() and sha(copy) == sha(original), 'template differs: ' + filename)

        # Independently recompute the frozen bundle before trusting its published
        # macros and tables; this imports no generation implementation.
        from check_paper_results import audit_bundle
        numeric = audit_bundle(root / BUNDLE, root)
        audit.require(numeric['status'] == 'PASS', 'frozen numeric verification failed: ' + '; '.join(numeric['issues']))
        details['frozen_numeric_checks'] = numeric['checks']
        rows = read_rows(root / BUNDLE / 'numbers.tsv')
        expected_macros = {r['id']: r['display_value'] for r in rows}
        # The supplementary export is independently reconstructed from its saved
        # evidence, never by importing prepare_article_support.py.
        support = root / 'results/article_step5/20261118'
        validation = root / 'results/paper_step4/validation_20261118'
        check_evidence = json.loads((validation/'numeric_verification.json').read_text())
        test_evidence = json.loads((validation/'numeric_tests.json').read_text())
        pair = [r for r in read_rows(root/BUNDLE/'comparison_pairs.tsv')
                if r['model_a']=='sklearn_silva138' and r['model_b']=='sklearn_gtdb' and r['rank']=='семейство']
        audit.require(len(pair)==1, 'supplementary range selector is not unique')
        extra_expected = {
            'ArticleTransferChecks': str(check_evidence['checks']),
            'ArticleTransferCases': str(len(test_evidence['cases'])),
            'ArticleDetectedMutations': str(sum(c['expected_exit_code']==c['actual_exit_code']==1 for c in test_evidence['cases'])),
            'ArticleSilvaGtdbFamilyDistanceMin': format(Decimal(pair[0]['distance_min']).quantize(Decimal('.000001'),rounding=ROUND_HALF_EVEN),'f'),
            'ArticleSilvaGtdbFamilyDistanceMax': format(Decimal(pair[0]['distance_max']).quantize(Decimal('.000001'),rounding=ROUND_HALF_EVEN),'f'),
        }
        audit.require(test_evidence['case_count']==len(test_evidence['cases']), 'saved test case count differs from actual cases')
        audit.need((support/'extra_numbers.tsv').is_file(), 'supplementary scalar map missing')
        if (support/'extra_numbers.tsv').is_file():
            extra_rows = read_rows(support/'extra_numbers.tsv')
            audit.require({r['macro']:r['value'] for r in extra_rows} == extra_expected and len(extra_rows)==len(extra_expected), 'supplementary scalar values differ from independent selectors')
            expected_sources = {
                'ArticleTransferChecks': 'results/paper_step4/validation_20261118/numeric_verification.json',
                'ArticleTransferCases': 'results/paper_step4/validation_20261118/numeric_tests.json',
                'ArticleDetectedMutations': 'results/paper_step4/validation_20261118/numeric_tests.json',
                'ArticleSilvaGtdbFamilyDistanceMin': str(BUNDLE/'comparison_pairs.tsv'),
                'ArticleSilvaGtdbFamilyDistanceMax': str(BUNDLE/'comparison_pairs.tsv'),
            }
            for r in extra_rows:
                audit.require(r.get('source')==expected_sources.get(r['macro']), 'supplementary scalar source differs: '+r['macro'])
        expected_macros.update(extra_expected)
        definitions = {}
        for _, _, args, _ in commands(preamble, 'newcommand', 2, optional=True):
            name = args[0].lstrip('\\')
            definitions.setdefault(name, []).append(args[1])
        for name, expected in expected_macros.items():
            audit.require(definitions.get(name) == [expected], f'macro {name}: expected sole definition {expected!r}, found {definitions.get(name)}')
        for name in expected_macros:
            audit.require(not re.search(r'\\(?:renewcommand|providecommand|def|edef|gdef|xdef)\s*\{?\\' + name + r'\b', text), 'numeric macro overridden: ' + name)
        frozen_used = set(re.findall(r'\\((?:Project|Cmp|Selected|Article)[A-Za-z]+)\b', body))
        audit.require(frozen_used <= expected_macros.keys(), 'unknown numeric macros: ' + ', '.join(sorted(frozen_used - expected_macros.keys())))
        details['frozen_macros_used'] = sorted(frozen_used)

        sections = [(start, args[0]) for start, _, args, _ in commands(body, 'section')]
        # The supplied .sty creates this heading inside the environment, so it
        # is a genuine section context even without a literal \section in TeX.
        ai_start = body.find('\\begin{aiusage}')
        if ai_start >= 0:
            sections.append((ai_start, 'Использование ИИ-ассистентов'))
            sections.sort()
        ordered_main = [name for _, name in sections if name in SECTIONS]
        audit.need(ordered_main == list(SECTIONS), 'eight required sections are missing or out of order: ' + repr(ordered_main))
        for name in SECTIONS:
            audit.require(sum(x == name for _, x in sections) <= 1, 'duplicate section: ' + name)
        positions = {name: pos for pos, name in sections}
        main_end = body.find('\\label{main:end}')
        audit.need('\\label{main:start}' in body and main_end >= 0, 'main page-boundary labels are missing')
        if main_end >= 0 and 'Заключение' in positions:
            audit.require(main_end > positions['Заключение'], 'main:end precedes conclusion')
        if 'Ограничения' in positions and 'Эксперимент' in positions:
            limitation = body[positions['Ограничения']:positions['Эксперимент']]
            limit_claims = list(commands(limitation, 'Claim', 2))
            # One or several claims may each contain complete sentences.
            content = ' '.join(args[1] for _, _, args, _ in limit_claims)
            content = re.sub(r'\\(?:cite\w*|ref|label)\{[^}]*\}', '', content)
            content = re.sub(r'\\[A-Za-z]+', ' ', content)
            content = content.translate(str.maketrans('', '', '{}'))
            sentences = re.split(r'(?<=[.!?])\s+(?=[А-ЯЁA-Z])', content.strip()) if content.strip() else []
            audit.need(5 <= len(sentences) <= 8, f'limitations must contain 5–8 sentences; found {len(sentences)}')
            details['limitation_sentences'] = len(sentences)
            for sentence in sentences:
                audit.require(bool(re.search(r'\b(?:Если|При|Без|Когда|Поскольку|Пока|В отсутствие|В условиях)\b', sentence)), 'limitation lacks an explicit applicability condition: ' + sentence[:100])

        bib_path = directory / 'references.bib'
        audit.need(bib_path.is_file(), 'references.bib is missing')
        if bib_path.is_file():
            audit.require(bib_path.read_bytes() == (root / 'lit/references.bib').read_bytes(), 'bibliography differs from approved registry')
        keys = set(bib_keys((root / 'lit/references.bib').read_text()))
        used_citations = citation_keys(body)
        audit.require(used_citations <= keys, 'unknown citations: ' + ', '.join(sorted(used_citations - keys)))
        audit.require(not re.search(r'\\nocite\s*\{\*\}', body), 'nocite{*} would hide uncited records')
        audit.need(bool(list(commands(body, 'bibliography'))), 'BibTeX bibliography command is missing')
        for _, _, args, _ in commands(body, 'bibliography'):
            audit.require(args == ['references'], 'unexpected bibliography input')
        styles = [args[0] for _,_,args,_ in commands(body,'bibliographystyle')]
        audit.need(styles == ['ai4science_ysda2026'], 'required BibTeX style command missing/different')
        details['citation_keys'] = sorted(used_citations)

        claim_commands = list(commands(body, 'Claim', 2))
        claim_ids = [args[0] for _, _, args, _ in claim_commands]
        audit.require(len(claim_ids) == len(set(claim_ids)), 'duplicate Claim identifiers')
        claims_path = directory / 'claims.tsv'
        audit.need(claims_path.is_file(), 'claims.tsv is missing')
        claim_rows = read_rows(claims_path) if claims_path.is_file() else []
        if claim_rows:
            audit.require(set(FIELDS) <= claim_rows[0].keys(), 'claims.tsv schema does not contain all required fields')
        claims = {r.get('claim_id'): r for r in claim_rows}
        audit.require(len(claims) == len(claim_rows), 'duplicate claims.tsv identifiers')
        audit.need(bool(claim_commands), 'no traceable content claims')
        audit.require(set(claim_ids) == claims.keys(), 'TeX/map claim IDs differ: ' + repr(sorted(set(claim_ids) ^ claims.keys())))
        generated = list(commands(body, 'porozhdeno', optional=True))
        audit.need(bool(generated), 'generated text is not marked')
        audit.need('\\begin{aiusage}' in body and '\\end{aiusage}' in body, 'AI usage table is missing')
        # Generated prose is not silently reclassified as human-authored.
        for _, _, args, _ in commands(body, 'napisano', optional=True):
            audit.require(not args[0].strip(), 'human-authored text requires separate documented human evidence; none registered')
        for _, _, args, _ in commands(body, 'istochnik', optional=True):
            if args[0].strip():
                audit.require(bool(citation_keys(args[0])), 'source text lacks a citation inside istochnik')
        for start, finish_pos, args, _ in claim_commands:
            claim_id, content = args
            if claim_id not in claims:
                continue
            row = claims[claim_id]
            audit.require(all(row.get(field, '').strip() for field in ('section','summary','basis','source','locator','verification')), claim_id + ': empty provenance field')
            audit.require(row.get('basis') in KINDS, claim_id + ': unsupported basis')
            audit.require(any(a <= start and b >= finish_pos for a, b, _, _ in generated), claim_id + ': not enclosed in generated-text mark')
            active_section = next((name for pos, name in reversed(sections) if pos < start), 'Аннотация')
            audit.require(row.get('section') == active_section, f'{claim_id}: section mismatch ({active_section!r})')
            for source in separated(row.get('source','')):
                audit.require(not Path(source).is_absolute() and (root/source).is_file(), claim_id + ': unsupported/missing source ' + source)
            actual_macros = set(re.findall(r'\\((?:Project|Cmp|Selected|Article)[A-Za-z]+)\b', content))
            audit.require(actual_macros == separated(row.get('macros','')), claim_id + ': macro list does not match actual claim')
            actual_citations = citation_keys(content)
            audit.require(actual_citations == separated(row.get('citations','')), claim_id + ': citation list does not match actual claim')
            if row.get('basis') == 'literature':
                audit.require(bool(actual_citations), claim_id + ': literature claim without citation')
            # A literal empirical number cannot bypass the verified macros. Zero
            # and unity may occur as definitional boundary values, not estimates.
            if row.get('basis') in {'result', 'verification'}:
                scrub = re.sub(r'\\[A-Za-z]+', '', content)
                scrub = re.sub(r'(?:Kr\d+|LUB\d+)(?:\\?_S\d+)?|K\\?_S63|(?:QIIME\s*2|DADA2|R[12]|GG2|Greengenes2)|\d{4}-\d{2}-\d{2}', '', scrub)
                literals = re.findall(r'(?<![\w])\d+(?:[.,]\d+)?(?![\w])', scrub)
                audit.require(not [n for n in literals if n not in {'0','1','2'}], claim_id + ': literal result numbers outside verified macros: ' + repr(literals))
        details['claims_checked'] = len(claim_commands)

        # Every paragraph carrying prose in main sections must have a claim. Table
        # payloads are checked below as matrices, and title/labels are structural.
        main = body[:main_end] if main_end >= 0 else body
        spans = [(a,b) for a,b,_,_ in commands(main,'Claim',2)]
        for a,b in reversed(spans):
            main = main[:a] + ' CLAIM ' + main[b:]
        main = re.sub(r'\\begin\{tabular\}.*?\\end\{tabular\}', '', main, flags=re.S)
        for paragraph in re.split(r'\n\s*\n', main):
            if 'CLAIM' in paragraph:
                continue
            cleaned = paragraph
            for command in ('section','label','icmltitle','icmlkeywords','caption','icmltitlerunning'):
                for a,b,_,_ in reversed(list(commands(cleaned,command,optional=True))):
                    cleaned = cleaned[:a] + cleaned[b:]
            cleaned = re.sub(r'\\[A-Za-z]+', '', cleaned)
            words = re.findall(r'[А-Яа-яЁё]{3,}', cleaned)
            audit.require(len(words) < 3, 'unmapped prose paragraph: ' + ' '.join(words[:14]))

        reviews = directory / 'section_review.md'
        audit.need(reviews.is_file(), 'section_review.md is missing')
        if reviews.is_file():
            review_text = reviews.read_text()
            for section in SECTIONS:
                audit.need(section in review_text, 'section review not recorded: ' + section)
            audit.need('вслух' in review_text.lower(), 'review does not yet disclose oral-reading status')

        # Check the actually included table matrices, preserving row and model
        # order. Header translation is allowed; data-cell translation is not.
        tables = [(path, table_cells(value)) for path, value in all_files.items() if table_cells(value)]
        expected_tables = {name: table_cells((root/BUNDLE/name).read_text()) for name in ('comparison_main.tex','family_main.tex')}
        for name, expected_table in expected_tables.items():
            matching = [(path, matrix) for path,matrix in tables if len(matrix) == len(expected_table) and len(matrix[0]) == len(expected_table[0])]
            audit.need(len(matching) == 1, 'expected one actual table of shape ' + name)
            if len(matching) != 1:
                continue
            path, matrix = matching[0]
            audit.require([r[1:] for r in matrix] == [r[1:] for r in expected_table], 'actual numeric table differs: ' + str(path))
            if name == 'family_main.tex':
                audit.require([r[0] for r in matrix] == [r[0] for r in expected_table], 'sample row order differs: ' + str(path))
            else:
                aliases = ({'Phylum','Тип'}, {'Class','Класс'}, {'Order','Порядок'}, {'Family','Семейство'}, {'Genus','Род'})
                audit.require(all(row[0] in allowed for row,allowed in zip(matrix,aliases)), 'rank row order differs: ' + str(path))
            header = all_files[path].split('\\hline',2)[1] if '\\hline' in all_files[path] else all_files[path]
            audit.require(bool(re.search(r'SILVA.*?GTDB.*?(?:GG2|Greengenes)', header, re.S)), 'model header order is not SILVA/GTDB/GG2')
        audit.need(len(tables) >= 2, 'numeric tables missing from actual inputs')
        audit.require(len(tables) <= 2, 'unexpected numeric tables in actual inputs')

        figures = read_rows(root/BUNDLE/'figure_sources.tsv')
        expected_paths = {(root/r['path']).resolve():r for r in figures}
        included, placements = [], {}
        for start, _, args, starts in commands(body, 'includegraphics', optional=True):
            path = resolve(args[0], '.pdf'); included.append(path)
            placements.setdefault(path, []).append(body[start:starts[0]-1])
            audit.need(path.is_file(), 'missing figure: ' + str(path))
            audit.require(path in expected_paths, 'unregistered figure: ' + str(path))
            if path in expected_paths and path.is_file():
                audit.require(sha(path) == expected_paths[path]['sha256'], 'figure hash differs: ' + str(path))
        audit.need(set(included) == expected_paths.keys(), 'all four registered appendix figure assets are required')
        for path, options in placements.items():
            if path not in expected_paths:
                continue
            is_overview = expected_paths[path]['figure_id'] == 'prokaryotic_family_overview'
            if not is_overview or len(options) == 1:
                audit.require(len(options) == 1, 'model PDF must appear exactly once: ' + path.name)
                audit.require(not any('trim' in o or 'viewport' in o for o in options), 'single-page figure must preserve its full source: ' + path.name)
            else:
                # A tall overview may be split into three adjacent bands. Verify
                # the complete vertical partition of the unchanged registered PDF,
                # rather than accepting arbitrary cropping or duplicate pictures.
                audit.require(len(options) == 3, 'overview split requires exactly three complete adjacent bands')
                info = subprocess.run(['pdfinfo',str(path)],capture_output=True,text=True)
                size = re.search(r'Page size:\s*([\d.]+) x ([\d.]+) pts',info.stdout)
                audit.require(info.returncode == 0 and bool(size), 'cannot determine original overview page size')
                intervals=[]
                if size:
                    height=Decimal(size.group(2))
                    for option in options:
                        trim=re.search(r'trim\s*=\s*\{?([^},\]]+)',option)
                        values=re.findall(r'-?\d+(?:\.\d+)?(?:bp)?',trim.group(1)) if trim else []
                        clipping=bool(re.search(r'(?:,|\[)\s*clip(?:\s*=\s*true)?\s*(?:,|\])',option))
                        audit.require(len(values)==4 and clipping, 'overview band must have four explicit bp trim values and clipping')
                        if trim:
                            tokens=trim.group(1).split()
                            audit.require(len(tokens)==4 and all(re.fullmatch(r'(?:0(?:\.0+)?|\d+(?:\.\d+)?bp)',v) for v in tokens), 'overview trim must use PDF bp units (or unitless zero)')
                        audit.require('width=\\textwidth' in option.replace(' ',''), 'overview bands must retain full text width')
                        if len(values)==4:
                            left,bottom,right,top = [Decimal(v.removesuffix('bp')) for v in values]
                            audit.require(left==right==0 and 0<=top<height-bottom<=height, 'invalid overview band or horizontal crop')
                            intervals.append((top,height-bottom))
                    audit.require(len(intervals)==3 and intervals[0][0]==0 and intervals[-1][1]==height and all(abs(a[1]-b[0])<=Decimal('.01') for a,b in zip(intervals,intervals[1:])), 'overview bands leave a gap, overlap, or omit source content')
                    details['overview_bands_bp']=[[str(a),str(b)] for a,b in intervals]
        details['registered_figure_assets']=len(set(included))
        details['figure_placements']=len(included)
        for _, _, args, _ in commands(body,'ref'):
            labels = {a[0] for _,_,a,_ in commands(body,'label')}
            (audit.need if draft else audit.require)(args[0] in labels, 'unresolved internal ref: ' + args[0])

        if not draft:
            check_pdf(audit, details, directory, all_files, body, keys)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        audit.errors.append(f'{type(exc).__name__}: {exc}')
    return finish(audit, details)


def check_pdf(audit, details, directory, all_files, body, keys):
    pdf, aux, log = (directory / ('paper.'+ext) for ext in ('pdf','aux','log'))
    for path in (pdf,aux,log,directory/'paper.bbl'):
        audit.need(path.is_file(), 'missing build output: ' + path.name)
    if not all(path.is_file() for path in (pdf,aux,log,directory/'paper.bbl')):
        return
    for executable in ('pdftotext','pdfinfo'):
        audit.need(bool(shutil.which(executable)), 'PDF inspection tool missing: ' + executable)
    if not shutil.which('pdftotext') or not shutil.which('pdfinfo'):
        return
    extract = subprocess.run(['pdftotext','-layout',str(pdf),'-'],capture_output=True,text=True)
    info = subprocess.run(['pdfinfo',str(pdf)],capture_output=True,text=True)
    audit.require(extract.returncode == info.returncode == 0, 'PDF extraction/metadata command failed')
    pdf_text = extract.stdout
    normalized_pdf = ' '.join(pdf_text.split())
    pages = pdf_text.split('\f')
    if pages and not pages[-1].strip():
        pages.pop()
    details['pdf_total_pages'] = len(pages)
    audit.require(not re.search(r'/Users/|/home/|/private/tmp/|[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}',pdf_text+'\n'+info.stdout), 'PDF text/metadata contains private path or email')
    author = re.search(r'^Author:[ \t]*(.*)$',info.stdout,re.M)
    audit.require(not author or not author.group(1).strip() or author.group(1).strip() in {'Anonymous','Анонимные авторы'}, 'PDF Author metadata is not anonymous')
    audit.require('Анонимные авторы' in pdf_text, 'anonymous author label missing from PDF')
    audit.require('title suppressed due to excessive size' not in normalized_pdf.lower(), 'running title suppressed in PDF')
    audit.require(not re.search(r'\bРабота\s+подана\b|\bSubmitted\s+to\b', normalized_pdf, re.I), 'PDF incorrectly claims the draft was submitted')
    audit.require(bool(re.search(r'Локальный\s+обезличенный\s+черновик', normalized_pdf, re.I)) and bool(re.search(r'не\s+подана', normalized_pdf, re.I)), 'PDF must disclose local anonymous draft and no submission')
    audit.require(not re.search(r'\?\?|\[\?\]',pdf_text), 'unresolved PDF reference/citation')
    labels = {}
    for _,_,args,_ in commands(aux.read_text(),'newlabel',2):
        try:
            _,_,end = group(args[1],0)
            number,_,_ = group(args[1],end)
            if number.isdigit():
                labels[args[0]] = int(number)
        except ValueError:
            pass
    start,text_end = labels.get('main:start'),labels.get('main:end')
    audit.need(start is not None and text_end is not None, 'built PDF main boundary labels missing')
    if start is not None and text_end is not None:
        main_body=body[:body.find('\\label{main:end}')]
        table_labels=[]
        for match in re.finditer(r'\\begin\{table\*?\}(.*?)\\end\{table\*?\}',main_body,re.S):
            table_labels.extend(args[0] for _,_,args,_ in commands(match.group(1),'label'))
        audit.require(len(table_labels)==2, 'both main tables need page labels')
        for label in table_labels:
            audit.require(label in labels, 'main table has no resolved PDF page: '+label)
        # Full-width floats may follow the conclusion on their own page. They
        # count toward main-text volume; the conclusion label alone is not a
        # truthful end boundary in a two-column document.
        end=max([text_end]+[labels[label] for label in table_labels if label in labels])
        details['main_text_end']=text_end
        details['main_table_pages']={label:labels.get(label) for label in table_labels}
        details['main_page_start'],details['main_page_end'] = start,end
        details['main_page_count'] = end-start+1
        audit.require(1 <= start <= end <= len(pages) and 2 <= end-start+1 <= 4, 'main text must occupy 2–4 actual PDF pages')
        main_text = '\n'.join(pages[start-1:end])
        for section in SECTIONS:
            audit.require(section in main_text, 'main section absent from marked PDF pages: ' + section)
        audit.require('Литература' not in main_text and 'Использование ИИ-ассистентов' not in main_text, 'excluded sections overlap marked main pages')
        for label in table_labels:
            audit.require(label in labels and start <= labels[label] <= end, 'main table floats outside marked main pages: '+label)
    warnings = log.read_text(errors='replace')
    audit.require('Title exceeds size limitations' not in warnings, 'running title exceeds template size limit')
    audit.require(not re.search(r'(?m)^!|undefined references|Citation .* undefined|Reference .* undefined|Label\(s\) may have changed',warnings), 'unresolved LaTeX error/reference/rerun in final log')
    overfull = re.findall(r'Overfull \\[hv]box \(([^)]*)\)',warnings)
    details['overfull_boxes'] = overfull
    audit.require(not overfull, 'overfull boxes remain in final log')
    latest_source = max([p.stat().st_mtime for p in all_files] + [(directory/'references.bib').stat().st_mtime])
    audit.require(pdf.stat().st_mtime >= latest_source, 'PDF is older than an included TeX/BibTeX source')
    cited = citation_keys(body)
    bbl_keys = {args[0] for _,_,args,_ in commands((directory/'paper.bbl').read_text(),'bibitem',optional=True)}
    audit.require(cited <= bbl_keys <= keys, 'built bibliography does not match cited/approved keys')


def finish(audit, details):
    return {'status':'FAIL' if audit.errors else 'INCOMPLETE' if audit.incomplete else 'PASS',
            'checks':audit.count,'errors':audit.errors,'incomplete':audit.incomplete,**details}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--article-dir',type=Path)
    parser.add_argument('--draft',action='store_true',help='Inspect manuscript without requiring a built PDF; incomplete manuscript is still INCOMPLETE/2.')
    parser.add_argument('--report',type=Path,help='Write a new JSON report; never overwrite.')
    args = parser.parse_args()
    result = inspect(args.root,args.article_dir,args.draft)
    encoded = json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.report:
        with args.report.open('x',encoding='utf-8') as f:
            f.write(encoded)
    print(encoded,end='')
    return {'PASS':0,'FAIL':1,'INCOMPLETE':2}[result['status']]


if __name__ == '__main__':
    sys.exit(main())

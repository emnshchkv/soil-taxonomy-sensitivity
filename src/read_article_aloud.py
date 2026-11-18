"""Synthesize local speech from the current article, without speaker playback.

This is an accessibility/editorial aid, not a scientific verifier or a record
of human listening. Run through the project's existing conda environment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
from datetime import datetime, timezone

from check_article import commands, comments

ROOT = Path(__file__).resolve().parents[1]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def unwrap(text: str, name: str, replacement=None, optional=False) -> str:
    found = list(commands(text, name, 1, optional=optional))
    for start, end, args, _ in reversed(found):
        value = args[0] if replacement is None else replacement(args[0])
        text = text[:start] + value + text[end:]
    return text


def spoken_text(text: str, macros: dict[str, str]) -> str:
    """Preserve Claim wording; normalize only TeX presentation for speech."""
    text = comments(text)
    for name, value in sorted(macros.items(), key=lambda item: -len(item[0])):
        text = re.sub(r'\\' + re.escape(name) + r'(?![A-Za-z])(?:\{\})?',
                      lambda _: value, text)
    for name in ('istochnik', 'porozhdeno', 'napisano', 'texttt', 'textit',
                 'textbf', 'emph', 'mbox', 'ensuremath', 'url', 'path'):
        text = unwrap(text, name, optional=name in {'istochnik', 'porozhdeno', 'napisano'})
    for name in ('cite', 'citep', 'citet'):
        text = unwrap(text, name, lambda _: '', optional=True)
    references = {
        'tab:comparison': 'сравнения', 'tab:family': 'семейственных профилей',
        'fig:overview': 'обзора', 'fig:ggtwo': 'модели Greengenes два',
        'fig:silva': 'модели SILVA', 'fig:gtdb': 'модели GTDB',
    }
    text = unwrap(text, 'ref', lambda key: references.get(key, key))
    text = text.replace(r'\frac{1}{2}\sum_t |p_t-q_t|',
                        'половина суммы по таксонам модулей разностей долей пэ тэ и ку тэ')
    text = text.replace(r'\[', ' ').replace(r'\]', ' ')
    text = text.replace(r'\%', ' процентов').replace(r'\_', '_')
    text = text.replace(r'\&', ' и ').replace('~', ' ').replace('$', '')
    text = text.replace('p_t', 'пэ тэ').replace('q_t', 'ку тэ')
    text = text.replace('---', ' — ').replace('--', ' — ')
    text = text.replace('{', '').replace('}', '')
    text = re.sub(r'(?<=\d)\.(?=\d)', ',', text)
    text = re.sub(r'\s+([,.;:])', r'\1', text)
    text = re.sub(r'\s+', ' ', text).strip()
    if re.search(r'\\[A-Za-z]', text):
        raise ValueError('Unrendered TeX command: ' + repr(text))
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--article-dir', type=Path,
                        default=ROOT / 'paper/ai4science_ysda2026')
    parser.add_argument('--out-dir', type=Path,
                        default=ROOT / 'results/article_step5/reading_20261118')
    parser.add_argument('--voice', default='Milena')
    parser.add_argument('--rate', type=int, default=220)
    args = parser.parse_args()
    article = args.article_dir.resolve()
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f'Refusing to overwrite existing output: {out}')
    paper = article / 'paper.tex'
    inputs = [args_[0] for _, _, args_, _ in commands(comments(paper.read_text()), 'input')]
    macro_paths = [(article / value).resolve() for value in inputs
                   if value.endswith('numbers.tex')]
    sections = [(article / value).resolve() for value in inputs
                if value.startswith('sections/')]
    if len(sections) != 9:
        raise ValueError('Expected abstract and eight section files')
    source_paths = [paper, *macro_paths, *sections, Path(__file__).resolve()]
    source_bytes = {str(p): p.read_bytes() for p in source_paths}
    macros = {}
    for p in macro_paths:
        for _, _, values, _ in commands(comments(source_bytes[str(p)].decode()),
                                       'newcommand', 2):
            name = values[0].lstrip('\\')
            if name in macros:
                raise ValueError('Duplicate macro: ' + name)
            macros[name] = values[1]
    prepared = []
    for index, p in enumerate(sections, 1):
        content = source_bytes[str(p)].decode()
        titles = list(commands(content, 'section'))
        title = titles[0][2][0] if titles else 'Аннотация'
        claims = [values for _, _, values, _ in commands(content, 'Claim', 2)]
        if not claims:
            raise ValueError('No claims in ' + str(p))
        text = title + '.\n\n' + '\n\n'.join(spoken_text(c[1], macros) for c in claims) + '\n'
        prepared.append((index, p, title, claims, text))
    out.mkdir(parents=True)
    audio_dir = Path(tempfile.mkdtemp(prefix='article-step5-readings-', dir='/private/tmp'))
    manifest = {
        'status': 'RUNNING', 'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Local speech synthesis from actual Claim text; no speaker playback or human listening.',
        'voice': args.voice, 'rate_words_per_minute': args.rate,
        'macro_count': len(macros), 'macros': macros,
        'audio_directory': str(audio_dir),
        'sources': [{'path': name, 'sha256': digest(data)} for name, data in source_bytes.items()],
        'transformations': [
            'Numeric macros expanded from the actual included TeX files.',
            'Authorship/style wrappers unwrapped; citation commands omitted for speech.',
            'Cross-reference labels rendered as descriptions; decimal separator localized.',
            'Distance equation normalized to a spoken description without changing its meaning.'
        ],
        'human_listening_performed': False, 'speaker_playback_performed': False,
        'sections': [],
    }
    manifest_path = out / 'manifest.json'

    def save():
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')

    save()
    for index, p, title, claims, text in prepared:
        stem = f'{index:02d}_{p.stem}'
        text_path = out / (stem + '.txt')
        audio_path = audio_dir / (stem + '.aiff')
        text_path.write_text(text)
        command = ['/usr/bin/say', '-v', args.voice, '-r', str(args.rate),
                   '-f', str(text_path), '-o', str(audio_path), '--file-format=AIFF']
        started = datetime.now(timezone.utc).isoformat()
        proc = subprocess.run(command, capture_output=True, text=True, timeout=300)
        row = {'title': title, 'source': str(p), 'claims': [c[0] for c in claims],
               'text_path': str(text_path), 'text_sha256': digest(text_path.read_bytes()),
               'audio_path': str(audio_path), 'command': command,
               'started_at_utc': started, 'finished_at_utc': datetime.now(timezone.utc).isoformat(),
               'exit_code': proc.returncode, 'stdout': proc.stdout, 'stderr': proc.stderr,
               'audio_sha256': digest(audio_path.read_bytes()) if audio_path.exists() else None,
               'audio_bytes': audio_path.stat().st_size if audio_path.exists() else 0}
        manifest['sections'].append(row)
        save()
        print(f'{title}: exit {proc.returncode}; {row["audio_bytes"]} bytes', flush=True)
        if proc.returncode or not row['audio_bytes']:
            manifest['status'] = 'FAIL'
            save()
            return 1
    changed = [name for name, data in source_bytes.items()
               if Path(name).read_bytes() != data]
    manifest['changed_sources_during_synthesis'] = changed
    manifest['status'] = 'STALE' if changed else 'PASS'
    manifest['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    save()
    print(f'{manifest["status"]}: {manifest_path}', flush=True)
    return 1 if changed else 0


if __name__ == '__main__':
    raise SystemExit(main())

"""Draw already prepared bacterial taxonomy profiles; AI4Science-агент.

This module makes no domain selection, taxon selection, aggregation, or fraction
calculation. The caller supplies all rows, labels, colours and stacking orders.
Use render_barplots from the project's existing conda environment.
"""
from collections import defaultdict
import math
from pathlib import Path
import re
import textwrap


RANK_FILES = {
    'тип': 'phylum', 'phylum': 'phylum',
    'класс': 'class', 'class': 'class',
    'порядок': 'order', 'order': 'order',
    'семейство': 'family', 'family': 'family',
    'род': 'genus', 'genus': 'genus',
}
RANK_TITLES = {
    'phylum': 'Тип', 'class': 'Класс', 'order': 'Порядок',
    'family': 'Семейство', 'genus': 'Род',
}
MODEL_TITLES = {
    'sklearn_silva138': 'SILVA',
    'sklearn_gtdb': 'GTDB',
    'sklearn_gg2': 'Greengenes2',
}
KINDS = {'selected', 'other_assigned', 'unclassified'}
HATCHES = {'selected': None, 'other_assigned': '///', 'unclassified': 'xxx'}
REQUIRED = {'model', 'sample', 'rank', 'category', 'category_kind', 'count',
            'total_reads', 'fraction', 'order', 'label', 'color'}


def _integer(value, field):
    if isinstance(value, bool):
        raise ValueError(f'{field} must be a nonnegative integer')
    number = int(value)
    if number < 0 or float(value) != number:
        raise ValueError(f'{field} must be a nonnegative integer')
    return number


def _wrap_label(label):
    # Preserve supplied line breaks; wrapping never changes the key or label text.
    return '\n'.join(textwrap.fill(line, width=33, break_long_words=True,
                                   break_on_hyphens=False)
                     for line in str(label).splitlines())


def _prepare(rows, samples, models, ranks):
    """Validate drawing inputs without selecting taxa or recomputing fractions."""
    if len(models) != 3 or len(set(models)) != 3:
        raise ValueError('exactly three distinct models are required')
    if len(samples) != 10 or len(set(samples)) != 10:
        raise ValueError('exactly ten distinct soil samples are required')
    if not 1 <= len(ranks) <= 5 or len({RANK_FILES.get(rank) for rank in ranks}) != len(ranks):
        raise ValueError('supply one to five distinct supported ranks')
    if any(rank not in RANK_FILES for rank in ranks):
        raise ValueError('unsupported rank')
    groups = defaultdict(dict)
    for source in rows:
        if not REQUIRED.issubset(source):
            raise ValueError('row lacks required drawing fields: ' +
                             ', '.join(sorted(REQUIRED - source.keys())))
        row = dict(source)
        model, sample, rank = row['model'], row['sample'], row['rank']
        if model not in models or sample not in samples or rank not in ranks:
            raise ValueError('row refers to a model, sample or rank outside the requested order')
        kind = row['category_kind']
        if kind not in KINDS:
            raise ValueError(f'unknown category_kind: {kind}')
        row['count'] = _integer(row['count'], 'count')
        row['total_reads'] = _integer(row['total_reads'], 'total_reads')
        row['order'] = _integer(row['order'], 'order')
        if row['total_reads'] == 0:
            if row['fraction'] != 'NA' or row['count'] != 0:
                raise ValueError('zero-denominator rows must have fraction NA and count zero')
        else:
            value = float(row['fraction'])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError('positive-denominator fraction must be finite and in [0, 1]')
            row['fraction'] = value
        if not re.fullmatch(r'#[0-9A-Fa-f]{6}(?:[0-9A-Fa-f]{2})?', str(row['color'])):
            raise ValueError('color must be a supplied hexadecimal colour')
        row['label'] = str(row['label'])
        if not row['label'].strip():
            raise ValueError('legend label must not be empty')
        key = (rank, model, sample)
        if row['category'] in groups[key]:
            raise ValueError('duplicate category in one model/sample/rank bar')
        groups[key][row['category']] = row
    categories = {}
    for rank in ranks:
        for model in models:
            first = groups.get((rank, model, samples[0]))
            if not first:
                raise ValueError('missing model/sample/rank bar')
            ordered = sorted(first.values(), key=lambda row: row['order'])
            if len({row['order'] for row in ordered}) != len(ordered):
                raise ValueError('stack orders must be distinct within each model/rank')
            kinds = [row['category_kind'] for row in ordered]
            if kinds.count('other_assigned') != 1 or kinds.count('unclassified') != 1:
                raise ValueError('each model/rank needs separate other_assigned and unclassified categories')
            for sample in samples:
                current = groups.get((rank, model, sample))
                if not current or current.keys() != first.keys():
                    raise ValueError('all samples need all supplied categories, including zero rows')
                totals = {row['total_reads'] for row in current.values()}
                if len(totals) != 1:
                    raise ValueError('one bar has inconsistent bacterial denominators')
                total = next(iter(totals))
                if sum(row['count'] for row in current.values()) != total:
                    raise ValueError('supplied category counts do not sum to the bar denominator')
                if total and not math.isclose(sum(row['fraction'] for row in current.values()),
                                               1.0, rel_tol=0, abs_tol=1e-10):
                    raise ValueError('supplied fractions must sum to one; renderer does not renormalize')
                for category, row in current.items():
                    if any(row[field] != first[category][field]
                           for field in ('order', 'label', 'color', 'category_kind')):
                        raise ValueError('category display metadata changes between samples')
            categories[(rank, model)] = ordered
    return groups, categories


def render_barplots(rows, sample_order, model_order, rank_order, out_dir):
    """Render requested ranks as three model panels, returning three Paths per rank.

    ``out_dir`` is the result directory; outputs are written under ``figures/``.
    Each input dict has model, sample, rank, category, category_kind, count,
    total_reads, fraction, order, label and color. Fractions are supplied numbers
    in [0, 1], or the string 'NA' when the bacterial denominator is zero.
    All categories, including zero rows, must exist in every sample of a panel.

    Counts under bars are the supplied bacterial denominators. A zero denominator
    is shown by a short dash and NA, never by an invented percentage. Figure width
    is 180 mm; height grows only when the supplied legends require more room.
    """
    samples, models, ranks = list(sample_order), list(model_order), list(rank_order)
    groups, categories = _prepare(rows, samples, models, ranks)
    # Imports stay local so importing this helper neither chooses a backend nor
    # changes global plot styles used elsewhere in the project.
    import matplotlib as mpl
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    figures = Path(out_dir) / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    outputs = []
    style = {'font.family': 'DejaVu Sans', 'font.size': 8,
             'axes.labelsize': 8, 'axes.titlesize': 9,
             'xtick.labelsize': 8, 'ytick.labelsize': 8,
             'legend.fontsize': 8, 'pdf.fonttype': 42,
             'svg.fonttype': 'none', 'svg.hashsalt': 'taxonomy-barplots',
             'hatch.linewidth': 0.35}
    with mpl.rc_context(style):
        for rank in ranks:
            rank_file = RANK_FILES[rank]
            wrapped = {(model, row['category']): _wrap_label(row['label'])
                       for model in models for row in categories[(rank, model)]}
            # Legend height accounts for every wrapped line at the actual 8 pt.
            legend_points = max(
                sum(9.2 * (wrapped[(model, row['category'])].count('\n') + 1)
                    for row in categories[(rank, model)]) +
                1.0 * (len(categories[(rank, model)]) - 1)
                for model in models)
            panel_mm = max(33.0, legend_points * 25.4 / 72.0 + 1.0)
            gap_mm, bottom_mm, top_mm = 11.5, 30.0, 9.0
            height_mm = 3 * panel_mm + 2 * gap_mm + bottom_mm + top_mm
            figure = Figure(figsize=(180 / 25.4, height_mm / 25.4))
            FigureCanvasAgg(figure)
            axes = figure.subplots(3, 1, sharex=True)
            figure.subplots_adjust(left=0.115, right=0.655,
                                   bottom=bottom_mm / height_mm,
                                   top=1 - top_mm / height_mm,
                                   hspace=gap_mm / panel_mm)
            figure.supylabel('Доля бактериальных ридов, %', x=0.024, fontsize=8)
            figure.text(0.115, 2.5 / height_mm,
                        'N — число ридов после отбора; NA — нет отобранных ридов.',
                        ha='left', va='bottom', fontsize=8)
            positions = list(range(len(samples)))
            for panel_index, (model, axis) in enumerate(zip(models, axes)):
                axis.set_ylim(0, 100)
                axis.set_xlim(-0.6, len(samples) - 0.4)
                axis.set_yticks([0, 25, 50, 75, 100])
                axis.grid(axis='y', color='#DDDDDD', linewidth=0.4)
                axis.set_axisbelow(True)
                for side in ('top', 'right'):
                    axis.spines[side].set_visible(False)
                axis.spines['left'].set_linewidth(0.5)
                axis.spines['bottom'].set_linewidth(0.5)
                axis.tick_params(width=0.5, length=2.5)
                axis.set_title(MODEL_TITLES.get(model, str(model)), loc='left', pad=5)
                if panel_index == 0:
                    axis.set_title(RANK_TITLES[rank_file], loc='right', pad=5)
                totals = [next(iter(groups[(rank, model, sample)].values()))['total_reads']
                          for sample in samples]
                bottom = [0.0] * len(samples)
                handles = []
                for display in categories[(rank, model)]:
                    category, kind = display['category'], display['category_kind']
                    values = [100.0 * groups[(rank, model, sample)][category]['fraction']
                              if total else 0.0 for sample, total in zip(samples, totals)]
                    edge = '#FFFFFF' if kind == 'selected' else '#555555'
                    axis.bar(positions, values, width=0.76, bottom=bottom,
                             color=display['color'], edgecolor=edge, linewidth=0.25,
                             hatch=HATCHES[kind], zorder=3)
                    bottom = [old + value for old, value in zip(bottom, values)]
                    handles.append(Patch(facecolor=display['color'], edgecolor=edge,
                                         linewidth=0.4, hatch=HATCHES[kind],
                                         label=wrapped[(model, category)]))
                for position, total in zip(positions, totals):
                    axis.text(position, -0.075, str(total),
                              transform=axis.get_xaxis_transform(), ha='center',
                              va='top', fontsize=8, clip_on=False)
                    if total == 0:
                        axis.plot([position - 0.22, position + 0.22], [1.2, 1.2],
                                  color='#555555', linewidth=0.9, zorder=5)
                        axis.text(position, 4.0, 'NA', ha='center', va='bottom',
                                  fontsize=8, color='#555555', zorder=5)
                axis.text(-0.08, -0.075, 'N', transform=axis.transAxes,
                          ha='right', va='top', fontsize=8, clip_on=False)
                axis.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.025, 1.0),
                            borderaxespad=0, frameon=False, handlelength=1.15,
                            handletextpad=0.45, labelspacing=0.12,
                            borderpad=0, fontsize=8)
                axis.set_xticks(positions)
                if panel_index == len(models) - 1:
                    axis.set_xticklabels(samples, rotation=55, ha='right',
                                         rotation_mode='anchor')
                    axis.tick_params(axis='x', pad=16)
                else:
                    axis.tick_params(axis='x', labelbottom=False, bottom=False)
            for extension in ('pdf', 'png', 'svg'):
                path = figures / f'taxonomy_{rank_file}.{extension}'
                if extension == 'pdf':
                    metadata = {'Creator': 'AI4Science-агент', 'Author': '',
                                'CreationDate': None, 'ModDate': None}
                elif extension == 'svg':
                    metadata = {'Creator': 'AI4Science-агент', 'Date': None}
                else:
                    metadata = {'Software': 'AI4Science-агент'}
                figure.savefig(path, format=extension, dpi=300,
                               facecolor='white', metadata=metadata)
                outputs.append(path)
            figure.clear()
    return outputs

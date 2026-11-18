"""Draw prepared per-sample family segments; AI4Science-агент.

The caller supplies the mask, thresholds, counts, fractions, labels and colours.
This helper does not select ASVs, infer taxonomy, or aggregate rare families.
"""
from collections import defaultdict
from fractions import Fraction
import math
from pathlib import Path
import re
import textwrap


MODELS = {'sklearn_silva138': ('SILVA', 'silva'),
          'sklearn_gtdb': ('GTDB', 'gtdb'), 'sklearn_gg2': ('Greengenes2', 'gg2')}
KINDS = {'selected', 'other_assigned', 'unclassified'}


def prepare(rows, samples, models):
    if len(samples) != 10 or len(set(samples)) != 10 or len(models) != 3 or set(models) != set(MODELS):
        raise ValueError('Expected ten distinct soil samples and the three project models')
    bars = defaultdict(dict)
    for source in rows:
        row = dict(source)
        if row['model'] not in models or row['sample'] not in samples or row['category_kind'] not in KINDS:
            raise ValueError('Unexpected model, sample or display category kind')
        for field in ('count', 'total_reads', 'order'):
            number = int(row[field])
            if number < 0 or float(row[field]) != number:
                raise ValueError('Expected nonnegative integer ' + field)
            row[field] = number
        if row['total_reads'] == 0:
            if row['fraction'] != 'NA' or row['count']:
                raise ValueError('N=0 requires NA fractions and zero segment counts')
        else:
            row['fraction'] = float(row['fraction'])
            if not math.isfinite(row['fraction']) or not 0 <= row['fraction'] <= 1:
                raise ValueError('Invalid supplied segment fraction')
        if not re.fullmatch(r'#[0-9a-fA-F]{6}', row['color']) or not row['label'].strip():
            raise ValueError('Expected a hex colour and a nonempty supplied label')
        key = row['model'], row['sample']
        if row['category'] in bars[key]:
            raise ValueError('Duplicate segment')
        bars[key][row['category']] = row
    categories = {}
    for model in models:
        first = bars.get((model, samples[0]), {})
        if not first:
            raise ValueError('Missing model/sample bar')
        ordered = sorted(first.values(), key=lambda row: row['order'])
        if len({row['order'] for row in ordered}) != len(ordered):
            raise ValueError('Repeated display order')
        kinds = [row['category_kind'] for row in ordered]
        if kinds.count('other_assigned') != 1 or kinds.count('unclassified') != 1:
            raise ValueError('Separate rare-family and unknown-family categories are required')
        for sample in samples:
            current = bars.get((model, sample), {})
            if current.keys() != first.keys():
                raise ValueError('Every sample needs every panel-legend category, including zeros')
            totals = {row['total_reads'] for row in current.values()}
            if len(totals) != 1 or sum(row['count'] for row in current.values()) != next(iter(totals)):
                raise ValueError('Supplied counts and per-bar denominator disagree')
            if next(iter(totals)) and not math.isclose(sum(row['fraction'] for row in current.values()), 1, abs_tol=1e-10, rel_tol=0):
                raise ValueError('Supplied fractions do not sum to one; no renormalization is performed')
            for category, row in current.items():
                if any(row[field] != first[category][field] for field in ('order', 'label', 'color', 'category_kind', 'full_key')):
                    raise ValueError('Display metadata changes between samples')
        categories[model] = ordered
    return bars, categories


def render_family_profiles(rows, sample_order, model_order, out_dir, threshold=0.02):
    """Return twelve Paths: overview and three model figures, PDF/PNG/SVG each."""
    samples, models = list(sample_order), list(model_order)
    bars, categories = prepare(rows, samples, models)
    threshold = Fraction(str(threshold))
    if not 0 < threshold <= 1:
        raise ValueError('Display threshold must be in (0,1]')
    percent = format(float(threshold * 100), 'g').replace('.', ',')
    import matplotlib as mpl
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.patches import Patch
    destination = Path(out_dir) / 'figures'
    destination.mkdir(parents=True, exist_ok=True)
    outputs = []
    settings = {'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9,
                'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'prokaryotic-families',
                'hatch.linewidth': 0.35}
    with mpl.rc_context(settings):
        for group, filename in [(models, 'overview')] + [([model], MODELS[model][1]) for model in models]:
            width_mm = 180 if len(group) > 1 else 220
            wrapped = {(model, row['category']): '\n'.join(textwrap.fill(line, width=33,
                       break_long_words=True, break_on_hyphens=False) for line in row['label'].splitlines())
                       for model in group for row in categories[model]}
            legend_pt = max(sum(9.2 * (wrapped[model, row['category']].count('\n') + 1)
                               for row in categories[model]) + len(categories[model]) - 1 for model in group)
            panel_mm = max(36, legend_pt * 25.4 / 72 + 1)
            gap_mm, bottom_mm, top_mm = 11.5, 35, 9
            height_mm = len(group) * panel_mm + (len(group) - 1) * gap_mm + bottom_mm + top_mm
            figure = Figure(figsize=(width_mm / 25.4, height_mm / 25.4))
            FigureCanvasAgg(figure)
            axes = figure.subplots(len(group), 1, sharex=True, squeeze=False)[:, 0]
            figure.subplots_adjust(left=0.115, right=0.655, bottom=bottom_mm / height_mm,
                                   top=1 - top_mm / height_mm, hspace=gap_mm / panel_mm)
            figure.supylabel('Доля оставшихся ридов Bacteria + Archaea, %', x=0.025, fontsize=8)
            figure.text(0.115, 7 / height_mm, f'Порог {percent}% применяется отдельно к каждому образцу.', fontsize=8)
            figure.text(0.115, 2.5 / height_mm, 'N — число ридов после отбора; NA — нет отобранных ридов.', fontsize=8)
            positions = list(range(len(samples)))
            for index, (model, axis) in enumerate(zip(group, axes)):
                axis.set(xlim=(-0.6, len(samples) - 0.4), ylim=(0, 100), yticks=[0, 25, 50, 75, 100])
                axis.grid(axis='y', color='#DDDDDD', linewidth=0.4)
                axis.set_axisbelow(True)
                for side in ('top', 'right'):
                    axis.spines[side].set_visible(False)
                for side in ('left', 'bottom'):
                    axis.spines[side].set_linewidth(0.5)
                axis.tick_params(width=0.5, length=2.5, labelsize=8)
                axis.set_title(MODELS[model][0], loc='left', pad=5)
                if index == 0:
                    axis.set_title('Семейства', loc='right', pad=5)
                totals = [next(iter(bars[model, sample].values()))['total_reads'] for sample in samples]
                bottom = [0.0] * len(samples)
                handles = []
                for display in categories[model]:
                    category, kind = display['category'], display['category_kind']
                    values = [100 * bars[model, sample][category]['fraction'] if total else 0
                              for sample, total in zip(samples, totals)]
                    hatch = 'xxx' if kind == 'unclassified' else None
                    edge = '#555555' if kind == 'unclassified' else '#FFFFFF'
                    axis.bar(positions, values, width=0.76, bottom=bottom, color=display['color'],
                             edgecolor=edge, linewidth=0.25, hatch=hatch, zorder=3)
                    bottom = [old + value for old, value in zip(bottom, values)]
                    handles.append(Patch(facecolor=display['color'], edgecolor=edge, hatch=hatch,
                                         linewidth=0.4, label=wrapped[model, category]))
                for position, total in zip(positions, totals):
                    axis.annotate(str(total), xy=(position, 0), xycoords=axis.get_xaxis_transform(),
                                  xytext=(0, -6), textcoords='offset points', ha='center',
                                  va='top', fontsize=8, annotation_clip=False)
                    if total == 0:
                        axis.plot([position - 0.22, position + 0.22], [1.2, 1.2], color='#555555', linewidth=0.9)
                        axis.text(position, 4, 'NA', ha='center', va='bottom', fontsize=8, color='#555555')
                axis.annotate('N', xy=(-0.08, 0), xycoords='axes fraction',
                              xytext=(0, -6), textcoords='offset points', ha='right',
                              va='top', fontsize=8, annotation_clip=False)
                axis.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.025, 1), borderaxespad=0,
                            frameon=False, handlelength=1.15, handletextpad=0.45, labelspacing=0.12,
                            borderpad=0, fontsize=8)
                axis.set_xticks(positions)
                if index == len(group) - 1:
                    axis.set_xticklabels(samples, rotation=55, ha='right', rotation_mode='anchor')
                    axis.tick_params(axis='x', pad=20)
                else:
                    axis.tick_params(axis='x', labelbottom=False, bottom=False)
            for extension in ('pdf', 'png', 'svg'):
                path = destination / f'prokaryotic_family_{filename}.{extension}'
                metadata = ({'Creator': 'AI4Science-агент', 'Author': '', 'CreationDate': None, 'ModDate': None}
                            if extension == 'pdf' else {'Creator': 'AI4Science-агент', 'Date': None}
                            if extension == 'svg' else {'Software': 'AI4Science-агент'})
                figure.savefig(path, format=extension, dpi=300, facecolor='white', metadata=metadata)
                outputs.append(path)
            figure.clear()
    return outputs

#!/usr/bin/env python3
"""
Plot u errors from ADRSolver parameter sweeps.

Examples:
  # u_L2 vs time, one line per NumRuns value (fix h and p):
  python plot_errors.py --xaxis time --color n --fix-h 0.2 --fix-p 5

  # Final u_L2 vs NumRuns on log-log axes, one line per AdaptBL_h value (fix p):
  python plot_errors.py --xaxis n --color h --fix-p 5

  # Final u_L2 vs NUMMODES, one line per h value:
  python plot_errors.py --xaxis p --color h --fix-n 32
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

METRIC_COL = {'u_L2': 1, 'u_Linf': 2, 'u_H1': 3}
FILE_RE = re.compile(r'ErrorFile_n_(\d+)_h_(.+)_p_(\d+)(?:_m_(.+))?\.err$')

AXIS_LABELS = {
    'time': 'Time',
    'n':    'Number of runs (n)',
    'h':    'Maximum layer height (h)',
    'p':    'Polynomial order (p)',
}

METRIC_LABELS = {
    'u_L2':  'L2 error',
    'u_Linf': 'Linf error',
    'u_H1':  'H1 error',
}

# Colormaps:
#   primary   = overall L2 from Projection runs  (--combine-both)
#   secondary = Projection solution transfer L2   (--combine-both)
#   tertiary  = ALE solution transfer L2          (--combine-both)
_CMAPS  = {'n': 'Blues',   'h': 'Greens',  'p': 'Oranges'}
_CMAPS2 = {'n': 'Reds',    'h': 'Purples', 'p': 'Blues'}
_CMAPS3 = {'n': 'Greens',  'h': 'Reds',    'p': 'Purples'}
_CMAPS4 = {'n': 'Purples', 'h': 'Oranges', 'p': 'Greens'}

# Per-method line style, marker, and legend-column order in non-combine mode
_METHOD_LS     = {'Projection': '-',    'ALE': '-',           'FullInterp': '-'}
_METHOD_MARKER = {'Projection': 'o',    'ALE': 's',           'FullInterp': '^'}
_METHOD_LABEL  = {'Projection': 'Proj', 'ALE': 'ALE',         'FullInterp': 'FullInterp'}
_METHOD_ORDER  = {'Projection': 0,      'ALE': 1,             'FullInterp': 2}


def parse_filename(fname):
    m = FILE_RE.match(os.path.basename(fname))
    if not m:
        return None
    n = int(m.group(1))
    h = float(m.group(2).replace('dot', '.'))
    p = int(m.group(3))
    method = m.group(4) if m.group(4) else 'Projection'
    return n, h, p, method


def get_var(run, key):
    n, h, p, method, _ = run
    return {'n': n, 'h': h, 'p': p}[key]


def main():
    parser = argparse.ArgumentParser(
        description='Plot ADRSolver error files from a parameter sweep.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--xaxis', choices=['time', 'n', 'h', 'p'], default='time',
                        help='Variable on the x-axis (default: time)')
    parser.add_argument('--color', choices=['n', 'h', 'p'], required=True,
                        help='Variable to distinguish lines by colour')
    parser.add_argument('--fix-n', type=int,   default=None, metavar='N',
                        help='Restrict to runs with this NumRuns value')
    parser.add_argument('--rename-n', action='store_true',
                        help='Relabel n as cycles (= n − 1) on axes and in legends')
    parser.add_argument('--fix-h', type=float, default=None, metavar='H',
                        help='Restrict to runs with this AdaptBL_h value')
    parser.add_argument('--fix-p', type=int,   default=None, metavar='P',
                        help='Restrict to runs with this NUMMODES value')
    parser.add_argument('--metric', choices=['u_L2', 'u_Linf', 'u_H1'], default='u_L2',
                        help='Error metric to plot (default: u_L2)')
    parser.add_argument('--results-dir', default='results', metavar='DIR',
                        help='Directory containing .err files (default: results)')
    parser.add_argument('--sub-n1', action='store_true',
                        help='Subtract the equivalent direct projection error from each value')
    parser.add_argument('--combine-both', action='store_true',
                        help='Overlay overall (Proj), Projection solution transfer, and ALE solution transfer '
                             'on one plot with three distinct colour progressions. '
                             'Implies --sub-n1 for the solution transfer series.')
    parser.add_argument('--exclude', action='append', default=[], metavar='KEY=VALUE',
                        help='Exclude runs where KEY matches VALUE. '
                             'Keys: n (int), h (float), p (int), method (str). '
                             'Repeatable, e.g. --exclude method=ALE --exclude p=3')
    parser.add_argument('--logy', action='store_true',
                        help='Logarithmic y axis')
    parser.add_argument('--logx', action='store_true',
                        help='Logarithmic x axis')
    parser.add_argument('--save', metavar='FILE',
                        help='Save figure to FILE instead of displaying it')
    parser.add_argument('--subslides', action='store_true',
                        help='Generate a sequence of PNGs revealing one colour-group at a time '
                             '(stem taken from --save, or "plot" by default)')
    args = parser.parse_args()

    if args.rename_n:
        AXIS_LABELS['n'] = 'Number of adaptive cycles'
    x_off = -1 if (args.rename_n and args.xaxis == 'n') else 0
    color_key = 'cycles' if (args.rename_n and args.color == 'n') else args.color
    c_disp = (lambda v: v - 1) if (args.rename_n and args.color == 'n') else (lambda v: v)

    # Parse --exclude specs into typed predicates
    _EXCL_TYPES = {'n': int, 'h': float, 'p': int, 'method': str}
    exclude_filters = []
    for spec in args.exclude:
        if '=' not in spec:
            sys.exit(f'--exclude: expected KEY=VALUE, got "{spec}"')
        key, _, val = spec.partition('=')
        if key not in _EXCL_TYPES:
            sys.exit(f'--exclude: unknown key "{key}". Valid keys: {list(_EXCL_TYPES)}')
        for v in val.split(','):
            try:
                typed_val = _EXCL_TYPES[key](v)
            except ValueError:
                sys.exit(f'--exclude: cannot parse "{v}" as {_EXCL_TYPES[key].__name__} for key "{key}"')
            exclude_filters.append((key, typed_val))

    def is_excluded(n, h, p, method):
        for key, val in exclude_filters:
            run_val = {'n': n, 'h': h, 'p': p, 'method': method}[key]
            if key == 'h':
                if np.isclose(run_val, val):
                    return True
            elif run_val == val:
                return True
        return False

    combine = args.combine_both

    if args.xaxis != 'time' and args.xaxis == args.color:
        sys.exit(f'Error: --xaxis and --color cannot both be "{args.xaxis}"')

    mcol = METRIC_COL[args.metric]

    # Discover and filter files
    pattern = os.path.join(args.results_dir, 'ErrorFile_n_*_h_*_p_*.err')
    all_files = sorted(glob.glob(pattern))
    if not all_files:
        sys.exit(f'No error files found matching {pattern}')

    runs = []
    for f in all_files:
        params = parse_filename(f)
        if params is None:
            continue
        n, h, p, method = params
        if args.fix_n is not None and n != args.fix_n:
            continue
        if args.fix_h is not None and not np.isclose(h, args.fix_h):
            continue
        if args.fix_p is not None and p != args.fix_p:
            continue
        if is_excluded(n, h, p, method):
            continue
        runs.append((n, h, p, method, f))

    # In single-mode with --sub-n1, drop n=1 runs (they are the baseline).
    # In combine mode the same filtering applies for the solution transfer series,
    # but we keep runs as-is here and handle it in the plot loop.
    if args.sub_n1 and not combine:
        runs = [(n, h, p, method, f) for n, h, p, method, f in runs if n != 1]

    if not runs:
        sys.exit('No runs matched the given filters.')

    # Build n=1 baseline lookup: (h, p) -> final metric value
    baseline = {}
    if args.sub_n1 or combine:
        for f in all_files:
            params = parse_filename(f)
            if params is None:
                continue
            n, h, p, method = params
            if n == 1:
                data = np.loadtxt(f, skiprows=1)
                baseline[(h, p)] = data[-1, mcol]

    def get_baseline(h, p):
        key = next((k for k in baseline if np.isclose(k[0], h) and k[1] == p), None)
        if key is None:
            sys.exit(f'--sub-n1: no n=1 baseline found for h={h}, p={p}')
        return baseline[key]

    # Split runs by method
    km_runs  = [(n, h, p, m, f) for (n, h, p, m, f) in runs if m == 'Projection']
    ale_runs = [(n, h, p, m, f) for (n, h, p, m, f) in runs if m == 'ALE']

    # Build colour maps (colour_vals from all runs so both methods share the same palette)
    color_vals = sorted(set(get_var(run, args.color) for run in runs))
    n_cols = len(color_vals)
    # Sample 0.35–0.9 to avoid near-white at the light end
    _samples = [0.35 + 0.55 * i / max(n_cols - 1, 1) for i in range(n_cols)]

    cmap = plt.get_cmap(_CMAPS[args.color])
    color_map = {v: cmap(s) for v, s in zip(color_vals, _samples)}

    cmap2 = plt.get_cmap(_CMAPS2[args.color])
    color_map2 = {v: cmap2(s) for v, s in zip(color_vals, _samples)}
    cmap3 = plt.get_cmap(_CMAPS3[args.color])
    color_map3 = {v: cmap3(s) for v, s in zip(color_vals, _samples)}
    cmap4 = plt.get_cmap(_CMAPS4[args.color])
    color_map4 = {v: cmap4(s) for v, s in zip(color_vals, _samples)}

    _METHOD_CMAP = {'Projection': color_map2, 'ALE': color_map3, 'FullInterp': color_map4}

    fig, ax = plt.subplots(figsize=(8, 5))
    seen_labels = set()

    print(f"{len(runs)} runs:")
    for run in runs:
        print(run)

    if args.xaxis == 'time':
        if combine:
            # Overall + Projection solution transfer from km_runs
            for run in km_runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                cv = get_var(run, args.color)

                label = f'{color_key} = {c_disp(cv)} (overall)'
                ax.plot(data[:, 0], data[:, mcol],
                        color=color_map[cv],
                        label=label if label not in seen_labels else '_nolegend_')
                seen_labels.add(label)

                label2 = f'{color_key} = {c_disp(cv)} (Proj sol. transfer)'
                ax.plot(data[:, 0], data[:, mcol] - get_baseline(h, p),
                        color=color_map2[cv],
                        label=label2 if label2 not in seen_labels else '_nolegend_')
                seen_labels.add(label2)

            # ALE solution transfer from ale_runs
            for run in ale_runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                cv = get_var(run, args.color)
                label3 = f'{color_key} = {c_disp(cv)} (ALE sol. transfer)'
                ax.plot(data[:, 0], data[:, mcol] - get_baseline(h, p),
                        color=color_map3[cv],
                        label=label3 if label3 not in seen_labels else '_nolegend_')
                seen_labels.add(label3)

        else:
            # One line per (method, cv) combination; method → linestyle + colour palette
            for run in runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                cv = get_var(run, args.color)
                yvals = data[:, mcol]
                if args.sub_n1:
                    yvals = yvals - get_baseline(h, p)
                mlabel = _METHOD_LABEL.get(method, method)
                label = f'{color_key} = {c_disp(cv)} ({mlabel})'
                cmap_m = _METHOD_CMAP.get(method, color_map)
                ax.plot(data[:, 0], yvals,
                        color=cmap_m[cv],
                        linestyle=_METHOD_LS.get(method, '-'),
                        label=label if label not in seen_labels else '_nolegend_')
                seen_labels.add(label)

    else:
        # One point per run: final-timestep error (last row = post-adaptation value)
        if combine:
            groups  = defaultdict(list)  # overall KM
            groups2 = defaultdict(list)  # Projection solution transfer
            groups3 = defaultdict(list)  # ALE solution transfer

            for run in km_runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                final_err = data[-1, mcol]
                xv = get_var(run, args.xaxis) + x_off
                cv = get_var(run, args.color)
                groups[cv].append((xv, final_err))
                groups2[cv].append((xv, final_err - get_baseline(h, p)))

            for run in ale_runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                final_err = data[-1, mcol]
                xv = get_var(run, args.xaxis) + x_off
                cv = get_var(run, args.color)
                groups3[cv].append((xv, final_err - get_baseline(h, p)))

            for cv in sorted(groups.keys()):
                pts = sorted(groups[cv])
                xs, ys = zip(*pts)
                print(ys)
                ax.plot(xs, ys, marker='o', color=color_map[cv],
                        label=f'{color_key} = {c_disp(cv)} (overall)')

            for cv in sorted(groups2.keys()):
                pts = sorted(groups2[cv])
                xs, ys = zip(*pts)
                ax.plot(xs, ys, marker='s', color=color_map2[cv],
                        label=f'{color_key} = {c_disp(cv)} (Proj sol. transfer)')

            for cv in sorted(groups3.keys()):
                pts = sorted(groups3[cv])
                xs, ys = zip(*pts)
                ax.plot(xs, ys, marker='^', color=color_map3[cv],
                        label=f'{color_key} = {c_disp(cv)} (ALE sol. transfer)')

        else:
            # Key by (method, cv) so KM and ALE are separate lines
            groups = defaultdict(list)

            for run in runs:
                n, h, p, method, f = run
                data = np.loadtxt(f, skiprows=1)
                final_err = data[-1, mcol]
                if args.sub_n1:
                    final_err -= get_baseline(h, p)
                xv = get_var(run, args.xaxis) + x_off
                cv = get_var(run, args.color)
                groups[(method, cv)].append((xv, final_err))

            for (method, cv) in sorted(groups.keys(), key=lambda k: (_METHOD_ORDER.get(k[0], 99), k[1])):
                pts = sorted(groups[(method, cv)])
                xs, ys = zip(*pts)
                print(ys)
                mlabel = _METHOD_LABEL.get(method, method)
                cmap_m = _METHOD_CMAP.get(method, color_map)
                ax.plot(xs, ys,
                        marker=_METHOD_MARKER.get(method, 'o'),
                        linestyle=_METHOD_LS.get(method, '-'),
                        color=cmap_m[cv],
                        label=f'{color_key} = {c_disp(cv)} ({mlabel})')

            if args.sub_n1 and args.xaxis == 'n':
                # Draw a horizontal orange reference line at the n=1 baseline error
                # for each unique (h, p) combination, labelled by polynomial order.
                seen_p = {}  # p -> baseline value (use first h encountered per p)
                for n, h, p, method, f in runs:
                    if p not in seen_p:
                        seen_p[p] = get_baseline(h, p)
                for p in sorted(seen_p):
                    bval = seen_p[p]
                    col = color_map[p] if args.color == 'p' else 'orange'
                    ref_label = f'cycles=0 error (p={p})' if args.rename_n else f'n=1 error (p={p})'
                    ax.axhline(bval, color=col, linestyle='--', linewidth=1.5,
                               label=ref_label, zorder=1)

        if args.logy:
            ax.set_yscale('log')
        if args.logx:
            ax.set_xscale('log')

    if combine:
        base_label = METRIC_LABELS[args.metric]
        ylabel = f'Overall & Solution Transfer {base_label}'
    else:
        base_label = METRIC_LABELS[args.metric]
        ylabel = f'Solution Transfer {base_label}' if args.sub_n1 else f'Overall {base_label}'

    ax.set_xlabel(AXIS_LABELS[args.xaxis])
    ax.set_ylabel(ylabel)
    fixed_parts = []
    if args.fix_n is not None:
        fixed_parts.append(f'cycles={args.fix_n - 1}' if args.rename_n else f'n={args.fix_n}')
    if args.fix_h is not None: fixed_parts.append(f'h={args.fix_h}')
    if args.fix_p is not None: fixed_parts.append(f'p={args.fix_p}')
    fixed_str = ', '.join(fixed_parts)
    ax.set_title(f'{ylabel} vs {AXIS_LABELS[args.xaxis].lower()}  ({fixed_str})')
    # ax.set_title(f'{ylabel} vs {AXIS_LABELS[args.xaxis].lower()}')
    n_legend_items = len(ax.get_legend_handles_labels()[0])
    ncol = max(1, n_legend_items // 4)
    if combine:
        tmp_legend = ax.legend(ncol=ncol, fontsize='small', framealpha=0.5, loc='best')#,
                #   loc='upper right', bbox_to_anchor=(0.99, 0.805))
    else:
        tmp_legend = ax.legend(ncol=ncol, fontsize='small', framealpha=0.5, loc='best')
    best_loc = tmp_legend._loc
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    plt.tight_layout()

    if args.subslides:
        handles, labels = ax.get_legend_handles_labels()

        # Group handles by series name (text in trailing parentheses), preserving order
        series_re = re.compile(r'\((.+)\)$')
        series_order = []
        series_groups = {}
        extra = []  # handles without a series group, e.g. axhline reference lines
        for h, l in zip(handles, labels):
            m = series_re.search(l)
            if m:
                key = m.group(1)
                if key not in series_groups:
                    series_groups[key] = []
                    series_order.append(key)
                series_groups[key].append((h, l))
            else:
                extra.append((h, l))

        base = args.save if args.save else 'plot'
        stem, ext = os.path.splitext(base)
        ext = ext or '.png'

        for i in range(1, len(series_order) + 1):
            visible = set(series_order[:i])
            for h, l in zip(handles, labels):
                m = series_re.search(l)
                h.set_visible(not m or m.group(1) in visible)

            vis_h = [h for key in series_order[:i] for h, _ in series_groups[key]]
            vis_l = [l for key in series_order[:i] for _, l in series_groups[key]]
            for h, l in extra:
                vis_h.append(h)
                vis_l.append(l)

            slide_ncol = max(1, len(vis_h) // 4)
            ax.get_legend().remove()
            ax.legend(vis_h, vis_l, ncol=slide_ncol, fontsize='small',
                      framealpha=0.5, loc=best_loc)

            out = f'{stem}_slide{i}{ext}'
            plt.savefig(out, dpi=300)
            print(f'Saved {out}')

        for h in handles:
            h.set_visible(True)
    elif args.save:
        plt.savefig(args.save, dpi=300)
        print(f'Saved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()

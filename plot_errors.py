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
FILE_RE = re.compile(r'ErrorFile_n_(\d+)_h_(.+)_p_(\d+)\.err$')

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

# ── Combined-plot toggle ─────────────────────────────────────────────────────
# Set True to overlay overall L2 and reprojection L2 (--sub-n1) on one axes
# with distinct hue progressions.  The --sub-n1 flag is ignored in this mode.
COMBINE_BOTH = True

# Colormaps: primary = overall L2, secondary = reprojection L2 (COMBINE_BOTH)
_CMAPS  = {'n': 'Blues',   'h': 'Greens',  'p': 'Oranges'}
_CMAPS2 = {'n': 'Reds',    'h': 'Purples', 'p': 'Blues'}
# ─────────────────────────────────────────────────────────────────────────────


def parse_filename(fname):
    m = FILE_RE.match(os.path.basename(fname))
    if not m:
        return None
    n = int(m.group(1))
    h = float(m.group(2).replace('dot', '.'))
    p = int(m.group(3))
    return n, h, p


def get_var(run, key):
    n, h, p, _ = run
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
    parser.add_argument('--logy', action='store_true',
                        help='Logarithmic y axis')
    parser.add_argument('--save', metavar='FILE',
                        help='Save figure to FILE instead of displaying it')
    args = parser.parse_args()

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
        n, h, p = params
        if args.fix_n is not None and n != args.fix_n:
            continue
        if args.fix_h is not None and not np.isclose(h, args.fix_h):
            continue
        if args.fix_p is not None and p != args.fix_p:
            continue
        runs.append((n, h, p, f))

    # In single-mode with --sub-n1, drop n=1 runs (they are the baseline).
    # In COMBINE_BOTH mode the same filtering applies for the reprojection
    # series, but we keep runs as-is here and handle it in the plot loop.
    if args.sub_n1 and not COMBINE_BOTH:
        runs = [(n, h, p, f) for n, h, p, f in runs if n != 1]

    if not runs:
        sys.exit('No runs matched the given filters.')

    # Build n=1 baseline lookup: (h, p) -> final metric value
    baseline = {}
    if args.sub_n1 or COMBINE_BOTH:
        for f in all_files:
            params = parse_filename(f)
            if params is None:
                continue
            n, h, p = params
            if n == 1:
                data = np.loadtxt(f, skiprows=1)
                baseline[(h, p)] = data[-1, mcol]

    def get_baseline(h, p):
        key = next((k for k in baseline if np.isclose(k[0], h) and k[1] == p), None)
        if key is None:
            sys.exit(f'--sub-n1: no n=1 baseline found for h={h}, p={p}')
        return baseline[key]

    # Build colour maps
    color_vals = sorted(set(get_var(run, args.color) for run in runs))
    n_cols = len(color_vals)
    # Sample 0.35–0.9 to avoid near-white at the light end
    _samples = [0.35 + 0.55 * i / max(n_cols - 1, 1) for i in range(n_cols)]

    cmap = plt.get_cmap(_CMAPS[args.color])
    color_map = {v: cmap(s) for v, s in zip(color_vals, _samples)}

    if COMBINE_BOTH:
        cmap2 = plt.get_cmap(_CMAPS2[args.color])
        color_map2 = {v: cmap2(s) for v, s in zip(color_vals, _samples)}

    fig, ax = plt.subplots(figsize=(8, 5))
    seen_labels = set()

    print(f"{len(runs)} runs:")
    for run in runs:
        print(run)

    if args.xaxis == 'time':
        for run in runs:
            n, h, p, f = run
            data = np.loadtxt(f, skiprows=1)
            cv = get_var(run, args.color)

            if COMBINE_BOTH:
                # Overall series
                label = f'{args.color} = {cv} (overall)'
                ax.plot(data[:, 0], data[:, mcol],
                        color=color_map[cv],
                        label=label if label not in seen_labels else '_nolegend_')
                seen_labels.add(label)
                # Reprojection series
                label2 = f'{args.color} = {cv} (reprojection)'
                ax.plot(data[:, 0], data[:, mcol] - get_baseline(h, p),
                        color=color_map2[cv],
                        label=label2 if label2 not in seen_labels else '_nolegend_')
                seen_labels.add(label2)
            else:
                yvals = data[:, mcol]
                if args.sub_n1:
                    yvals = yvals - get_baseline(h, p)
                label = f'{args.color} = {cv}'
                ax.plot(data[:, 0], yvals,
                        color=color_map[cv],
                        label=label if label not in seen_labels else '_nolegend_')
                seen_labels.add(label)

    else:
        # One point per run: final-timestep error (last row = post-adaptation value)
        groups = defaultdict(list)
        groups2 = defaultdict(list)  # reprojection series for COMBINE_BOTH

        for run in runs:
            n, h, p, f = run
            data = np.loadtxt(f, skiprows=1)
            final_err = data[-1, mcol]
            xv = get_var(run, args.xaxis)
            cv = get_var(run, args.color)

            if COMBINE_BOTH:
                groups[cv].append((xv, final_err))
                groups2[cv].append((xv, final_err - get_baseline(h, p)))
            else:
                if args.sub_n1:
                    final_err -= get_baseline(h, p)
                groups[cv].append((xv, final_err))

        for cv in sorted(groups.keys()):
            pts = sorted(groups[cv])
            xs, ys = zip(*pts)
            print(ys)
            label_suffix = ' (overall)' if COMBINE_BOTH else ''
            ax.plot(xs, ys, marker='o', color=color_map[cv],
                    label=f'{args.color} = {cv}{label_suffix}')

        if COMBINE_BOTH:
            for cv in sorted(groups2.keys()):
                pts = sorted(groups2[cv])
                xs, ys = zip(*pts)
                ax.plot(xs, ys, marker='s', color=color_map2[cv],
                        label=f'{args.color} = {cv} (reprojection)')

        if args.logy:
            ax.set_yscale('log')

    if COMBINE_BOTH:
        base_label = METRIC_LABELS[args.metric]
        ylabel = f'Overall & Reprojection {base_label}'
    else:
        base_label = METRIC_LABELS[args.metric]
        ylabel = f'Reprojection {base_label}' if args.sub_n1 else f'Overall {base_label}'

    ax.set_xlabel(AXIS_LABELS[args.xaxis])
    ax.set_ylabel(ylabel)
    fixed_parts = []
    if args.fix_n is not None: fixed_parts.append(f'n={args.fix_n}')
    if args.fix_h is not None: fixed_parts.append(f'h={args.fix_h}')
    if args.fix_p is not None: fixed_parts.append(f'p={args.fix_p}')
    fixed_str = ',  '.join(fixed_parts) if fixed_parts else 'none fixed'
    # ax.set_title(f'{ylabel} vs {AXIS_LABELS[args.xaxis].lower()}  ({fixed_str})')
    ax.set_title(f'{ylabel} vs {AXIS_LABELS[args.xaxis].lower()}')
    ax.legend(ncol=2 if COMBINE_BOTH else 1, fontsize='small', framealpha=0.5,
              loc='upper right', bbox_to_anchor=(0.99, 0.805) if COMBINE_BOTH else (0.01, 0.99))
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    plt.tight_layout()

    if args.save:
        plt.savefig(args.save, dpi=300)
        print(f'Saved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()

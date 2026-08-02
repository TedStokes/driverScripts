#!/usr/bin/env python3
"""
Heatmap of final u error over the (minimum layer height, advection velocity) plane.

Runs are produced by stability_steps.sh, which writes an .err file plus a .status
sidecar holding the solver exit code. A non-zero exit code means the solver hit
"NaN found during time integration" and aborted; those cells are drawn in a flat
critical colour and labelled "NaN" in the cell.

Instability does not always reach NaN within NumSteps - a diverging run often
completes cleanly with a final L2 error of 1e130 or so. Those keep their value and
simply clip to the top of the colour scale, which --vmax fixes at the O(1) scale of
the exact solution by default. A clipped cell therefore reads as "at or past total
loss of the solution".

Examples:
  # No adaptation, mesh built directly at the deformed final position
  python plot_stability.py --fix-n 1  --fix-p 5 --method Projection \
      --title "No adaptation (deformed mesh)" --save heat_noadapt.png

  # Projection solution transfer over 16 adaptive runs
  python plot_stability.py --fix-n 16 --fix-p 5 --method Projection \
      --title "Projection solution transfer" --save heat_proj.png

  # ALE solution transfer over 16 adaptive runs
  python plot_stability.py --fix-n 16 --fix-p 5 --method ALE \
      --title "ALE solution transfer" --save heat_ale.png

Pass a shared --vmin/--vmax to all three so the colour scales are comparable;
the script prints the data range it found to help you choose them.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

METRIC_COL = {'u_L2': 1, 'u_Linf': 2, 'u_H1': 3}
FILE_RE = re.compile(r'ErrorFile_v_(.+?)_h_(.+?)_n_(\d+)_p_(\d+)_m_(.+)\.err$')
LOG_RE = re.compile(r'log_v(.+?)_h(.+?)_n(\d+)_p(\d+)_m(.+)\.txt$')
# Printed once per adaptive cycle by DriverAdaptBL in ALE mode.
GRIDVEL_RE = re.compile(r'GridVel n=(\d+) vy_min=(\S+) vy_max=(\S+)')

METRIC_LABELS = {
    'u_L2':   'L2 error',
    'u_Linf': 'Linf error',
    'u_H1':   'H1 error',
}

EXPLODED_COLOR = '#8b1a1a'   # critical: solver aborted on NaN
MISSING_COLOR  = '#d9d9d9'   # neutral: no run at this grid point

# Default colour scale: decades from 1e-2 (converged) to 1e0, the O(1) scale of
# the exact solution. Anything larger clips to the top of the ramp.
DEFAULT_VMIN = 1e-2
DEFAULT_VMAX = 1e0

# Cell geometry in inches. The heatmap axes is placed at exactly
# (ncols*CELL_W) x (nrows*CELL_W*CELL_ASPECT) inches via a fixed Divider, so
# cell size and shape are identical on every figure regardless of what else is
# drawn - the grid-velocity top axis, longer tick labels, or different fonts.
# Nothing is left to tight_layout to negotiate.
CELL_W = 0.589
CELL_ASPECT = 0.856

# Margins reserved around the axes. Only need to be large enough that labels are
# not clipped before bbox_inches='tight' crops the figure down.
MARGIN_L, MARGIN_B = 1.0, 0.9
MARGIN_R, MARGIN_T = 1.6, 1.2

# Colourbar geometry, fixed so it is identical on every figure: an absolute
# width in inches, and a gap given as a fraction of the heatmap width.
CBAR_WIDTH_INCHES = 0.22
CBAR_PAD_FRAC = 0.025


def unsanitise(s):
    """Inverse of stability_steps.sh's sanitise(): 0dot05 -> 0.05, m2dot5 -> -2.5."""
    if s.startswith('m'):
        return -float(s[1:].replace('dot', '.'))
    return float(s.replace('dot', '.'))


def read_gridvel(logfile):
    """Most negative wall-normal grid velocity over all adaptive cycles.

    The boundary layers collapse towards the wall, so vy is negative and the
    peak mesh speed is the minimum. Returns None for a log with no GridVel
    lines, which is every non-ALE run.
    """
    try:
        with open(logfile) as f:
            vals = [float(m.group(2))
                    for m in (GRIDVEL_RE.search(line) for line in f) if m]
    except OSError:
        return None
    return min(vals) if vals else None


def parse_filename(fname):
    m = FILE_RE.match(os.path.basename(fname))
    if not m:
        return None
    advy = unsanitise(m.group(1))
    h    = unsanitise(m.group(2))
    n    = int(m.group(3))
    p    = int(m.group(4))
    method = m.group(5)
    return advy, h, n, p, method


def read_run(errfile, mcol, stat='peak'):
    """Return (value, exploded).

    stat selects the peak error over the whole run (the default, which captures
    the damage done by an adaptive step rather than how well it recovered) or
    the final-timestep error.

    exploded is True when the solver aborted on NaN, in which case there is no
    usable value. A run that diverged without reaching NaN is not special cased
    - it keeps its (huge) value and clips to the top of the colour scale.
    """
    status_file = errfile[:-len('.err')] + '.status'
    if os.path.exists(status_file):
        with open(status_file) as f:
            try:
                rc = int(f.read().strip())
            except ValueError:
                rc = 1
        if rc != 0:
            return None, True

    try:
        data = np.loadtxt(errfile, skiprows=1, ndmin=2)
    except (OSError, ValueError):
        return None, True
    if data.size == 0:
        return None, True

    val = data[:, mcol].max() if stat == 'peak' else data[-1, mcol]
    if not np.isfinite(val):
        return None, True
    return val, False


def main():
    parser = argparse.ArgumentParser(
        description='Heatmap of final error vs advection velocity and minimum layer height.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--fix-n', type=int, required=True, metavar='N',
                        help='NumRuns value to plot')
    parser.add_argument('--fix-p', type=int, required=True, metavar='P',
                        help='NUMMODES value to plot')
    parser.add_argument('--method', default='Projection',
                        help='AdaptBL_transfer method to plot (default: Projection). '
                             'Irrelevant to the solver when n=1, but still part of the filename.')
    parser.add_argument('--results-dir', default='results_stability', metavar='DIR',
                        help='Directory containing .err/.status files (default: results_stability)')
    parser.add_argument('--metric', choices=list(METRIC_COL), default='u_L2',
                        help='Error metric to plot (default: u_L2)')
    parser.add_argument('--vmin', type=float, default=DEFAULT_VMIN,
                        help=f'Colour scale minimum (default: {DEFAULT_VMIN:g})')
    parser.add_argument('--vmax', type=float, default=DEFAULT_VMAX,
                        help=f'Colour scale maximum (default: {DEFAULT_VMAX:g})')
    parser.add_argument('--stat', choices=['peak', 'final'], default='peak',
                        help='Peak error over the run (default) or final-timestep error. '
                             'Peak measures the damage an adaptive step does; final also '
                             'reflects how far the solution recovered afterwards.')
    parser.add_argument('--cmap', default='viridis',
                        help='Named matplotlib colormap (default: viridis)')
    parser.add_argument('--title', default=None,
                        help='Figure title (default: describes the fixed parameters)')
    parser.add_argument('--annotate', action='store_true',
                        help='Print the error value inside each cell')
    parser.add_argument('--paper', action='store_true',
                        help='Paper mode: smaller figure so elements scale up when embedded')
    parser.add_argument('--save', metavar='FILE',
                        help='Save figure to FILE instead of displaying it')
    args = parser.parse_args()

    mcol = METRIC_COL[args.metric]

    pattern = os.path.join(args.results_dir, 'ErrorFile_v_*_h_*_n_*_p_*_m_*.err')
    all_files = sorted(glob.glob(pattern))
    if not all_files:
        sys.exit(f'No error files found matching {pattern}')

    # cell[(advy, h)] -> (value, exploded)
    cells = {}
    for f in all_files:
        params = parse_filename(f)
        if params is None:
            continue
        advy, h, n, p, method = params
        if n != args.fix_n or p != args.fix_p or method != args.method:
            continue
        cells[(advy, h)] = read_run(f, mcol, args.stat)

    if not cells:
        sys.exit(f'No runs matched n={args.fix_n}, p={args.fix_p}, method={args.method}')

    h_vals    = sorted({h for _, h in cells})
    advy_vals = sorted({v for v, _ in cells})

    grid     = np.full((len(advy_vals), len(h_vals)), np.nan)
    exploded = np.zeros_like(grid, dtype=bool)
    present  = np.zeros_like(grid, dtype=bool)

    for (advy, h), (val, boom) in cells.items():
        i, j = advy_vals.index(advy), h_vals.index(h)
        present[i, j] = True
        exploded[i, j] = boom
        if val is not None:
            grid[i, j] = val

    # Peak wall-normal mesh speed per h column. It is set purely by the r-ramp
    # geometry and so does not depend on advy - but a run that aborted early
    # logged fewer cycles, so take the extreme over every run in the column.
    gridvel = {}
    for logfile in glob.glob(os.path.join(args.results_dir, 'log_v*.txt')):
        m = LOG_RE.match(os.path.basename(logfile))
        if not m:
            continue
        if (int(m.group(3)) != args.fix_n or int(m.group(4)) != args.fix_p
                or m.group(5) != args.method):
            continue
        v = read_gridvel(logfile)
        if v is None:
            continue
        h = unsanitise(m.group(2))
        gridvel[h] = min(gridvel.get(h, v), v)

    n_ok = int(np.isfinite(grid).sum())
    n_boom = int(exploded.sum())
    n_gap = int((~present).sum())
    finite = grid[np.isfinite(grid)]
    print(f'{len(cells)} runs: {n_ok} completed, {n_boom} exploded (NaN abort), '
          f'{n_gap} grid point(s) with no run')
    if finite.size:
        print(f'{args.metric} range: {finite.min():.6g} to {finite.max():.6g}')
    else:
        sys.exit('Every matching run exploded - nothing to colour scale against.')

    vmin, vmax = args.vmin, args.vmax
    if vmin <= 0:
        sys.exit(f'--vmin must be positive for a log colour scale (got {vmin})')
    if vmax <= vmin:
        sys.exit(f'--vmax ({vmax}) must exceed --vmin ({vmin})')

    # Fixed axes box: ncols x nrows cells of an exact size in inches.
    cell_w = CELL_W * (0.72 if args.paper else 1.0)
    ax_w = len(h_vals) * cell_w
    ax_h = len(advy_vals) * cell_w * CELL_ASPECT
    fig = plt.figure(figsize=(MARGIN_L + ax_w + MARGIN_R,
                              MARGIN_B + ax_h + MARGIN_T))
    divider = Divider(fig, (0, 0, 1, 1),
                      [Size.Fixed(MARGIN_L), Size.Fixed(ax_w)],
                      [Size.Fixed(MARGIN_B), Size.Fixed(ax_h)], aspect=False)
    ax = fig.add_axes(divider.get_position(),
                      axes_locator=divider.new_locator(nx=1, ny=1))

    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad(MISSING_COLOR)

    xs = np.arange(len(h_vals))
    ys = np.arange(len(advy_vals))
    mesh = ax.pcolormesh(xs, ys, np.ma.masked_invalid(grid),
                         cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                         shading='nearest', edgecolors='white', linewidth=0.5)

    # Exploded cells get the reserved critical colour, drawn over the masked
    # grid, and label themselves so no legend is needed to read them.
    for i, j in zip(*np.nonzero(exploded)):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                               facecolor=EXPLODED_COLOR, edgecolor='white',
                               linewidth=0.5, zorder=2))
        ax.text(j, i, 'NaN', ha='center', va='center', color='white',
                fontsize='x-small', zorder=3)

    if args.annotate:
        for i in range(len(advy_vals)):
            for j in range(len(h_vals)):
                if exploded[i, j] or not np.isfinite(grid[i, j]):
                    continue
                # Contrast against whatever the colormap actually put in the cell.
                r, g, b, _ = cmap(LogNorm(vmin=vmin, vmax=vmax)(grid[i, j]))
                lum = 0.299 * r + 0.587 * g + 0.114 * b
                ax.text(j, i, f'{grid[i, j]:.2g}', ha='center', va='center',
                        fontsize='x-small', zorder=3,
                        color='#222222' if lum > 0.55 else 'white')

    ax.set_xticks(xs)
    ax.set_xticklabels([f'{h:g}' for h in h_vals])
    ax.set_yticks(ys)
    ax.set_yticklabels([f'{v:g}' for v in advy_vals])
    ax.set_xlabel('Minimum layer height (h)')
    ax.set_ylabel('Advection velocity (advy)')

    # ALE runs log their mesh motion, so label each column with it. Columns are
    # ordinal, so this is a second row of labels rather than a rescaled axis.
    if gridvel:
        secax = ax.secondary_xaxis('top')
        secax.set_xticks(xs)
        secax.set_xticklabels(
            [f'{gridvel[h]:.3g}' if h in gridvel else '-' for h in h_vals])
        secax.set_xlabel('Peak grid velocity ($v_y$, minimum over cycles)')
        missing = [h for h in h_vals if h not in gridvel]
        if missing:
            print(f'no grid velocity logged for h = '
                  f'{", ".join(f"{h:g}" for h in missing)}')

    # Clear the secondary axis and its label when one is present.
    title_pad = 32 if gridvel else None
    if args.title is not None:
        ax.set_title(args.title, pad=title_pad)
    else:
        ax.set_title(f'{METRIC_LABELS[args.metric]}  '
                     f'(n={args.fix_n}, p={args.fix_p}, {args.method})',
                     pad=title_pad)

    # Values can sit outside the fixed scale, so flag which ends are clipped.
    below, above = bool((finite < vmin).any()), bool((finite > vmax).any())
    extend = ('both' if below and above else
              'min' if below else 'max' if above else 'neither')

    # Anchor the colourbar to the axes as actually drawn, at a fixed width in
    # inches. fig.colorbar(ax=ax) and axes_grid1 dividers both size themselves
    # from the axes' pre-aspect slot, which differs between the taller ALE
    # figure and the others; an inset anchored in ax coordinates tracks the real
    # box, so the bar is identical in size and position on every plot.
    cax = inset_axes(ax, width=CBAR_WIDTH_INCHES, height='100%',
                     loc='lower left', bbox_to_anchor=(1.0 + CBAR_PAD_FRAC, 0., 1., 1.),
                     bbox_transform=ax.transAxes, borderpad=0)
    cbar = fig.colorbar(mesh, cax=cax, extend=extend)
    cbar.set_label(f'{args.stat.capitalize()} {METRIC_LABELS[args.metric]}')

    ax.set_xlim(-0.5, len(h_vals) - 0.5)
    ax.set_ylim(-0.5, len(advy_vals) - 0.5)

    if args.save:
        plt.savefig(args.save, dpi=300, bbox_inches='tight')
        print(f'Saved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()

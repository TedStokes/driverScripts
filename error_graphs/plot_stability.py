#!/usr/bin/env python3
"""
Heatmap of final u error against advection velocity and one of the two sweep
parameters: minimum layer height (fix n, --fix-n) or number of adaptive runs
(fix h, --fix-h). Exactly one of the two must be fixed; the other becomes the
x-axis.

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

  # Number of adaptive runs along the x-axis, at a fixed layer height
  python plot_stability.py --fix-h 0.1 --fix-p 5 --method ALE \
      --title "ALE solution transfer, h=0.1" --save heat_h0dot1_ale.png

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
# Printed every IO_CFLSteps by UnsteadyAdvection. Normalised so that CFL = 1 is
# the predicted stability limit; the startup summary line has no bare value
# after the colon, so the number is what distinguishes it.
CFL_RE = re.compile(r'CFL_GLL:\s+([0-9eE.+-]+)\s+\(in elmt')

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

# For --stat growth the metric column is an amplification factor, not an error:
# 1 is "unchanged", so the scale straddles it instead.
GROWTH_VMIN = 1e-1
GROWTH_VMAX = 1e3

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


def read_cfl(logfile):
    """Largest CFL_GLL the solver reported over the whole run.

    This is the predicted stability parameter: blow-up is expected above 1. A
    run that aborted early logged fewer steps, so this is the peak over what it
    managed, not over the schedule it was asked for.
    """
    try:
        with open(logfile) as f:
            vals = [float(m.group(1))
                    for m in (CFL_RE.search(line) for line in f) if m]
    except OSError:
        return None
    return max(vals) if vals else None


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
    the damage done by an adaptive step rather than how well it recovered), the
    final-timestep error, or 'growth' - the ratio of the last value to the
    first. Under -H (zero exact solution, noise initial condition) the metric
    column is ||u_h|| itself, so growth is the amplification the scheme applied
    over the run and is the quantity the CFL prediction is about.

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

    if stat == 'peak':
        val = data[:, mcol].max()
    elif stat == 'growth':
        # Duplicate rows share a timestamp (the solver writes one per solve and
        # one per transfer), which does not matter for a first/last ratio.
        if data[0, mcol] == 0.0:
            return None, True
        val = data[-1, mcol] / data[0, mcol]
    else:
        val = data[-1, mcol]
    if not np.isfinite(val):
        return None, True
    return val, False


def main():
    parser = argparse.ArgumentParser(
        description='Heatmap of final error vs advection velocity and either '
                    'minimum layer height (--fix-n) or number of adaptive runs (--fix-h).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--fix-n', type=int, default=None, metavar='N',
                        help='NumRuns value to plot; layer height goes on the x-axis')
    parser.add_argument('--fix-h', type=float, default=None, metavar='H',
                        help='AdaptBL_h_init value to plot; NumRuns goes on the x-axis')
    parser.add_argument('--fix-p', type=int, required=True, metavar='P',
                        help='NUMMODES value to plot')
    parser.add_argument('--method', default='Projection',
                        help='AdaptBL_transfer method to plot (default: Projection). '
                             'Irrelevant to the solver when n=1, but still part of the filename.')
    parser.add_argument('--x-values', default=None, metavar='LIST',
                        help='Comma-separated subset of x-axis values to plot, e.g. '
                             '"1,2,3,4,5". Default: every value found on disk.')
    parser.add_argument('--advy-range', default=None, metavar='MIN,MAX',
                        help='Restrict the y-axis to advection velocities in [MIN, MAX], '
                             'e.g. "=-1.5,1.5". Default: every value found on disk. '
                             'Note the "=" - argparse reads a leading minus as a flag.')
    parser.add_argument('--advy-values', default=None, metavar='LIST',
                        help='Comma-separated subset of advection velocities to plot. '
                             'Combines with --advy-range; a row must satisfy both.')
    parser.add_argument('--results-dir', default='results_stability', metavar='DIR',
                        help='Directory containing .err/.status files (default: results_stability)')
    parser.add_argument('--metric', choices=list(METRIC_COL), default='u_L2',
                        help='Error metric to plot (default: u_L2)')
    parser.add_argument('--vmin', type=float, default=None,
                        help=f'Colour scale minimum (default: {DEFAULT_VMIN:g}, '
                             f'or {GROWTH_VMIN:g} with --stat growth)')
    parser.add_argument('--vmax', type=float, default=None,
                        help=f'Colour scale maximum (default: {DEFAULT_VMAX:g}, '
                             f'or {GROWTH_VMAX:g} with --stat growth)')
    parser.add_argument('--stat', choices=['peak', 'final', 'growth'], default='peak',
                        help='Peak error over the run (default), final-timestep error, or '
                             'growth (last/first). Peak measures the damage an adaptive '
                             'step does; final also reflects how far the solution '
                             'recovered afterwards; growth is the amplification factor, '
                             'and is what the CFL prediction refers to - use it with the '
                             'homogeneous runs from stability_steps.sh -H.')
    parser.add_argument('--cmap', default='viridis',
                        help='Named matplotlib colormap (default: viridis)')
    parser.add_argument('--title', default=None,
                        help='Figure title (default: describes the fixed parameters)')
    parser.add_argument('--annotate', nargs='?', const='value', default=None,
                        choices=['value', 'cfl', 'both'],
                        help='Print a number inside each cell: the plotted value '
                             '(default), the predicted CFL_GLL, or both.')
    parser.add_argument('--cfl-contour', action='store_true',
                        help='Draw the CFL_GLL = 1 contour, i.e. the predicted '
                             'stability boundary, over the heatmap.')
    parser.add_argument('--paper', action='store_true',
                        help='Paper mode: smaller figure so elements scale up when embedded')
    parser.add_argument('--save', metavar='FILE',
                        help='Save figure to FILE instead of displaying it')
    args = parser.parse_args()

    if (args.fix_n is None) == (args.fix_h is None):
        sys.exit('Give exactly one of --fix-n (layer height on the x-axis) '
                 'or --fix-h (number of adaptive runs on the x-axis)')
    # Whichever of h/n is not pinned becomes the x-axis.
    x_is_n = args.fix_n is None

    x_filter = None
    if args.x_values is not None:
        try:
            x_filter = [int(s) if x_is_n else float(s)
                        for s in args.x_values.split(',') if s.strip()]
        except ValueError:
            sys.exit(f'--x-values must be a comma-separated list of '
                     f'{"integers" if x_is_n else "numbers"} (got {args.x_values!r})')

    def keep_x(x):
        return x_filter is None or any(np.isclose(x, t) for t in x_filter)

    advy_lo, advy_hi = -np.inf, np.inf
    if args.advy_range is not None:
        try:
            advy_lo, advy_hi = (float(s) for s in args.advy_range.split(','))
        except ValueError:
            sys.exit(f'--advy-range must be MIN,MAX (got {args.advy_range!r})')
        if advy_hi < advy_lo:
            sys.exit(f'--advy-range MAX ({advy_hi}) must not be below MIN ({advy_lo})')

    advy_filter = None
    if args.advy_values is not None:
        try:
            advy_filter = [float(s) for s in args.advy_values.split(',') if s.strip()]
        except ValueError:
            sys.exit(f'--advy-values must be a comma-separated list of numbers '
                     f'(got {args.advy_values!r})')

    def keep_advy(v):
        if not advy_lo - 1e-9 <= v <= advy_hi + 1e-9:
            return False
        return advy_filter is None or any(np.isclose(v, t) for t in advy_filter)

    mcol = METRIC_COL[args.metric]

    growth = args.stat == 'growth'
    if args.vmin is None:
        args.vmin = GROWTH_VMIN if growth else DEFAULT_VMIN
    if args.vmax is None:
        args.vmax = GROWTH_VMAX if growth else DEFAULT_VMAX
    # "Peak L2 error" reads fine; "Growth L2 error" does not.
    quantity = (f'{METRIC_LABELS[args.metric]} growth' if growth
                else f'{args.stat.capitalize()} {METRIC_LABELS[args.metric]}')

    pattern = os.path.join(args.results_dir, 'ErrorFile_v_*_h_*_n_*_p_*_m_*.err')
    all_files = sorted(glob.glob(pattern))
    if not all_files:
        sys.exit(f'No error files found matching {pattern}')

    # cell[(advy, x)] -> (value, exploded), where x is h or n per x_is_n.
    cells = {}
    for f in all_files:
        params = parse_filename(f)
        if params is None:
            continue
        advy, h, n, p, method = params
        if p != args.fix_p or method != args.method or not keep_advy(advy):
            continue
        if x_is_n:
            if not np.isclose(h, args.fix_h):
                continue
            x = n
        else:
            if n != args.fix_n:
                continue
            x = h
        if not keep_x(x):
            continue
        cells[(advy, x)] = read_run(f, mcol, args.stat)

    fixed_desc = (f'h={args.fix_h:g}' if x_is_n else f'n={args.fix_n}')
    if not cells:
        sys.exit(f'No runs matched {fixed_desc}, p={args.fix_p}, method={args.method}')

    x_vals    = sorted({x for _, x in cells})
    advy_vals = sorted({v for v, _ in cells})

    grid     = np.full((len(advy_vals), len(x_vals)), np.nan)
    exploded = np.zeros_like(grid, dtype=bool)
    present  = np.zeros_like(grid, dtype=bool)

    for (advy, x), (val, boom) in cells.items():
        i, j = advy_vals.index(advy), x_vals.index(x)
        present[i, j] = True
        exploded[i, j] = boom
        if val is not None:
            grid[i, j] = val

    # Peak wall-normal mesh speed per column. It is set purely by the r-ramp
    # geometry and so does not depend on advy - but a run that aborted early
    # logged fewer cycles, so take the extreme over every run in the column.
    gridvel = {}
    cfl_cells = {}
    for logfile in glob.glob(os.path.join(args.results_dir, 'log_v*.txt')):
        m = LOG_RE.match(os.path.basename(logfile))
        if not m:
            continue
        if int(m.group(4)) != args.fix_p or m.group(5) != args.method:
            continue
        h, n = unsanitise(m.group(2)), int(m.group(3))
        if x_is_n:
            if not np.isclose(h, args.fix_h):
                continue
            x = n
        else:
            if n != args.fix_n:
                continue
            x = h
        if not keep_x(x):
            continue
        advy = unsanitise(m.group(1))
        if keep_advy(advy):
            c = read_cfl(logfile)
            if c is not None:
                cfl_cells[(advy, x)] = c

        v = read_gridvel(logfile)
        if v is None:
            continue
        gridvel[x] = min(gridvel.get(x, v), v)

    # Unlike the grid velocity, CFL depends on advy as well, so it is a full
    # grid rather than one value per column.
    cflgrid = np.full_like(grid, np.nan)
    for (advy, x), c in cfl_cells.items():
        if advy in advy_vals and x in x_vals:
            cflgrid[advy_vals.index(advy), x_vals.index(x)] = c

    n_ok = int(np.isfinite(grid).sum())
    n_boom = int(exploded.sum())
    n_gap = int((~present).sum())
    finite = grid[np.isfinite(grid)]
    print(f'{len(cells)} runs: {n_ok} completed, {n_boom} exploded (NaN abort), '
          f'{n_gap} grid point(s) with no run')
    n_cfl = int(np.isfinite(cflgrid).sum())
    print(f'CFL_GLL logged for {n_cfl} of {present.sum()} cell(s)')
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
    ax_w = len(x_vals) * cell_w
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

    xs = np.arange(len(x_vals))
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
            for j in range(len(x_vals)):
                if not present[i, j]:
                    continue
                parts = []
                if args.annotate in ('value', 'both') and np.isfinite(grid[i, j]):
                    parts.append(f'{grid[i, j]:.2g}')
                if args.annotate in ('cfl', 'both') and np.isfinite(cflgrid[i, j]):
                    parts.append(f'CFL {cflgrid[i, j]:.2g}')
                if not parts:
                    continue

                # Contrast against whatever was actually drawn in the cell:
                # the reserved colour for an exploded run, the colormap
                # otherwise.
                if exploded[i, j] or not np.isfinite(grid[i, j]):
                    lum = 0.0
                    dy = 0.22 if args.annotate == 'both' else 0.0
                else:
                    r, g, b, _ = cmap(LogNorm(vmin=vmin, vmax=vmax)(grid[i, j]))
                    lum = 0.299 * r + 0.587 * g + 0.114 * b
                    dy = 0.0
                colour = '#222222' if lum > 0.55 else 'white'

                ax.text(j, i + dy, '\n'.join(parts), ha='center', va='center',
                        fontsize='xx-small' if len(parts) > 1 else 'x-small',
                        zorder=3, color=colour, linespacing=1.1)

    # Predicted stability boundary. Drawn from the CFL grid rather than the
    # measured errors, so where it sits relative to the blow-up region is the
    # comparison the plot exists to make.
    if args.cfl_contour:
        if np.isfinite(cflgrid).sum() < 4:
            print('--cfl-contour: too few CFL_GLL values logged to contour '
                  '(is IO_CFLSteps set in the sessions?)')
        else:
            ax.contour(xs, ys, np.ma.masked_invalid(cflgrid), levels=[1.0],
                       colors='white', linewidths=2.0, zorder=4)
            ax.contour(xs, ys, np.ma.masked_invalid(cflgrid), levels=[1.0],
                       colors='#8b1a1a', linewidths=1.0, zorder=5)

    ax.set_xticks(xs)
    ax.set_xticklabels([f'{x:g}' for x in x_vals])
    ax.set_yticks(ys)
    ax.set_yticklabels([f'{v:g}' for v in advy_vals])
    ax.set_xlabel('Number of adaptive runs (n)' if x_is_n
                  else 'Minimum layer height (h)')
    ax.set_ylabel('Advection velocity (advy)')

    # ALE runs log their mesh motion, so label each column with it. Columns are
    # ordinal, so this is a second row of labels rather than a rescaled axis.
    if gridvel:
        secax = ax.secondary_xaxis('top')
        secax.set_xticks(xs)
        secax.set_xticklabels(
            [f'{gridvel[x]:.3g}' if x in gridvel else '-' for x in x_vals])
        secax.set_xlabel('Peak grid velocity ($v_y$, minimum over cycles)')
        missing = [x for x in x_vals if x not in gridvel]
        if missing:
            print(f'no grid velocity logged for {"n" if x_is_n else "h"} = '
                  f'{", ".join(f"{x:g}" for x in missing)}')

    # Clear the secondary axis and its label when one is present.
    title_pad = 32 if gridvel else None
    if args.title is not None:
        ax.set_title(args.title, pad=title_pad)
    else:
        ax.set_title(f'{quantity}  '
                     f'({fixed_desc}, p={args.fix_p}, {args.method})',
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
    cbar.set_label(quantity)

    ax.set_xlim(-0.5, len(x_vals) - 0.5)
    ax.set_ylim(-0.5, len(advy_vals) - 0.5)

    if args.save:
        plt.savefig(args.save, dpi=300, bbox_inches='tight')
        print(f'Saved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()

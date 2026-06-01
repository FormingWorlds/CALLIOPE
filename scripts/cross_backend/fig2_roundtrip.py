"""Figure 2: Each backend round-trips internally across the magma-ocean
temperature and redox range.

We sweep T_magma because the chemistry is genuinely T-dependent. The
modified equilibrium constants are evaluated at T; the Dasgupta (2022)
nitrogen solubility has explicit T and fO2 dependences; the Gaillard
(2022) sulfur solubility has an explicit fO2 dependence. A round-trip
that only worked at one T would not be evidence of internal
consistency. Each backend must invert cleanly across the full range
where the calibrated chemistry is valid.

Each (T, dIW_input) combination produces one residual dIW_recovered −
dIW_input. The figure plots this residual vs T_magma for each input
dIW, separately for the two backends. If the chemistry path is
internally consistent the residual should be sub-tolerance everywhere.

Reused output: writes `data/fig2_roundtrip.csv` so the docs page can
quote the worst-case residual without re-running the harness.
"""

from __future__ import annotations

import csv
import logging
import time

import matplotlib.pyplot as plt
import numpy as np

from .plot_style import DATA_DIR, apply_style, panel_label, save
from .verification import round_trip_atmodeller, round_trip_calliope

log = logging.getLogger('cross_backend.fig2')


T_GRID = [1500.0, 2000.0, 2500.0, 3000.0]
DIW_GRID = [-2.0, 0.0, 2.0, 4.0]

# Discrete high-contrast palette, one colour per T_magma. We sweep T
# because the chemistry is genuinely T-dependent: the equilibrium
# constants, the Dasgupta nitrogen solubility, and the Gaillard sulfur
# solubility all carry explicit T (and fO2) terms. A round-trip that
# only worked at one T would not be evidence of internal consistency.
# Cool -> warm matches the colour scale to the temperature.
T_COLORS = {
    1500.0: '#3949ab',  # indigo (coolest)
    2000.0: '#00897b',  # teal
    2500.0: '#fb8c00',  # orange
    3000.0: '#d81b60',  # magenta (hottest)
}


def collect() -> dict:
    """Return {backend: [(T, dIW_in, dIW_out), ...]} for both backends."""
    out = {'calliope': [], 'atmodeller': []}
    for backend, rt in (
        ('calliope', round_trip_calliope),
        ('atmodeller', round_trip_atmodeller),
    ):
        log.info('Round-trip: %s', backend)
        for T in T_GRID:
            for dIW in DIW_GRID:
                t0 = time.time()
                try:
                    _, recov = rt(T, dIW)
                except Exception as exc:  # noqa: BLE001
                    log.warning('  %s T=%.0f dIW=%.1f raised %s', backend, T, dIW, exc)
                    recov = float('nan')
                dt = time.time() - t0
                log.info('  T=%.0f dIW=%+.1f -> %+.3f  (%.1fs)', T, dIW, recov, dt)
                out[backend].append((T, dIW, recov))
    return out


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'fig2_roundtrip.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['backend', 'T_K', 'dIW_in', 'dIW_recovered', 'residual_dex'])
        for backend, rows in data.items():
            for T, dIW, recov in rows:
                resid = recov - dIW if np.isfinite(recov) else float('nan')
                w.writerow([backend, T, dIW, recov, resid])
    log.info('Wrote %s', csv_path)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), sharey=True)

    tol_band = 0.01
    y_window = 0.5  # dex on each side; off-scale points get a triangle

    for ax, backend, label_short in (
        (axes[0], 'calliope', 'CALLIOPE'),
        (axes[1], 'atmodeller', 'atmodeller'),
    ):
        ax.axhspan(
            -tol_band,
            tol_band,
            color='k',
            alpha=0.07,
            linewidth=0,
            label=rf'$\pm {tol_band:g}$ dex band',
        )
        ax.axhline(0.0, color='k', alpha=0.4, linewidth=0.7)

        # Small x-jitter per T so the four T markers fan out
        # horizontally at each input Delta-IW rather than stacking on
        # one another. Each T sits at the same y-residual; the offset
        # only affects horizontal placement so the colours are
        # individually visible. Jitter is symmetric around the integer
        # dIW tick to keep the eye on the underlying Delta-IW value.
        n_T = len(T_GRID)
        jitter_step = 0.16
        for i_T, T in enumerate(T_GRID):
            dx = (i_T - (n_T - 1) / 2.0) * jitter_step
            xs_ok = []
            ys_ok = []
            xs_off = []
            for Tx, dIWx, recov in data[backend]:
                if Tx != T:
                    continue
                if not np.isfinite(recov):
                    xs_off.append(dIWx)
                    continue
                resid = recov - dIWx
                if abs(resid) > y_window:
                    xs_off.append(dIWx)
                else:
                    xs_ok.append(dIWx)
                    ys_ok.append(resid)
            order = np.argsort(xs_ok) if xs_ok else []
            if xs_ok:
                xs_ok = np.array(xs_ok)[order] + dx
                ys_ok = np.array(ys_ok)[order]
                ax.plot(
                    xs_ok,
                    ys_ok,
                    marker='o',
                    markersize=7.0,
                    linewidth=0,
                    color=T_COLORS[T],
                    markeredgecolor='k',
                    markeredgewidth=0.5,
                    label=rf'$T_\mathrm{{magma}} = {int(T)}$ K',
                )
            for dIWx in xs_off:
                ax.plot(
                    dIWx + dx,
                    -y_window * 0.9,
                    marker='v',
                    markersize=12,
                    linewidth=0,
                    color=T_COLORS[T],
                    markeredgecolor='k',
                    markeredgewidth=0.6,
                    clip_on=False,
                )

        ax.set_xlabel(r'input $\Delta\mathrm{IW}$ [dex]')
        ax.set_title(label_short)
        ax.set_xlim(min(DIW_GRID) - 0.5, max(DIW_GRID) + 0.5)
        ax.set_xticks(DIW_GRID)
        ax.set_ylim(-y_window, y_window)

    axes[0].set_ylabel(r'residual: recovered $-$ input $\Delta\mathrm{IW}$ [dex]')

    # One legend total, below the two panels. Reorder so input dIW
    # entries come first and the tolerance-band entry last.
    handles, labels = axes[0].get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: ('band' in labels[i], labels[i]))
    handles = [handles[i] for i in order]
    labels = [labels[i] for i in order]
    fig.legend(
        handles,
        labels,
        loc='lower center',
        ncol=len(labels),
        bbox_to_anchor=(0.5, -0.03),
        frameon=False,
        fontsize=9.5,
    )

    panel_label(axes[0], '(a)')
    panel_label(axes[1], '(b)')

    fig.suptitle(
        'Internal round-trip: buffered mode $\\to$ authoritative-O recovers the input $\\Delta$IW',
        fontsize=11.5,
        y=0.99,
    )

    # Compact the bottom margin to make room for the bottom legend.
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.96))

    paths = save(fig, 'fig2_roundtrip')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s %(message)s'
    )
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

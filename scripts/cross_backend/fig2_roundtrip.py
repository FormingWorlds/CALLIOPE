"""Figure 2: Each backend round-trips internally.

For each backend independently: starting from a buffered-mode call at
known Delta-IW, extract the resulting O budget, feed it back into the
authoritative-O entry point, and verify the recovered Delta-IW matches
the input. Two panels (one per backend), identity line, points coloured
by T_magma.

Reused output: writes `data/fig2_roundtrip.csv` so the docs page can
quote the worst-case residual without re-running the harness.
"""

from __future__ import annotations

import csv
import logging
import time

import matplotlib.pyplot as plt
import numpy as np

from .plot_style import COLOR_ATM, COLOR_CAL, DATA_DIR, apply_style, panel_label, save
from .verification import round_trip_atmodeller, round_trip_calliope

log = logging.getLogger('cross_backend.fig2')


T_GRID = [1500.0, 2000.0, 2500.0, 3000.0]
DIW_GRID = [-2.0, 0.0, 2.0, 4.0]


def collect() -> dict:
    """Return {backend: [(T, dIW_in, dIW_out), ...]} for both backends."""
    out = {'calliope': [], 'atmodeller': []}
    for backend, rt in (('calliope', round_trip_calliope), ('atmodeller', round_trip_atmodeller)):
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

    # Write CSV alongside.
    csv_path = DATA_DIR / 'fig2_roundtrip.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['backend', 'T_K', 'dIW_in', 'dIW_recovered', 'residual_dex'])
        for backend, rows in data.items():
            for T, dIW, recov in rows:
                resid = recov - dIW if np.isfinite(recov) else float('nan')
                w.writerow([backend, T, dIW, recov, resid])
    log.info('Wrote %s', csv_path)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
    colour_cycle = plt.colormaps['viridis'](np.linspace(0.0, 0.85, len(T_GRID)))

    for ax, backend, label_short in (
        (axes[0], 'calliope', 'CALLIOPE'),
        (axes[1], 'atmodeller', 'atmodeller'),
    ):
        for i_T, T in enumerate(T_GRID):
            xs = [dIW for (Tx, dIW, _) in data[backend] if Tx == T and np.isfinite(_)]
            ys = [recov for (Tx, dIW, recov) in data[backend] if Tx == T and np.isfinite(recov)]
            ax.scatter(xs, ys, s=50, color=colour_cycle[i_T],
                       edgecolor='k', linewidth=0.5,
                       label=f'$T = {T:.0f}$ K')
        identity = np.array([min(DIW_GRID) - 0.5, max(DIW_GRID) + 0.5])
        ax.plot(identity, identity, color='k', linewidth=0.8, alpha=0.5)
        ax.set_xlabel(r'input $\Delta\mathrm{IW}$ [dex]')
        ax.set_title(label_short)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(identity[0], identity[1])
        ax.set_ylim(identity[0], identity[1])

    axes[0].set_ylabel(r'recovered $\Delta\mathrm{IW}$ [dex]')
    axes[0].legend(loc='lower right', fontsize=8.5)
    panel_label(axes[0], '(a)')
    panel_label(axes[1], '(b)')

    fig.suptitle('Internal round-trip: buffered mode → authoritative-O → recovered $\\Delta$IW',
                 fontsize=11.5, y=0.99)

    paths = save(fig, 'fig2_roundtrip')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

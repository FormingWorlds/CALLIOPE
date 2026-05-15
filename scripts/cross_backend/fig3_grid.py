"""Figure 3: Cross-backend Delta-IW disagreement as a function of T_magma.

Both backends called at the canonical Earth-BSE Krijt+2023 volatile
budget with H/C/N/S fixed and the volatile O reference set by a
buffered-mode call at Delta-IW = +3.5 (Sossi 2020). For each T_magma
in the grid the two backends produce a converged Delta-IW; the figure
shows them side by side together with the analytical buffer offset
(Hirschmann minus O'Neill at that T).

A 2D (T, O-budget) heatmap was attempted first and abandoned: the
authoritative-O entry point has a known non-monotonic regime at
sub-trace O budgets which makes a wider O sweep dominated by basin-
selection effects rather than backend-physics differences. The single-
axis T sweep at fixed Earth-like O reports the backend-physics
contrast cleanly.
"""

from __future__ import annotations

import csv
import logging
import time

import matplotlib.pyplot as plt
import numpy as np

from . import buffers
from .inventories import EARTH_BSE_KRIJT23
from .plot_style import COLOR_ATM, COLOR_CAL, DATA_DIR, apply_style, panel_label, save
from .runners import run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.fig3')


T_GRID = np.array([1800.0, 2000.0, 2400.0, 2800.0, 3000.0])
FO2_HINT = 3.5


def collect() -> dict:
    nT = len(T_GRID)
    dIW_cal = np.full(nT, np.nan)
    dIW_atm = np.full(nT, np.nan)
    P_cal = np.full(nT, np.nan)
    P_atm = np.full(nT, np.nan)
    for i, T in enumerate(T_GRID):
        t0 = time.time()
        r_cal = run_calliope(EARTH_BSE_KRIJT23, T_magma=float(T), fO2_hint=FO2_HINT)
        r_atm = run_atmodeller(EARTH_BSE_KRIJT23, T_magma=float(T))
        dt = time.time() - t0
        log.info('T=%4.0f cal=%+6.3f (%s)  atm=%+6.3f (%s)  %.1fs',
                 T,
                 r_cal.fO2_shift_derived, 'ok' if r_cal.converged else 'fail',
                 r_atm.fO2_shift_derived, 'ok' if r_atm.converged else 'fail',
                 dt)
        if r_cal.converged:
            dIW_cal[i] = r_cal.fO2_shift_derived
            P_cal[i] = r_cal.total_P_bar
        if r_atm.converged:
            dIW_atm[i] = r_atm.fO2_shift_derived
            P_atm[i] = r_atm.total_P_bar
    return dict(dIW_cal=dIW_cal, dIW_atm=dIW_atm, P_cal=P_cal, P_atm=P_atm)


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()
    dIW_cal = data['dIW_cal']
    dIW_atm = data['dIW_atm']

    buf_offset = np.array([
        buffers.hirschmann_minus_oneill_offset(np.array([T]))[0] for T in T_GRID
    ])
    # Predicted cross-backend gap (atm - cal) under identical chemistry
    # is -(hirschmann - oneill) = -buf_offset. So the buffer-predicted
    # atmodeller curve is dIW_cal - buf_offset.
    dIW_atm_predicted = dIW_cal - buf_offset
    raw_gap = dIW_atm - dIW_cal
    corrected = raw_gap + buf_offset  # see fig4 for sign derivation

    csv_path = DATA_DIR / 'fig3_grid.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['T_K', 'dIW_calliope', 'dIW_atmodeller',
                    'buffer_offset_dex', 'raw_delta_atm_minus_cal',
                    'residual_after_buffer_correction',
                    'P_total_cal_bar', 'P_total_atm_bar'])
        for i, T in enumerate(T_GRID):
            w.writerow([T, dIW_cal[i], dIW_atm[i],
                        buf_offset[i], raw_gap[i], corrected[i],
                        data['P_cal'][i], data['P_atm'][i]])
    log.info('Wrote %s', csv_path)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7.2, 6.0), sharex=True,
        gridspec_kw={'height_ratios': [1.6, 1.0], 'hspace': 0.10},
    )

    ax_top.plot(T_GRID, dIW_cal, marker='o', color=COLOR_CAL, label='CALLIOPE')
    ax_top.plot(T_GRID, dIW_atm, marker='s', color=COLOR_ATM, label='atmodeller (default)')
    ax_top.plot(T_GRID, dIW_atm_predicted, marker='x', linestyle=':',
                color=COLOR_ATM, alpha=0.6,
                label="atmodeller predicted from buffer alone\n(= CALLIOPE − (Hirschmann − O'Neill))")
    ax_top.set_ylabel(r'$\Delta\mathrm{IW}$ [dex]')
    ax_top.set_title('Cross-backend $\\Delta$IW at Earth-BSE volatile inventory, $\\Phi = 1$')
    # Pad the y-axis so the legend fits below the CALLIOPE line at the
    # cold end without clipping the converged data, and so the (a)
    # panel label in the upper-left has clear headroom above the
    # CALLIOPE data point at T = 1800 K.
    y_top_min = min(float(np.nanmin(dIW_atm_predicted)), float(np.nanmin(dIW_atm))) - 0.45
    y_top_max = float(np.nanmax(dIW_cal)) + 0.30
    ax_top.set_ylim(y_top_min, y_top_max)
    ax_top.legend(loc='lower left', fontsize=9.0,
                  framealpha=0.92, facecolor='white', edgecolor='none')
    panel_label(ax_top, '(a)')

    ax_bot.axhline(0.0, color='k', alpha=0.4, linewidth=0.7)
    ax_bot.axhline(0.1, color='k', alpha=0.25, linestyle='--', linewidth=0.7)
    ax_bot.axhline(-0.1, color='k', alpha=0.25, linestyle='--', linewidth=0.7)
    ax_bot.plot(T_GRID, raw_gap, marker='o', color='#7a7a7a', label='raw $\\Delta$IW$_\\mathrm{atm}-\\Delta$IW$_\\mathrm{cal}$')
    ax_bot.plot(T_GRID, corrected, marker='D', color=COLOR_CAL,
                label='after buffer correction (residual chemistry gap)')
    ax_bot.text(T_GRID[-1], 0.115, r'$\pm 0.1$ dex solver tolerance',
                fontsize=8.5, va='bottom', ha='right', alpha=0.6)
    ax_bot.set_xlabel(r'$T_\mathrm{magma}$ [K]')
    ax_bot.set_ylabel('disagreement [dex]')
    # Pad the y-axis: bottom margin makes room for the legend, top
    # margin makes room for the (b) panel label and the tolerance text
    # so neither overlaps the corrected-residual diamond at the cold
    # end where the corrected residual sits at about +0.08 dex.
    y_bot_min = float(np.nanmin(raw_gap)) - 0.55
    y_bot_max = max(0.75, float(np.nanmax(corrected)) + 0.55)
    ax_bot.set_ylim(y_bot_min, y_bot_max)
    ax_bot.legend(loc='lower left', fontsize=9.0,
                  framealpha=0.92, facecolor='white', edgecolor='none')
    panel_label(ax_bot, '(b)')

    paths = save(fig, 'fig3_grid')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

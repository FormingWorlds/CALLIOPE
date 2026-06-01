"""Figure 3: Cross-backend Delta-IW disagreement as a function of T_magma.

Both backends called at the canonical Earth-BSE Krijt+2023 volatile
budget with H/C/N/S fixed and the volatile O reference set by a
buffered-mode call at Delta-IW = +3.5 (Sossi 2020). For each T_magma
in the grid the two backends produce a converged Delta-IW; the figure
shows them side by side together with the analytical buffer offsets
(Hirschmann minus Fischer for the current default; Hirschmann minus
O'Neill for the legacy buffer).

CALLIOPE is run twice at every grid point: once with the current
default Fischer 2011 buffer and once with the legacy O'Neill 2002
buffer. The Fischer trace is much closer to atmodeller because
Fischer 2011 is closer to Hirschmann than O'Neill 2002 is across the
magma-ocean range (~0.1 dex residual at 2000 K vs ~0.95 dex).

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
from .plot_style import (
    COLOR_ATM,
    COLOR_CAL,
    COLOR_FIS,
    DATA_DIR,
    apply_style,
    panel_label,
    save,
)
from .runners import run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.fig3')


T_GRID = np.array([1800.0, 2000.0, 2400.0, 2800.0, 3000.0])
FO2_HINT = 3.5


def collect() -> dict:
    nT = len(T_GRID)
    dIW_cal_fis = np.full(nT, np.nan)
    dIW_cal_one = np.full(nT, np.nan)
    dIW_atm = np.full(nT, np.nan)
    P_cal_fis = np.full(nT, np.nan)
    P_cal_one = np.full(nT, np.nan)
    P_atm = np.full(nT, np.nan)
    for i, T in enumerate(T_GRID):
        t0 = time.time()
        r_cal_fis = run_calliope(
            EARTH_BSE_KRIJT23,
            T_magma=float(T),
            fO2_hint=FO2_HINT,
            buffer='fischer',
        )
        r_cal_one = run_calliope(
            EARTH_BSE_KRIJT23,
            T_magma=float(T),
            fO2_hint=FO2_HINT,
            buffer='oneill',
        )
        r_atm = run_atmodeller(EARTH_BSE_KRIJT23, T_magma=float(T))
        dt = time.time() - t0
        log.info(
            'T=%4.0f cal_F=%+6.3f cal_O=%+6.3f atm=%+6.3f  %.1fs',
            T,
            r_cal_fis.fO2_shift_derived,
            r_cal_one.fO2_shift_derived,
            r_atm.fO2_shift_derived,
            dt,
        )
        if r_cal_fis.converged:
            dIW_cal_fis[i] = r_cal_fis.fO2_shift_derived
            P_cal_fis[i] = r_cal_fis.total_P_bar
        if r_cal_one.converged:
            dIW_cal_one[i] = r_cal_one.fO2_shift_derived
            P_cal_one[i] = r_cal_one.total_P_bar
        if r_atm.converged:
            dIW_atm[i] = r_atm.fO2_shift_derived
            P_atm[i] = r_atm.total_P_bar
    return dict(
        dIW_cal_fis=dIW_cal_fis,
        dIW_cal_one=dIW_cal_one,
        dIW_atm=dIW_atm,
        P_cal_fis=P_cal_fis,
        P_cal_one=P_cal_one,
        P_atm=P_atm,
    )


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()
    dIW_cal_fis = data['dIW_cal_fis']
    dIW_cal_one = data['dIW_cal_one']
    dIW_atm = data['dIW_atm']

    buf_offset_one = np.array(
        [buffers.hirschmann_minus_oneill_offset(np.array([T]))[0] for T in T_GRID]
    )
    buf_offset_fis = np.array(
        [buffers.hirschmann_minus_fischer_offset(np.array([T]))[0] for T in T_GRID]
    )
    # Predicted cross-backend gap (atm - cal) under identical chemistry
    # is -(hirschmann - cal_buffer). The buffer-predicted atmodeller
    # curve is therefore (cal − buffer_offset).
    dIW_atm_pred_from_fis = dIW_cal_fis - buf_offset_fis
    raw_gap_one = dIW_atm - dIW_cal_one
    raw_gap_fis = dIW_atm - dIW_cal_fis
    corrected_one = raw_gap_one + buf_offset_one
    corrected_fis = raw_gap_fis + buf_offset_fis

    csv_path = DATA_DIR / 'fig3_grid.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                'T_K',
                'dIW_calliope_fischer',
                'dIW_calliope_oneill',
                'dIW_atmodeller',
                'buffer_offset_HminusF_dex',
                'buffer_offset_HminusO_dex',
                'raw_gap_F_dex',
                'raw_gap_O_dex',
                'residual_after_buffer_F_dex',
                'residual_after_buffer_O_dex',
                'P_total_cal_fischer_bar',
                'P_total_cal_oneill_bar',
                'P_total_atm_bar',
            ]
        )
        for i, T in enumerate(T_GRID):
            w.writerow(
                [
                    T,
                    dIW_cal_fis[i],
                    dIW_cal_one[i],
                    dIW_atm[i],
                    buf_offset_fis[i],
                    buf_offset_one[i],
                    raw_gap_fis[i],
                    raw_gap_one[i],
                    corrected_fis[i],
                    corrected_one[i],
                    data['P_cal_fis'][i],
                    data['P_cal_one'][i],
                    data['P_atm'][i],
                ]
            )
    log.info('Wrote %s', csv_path)

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(7.6, 6.4),
        sharex=True,
        gridspec_kw={'height_ratios': [1.6, 1.0], 'hspace': 0.10},
    )

    ax_top.plot(
        T_GRID,
        dIW_cal_fis,
        marker='o',
        color=COLOR_CAL,
        label='CALLIOPE (Fischer 2011, default)',
    )
    ax_top.plot(
        T_GRID,
        dIW_cal_one,
        marker='o',
        linestyle='--',
        color=COLOR_FIS,
        alpha=0.85,
        label="CALLIOPE (O'Neill 2002, legacy)",
    )
    ax_top.plot(
        T_GRID, dIW_atm, marker='s', color=COLOR_ATM, label='atmodeller (Hirschmann composite)'
    )
    ax_top.plot(
        T_GRID,
        dIW_atm_pred_from_fis,
        marker='x',
        linestyle=':',
        color=COLOR_ATM,
        alpha=0.6,
        label='atmodeller predicted from buffer alone\n(= CALLIOPE-F − (Hirschmann − Fischer))',
    )
    ax_top.set_ylabel(r'$\Delta\mathrm{IW}$ [dex]')
    ax_top.set_title('Cross-backend $\\Delta$IW at Earth-BSE volatile inventory, $\\Phi = 1$')
    ymins = [
        float(np.nanmin(dIW_atm_pred_from_fis)),
        float(np.nanmin(dIW_atm)),
        float(np.nanmin(dIW_cal_one)),
    ]
    y_top_min = min(ymins) - 0.55
    y_top_max = (
        max(
            float(np.nanmax(dIW_cal_fis)),
            float(np.nanmax(dIW_cal_one)),
        )
        + 0.40
    )
    ax_top.set_ylim(y_top_min, y_top_max)
    ax_top.legend(
        loc='lower left', fontsize=8.5, framealpha=0.92, facecolor='white', edgecolor='none'
    )
    panel_label(ax_top, '(a)')

    ax_bot.axhline(0.0, color='k', alpha=0.4, linewidth=0.7)
    ax_bot.axhline(0.1, color='k', alpha=0.25, linestyle='--', linewidth=0.7)
    ax_bot.axhline(-0.1, color='k', alpha=0.25, linestyle='--', linewidth=0.7)
    ax_bot.plot(
        T_GRID,
        raw_gap_fis,
        marker='o',
        color='#7a7a7a',
        label='raw, Fischer default ($\\Delta$IW$_\\mathrm{atm}-\\Delta$IW$_\\mathrm{cal,F}$)',
    )
    ax_bot.plot(
        T_GRID,
        raw_gap_one,
        marker='o',
        linestyle='--',
        color='#444444',
        alpha=0.7,
        label="raw, O'Neill legacy ($\\Delta$IW$_\\mathrm{atm}-\\Delta$IW$_\\mathrm{cal,O}$)",
    )
    ax_bot.plot(
        T_GRID,
        corrected_fis,
        marker='D',
        color=COLOR_CAL,
        label='after buffer correction\n(residual chemistry gap, either buffer)',
    )
    ax_bot.text(
        T_GRID[-1],
        0.115,
        r'$\pm 0.1$ dex solver tolerance',
        fontsize=8.5,
        va='bottom',
        ha='right',
        alpha=0.6,
    )
    ax_bot.set_xlabel(r'$T_\mathrm{magma}$ [K]')
    ax_bot.set_ylabel('disagreement [dex]')
    y_bot_min = (
        min(
            float(np.nanmin(raw_gap_fis)),
            float(np.nanmin(raw_gap_one)),
        )
        - 0.60
    )
    y_bot_max = max(
        0.75,
        float(np.nanmax(raw_gap_fis)) + 0.55,
        float(np.nanmax(corrected_fis)) + 0.55,
    )
    ax_bot.set_ylim(y_bot_min, y_bot_max)
    ax_bot.legend(
        loc='lower left', fontsize=8.5, framealpha=0.92, facecolor='white', edgecolor='none'
    )
    panel_label(ax_bot, '(b)')

    paths = save(fig, 'fig3_grid')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s %(message)s'
    )
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

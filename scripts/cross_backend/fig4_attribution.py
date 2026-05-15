"""Figure 4: Attribution of cross-backend Delta-IW disagreement at a
single canonical Earth scenario.

Bar chart. From left to right, each bar shows the |Delta-IW|
disagreement at the canonical Earth fiducial:

1. Raw, legacy: CALLIOPE with the O'Neill 2002 IW buffer (the legacy
   choice) vs atmodeller default.
2. Raw, current default: CALLIOPE with the Fischer 2011 IW buffer (the
   library default since the buffer audit) vs atmodeller default.
3. After analytical buffer correction (Hirschmann - Fischer
   subtracted), residual chemistry gap with default solubility maps.
4. After also disabling H2 / CO / CH4 solubility in atmodeller
   (matching CALLIOPE's Bower 2022 §2.2.3 convention).

The drop from bar 1 to bar 2 is the buffer-default change alone;
bars 2-4 are the residuals after successive alignment moves. A
horizontal line marks 0.1 dex, the per-element solver tolerance.

Reusable: pass `T_magma` and inventory factor to reuse the script for a
different fiducial scenario.
"""

from __future__ import annotations

import csv
import logging

import matplotlib.pyplot as plt
import numpy as np

from . import buffers
from .inventories import EARTH_BSE_KRIJT23
from .plot_style import COLOR_ATM, COLOR_CAL, COLOR_FIS, DATA_DIR, apply_style, save
from .runners import _CALLIOPE_ALIGNED_ATM_SOL, _DEFAULT_ATM_SOL, run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.fig4')


def collect(T_magma: float = 2000.0) -> dict:
    """Compute the attribution stages with both CALLIOPE buffer
    choices.
    """
    inv = EARTH_BSE_KRIJT23
    log.info('Step 1: as-shipped defaults (Fischer 2011 buffer)')
    cal_fis = run_calliope(inv, T_magma=T_magma, fO2_hint=3.5, buffer='fischer')
    cal_one = run_calliope(inv, T_magma=T_magma, fO2_hint=3.5, buffer='oneill')
    atm_def = run_atmodeller(inv, T_magma=T_magma)
    log.info('  CALLIOPE dIW (Fischer)  = %+.3f', cal_fis.fO2_shift_derived)
    log.info("  CALLIOPE dIW (O'Neill)  = %+.3f", cal_one.fO2_shift_derived)
    log.info('  atmodeller dIW (default) = %+.3f', atm_def.fO2_shift_derived)

    raw_gap_fis = atm_def.fO2_shift_derived - cal_fis.fO2_shift_derived
    raw_gap_one = atm_def.fO2_shift_derived - cal_one.fO2_shift_derived
    buf_offset_fis = float(buffers.hirschmann_minus_fischer_offset(np.array([T_magma]))[0])
    buf_offset_one = float(buffers.hirschmann_minus_oneill_offset(np.array([T_magma]))[0])
    # Buffer contribution to (dIW_atm - dIW_cal) under identical
    # chemistry: -(hirschmann - cal_buffer) = -buf_offset. Residual is
    # raw_gap - (-buf_offset) = raw_gap + buf_offset.
    after_buffer_fis = raw_gap_fis + buf_offset_fis

    log.info('Step 2: align atmodeller solubility to CALLIOPE convention')
    atm_align = run_atmodeller(
        inv, T_magma=T_magma, solubility_map=_CALLIOPE_ALIGNED_ATM_SOL,
    )
    log.info('  atmodeller dIW (aligned) = %+.3f', atm_align.fO2_shift_derived)
    after_solubility = atm_align.fO2_shift_derived - cal_fis.fO2_shift_derived + buf_offset_fis

    return dict(
        T_magma=T_magma,
        cal_fischer_dIW=cal_fis.fO2_shift_derived,
        cal_oneill_dIW=cal_one.fO2_shift_derived,
        atm_default_dIW=atm_def.fO2_shift_derived,
        atm_aligned_dIW=atm_align.fO2_shift_derived,
        raw_gap_fischer=raw_gap_fis,
        raw_gap_oneill=raw_gap_one,
        buffer_offset_fischer=buf_offset_fis,
        buffer_offset_oneill=buf_offset_one,
        after_buffer_fischer=after_buffer_fis,
        after_solubility=after_solubility,
    )


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'fig4_attribution.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        for k, v in data.items():
            w.writerow([k, v])
    log.info('Wrote %s', csv_path)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))

    labels = [
        "Raw, O'Neill\n(legacy buffer)",
        'Raw, Fischer\n(default buffer)',
        'After buffer\ncorrection',
        'After also matching\nsolubility selection',
    ]
    values = [
        abs(data['raw_gap_oneill']),
        abs(data['raw_gap_fischer']),
        abs(data['after_buffer_fischer']),
        abs(data['after_solubility']),
    ]
    colors = [COLOR_FIS, COLOR_ATM, '#c79f3a', COLOR_CAL]

    bars = ax.bar(labels, values, color=colors, edgecolor='k', linewidth=0.6, width=0.55)
    tol = 0.10
    # Always label above the bar in dark text so labels stay readable
    # regardless of the bar fill colour. Small bars whose top sits
    # within +/-0.03 dex of the tolerance line get pushed higher so
    # their label clears the dashed line and the per-element-solver
    # annotation drawn just above it.
    for bar, val in zip(bars, values):
        x = bar.get_x() + bar.get_width() / 2
        if abs(val - tol) <= 0.04:
            y = tol + 0.045
        else:
            y = val + 0.018
        ax.text(x, y, f'{val:.2f} dex',
                ha='center', va='bottom', fontsize=10, color='#1a1a1a')

    ax.axhline(tol, color='k', alpha=0.4, linestyle='--', linewidth=0.8)
    # Tolerance annotation parked in the gap between bars 1 and 2,
    # well clear of every bar label and the dashed line itself.
    ax.text(0.5, tol + 0.005, 'per-element solver tolerance',
            fontsize=8.5, va='bottom', ha='center', alpha=0.6)

    ax.set_ylabel(r'$|\Delta\mathrm{IW}_{\mathrm{atm}}-\Delta\mathrm{IW}_{\mathrm{cal}}|$ [dex]')
    ax.set_title(
        f'Attribution of cross-backend $\\Delta$IW disagreement '
        f'at Earth-BSE, $T={data["T_magma"]:.0f}$ K, $\\Phi=1$',
        fontsize=11.0,
    )
    ax.set_ylim(0, max(values) * 1.25 + 0.1)
    ax.grid(axis='y', alpha=0.25)

    paths = save(fig, 'fig4_attribution')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

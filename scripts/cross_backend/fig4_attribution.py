"""Figure 4: Attribution of cross-backend Delta-IW disagreement at a
single canonical Earth scenario.

Bar chart. From left to right, each bar shows the |Delta-IW| disagreement
under a successively more aligned configuration:

1. Raw: both backends in their as-shipped default configuration.
2. After analytical buffer correction (Hirschmann - O'Neill subtracted).
3. After also disabling H2 / CO / CH4 solubility in atmodeller (matching
   CALLIOPE's Bower 2022 §2.2.3 convention).
4. Residual (interpretation): solubility-law differences for S2 that
   cannot be aligned at the wrapper level (Gaillard 2022 in CALLIOPE
   vs Boulliung 2023 in atmodeller) plus equilibrium-constant fits and
   numerical solver tolerance.

A horizontal line marks 0.1 dex, the per-element solver tolerance.

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
from .plot_style import COLOR_ATM, COLOR_CAL, DATA_DIR, apply_style, save
from .runners import _CALLIOPE_ALIGNED_ATM_SOL, _DEFAULT_ATM_SOL, run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.fig4')


def collect(T_magma: float = 2000.0) -> dict:
    """Compute the four attribution stages.

    Returns
    -------
    dict with keys 'raw_gap', 'after_buffer', 'after_solubility',
    'buffer_offset', and the underlying Delta-IW values.
    """
    inv = EARTH_BSE_KRIJT23
    log.info('Step 1: as-shipped defaults')
    cal = run_calliope(inv, T_magma=T_magma, fO2_hint=3.5)
    atm_def = run_atmodeller(inv, T_magma=T_magma)
    log.info('  CALLIOPE dIW = %+.3f', cal.fO2_shift_derived)
    log.info('  atmodeller dIW (default) = %+.3f', atm_def.fO2_shift_derived)

    raw_gap = atm_def.fO2_shift_derived - cal.fO2_shift_derived
    buf_offset = float(buffers.hirschmann_minus_oneill_offset(np.array([T_magma]))[0])
    # Buffer contribution to (dIW_atm - dIW_cal) under identical
    # chemistry: -(hirschmann - oneill) = -buf_offset. Residual after
    # removing that contribution: raw_gap - (-buf_offset) = raw_gap +
    # buf_offset. Missing the sign flip doubles the apparent
    # disagreement instead of cancelling it.
    after_buffer = raw_gap + buf_offset

    log.info('Step 2: align atmodeller solubility to CALLIOPE convention')
    atm_align = run_atmodeller(
        inv, T_magma=T_magma, solubility_map=_CALLIOPE_ALIGNED_ATM_SOL,
    )
    log.info('  atmodeller dIW (aligned) = %+.3f', atm_align.fO2_shift_derived)
    after_solubility = atm_align.fO2_shift_derived - cal.fO2_shift_derived + buf_offset

    return dict(
        T_magma=T_magma,
        cal_dIW=cal.fO2_shift_derived,
        atm_default_dIW=atm_def.fO2_shift_derived,
        atm_aligned_dIW=atm_align.fO2_shift_derived,
        raw_gap=raw_gap,
        buffer_offset=buf_offset,
        after_buffer=after_buffer,
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

    fig, ax = plt.subplots(figsize=(7.0, 4.2))

    labels = [
        'Raw\n(as-shipped)',
        'After buffer\ncorrection',
        'After also matching\nsolubility selection',
    ]
    values = [
        abs(data['raw_gap']),
        abs(data['after_buffer']),
        abs(data['after_solubility']),
    ]
    colors = [COLOR_ATM, '#c79f3a', COLOR_CAL]

    bars = ax.bar(labels, values, color=colors, edgecolor='k', linewidth=0.6, width=0.55)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                f'{val:.2f} dex',
                ha='center', va='bottom', fontsize=10)

    ax.axhline(0.10, color='k', alpha=0.4, linestyle='--', linewidth=0.8)
    ax.text(2.4, 0.115, 'per-element solver tolerance',
            fontsize=8.5, va='bottom', ha='right', alpha=0.6)

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

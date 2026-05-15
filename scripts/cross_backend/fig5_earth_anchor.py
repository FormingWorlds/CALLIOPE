"""Figure 5: Cross-backend Delta-IW at Earth-BSE against the empirical
Sossi (2020) anchor.

A single horizontal axis is Delta-IW. Two vertical markers show the
two backends' converged Delta-IW at the canonical Earth fiducial
(T_magma = 2000 K, Phi = 1, Krijt+2023 BSE H/C/N/S, volatile O derived
self-consistently). A shaded band is the Frost & McCammon (2008)
"Earth's mantle redox state" range (Delta-IW = +1 to +5, corresponding
to FMQ-3 to FMQ+1, with Sossi 2020 placing modern upper mantle at the
+3.5 centre).

This figure stress-tests the cross-backend disagreement against an
empirical anchor: if both backends fall inside the empirical range,
neither parameterisation is in tension with petrology at this single
fiducial. If one falls outside, that backend is in tension and should
be flagged to the reader.
"""

from __future__ import annotations

import csv
import logging

import matplotlib.pyplot as plt

from .inventories import EARTH_BSE_KRIJT23
from .plot_style import COLOR_ATM, COLOR_BG, COLOR_CAL, DATA_DIR, apply_style, save
from .runners import run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.fig5')


# Frost & McCammon (2008) Earth-mantle redox-state range, IW reference.
EARTH_MANTLE_DIW_LOW = 1.0
EARTH_MANTLE_DIW_HIGH = 5.0
SOSSI_2020_CENTER = 3.5


def collect(T_magma: float = 2000.0) -> dict:
    inv = EARTH_BSE_KRIJT23
    cal = run_calliope(inv, T_magma=T_magma, fO2_hint=3.5)
    atm = run_atmodeller(inv, T_magma=T_magma)
    log.info('CALLIOPE dIW = %+.3f', cal.fO2_shift_derived)
    log.info('atmodeller dIW = %+.3f', atm.fO2_shift_derived)
    return dict(
        T_magma=T_magma,
        cal_dIW=cal.fO2_shift_derived,
        atm_dIW=atm.fO2_shift_derived,
        cal_P_bar=cal.total_P_bar,
        atm_P_bar=atm.total_P_bar,
        cal_p_bar=cal.p_bar,
        atm_p_bar=atm.p_bar,
    )


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'fig5_earth_anchor.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['quantity', 'value'])
        w.writerow(['T_magma_K', data['T_magma']])
        w.writerow(['cal_dIW', data['cal_dIW']])
        w.writerow(['atm_dIW', data['atm_dIW']])
        w.writerow(['cal_P_total_bar', data['cal_P_bar']])
        w.writerow(['atm_P_total_bar', data['atm_P_bar']])
        w.writerow(['empirical_low_dIW', EARTH_MANTLE_DIW_LOW])
        w.writerow(['empirical_high_dIW', EARTH_MANTLE_DIW_HIGH])
        w.writerow(['sossi_2020_center', SOSSI_2020_CENTER])
    log.info('Wrote %s', csv_path)

    fig, ax = plt.subplots(figsize=(7.4, 3.6))

    ax.axhspan(0.4, 0.6, xmin=0.0, xmax=1.0, color=COLOR_BG, alpha=0.0)

    ax.axvspan(EARTH_MANTLE_DIW_LOW, EARTH_MANTLE_DIW_HIGH,
               color=COLOR_BG, alpha=0.6,
               label='Frost & McCammon 2008 Earth-mantle range')
    # Sossi 2020 anchor line: thicker dotted, raised z-order so it stays
    # visible even if a backend lands at the same dIW (e.g. CALLIOPE at
    # +3.50 here would otherwise hide the dotted line entirely).
    ax.axvline(SOSSI_2020_CENTER, color='k', alpha=0.7, linestyle=':',
               linewidth=1.8, zorder=3,
               label=fr'Sossi 2020 upper-mantle anchor: $\Delta\mathrm{{IW}} = {SOSSI_2020_CENTER:+.2f}$')

    ax.axvline(data['cal_dIW'], color=COLOR_CAL, linewidth=2.0, zorder=2,
               label=fr'CALLIOPE: $\Delta\mathrm{{IW}} = {data["cal_dIW"]:+.2f}$')
    ax.axvline(data['atm_dIW'], color=COLOR_ATM, linewidth=2.0, zorder=2,
               label=fr'atmodeller: $\Delta\mathrm{{IW}} = {data["atm_dIW"]:+.2f}$')

    # Suppress y-axis: this is a 1D figure.
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(-1.0, 7.0)
    ax.grid(axis='x', alpha=0.25)

    ax.set_xlabel(r'$\Delta\mathrm{IW}$ at $T_\mathrm{magma} = $' + f' {data["T_magma"]:.0f} K')
    ax.set_title(
        'Earth fiducial: both backends vs. the empirical mantle-fO$_2$ range',
        fontsize=11.0,
    )
    # Legend below the plot so the four entries do not crowd into the
    # data region where the three vertical lines already sit close
    # together at dIW = +3 to +3.5.
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2,
              frameon=False, fontsize=9.0, columnspacing=1.6,
              handlelength=2.4)

    paths = save(fig, 'fig5_earth_anchor')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

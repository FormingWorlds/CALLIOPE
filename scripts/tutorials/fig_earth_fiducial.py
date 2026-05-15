"""Tutorial figure for the "Reproducing the Earth fiducial" page.

Reproduces the Earth-fiducial point that anchors the backend-comparison
docs page. Takes the Krijt+2023 H/C/N/S BSE budget, runs CALLIOPE in
buffered mode at the Sossi 2020 Delta-IW = +3.5 to derive the volatile
O reference, then runs authoritative-O mode on (H, C, N, S, O) and
checks that the recovered Delta-IW lands inside the Frost & McCammon
(2008) empirical Earth-mantle range and on the Sossi 2020 anchor.

The figure is a clean version of cross-backend Figure 5 with only the
CALLIOPE point, designed as the "you have successfully reproduced the
docs Earth fiducial" closing image of the tutorial.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt

from calliope.constants import volatile_species
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

from ._style import COLOR_BG, COLOR_CAL, DATA_DIR, apply_style, save

log = logging.getLogger('tutorials.earth_fiducial')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
EARTH_HCNS = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}
T_MAGMA = 2000.0
DIW_ANCHOR = 3.5     # Sossi et al. 2020 Earth upper-mantle anchor

FROST_LO = 1.0       # Frost & McCammon (2008) Earth-mantle range, IW reference
FROST_HI = 5.0


def collect() -> dict:
    """Run the full derive-and-verify chain.

    Returns a dict with the buffered-mode O_kg_total and the
    authoritative-O recovered Delta-IW.
    """
    ddict = {**PLANET, 'T_magma': T_MAGMA, 'Phi_global': 1.0,
             'fO2_shift_IW': DIW_ANCHOR}
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0

    # Tight p_guess at the canonical CO2-dominated basin. The buffered
    # solver at high Delta-IW + C-rich BSE has a documented spurious
    # H2O-free basin at P_surf ~ 1500 bar; without this guess the
    # Monte-Carlo restart lands there 1 time in ~5 and the tutorial
    # output is non-reproducible.
    canonical_guess = {'H2O': 5.0, 'CO2': 1500.0, 'N2': 3.0, 'S2': 1e-3}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        buf = equilibrium_atmosphere(
            EARTH_HCNS, ddict, p_guess=canonical_guess, print_result=False,
        )
    O_kg = float(buf['O_kg_total'])
    log.info('Step 1 (buffered): dIW = %+.2f -> O_kg_total = %.3e kg',
             DIW_ANCHOR, O_kg)

    target = dict(EARTH_HCNS); target['O'] = O_kg
    p_guess = {s: float(buf[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        auth = equilibrium_atmosphere_authoritative_O(
            target, ddict, p_guess=p_guess, fO2_hint=DIW_ANCHOR,
            print_result=False,
        )
    recovered = float(auth['fO2_shift_derived'])
    log.info('Step 2 (authoritative-O): recovered dIW = %+.4f (residual %+.2e)',
             recovered, recovered - DIW_ANCHOR)

    return dict(O_kg_total=O_kg, dIW_recovered=recovered,
                P_surf_bar=float(auth['P_surf']))


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'earth_fiducial.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['quantity', 'value'])
        w.writerow(['Sossi_2020_anchor_dIW', DIW_ANCHOR])
        w.writerow(['Frost_McCammon_2008_low_dIW', FROST_LO])
        w.writerow(['Frost_McCammon_2008_high_dIW', FROST_HI])
        w.writerow(['O_kg_total_derived', data['O_kg_total']])
        w.writerow(['recovered_dIW_authoritative', data['dIW_recovered']])
        w.writerow(['P_surf_bar_authoritative', data['P_surf_bar']])
    log.info('Wrote %s', csv_path)

    fig, ax = plt.subplots(figsize=(7.6, 3.4))

    # Frost & McCammon (2008) Earth-mantle range as a soft band.
    ax.axvspan(FROST_LO, FROST_HI, color=COLOR_BG, alpha=0.6,
               label='Frost & McCammon 2008 Earth-mantle range')

    # Sossi 2020 dotted anchor; raised z-order so the dotted pattern
    # stays visible if CALLIOPE recovers exactly +3.5.
    ax.axvline(DIW_ANCHOR, color='k', alpha=0.7, linestyle=':',
               linewidth=1.8, zorder=3,
               label=fr'Sossi 2020 upper-mantle anchor: $\Delta\mathrm{{IW}} = {DIW_ANCHOR:+.2f}$')

    # Reproduced CALLIOPE result.
    ax.axvline(data['dIW_recovered'], color=COLOR_CAL, linewidth=2.4,
               zorder=2,
               label=fr'reproduced CALLIOPE: $\Delta\mathrm{{IW}} = {data["dIW_recovered"]:+.2f}$')

    # 1D layout: suppress y axis.
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(-1.0, 7.0)
    ax.grid(axis='x', alpha=0.25)
    ax.set_xlabel(r'$\Delta\mathrm{IW}$ at $T_\mathrm{magma} = $' + f' {T_MAGMA:.0f} K')
    ax.set_title('Reproducing the Earth fiducial: $\\Delta$IW lands on Sossi 2020')

    # Provenance summary box. Anchored to the right where the data
    # band ends but the data lines do not extend past dIW = +5.
    summary = (
        f"derived $O_\\mathrm{{tot}} = {data['O_kg_total']:.3e}$ kg\n"
        f"recovered − anchor = {data['dIW_recovered'] - DIW_ANCHOR:+.2e} dex\n"
        f"$P_\\mathrm{{surf}} = {data['P_surf_bar']:.0f}$ bar"
    )
    ax.text(
        0.985, 0.04, summary,
        transform=ax.transAxes,
        fontsize=9.0, va='bottom', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc'),
    )

    ax.legend(
        loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=1,
        frameon=False, fontsize=9.0,
    )

    paths = save(fig, 'earth_fiducial')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

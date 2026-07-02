"""Validation figure for the noble gas explanation page.

Runs each noble gas through ``equilibrium_atmosphere`` over a sweep of
budgets against a fixed Earth-like C-H-O-N-S background, then shows two
things the reader cannot see from the CHNOS-only tutorial figures:

* panel (a): the surface partial pressure rises with the supplied
  budget, approximately linearly where the noble gas is a trace
  component, with the curves bending together as the gas comes to
  dominate the atmosphere and shift its mean molar mass.
* panel (b): the split between atmosphere and melt is set by the gas's
  volume-based Henry constant (cm3 STP per gram per bar) and the
  atmosphere's mean molar mass. The molar mass cancels out of the
  mass-based ppmw-per-bar constant, so retention does not simply track
  that constant; in this background neon is the most retained in the melt
  and xenon the least, so xenon sits almost entirely in the atmosphere.

The split is budget-independent because Henry's law is linear, so panel
(b) is evaluated at a single representative budget.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

from calliope.constants import dict_colors, noble_gases, volatile_species
from calliope.solve import equilibrium_atmosphere, get_target_from_params

from ._style import DATA_DIR, apply_style, panel_label, save, species_label

log = logging.getLogger('tutorials.noble_gases')

PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
STATE = {
    'T_magma': 1800.0,
    'Phi_global': 1.0,
    'fO2_shift_IW': 0.5,
}
COMPOSITION = {
    'hydrogen_earth_oceans': 1.0,
    'CH_ratio': 0.1,
    'nitrogen_ppmw': 2.0,
    'sulfur_ppmw': 200.0,
}

# Budget sweep, in ppmw relative to the mantle mass.
BUDGETS_PPMW = [1.0, 3.0, 10.0, 30.0, 100.0]
# Representative budget for the partitioning bar chart.
SPLIT_PPMW = 10.0


def _base_ddict() -> dict:
    ddict = {**PLANET, **STATE, **COMPOSITION}
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    for gas in noble_gases:
        ddict[f'{gas}_included'] = 0
        ddict[f'{gas}_ppmw'] = 0.0
    return ddict


def _solve(gas: str, ppmw: float) -> dict:
    ddict = _base_ddict()
    ddict[f'{gas}_included'] = 1
    ddict[f'{gas}_ppmw'] = ppmw
    target = get_target_from_params(ddict)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return equilibrium_atmosphere(target, ddict, print_result=False, nguess=3000)


def collect() -> dict:
    pressures = {gas: [] for gas in noble_gases}
    atm_fraction = {}
    for gas in noble_gases:
        for ppmw in BUDGETS_PPMW:
            result = _solve(gas, ppmw)
            pressures[gas].append(float(result[f'{gas}_bar']))
        split = _solve(gas, SPLIT_PPMW)
        atm = float(split[f'{gas}_kg_atm'])
        liq = float(split[f'{gas}_kg_liquid'])
        atm_fraction[gas] = atm / (atm + liq)
        log.info(
            '%s: atmospheric fraction %.3f at %.0f ppmw', gas, atm_fraction[gas], SPLIT_PPMW
        )
    return {'pressures': pressures, 'atm_fraction': atm_fraction}


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()
    pressures = data['pressures']
    atm_fraction = data['atm_fraction']

    csv_path = DATA_DIR / 'noble_gases.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['gas', 'budget_ppmw', 'partial_pressure_bar'])
        for gas in noble_gases:
            for ppmw, p in zip(BUDGETS_PPMW, pressures[gas]):
                w.writerow([gas, ppmw, p])
        w.writerow([])
        w.writerow(['gas', 'atmospheric_fraction_at_%g_ppmw' % SPLIT_PPMW])
        for gas in noble_gases:
            w.writerow([gas, atm_fraction[gas]])
    log.info('Wrote %s', csv_path)

    fig, (ax_p, ax_f) = plt.subplots(1, 2, figsize=(10, 4.2))

    for gas in noble_gases:
        ax_p.plot(
            BUDGETS_PPMW,
            pressures[gas],
            marker='o',
            color=dict_colors[gas],
            label=species_label(gas),
        )
    ax_p.set_xscale('log')
    ax_p.set_yscale('log')
    ax_p.set_xlabel('Noble gas budget [ppmw of mantle]')
    ax_p.set_ylabel('Surface partial pressure [bar]')
    ax_p.legend(frameon=False, ncol=2, fontsize=9)
    panel_label(ax_p, 'a')

    xs = np.arange(len(noble_gases))
    fracs = [atm_fraction[gas] for gas in noble_gases]
    colors = [dict_colors[gas] for gas in noble_gases]
    ax_f.bar(xs, fracs, color=colors)
    ax_f.bar(xs, [1.0 - f for f in fracs], bottom=fracs, color=colors, alpha=0.35)
    ax_f.set_xticks(xs)
    ax_f.set_xticklabels([species_label(gas) for gas in noble_gases])
    ax_f.set_ylim(0.0, 1.0)
    ax_f.set_ylabel('Mass fraction (solid: atmosphere, faded: melt)')
    ax_f.set_xlabel('Noble gas (increasing atomic number)')
    panel_label(ax_f, 'b')

    fig.tight_layout()
    paths = save(fig, 'noble_gases')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    out = make_figure()
    for ext, path in out.items():
        print(f'{ext}: {path}')

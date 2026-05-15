"""Reference figure for the "First run" tutorial.

Runs the exact inputs the tutorial walks through (1 ocean of H, C/H =
0.1, 2 ppmw N, 200 ppmw S, T_magma = 2500 K, Phi = 1, Delta-IW =
+0.5), then plots the resulting surface partial pressures. Saved into
the docs so a reader who has just finished the tutorial can compare
their output against the canonical answer.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

from calliope.constants import dict_colors, volatile_species
from calliope.solve import equilibrium_atmosphere, get_target_from_params

from ._style import DATA_DIR, apply_style, save

log = logging.getLogger('tutorials.firstrun_reference')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
STATE = {
    'T_magma': 2500.0,
    'Phi_global': 1.0,
    'fO2_shift_IW': 0.5,
}
COMPOSITION = {
    'hydrogen_earth_oceans': 1.0,
    'CH_ratio': 0.1,
    'nitrogen_ppmw': 2.0,
    'sulfur_ppmw': 200.0,
}
SPECIES_TO_PLOT = ['H2O', 'CO2', 'H2', 'CO', 'CH4', 'N2', 'NH3', 'S2', 'SO2', 'H2S']


def collect() -> dict:
    ddict = {**PLANET, **STATE, **COMPOSITION}
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0

    target = get_target_from_params(ddict)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result = equilibrium_atmosphere(target, ddict, print_result=False)

    pressures = {sp: float(result[f'{sp}_bar']) for sp in SPECIES_TO_PLOT}
    return {
        'P_surf_bar': float(result['P_surf']),
        'M_atm_kg': float(result['M_atm']),
        'mean_mol_mass_g_per_mol': float(result['atm_kg_per_mol']) * 1e3,
        'pressures': pressures,
    }


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'firstrun_reference.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['species', 'partial_pressure_bar'])
        for sp, p in data['pressures'].items():
            w.writerow([sp, p])
        w.writerow(['__P_surf_bar', data['P_surf_bar']])
        w.writerow(['__M_atm_kg', data['M_atm_kg']])
        w.writerow(['__mean_mol_mass_g_per_mol', data['mean_mol_mass_g_per_mol']])
    log.info('Wrote %s', csv_path)

    # Plot species in descending pressure order so the most-abundant
    # species sit at the top of the figure; horizontal layout means
    # every numeric label has a fixed amount of space to its right
    # regardless of how small or large the value is.
    order = np.argsort([data['pressures'][sp] for sp in SPECIES_TO_PLOT])[::-1]
    species_sorted = [SPECIES_TO_PLOT[i] for i in order]
    pressures = np.array([data['pressures'][sp] for sp in species_sorted])
    colors = [dict_colors[sp] for sp in species_sorted]

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    y_pos = np.arange(len(species_sorted))
    ax.barh(y_pos, pressures, color=colors,
            edgecolor='black', linewidth=0.5)

    ax.set_xscale('log')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(species_sorted)
    ax.invert_yaxis()  # largest pressure at the top
    ax.set_xlabel('Surface partial pressure (bar)')
    ax.set_title(
        rf'First-run reference: $T_\mathrm{{magma}} = {STATE["T_magma"]:.0f}$ K, '
        rf'$\Phi = {STATE["Phi_global"]:.0f}$, '
        rf'$\Delta\mathrm{{IW}} = {STATE["fO2_shift_IW"]:+.1f}$, '
        rf'1 Earth-ocean H'
    )
    x_lo = 1e-10
    x_hi = max(1e4, float(pressures.max()) * 5000.0)
    ax.set_xlim(x_lo, x_hi)
    ax.grid(axis='x', which='both', alpha=0.3)

    # Numeric label to the right of each bar, in scientific notation
    # so all values share one format the reader can scan column-wise.
    for yi, p in zip(y_pos, pressures):
        if p <= 0 or not np.isfinite(p):
            label = 'below floor'
        else:
            label = f'{p:.2e} bar'
        ax.text(
            p * 2.0 if (p > 0 and np.isfinite(p)) else x_lo * 2.0,
            yi, label,
            ha='left', va='center', fontsize=9.0, color='#333333',
        )

    # Summary box anchored bottom-right so it does not collide with the
    # H2O label (top bar, label extends to about 1 bar on the x axis).
    summary = (
        f"$P_\\mathrm{{surf}} = {data['P_surf_bar']:.2f}$ bar\n"
        f"$M_\\mathrm{{atm}} = {data['M_atm_kg']:.2e}$ kg\n"
        f"mean $M$ = {data['mean_mol_mass_g_per_mol']:.2f} g/mol"
    )
    ax.text(
        0.985, 0.04, summary,
        transform=ax.transAxes,
        fontsize=9.0, va='bottom', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc'),
    )

    paths = save(fig, 'firstrun_reference')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

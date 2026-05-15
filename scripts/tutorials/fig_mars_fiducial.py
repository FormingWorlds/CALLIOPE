"""Tutorial figure for the "Mars-like atmosphere" page.

Runs the same redox and temperature conditions as the first-run
tutorial on a Mars-scaled inventory and Mars planetary parameters,
then overlays the Mars result on the Earth-BSE reference so the
reader can read off which species change most.

The Mars inventory is the Krijt+2023 Earth BSE H/C/N/S mass-scaled by
Mars / Earth mass (0.107) for illustrative purposes; this is not a
Mars-petrology BSE estimate. The pedagogical goal is to demonstrate
the workflow generalises to non-Earth planets, not to claim a
specific Mars composition.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

from calliope.constants import dict_colors, volatile_species
from calliope.solve import equilibrium_atmosphere

from ._style import DATA_DIR, apply_style, sci_fmt_plain, save, species_label

log = logging.getLogger('tutorials.mars_fiducial')


# Earth setup (Krijt+2023 BSE H/C/N/S, terrestrial planet parameters).
EARTH = {
    'name': 'Earth',
    'planet': {'M_mantle': 4.03e24, 'gravity': 9.81, 'radius': 6.371e6},
    'hcns': {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21},
}

# Mars: planetary parameters from standard literature values; inventory
# is the Earth BSE scaled by mass ratio (Mars / Earth = 0.107) to keep
# the comparison driven by planet structure rather than by a separate
# (and less well-constrained) Mars BSE estimate.
MARS_MASS_RATIO = 0.107
MARS = {
    'name': 'Mars-scaled',
    'planet': {'M_mantle': 5.03e23, 'gravity': 3.71, 'radius': 3.39e6},
    'hcns': {k: v * MARS_MASS_RATIO for k, v in EARTH['hcns'].items()},
}

T_MAGMA = 2500.0
DIW = 0.5
PHI = 1.0

SPECIES_TO_PLOT = ['H2O', 'CO2', 'H2', 'CO', 'CH4',
                   'N2', 'NH3', 'S2', 'SO2', 'H2S']


def _solve(cfg: dict) -> dict:
    ddict = {**cfg['planet'], 'T_magma': T_MAGMA, 'Phi_global': PHI,
             'fO2_shift_IW': DIW}
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = equilibrium_atmosphere(cfg['hcns'], ddict, print_result=False)
    return {
        'pressures': {sp: float(r[f'{sp}_bar']) for sp in SPECIES_TO_PLOT},
        'P_surf_bar': float(r['P_surf']),
        'M_atm_kg': float(r['M_atm']),
        'mean_mol_mass': float(r['atm_kg_per_mol']) * 1e3,
    }


def collect() -> dict:
    earth_out = _solve(EARTH)
    mars_out = _solve(MARS)
    log.info('Earth: P_surf = %.0f bar, M_atm = %.2e kg', earth_out['P_surf_bar'], earth_out['M_atm_kg'])
    log.info('Mars : P_surf = %.0f bar, M_atm = %.2e kg', mars_out['P_surf_bar'], mars_out['M_atm_kg'])
    return {'Earth': earth_out, 'Mars': mars_out}


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'mars_fiducial.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['planet'] + SPECIES_TO_PLOT + ['P_surf_bar', 'M_atm_kg', 'mean_mol_mass_g_per_mol'])
        for planet in ('Earth', 'Mars'):
            row = [planet] + [data[planet]['pressures'][sp] for sp in SPECIES_TO_PLOT]
            row += [data[planet]['P_surf_bar'],
                    data[planet]['M_atm_kg'],
                    data[planet]['mean_mol_mass']]
            w.writerow(row)
    log.info('Wrote %s', csv_path)

    # Grouped horizontal bars: for each species, an Earth bar above
    # a Mars bar at the same y. Species sorted by Earth pressure
    # descending so the most-abundant species sit at the top.
    order = sorted(
        range(len(SPECIES_TO_PLOT)),
        key=lambda i: -data['Earth']['pressures'][SPECIES_TO_PLOT[i]],
    )
    species = [SPECIES_TO_PLOT[i] for i in order]
    earth_p = np.array([data['Earth']['pressures'][sp] for sp in species])
    mars_p = np.array([data['Mars']['pressures'][sp] for sp in species])

    fig, ax = plt.subplots(figsize=(8.0, 5.4))

    y = np.arange(len(species)) * 1.8       # extra space between species
    bar_h = 0.7
    ax.barh(y - bar_h / 2, earth_p, height=bar_h,
            color=[dict_colors[sp] for sp in species],
            edgecolor='black', linewidth=0.5, label='Earth-BSE')
    ax.barh(y + bar_h / 2, mars_p, height=bar_h,
            color=[dict_colors[sp] for sp in species], alpha=0.45,
            edgecolor='black', linewidth=0.5, hatch='///',
            label='Mars-scaled (0.107 x Earth inventory)')

    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels([species_label(sp) for sp in species])
    ax.invert_yaxis()
    ax.set_xlabel('Surface partial pressure (bar)')
    ax.set_title(
        rf'Earth vs Mars-scaled at $T_\mathrm{{magma}} = {T_MAGMA:.0f}$ K, '
        rf'$\Phi = {PHI:.0f}$, $\Delta\mathrm{{IW}} = {DIW:+.1f}$'
    )
    x_lo = 1e-10
    x_hi = max(earth_p.max(), mars_p.max()) * 5000.0
    ax.set_xlim(x_lo, x_hi)
    ax.grid(axis='x', which='both', alpha=0.3)
    # Legend below the plot so it does not crowd the top bars.
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2,
              frameon=False, fontsize=9.5)

    # Per-planet diagnostics in the upper-right of the data area.
    # Use plain-text scientific notation (Unicode superscripts) so the
    # values line up in monospace columns.
    def _fmt(val, unit):
        return sci_fmt_plain(val, unit=unit)
    e_p = _fmt(data["Earth"]["P_surf_bar"], 'bar')
    m_p = _fmt(data["Mars"]["P_surf_bar"], 'bar')
    e_m = _fmt(data["Earth"]["M_atm_kg"], 'kg')
    m_m = _fmt(data["Mars"]["M_atm_kg"], 'kg')
    e_w = f'{data["Earth"]["mean_mol_mass"]:.2f} g/mol'
    m_w = f'{data["Mars"]["mean_mol_mass"]:.2f} g/mol'
    summary = (
        f'             Earth                Mars\n'
        f'P_surf   {e_p:<18s}   {m_p}\n'
        f'M_atm    {e_m:<18s}   {m_m}\n'
        f'mean M   {e_w:<18s}   {m_w}'
    )
    ax.text(
        0.985, 0.97, summary,
        transform=ax.transAxes,
        fontsize=8.5, family='monospace',
        va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc'),
    )

    paths = save(fig, 'mars_fiducial')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

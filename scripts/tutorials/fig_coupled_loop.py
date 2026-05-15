"""Tutorial figure for the "Coupled-loop driver" page.

Simulates a magma ocean cooling from 3000 K to 1500 K in 25 steps,
calling CALLIOPE at each step with the previous result's partial
pressures as the warm-start guess. Plots partial pressures as a
function of decreasing T_magma so the reader can see how species
populate the cooling atmosphere.

This is the realistic warm-start chain pattern, not the unrelated
sweep done in firstrun.md Step 6.
"""

from __future__ import annotations

import csv
import logging
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

from calliope.constants import dict_colors, volatile_species
from calliope.solve import equilibrium_atmosphere

from ._style import DATA_DIR, apply_style, save

log = logging.getLogger('tutorials.coupled_loop')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
EARTH_HCNS = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}
DIW_FIXED = 0.5         # holding redox fixed: the focus is the cooling sequence
T_SEQUENCE = np.linspace(3000.0, 1500.0, 25)

SPECIES_TO_PLOT = ['H2O', 'CO2', 'H2', 'CO', 'CH4',
                   'N2', 'NH3', 'S2', 'SO2', 'H2S']


def collect() -> dict:
    """Run the cooling sequence and return per-step partial pressures
    plus wall-time."""
    n = T_SEQUENCE.size
    pressures = {sp: np.full(n, np.nan) for sp in SPECIES_TO_PLOT}
    P_total = np.full(n, np.nan)
    wall = np.zeros(n)

    p_guess = None  # cold start the first call
    for i, T in enumerate(T_SEQUENCE):
        ddict = {**PLANET, 'T_magma': float(T), 'Phi_global': 1.0,
                 'fO2_shift_IW': DIW_FIXED}
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 1
            ddict[f'{sp}_initial_bar'] = 0.0
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = equilibrium_atmosphere(
                EARTH_HCNS, ddict, p_guess=p_guess, print_result=False,
            )
        wall[i] = time.time() - t0
        for sp in SPECIES_TO_PLOT:
            pressures[sp][i] = float(res[f'{sp}_bar'])
        P_total[i] = float(res['P_surf'])
        # Carry the four primary pressures forward as the next guess.
        p_guess = {s: float(res[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
        log.info('  step %2d  T=%4.0f K  P_surf=%7.1f bar  %5.3f s',
                 i, T, P_total[i], wall[i])
    return dict(T=T_SEQUENCE.copy(), pressures=pressures,
                P_total=P_total, wall_s=wall)


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'coupled_loop.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['T_K', 'wall_s', 'P_total_bar'] + SPECIES_TO_PLOT)
        for i, T in enumerate(data['T']):
            row = [T, data['wall_s'][i], data['P_total'][i]] + [
                data['pressures'][sp][i] for sp in SPECIES_TO_PLOT
            ]
            w.writerow(row)
    log.info('Wrote %s', csv_path)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))

    # Plot every species; label only those that exceed 1e-4 bar at the
    # cold end (the others are clutter at the bottom of the log scale).
    visible_threshold = 1e-4
    for sp in SPECIES_TO_PLOT:
        ys = data['pressures'][sp]
        label = sp if np.nanmax(ys) > visible_threshold else None
        ax.plot(data['T'], ys, color=dict_colors[sp], linewidth=1.8,
                marker='o', markersize=3.5, markeredgecolor='none',
                label=label, alpha=0.95 if label else 0.55)

    ax.set_yscale('log')
    ax.set_xlabel(r'$T_\mathrm{magma}$ [K] (cooling $\rightarrow$)')
    ax.set_ylabel('Surface partial pressure (bar)')
    ax.set_title(
        f'Cooling sequence at fixed $\\Delta\\mathrm{{IW}} = {DIW_FIXED:+.1f}$, '
        f'Earth-BSE inventory, $\\Phi = 1$'
    )
    ax.invert_xaxis()  # cooling reads left-to-right
    ax.set_ylim(1e-6, 1e4)
    ax.grid(which='both', alpha=0.3)
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
              frameon=False, title='major\nspecies', title_fontsize=9.5)

    # Wall-time annotation: highlight the warm-start speed-up. First
    # call is cold; subsequent are warm.
    t_cold_ms = float(data['wall_s'][0]) * 1e3
    t_warm_ms = float(np.median(data['wall_s'][1:])) * 1e3
    t_total_s = float(data['wall_s'].sum())
    summary = (
        f'cold start (step 0): {t_cold_ms:6.1f} ms\n'
        f'warm steps median:   {t_warm_ms:6.1f} ms\n'
        f'total wall time:     {t_total_s:6.2f} s'
    )
    ax.text(
        0.02, 0.04, summary,
        transform=ax.transAxes,
        fontsize=9.0, family='monospace',
        va='bottom', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#cccccc'),
    )

    paths = save(fig, 'coupled_loop')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

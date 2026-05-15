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

from ._style import DATA_DIR, apply_style, save, species_label

log = logging.getLogger('tutorials.coupled_loop')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
EARTH_HCNS = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}
DIW_FIXED = 0.5         # holding redox fixed: the focus is the cooling sequence
T_SEQUENCE = np.linspace(3000.0, 1500.0, 25)
T_FREEZE = 1500.0       # temperature at which the magma-volume sweep runs
PHI_SEQUENCE = np.linspace(1.0, 0.5, 11)   # melt-fraction crystallisation step

SPECIES_TO_PLOT = ['H2O', 'CO2', 'H2', 'CO', 'CH4',
                   'N2', 'NH3', 'S2', 'SO2', 'H2S']


def _run(ddict_template: dict, schedule: list[tuple[float, float]],
         label: str) -> dict:
    """Run a sequence of (T_magma, Phi_global) steps with warm-start
    threading. Returns per-step partial pressures, surface pressure,
    and wall time.

    ``schedule`` is a list of ``(T, Phi)`` tuples. Both quantities can
    vary; only the chemistry knobs touched here change between steps.
    """
    n = len(schedule)
    pressures = {sp: np.full(n, np.nan) for sp in SPECIES_TO_PLOT}
    P_total = np.full(n, np.nan)
    wall = np.zeros(n)
    p_guess = None
    for i, (T, phi) in enumerate(schedule):
        ddict = {**ddict_template, 'T_magma': float(T), 'Phi_global': float(phi)}
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
        p_guess = {s: float(res[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
        log.info('  %s step %2d  T=%4.0f K  Phi=%.2f  P_surf=%7.1f bar  %5.3f s',
                 label, i, T, phi, P_total[i], wall[i])
    return dict(pressures=pressures, P_total=P_total, wall_s=wall)


def _base_ddict() -> dict:
    d = {**PLANET, 'fO2_shift_IW': DIW_FIXED}
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


def collect() -> dict:
    """Run a cooling sequence (Phase 1: T 3000 -> 1500 K at Phi = 1)
    followed by a magma-volume sweep (Phase 2: Phi 1.0 -> 0.5 at fixed
    T = 1500 K). Both phases use the same warm-start chain.
    """
    ddict = _base_ddict()

    phase1_schedule = [(float(T), 1.0) for T in T_SEQUENCE]
    log.info('Phase 1: cooling at Phi = 1')
    phase1 = _run(ddict, phase1_schedule, 'cool')

    phase2_schedule = [(T_FREEZE, float(phi)) for phi in PHI_SEQUENCE]
    log.info('Phase 2: crystallisation at T = %.0f K', T_FREEZE)
    phase2 = _run(ddict, phase2_schedule, 'cryst')

    return dict(
        T_cool=T_SEQUENCE.copy(),
        pressures_cool=phase1['pressures'],
        P_total_cool=phase1['P_total'],
        wall_cool=phase1['wall_s'],
        Phi_cryst=PHI_SEQUENCE.copy(),
        T_cryst=T_FREEZE,
        pressures_cryst=phase2['pressures'],
        P_total_cryst=phase2['P_total'],
        wall_cryst=phase2['wall_s'],
    )


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'coupled_loop.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['phase', 'step', 'T_K', 'Phi_global', 'wall_s',
                    'P_total_bar'] + SPECIES_TO_PLOT)
        for i, T in enumerate(data['T_cool']):
            row = ['cooling', i, T, 1.0, data['wall_cool'][i],
                   data['P_total_cool'][i]] + [
                data['pressures_cool'][sp][i] for sp in SPECIES_TO_PLOT
            ]
            w.writerow(row)
        for i, phi in enumerate(data['Phi_cryst']):
            row = ['crystallisation', i, data['T_cryst'], phi,
                   data['wall_cryst'][i], data['P_total_cryst'][i]] + [
                data['pressures_cryst'][sp][i] for sp in SPECIES_TO_PLOT
            ]
            w.writerow(row)
    log.info('Wrote %s', csv_path)

    fig, (ax_cool, ax_cryst) = plt.subplots(
        1, 2, figsize=(11.4, 5.0), sharey=True,
        gridspec_kw={'width_ratios': [1.6, 1.0], 'wspace': 0.08},
    )

    visible_threshold = 1e-4

    # Phase 1: cooling at Phi = 1
    for sp in SPECIES_TO_PLOT:
        ys = data['pressures_cool'][sp]
        ax_cool.plot(data['T_cool'], ys, color=dict_colors[sp], linewidth=1.8,
                     marker='o', markersize=3.5, markeredgecolor='none',
                     alpha=0.95 if np.nanmax(ys) > visible_threshold else 0.55)
    ax_cool.set_yscale('log')
    ax_cool.set_xlabel(r'$T_\mathrm{magma}$ [K] (cooling $\rightarrow$)')
    ax_cool.set_ylabel('Surface partial pressure (bar)')
    ax_cool.set_title(
        f'(a) cooling at $\\Phi = 1$, '
        f'$\\Delta\\mathrm{{IW}} = {DIW_FIXED:+.1f}$'
    )
    ax_cool.invert_xaxis()
    ax_cool.set_ylim(1e-6, 1e4)
    ax_cool.grid(which='both', alpha=0.3)

    # Phase 2: crystallisation at T = T_FREEZE
    for sp in SPECIES_TO_PLOT:
        ys = data['pressures_cryst'][sp]
        label = species_label(sp) if np.nanmax(ys) > visible_threshold else None
        ax_cryst.plot(data['Phi_cryst'], ys, color=dict_colors[sp], linewidth=1.8,
                      marker='s', markersize=3.5, markeredgecolor='none',
                      label=label, alpha=0.95 if label else 0.55)
    ax_cryst.set_xlabel(r'$\Phi_\mathrm{global}$ (crystallisation $\rightarrow$)')
    ax_cryst.set_title(
        f'(b) crystallisation at $T = {int(data["T_cryst"])}$ K'
    )
    ax_cryst.invert_xaxis()
    ax_cryst.grid(which='both', alpha=0.3)
    ax_cryst.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                    frameon=False, title='major\nspecies', title_fontsize=9.5)

    # Wall-time annotation for the cooling panel.
    t_cold_ms = float(data['wall_cool'][0]) * 1e3
    t_warm_ms = float(np.median(data['wall_cool'][1:])) * 1e3
    t_total_s = float(data['wall_cool'].sum() + data['wall_cryst'].sum())
    summary = (
        f'cold start (step 0): {t_cold_ms:6.1f} ms\n'
        f'warm steps median:   {t_warm_ms:6.1f} ms\n'
        f'total wall (both):   {t_total_s:6.2f} s'
    )
    ax_cool.text(
        0.02, 0.04, summary,
        transform=ax_cool.transAxes,
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

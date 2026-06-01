"""Tutorial figure for the "Two-mode round-trip" page.

Walks an input Delta-IW through CALLIOPE's buffered mode, takes the
total volatile O the chemistry produces, then feeds (H, C, N, S, O)
into the authoritative-O mode and recovers a Delta-IW. The figure
shows recovered vs input Delta-IW with the y = x reference line so the
closure within solver tolerance is visually trivial to verify.

Fixed T_magma = 2000 K, Phi = 1, Earth-BSE Krijt+2023 H/C/N/S.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

from calliope.constants import volatile_species
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

from ._style import COLOR_CAL, DATA_DIR, apply_style, save, sci_fmt

log = logging.getLogger('tutorials.two_modes')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
T_MAGMA = 2000.0
PHI_GLOBAL = 1.0

# Krijt+2023 PPVII BSE H/C/N/S budget in kg.
EARTH_HCNS = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}

DIW_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])


def _base_ddict() -> dict:
    ddict = {**PLANET}
    ddict['T_magma'] = T_MAGMA
    ddict['Phi_global'] = PHI_GLOBAL
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    return ddict


def collect() -> dict:
    """Run buffered -> authoritative-O for each input Delta-IW.

    The buffered leg threads p_guess forward across the Delta-IW grid
    so the Monte-Carlo restart never has to find the canonical basin
    from a cold start at high Delta-IW (where the carbon-rich BSE
    inventory has a documented spurious H2O-free basin that the cold
    solver lands in ~20% of the time).

    Returns
    -------
    dict
        Keys ``dIW_input``, ``dIW_recovered``, ``O_kg_total``.
    """
    recovered = np.full(DIW_GRID.size, np.nan)
    O_total = np.full(DIW_GRID.size, np.nan)
    p_guess_buf = None
    for i, diw in enumerate(DIW_GRID):
        # Step 1: buffered call at the input Delta-IW. Warm-start
        # from the previous grid point's converged pressures so the
        # solver stays in the canonical basin across the sweep.
        ddict = _base_ddict()
        ddict['fO2_shift_IW'] = float(diw)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            buf = equilibrium_atmosphere(
                EARTH_HCNS,
                ddict,
                p_guess=p_guess_buf,
                hide_warnings=True,
                print_result=False,
            )
        O_kg = float(buf['O_kg_total'])
        p_guess_buf = {s: float(buf[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}

        # Step 2: authoritative-O call with the buffered run's O budget.
        target = dict(EARTH_HCNS)
        target['O'] = O_kg
        ddict_auth = _base_ddict()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            auth = equilibrium_atmosphere_authoritative_O(
                target,
                ddict_auth,
                p_guess=p_guess_buf,
                fO2_hint=float(diw),
                hide_warnings=True,
                print_result=False,
            )
        recovered[i] = float(auth['fO2_shift_derived'])
        O_total[i] = O_kg
        log.info(
            'dIW_in = %+.2f, O_kg = %.3e, recovered = %+.4f, residual = %+.2e',
            diw,
            O_kg,
            recovered[i],
            recovered[i] - diw,
        )
    return dict(dIW_input=DIW_GRID.copy(), dIW_recovered=recovered, O_kg_total=O_total)


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'two_modes_round_trip.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['dIW_input', 'dIW_recovered', 'residual_dex', 'O_kg_total'])
        for d_in, d_rec, o_kg in zip(
            data['dIW_input'], data['dIW_recovered'], data['O_kg_total']
        ):
            w.writerow([d_in, d_rec, d_rec - d_in, o_kg])
    log.info('Wrote %s', csv_path)

    residuals = data['dIW_recovered'] - data['dIW_input']
    finite = np.isfinite(residuals)
    worst_abs = float(np.nanmax(np.abs(residuals[finite]))) if finite.any() else float('nan')

    fig, ax = plt.subplots(figsize=(6.4, 5.6))

    # y = x reference (i.e., perfect closure).
    lo = float(data['dIW_input'].min()) - 0.5
    hi = float(data['dIW_input'].max()) + 0.5
    ax.plot(
        [lo, hi], [lo, hi], color='k', alpha=0.4, linewidth=1.0, label='perfect closure (y = x)'
    )

    # The recovered points themselves.
    ax.scatter(
        data['dIW_input'],
        data['dIW_recovered'],
        s=85,
        color=COLOR_CAL,
        edgecolor='k',
        linewidth=0.6,
        zorder=3,
        label='CALLIOPE buffered $\\rightarrow$ authoritative-O',
    )

    # Annotate worst-case residual so the reader has a number to quote.
    ax.text(
        0.04,
        0.96,
        f'worst-case |recovered − input| = {sci_fmt(worst_abs, unit="dex")}',
        transform=ax.transAxes,
        fontsize=9.5,
        va='top',
        ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc'),
    )

    ax.set_xlabel(r'input $\Delta\mathrm{IW}$ (buffered mode) [dex]')
    ax.set_ylabel(r'recovered $\Delta\mathrm{IW}$ (authoritative-O mode) [dex]')
    ax.set_title(
        f'Two-mode round-trip at Earth-BSE, '
        f'$T_\\mathrm{{magma}} = {T_MAGMA:.0f}$ K, $\\Phi = {PHI_GLOBAL:.0f}$'
    )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal')
    ax.legend(loc='lower right', framealpha=0.9, edgecolor='none')

    paths = save(fig, 'two_modes_round_trip')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s %(message)s'
    )
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

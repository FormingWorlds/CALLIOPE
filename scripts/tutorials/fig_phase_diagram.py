"""Tutorial figure for the "Speciation phase diagram" page.

Sweeps (T_magma, Delta-IW) on a 15 x 12 grid at Earth-BSE Krijt+2023
H/C/N/S and identifies the single dominant volatile species at each
point. Output is a discrete-colour map of dominant species in
(T, Delta-IW) space.

Warm-start chain along Delta-IW at fixed T to keep the runtime under a
minute on a modern laptop.
"""

from __future__ import annotations

import csv
import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from calliope.constants import dict_colors, volatile_species
from calliope.solve import equilibrium_atmosphere

from ._style import DATA_DIR, apply_style, save

log = logging.getLogger('tutorials.phase_diagram')


PLANET = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
}
EARTH_HCNS = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}

T_GRID = np.linspace(1500.0, 3000.0, 15)
DIW_GRID = np.linspace(-4.0, 5.0, 12)

REPORTED_SPECIES = ['H2O', 'CO2', 'H2', 'CO', 'CH4', 'N2', 'NH3', 'S2', 'SO2', 'H2S']


def _base_ddict(T_magma: float, diw: float) -> dict:
    ddict = {**PLANET, 'T_magma': T_magma, 'Phi_global': 1.0,
             'fO2_shift_IW': diw}
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    return ddict


def collect() -> dict:
    """Return dominant-species index, max partial pressure, and total P
    on the (T, dIW) grid.
    """
    n_T = T_GRID.size
    n_d = DIW_GRID.size
    dominant_idx = np.full((n_T, n_d), -1, dtype=int)
    P_total = np.full((n_T, n_d), np.nan)
    p_max = np.full((n_T, n_d), np.nan)

    for iT, T in enumerate(T_GRID):
        # Warm-start along dIW. Sweep ascending dIW, reset p_guess on
        # the first call at this T, then thread the previous result
        # forward.
        p_guess = None
        for jd, diw in enumerate(DIW_GRID):
            ddict = _base_ddict(float(T), float(diw))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    res = equilibrium_atmosphere(
                        EARTH_HCNS, ddict, p_guess=p_guess, print_result=False,
                    )
                except Exception as exc:  # noqa: BLE001
                    log.warning('  T=%.0f dIW=%+.2f raised %s', T, diw, exc)
                    p_guess = None
                    continue
            ps = np.array([float(res[f'{sp}_bar']) for sp in REPORTED_SPECIES])
            dominant_idx[iT, jd] = int(np.argmax(ps))
            p_max[iT, jd] = float(ps[dominant_idx[iT, jd]])
            P_total[iT, jd] = float(res['P_surf'])
            p_guess = {s: float(res[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
            log.info('  T=%4.0f dIW=%+5.2f  dominant=%-4s  P_surf=%.2e bar',
                     T, diw, REPORTED_SPECIES[dominant_idx[iT, jd]], P_total[iT, jd])

    return dict(T=T_GRID, dIW=DIW_GRID,
                dominant_idx=dominant_idx,
                p_max=p_max,
                P_total=P_total)


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'phase_diagram.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['T_K', 'dIW_dex', 'dominant_species_index',
                    'dominant_species', 'p_max_bar', 'P_total_bar'])
        for iT, T in enumerate(data['T']):
            for jd, diw in enumerate(data['dIW']):
                idx = int(data['dominant_idx'][iT, jd])
                sp = REPORTED_SPECIES[idx] if idx >= 0 else 'failed'
                w.writerow([T, diw, idx, sp,
                            data['p_max'][iT, jd], data['P_total'][iT, jd]])
    log.info('Wrote %s', csv_path)

    # Map each dominant-species index seen in the data to its PROTEUS
    # colour. Build a discrete colormap over the seen indices so the
    # legend is compact (only the species that actually appear).
    seen = sorted({int(i) for i in data['dominant_idx'].ravel() if i >= 0})
    species_palette = [dict_colors[REPORTED_SPECIES[i]] for i in seen]
    cmap = ListedColormap(species_palette)
    remap = {orig: new for new, orig in enumerate(seen)}
    remapped = np.vectorize(lambda v: remap.get(int(v), -1))(data['dominant_idx'])

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    # pcolormesh with bin edges so each cell is a tile. Build edge
    # arrays by extending each axis by half the step on both sides.
    def edges(arr: np.ndarray) -> np.ndarray:
        step = arr[1] - arr[0]
        return np.concatenate([[arr[0] - 0.5 * step],
                               0.5 * (arr[:-1] + arr[1:]),
                               [arr[-1] + 0.5 * step]])

    t_edges = edges(data['T'])
    d_edges = edges(data['dIW'])
    # remapped is (n_T, n_d). pcolormesh wants (n_T+1, n_d+1) edges,
    # with the array oriented so that x = dIW (cols), y = T (rows).
    pcm = ax.pcolormesh(
        d_edges, t_edges, np.ma.masked_less(remapped, 0),
        cmap=cmap, vmin=-0.5, vmax=len(seen) - 0.5,
        shading='flat', edgecolors='white', linewidth=0.6,
    )
    pcm.set_clim(-0.5, len(seen) - 0.5)

    ax.set_xlabel(r'$\Delta\mathrm{IW}$ [dex]')
    ax.set_ylabel(r'$T_\mathrm{magma}$ [K]')
    ax.set_title('Speciation phase diagram at Earth-BSE, $\\Phi = 1$')

    # Discrete legend with one Patch per actually-occurring species.
    handles = [Patch(facecolor=species_palette[k],
                     edgecolor='k', linewidth=0.4,
                     label=REPORTED_SPECIES[seen[k]])
               for k in range(len(seen))]
    ax.legend(
        handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
        title='dominant\nspecies', frameon=False,
        title_fontsize=9.5,
    )

    paths = save(fig, 'phase_diagram')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

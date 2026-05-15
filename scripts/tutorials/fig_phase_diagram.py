"""Tutorial figure for the "Speciation phase diagram" page.

Sweeps (T_magma, Delta-IW) on a 15 x 12 grid at Earth-BSE Krijt+2023
H/C/N/S and identifies the four most abundant volatile species at
each point. Each cell of the figure is subdivided into a 2 x 2
quartet in reading order: top-left = rank 1 (dominant), top-right =
rank 2, bottom-left = rank 3, bottom-right = rank 4. A reader can
then see at a glance not just what is dominant but how the second
through fourth species shift across the (T, Delta-IW) plane.

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

from ._style import DATA_DIR, apply_style, save, species_label

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
    """Return per-cell ranked top-4 species, plus full partial-pressure
    array and total surface pressure on the (T, dIW) grid.
    """
    n_T = T_GRID.size
    n_d = DIW_GRID.size
    n_sp = len(REPORTED_SPECIES)
    pressures = np.full((n_T, n_d, n_sp), np.nan)
    P_total = np.full((n_T, n_d), np.nan)
    rank_idx = np.full((n_T, n_d, 4), -1, dtype=int)

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
            pressures[iT, jd] = ps
            order = np.argsort(ps)[::-1]
            rank_idx[iT, jd] = order[:4]
            P_total[iT, jd] = float(res['P_surf'])
            p_guess = {s: float(res[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
            top4 = ', '.join(REPORTED_SPECIES[i] for i in order[:4])
            log.info('  T=%4.0f dIW=%+5.2f  top-4=[%s]  P_surf=%.2e bar',
                     T, diw, top4, P_total[iT, jd])

    return dict(T=T_GRID, dIW=DIW_GRID,
                pressures=pressures,
                rank_idx=rank_idx,
                P_total=P_total)


def make_figure(data: dict | None = None) -> dict:
    apply_style()
    data = data or collect()

    csv_path = DATA_DIR / 'phase_diagram.csv'
    with csv_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['T_K', 'dIW_dex', 'rank1', 'rank2', 'rank3', 'rank4',
                    'P_total_bar'] +
                   [f'p_{sp}_bar' for sp in REPORTED_SPECIES])
        for iT, T in enumerate(data['T']):
            for jd, diw in enumerate(data['dIW']):
                idx4 = data['rank_idx'][iT, jd]
                names = [REPORTED_SPECIES[i] if i >= 0 else 'failed' for i in idx4]
                row = [T, diw] + names + [data['P_total'][iT, jd]]
                row += [data['pressures'][iT, jd, s] for s in range(len(REPORTED_SPECIES))]
                w.writerow(row)
    log.info('Wrote %s', csv_path)

    # Build the 2x2-expanded grid. For each original cell (iT, jd),
    # the four sub-cells in reading order map to the rank-1..rank-4
    # species index. Reading order in (data x, data y) where x = dIW
    # and y = T (with y increasing upward in the plot):
    #   rank 1 (top-left)     -> sub_T = 1 (upper),  sub_d = 0 (lower)
    #   rank 2 (top-right)    -> sub_T = 1 (upper),  sub_d = 1 (upper)
    #   rank 3 (bottom-left)  -> sub_T = 0 (lower),  sub_d = 0 (lower)
    #   rank 4 (bottom-right) -> sub_T = 0 (lower),  sub_d = 1 (upper)
    rank_idx = data['rank_idx']
    n_T, n_d = rank_idx.shape[:2]
    expanded = np.full((2 * n_T, 2 * n_d), -1, dtype=int)
    for iT in range(n_T):
        for jd in range(n_d):
            r1, r2, r3, r4 = rank_idx[iT, jd]
            expanded[2 * iT + 1, 2 * jd + 0] = r1
            expanded[2 * iT + 1, 2 * jd + 1] = r2
            expanded[2 * iT + 0, 2 * jd + 0] = r3
            expanded[2 * iT + 0, 2 * jd + 1] = r4

    seen = sorted({int(v) for v in expanded.ravel() if v >= 0})
    species_palette = [dict_colors[REPORTED_SPECIES[i]] for i in seen]
    cmap = ListedColormap(species_palette)
    remap = {orig: new for new, orig in enumerate(seen)}
    remapped = np.vectorize(lambda v: remap.get(int(v), -1))(expanded)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    fig.subplots_adjust(right=0.78)

    # Build sub-cell edges that quarter each original cell.
    def edges_doubled(arr: np.ndarray) -> np.ndarray:
        step = arr[1] - arr[0]
        outer = np.concatenate([[arr[0] - 0.5 * step],
                                0.5 * (arr[:-1] + arr[1:]),
                                [arr[-1] + 0.5 * step]])
        # Insert a midpoint between every consecutive pair of outer
        # edges so each original cell becomes two sub-cells.
        mids = 0.5 * (outer[:-1] + outer[1:])
        result = np.empty(2 * outer.size - 1)
        result[0::2] = outer
        result[1::2] = mids
        return result

    t_edges = edges_doubled(data['T'])
    d_edges = edges_doubled(data['dIW'])
    pcm = ax.pcolormesh(
        d_edges, t_edges, np.ma.masked_less(remapped, 0),
        cmap=cmap, vmin=-0.5, vmax=len(seen) - 0.5,
        shading='flat', edgecolors='white', linewidth=0.35,
    )
    pcm.set_clim(-0.5, len(seen) - 0.5)

    # Heavier gridlines between original cells (every other sub-edge).
    for x in d_edges[::2]:
        ax.axvline(x, color='white', linewidth=0.9, alpha=0.95, zorder=4)
    for y in t_edges[::2]:
        ax.axhline(y, color='white', linewidth=0.9, alpha=0.95, zorder=4)
    ax.set_xlim(d_edges[0], d_edges[-1])
    ax.set_ylim(t_edges[0], t_edges[-1])

    ax.set_xlabel(r'$\Delta\mathrm{IW}$ [dex]')
    ax.set_ylabel(r'$T_\mathrm{magma}$ [K]')
    ax.set_title('Top-4 species per cell at Earth-BSE, $\\Phi = 1$')

    # Discrete species legend, plus a small rank-layout inset legend
    # so the reader knows which corner is rank 1 vs rank 4.
    handles = [Patch(facecolor=species_palette[k],
                     edgecolor='k', linewidth=0.4,
                     label=species_label(REPORTED_SPECIES[seen[k]]))
               for k in range(len(seen))]
    sp_legend = ax.legend(
        handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0),
        title='species\n(any rank)', frameon=False,
        title_fontsize=9.5,
    )
    ax.add_artist(sp_legend)

    # Rank-layout inset: 2x2 mini-grid showing where rank 1 / 2 / 3 / 4
    # sit inside each cell. Anchored below the species legend in the
    # right-hand margin of the figure so it does not overlap the data.
    from matplotlib.patches import Rectangle
    inset = fig.add_axes([0.82, 0.18, 0.09, 0.12])
    inset.set_xlim(0, 2); inset.set_ylim(0, 2)
    inset.set_xticks([]); inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_edgecolor('#888888'); spine.set_linewidth(0.6)
    rank_positions = {1: (0, 1), 2: (1, 1), 3: (0, 0), 4: (1, 0)}
    for rank, (rx, ry) in rank_positions.items():
        inset.add_patch(Rectangle((rx, ry), 1, 1, facecolor='#f3f3f3',
                                  edgecolor='white', linewidth=1.0))
        inset.text(rx + 0.5, ry + 0.5, f'#{rank}',
                   ha='center', va='center', fontsize=10, color='#333333')
    inset.set_title('rank layout', fontsize=8.5, color='#555555', pad=2)

    paths = save(fig, 'phase_diagram')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

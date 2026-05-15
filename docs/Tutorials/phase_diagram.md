# Speciation phase diagram

The first-run tutorial swept one parameter ($\Delta\mathrm{IW}$) at fixed temperature. Real magma oceans cool and re-equilibrate, so the dominant atmospheric species can shift as both $T$ and $\Delta\mathrm{IW}$ evolve. This tutorial generalises the 1D sweep to a 2D grid and builds a *speciation phase diagram*: at every $(T, \Delta\mathrm{IW})$ point, which volatile species has the largest partial pressure?

By the end of it you will:

- have run CALLIOPE on a 2D parameter grid using a warm-started serpentine sweep so the wall-time stays manageable;
- have produced the dominant-species map on the front of this page;
- understand why Earth-BSE atmospheres are CO- or CO$_2$-dominated rather than H$_2$O-dominated in the magma-ocean regime.

You should already have completed the [First run](firstrun.md) tutorial so the `equilibrium_atmosphere` call signature is familiar.

## Step 1: set up the sweep

```python
import warnings

import numpy as np

from calliope.constants import volatile_species
from calliope.solve import equilibrium_atmosphere

planet = {'M_mantle': 4.03e24, 'gravity': 9.81, 'radius': 6.371e6}

# Earth-BSE H/C/N/S in kg (Krijt et al. 2023, PPVII Tables 1 and 2).
earth_hcns = {'H': 5.6e20, 'C': 3.1e21, 'N': 3.7e19, 'S': 1.0e21}

T_grid   = np.linspace(1500.0, 3000.0, 15)
diw_grid = np.linspace(  -4.0,    5.0, 12)

species_to_report = ['H2O', 'CO2', 'H2', 'CO', 'CH4',
                     'N2',  'NH3', 'S2', 'SO2', 'H2S']

dominant_idx = np.full((T_grid.size, diw_grid.size), -1, dtype=int)
P_total      = np.full((T_grid.size, diw_grid.size), np.nan)


def base_ddict(T, diw):
    d = {**planet, 'T_magma': T, 'Phi_global': 1.0, 'fO2_shift_IW': diw}
    for sp in volatile_species:
        d[f'{sp}_included']    = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d
```

## Step 2: warm-start along $\Delta\mathrm{IW}$ at each $T$

Each call from cold takes about a second because the Monte-Carlo restart has to find a basin; with a `p_guess` from the previous call, subsequent solves complete in milliseconds. The cleanest schedule is a serpentine sweep: walk $\Delta\mathrm{IW}$ at fixed $T$, reset between rows.

```python
for iT, T in enumerate(T_grid):
    p_guess = None                              # cold-start each row
    for jd, diw in enumerate(diw_grid):
        ddict = base_ddict(float(T), float(diw))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                res = equilibrium_atmosphere(
                    earth_hcns, ddict, p_guess=p_guess, print_result=False,
                )
            except Exception:
                p_guess = None                  # invalidate guess on failure
                continue
        pressures = np.array([res[f'{sp}_bar'] for sp in species_to_report])
        dominant_idx[iT, jd] = int(np.argmax(pressures))
        P_total[iT, jd]      = float(res['P_surf'])
        # carry the four primary partial pressures forward as the next guess
        p_guess = {s: float(res[f'{s}_bar']) for s in ('H2O', 'CO2', 'N2', 'S2')}
```

A 15 × 12 grid runs in about 30 to 60 seconds on a modern laptop. The grid resolution is a tradeoff: fewer points runs faster but smears the redox boundary; more points sharpens the boundary but adds linearly to the wall time.

!!! tip "If a cell fails to converge"
    The solver can land in a secondary basin at extreme conditions. The `except` clause above invalidates the warm start so the next call cold-starts from a fresh Monte-Carlo draw. Failed cells show up as masked (white) in the figure rather than corrupting the colour map.

## Step 3: plot the dominant-species map

```python
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from calliope.constants import dict_colors

# Map each "dominant species" index that actually appears in the grid
# to its PROTEUS colour. We build a compact palette so the legend
# only lists species that occur somewhere in the data.
seen = sorted({int(i) for i in dominant_idx.ravel() if i >= 0})
palette = [dict_colors[species_to_report[i]] for i in seen]
remap = {old: new for new, old in enumerate(seen)}
remapped = np.vectorize(lambda v: remap.get(int(v), -1))(dominant_idx)

# Cell-edge arrays for pcolormesh
def edges(arr):
    s = arr[1] - arr[0]
    return np.concatenate([[arr[0] - 0.5 * s],
                           0.5 * (arr[:-1] + arr[1:]),
                           [arr[-1] + 0.5 * s]])

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.pcolormesh(edges(diw_grid), edges(T_grid),
              np.ma.masked_less(remapped, 0),
              cmap=ListedColormap(palette),
              shading='flat', edgecolors='white', linewidth=0.6)

ax.set_xlabel(r'$\Delta$IW [dex]')
ax.set_ylabel(r'$T_\mathrm{magma}$ [K]')
ax.set_title('Speciation phase diagram at Earth-BSE')

handles = [Patch(facecolor=palette[k], edgecolor='k', linewidth=0.4,
                 label=species_to_report[seen[k]]) for k in range(len(seen))]
ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
          title='dominant\nspecies', frameon=False)
fig.tight_layout()
fig.savefig('phase_diagram.pdf')
```

## The goal of this tutorial

![Speciation phase diagram](../assets/figures/tutorials/phase_diagram.png)

*Dominant volatile species in $(T_\mathrm{magma}, \Delta\mathrm{IW})$ at the Earth-BSE Krijt+2023 H/C/N/S budget and $\Phi = 1$. Each cell is the species with the largest partial pressure at that $(T, \Delta\mathrm{IW})$ point. The redox boundary at $\Delta\mathrm{IW} \approx +2$ separates a CO-dominated reducing regime (left) from a CO$_2$-dominated oxidising regime (right). A small CH$_4$ patch appears at low $T$ and the most reducing conditions where methane synthesis becomes thermodynamically favourable.*

The result that may surprise a reader who thinks of magma-ocean atmospheres as "steam-dominated" is that on the *Earth* BSE inventory (with C / H mass ratio $\sim 5.5$), carbon, not hydrogen, sets the dominant species. H$_2$O is always there at the $\sim 0.1$ to $\sim 1$ bar level, but it is never the largest reservoir on this inventory. A water-dominated atmosphere needs either a much higher H budget (gas-giant-like) or a much lower C budget (volatile-poor / dehydrated body); compare the planetary case study tutorial for an example.

The phase boundary is largely temperature-independent above $\sim 1700$ K because both the equilibrium constants and the dominant solubility laws scale with $T$ in roughly the same way; the redox dependence dominates.

## Where to go next

- For a finer look at the speciation transitions across $\Delta\mathrm{IW}$ at one $T$, return to [First run](firstrun.md) Step 6 and use a denser $\Delta\mathrm{IW}$ grid with this tutorial's warm-start pattern.
- For how PROTEUS handles repeated CALLIOPE calls across a time-stepped simulation, read the next tutorial: [Coupled-loop driver](coupled_loop.md).
- For the chemistry behind the dominant-species shifts, read [Equilibrium chemistry](../Explanations/equilibrium_chemistry.md).

## Reproducing the figure

Generated by [`scripts/tutorials/fig_phase_diagram.py`](https://github.com/FormingWorlds/CALLIOPE/blob/main/scripts/tutorials/fig_phase_diagram.py). Re-run with `python -m scripts.tutorials.fig_phase_diagram` from the repository root.

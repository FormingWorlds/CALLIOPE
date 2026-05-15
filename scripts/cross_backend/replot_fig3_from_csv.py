"""Re-render Fig 3 from `data/fig3_grid.csv` without rerunning solvers.

Used when the only change is in the plot layer (e.g. a sign fix in the
buffer correction) and re-running the 15-point grid would be wasteful.
"""

from __future__ import annotations

import csv
import sys

import numpy as np

from .fig3_grid import make_figure
from .plot_style import DATA_DIR


def load_csv() -> dict:
    csv_path = DATA_DIR / 'fig3_grid.csv'
    if not csv_path.exists():
        print(f'fig3_grid.csv not found at {csv_path}; run fig3_grid first',
              file=sys.stderr)
        return None

    # Reconstruct the (nT, nO) arrays expected by make_figure.
    rows = []
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)

    T_vals = sorted({float(r['T_K']) for r in rows})
    O_vals = sorted({float(r['O_factor_x_Earth']) for r in rows})
    nT, nO = len(T_vals), len(O_vals)

    dIW_cal = np.full((nT, nO), np.nan)
    dIW_atm = np.full((nT, nO), np.nan)
    conv = np.zeros((nT, nO), dtype=bool)
    P_cal = np.full((nT, nO), np.nan)
    P_atm = np.full((nT, nO), np.nan)

    for r in rows:
        i = T_vals.index(float(r['T_K']))
        j = O_vals.index(float(r['O_factor_x_Earth']))
        for arr, key in (
            (dIW_cal, 'dIW_calliope'),
            (dIW_atm, 'dIW_atmodeller'),
            (P_cal, 'P_total_cal_bar'),
            (P_atm, 'P_total_atm_bar'),
        ):
            v = r[key]
            arr[i, j] = float(v) if v and v != 'nan' else np.nan
        conv[i, j] = (
            not np.isnan(dIW_cal[i, j]) and not np.isnan(dIW_atm[i, j])
        )

    return dict(
        dIW_cal=dIW_cal,
        dIW_atm=dIW_atm,
        conv_cal=~np.isnan(dIW_cal),
        conv_atm=~np.isnan(dIW_atm),
        P_cal=P_cal,
        P_atm=P_atm,
    )


def main() -> None:
    data = load_csv()
    if data is None:
        sys.exit(1)
    paths = make_figure(data=data)
    for ext, path in paths.items():
        print(f'  {ext}: {path}')


if __name__ == '__main__':
    main()

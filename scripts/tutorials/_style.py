"""Shared style for tutorial figures.

Imports `apply_style` from the cross-backend plotting scripts so the
Roboto fonts bundled there are registered as a side effect, then
exports a `save()` helper that writes into the tutorial-specific
output directory `docs/assets/figures/tutorials/`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from scripts.cross_backend.plot_style import (  # noqa: F401  re-exported for convenience
    COLOR_ATM,
    COLOR_BG,
    COLOR_CAL,
    COLOR_FIS,
    apply_style,
    panel_label,
)

FIGURE_DIR = (
    Path(__file__).resolve().parent.parent.parent
    / 'docs'
    / 'assets'
    / 'figures'
    / 'tutorials'
)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = Path(__file__).resolve().parent / 'data'
DATA_DIR.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, stem: str, *, formats=('pdf', 'png')) -> dict:
    """Save `fig` to `FIGURE_DIR/<stem>.<ext>` for each requested format."""
    paths = {}
    for ext in formats:
        path = FIGURE_DIR / f'{stem}.{ext}'
        fig.savefig(path, bbox_inches='tight')
        paths[ext] = path
    return paths

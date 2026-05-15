"""Shared matplotlib styling for publication-quality figures.

All figures use the same fonts, line widths, and colour palette so the
docs page reads as one coherent set of figures, not five unrelated
matplotlib outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

FIGURE_DIR = Path(__file__).resolve().parent.parent.parent / 'docs' / 'assets' / 'figures' / 'cross_backend'
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = Path(__file__).resolve().parent / 'data'
DATA_DIR.mkdir(parents=True, exist_ok=True)

COLOR_CAL = '#0f6e9d'    # CALLIOPE: deep blue
COLOR_ATM = '#d1471f'    # atmodeller: rust orange
COLOR_FIS = '#7a7a7a'    # Fischer alternative: grey
COLOR_BG = '#f3f0e7'     # subtle anchor-range background


def apply_style() -> None:
    """Install matplotlib rcParams. Idempotent."""
    mpl.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'font.family': 'serif',
        'font.size': 10.5,
        'axes.titlesize': 11.5,
        'axes.labelsize': 11,
        'axes.linewidth': 0.9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.25,
        'grid.linewidth': 0.6,
        'legend.fontsize': 9.5,
        'legend.frameon': False,
        'lines.linewidth': 1.7,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
    })


def save(fig: plt.Figure, stem: str, *, formats=('pdf', 'png')) -> dict:
    """Save `fig` to FIGURE_DIR/`stem`.<ext> for each format.

    Returns mapping {ext: path}.
    """
    paths = {}
    for ext in formats:
        path = FIGURE_DIR / f'{stem}.{ext}'
        fig.savefig(path, bbox_inches='tight')
        paths[ext] = path
    return paths


def panel_label(ax, text, *, x=0.02, y=0.95):
    """Add a bold (a)/(b)/(c) panel label in figure coordinates."""
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        fontsize=12, fontweight='bold',
        va='top', ha='left',
    )

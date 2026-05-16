"""Shared style for tutorial figures.

Imports `apply_style` from the cross-backend plotting scripts so the
Roboto fonts bundled there are registered as a side effect, then
exports a `save()` helper that writes into the tutorial-specific
output directory `docs/assets/figures/tutorials/`.
"""

from __future__ import annotations

import math
import re
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

_SPECIES_RE = re.compile(r'(\d+)')


def species_label(species: str) -> str:
    """Render a species name in matplotlib mathtext with the digits as
    subscripts: ``species_label('CH4') -> r'$\\mathrm{CH_{4}}$'``.

    Use this anywhere a plain ASCII species name (legend, tick label,
    annotation) would otherwise read as ``CH4`` instead of ``CH4``.
    """
    return r'$\mathrm{' + _SPECIES_RE.sub(r'_{\1}', species) + '}$'


def _decimal(x: float, prec: int) -> str:
    """Render ``x`` as a plain decimal string with ``prec`` digits
    beyond the leading non-zero digit. Never falls back to ``e+``.
    """
    exp = int(math.floor(math.log10(abs(x))))
    decimals = max(0, prec - exp)
    return f'{x:.{decimals}f}'


def sci_fmt(x: float, *, prec: int = 2, unit: str = '') -> str:
    """Format ``x`` for use as a plot annotation.

    Numbers in [0.01, 10000) print with ``prec + 1`` significant
    figures as plain decimals (e.g. ``1712`` stays as ``1712``, not
    ``1.71e+03``). Everything outside that range prints in
    matplotlib mathtext scientific notation
    (``$a \\times 10^{b}$``) so labels never show ``e+19`` or
    ``e-09``. The trailing ``unit`` string is appended with one space.
    """
    suffix = f' {unit}' if unit else ''
    if x == 0 or not math.isfinite(x):
        return f'0{suffix}'
    abs_x = abs(x)
    if 1e-2 <= abs_x < 1e4:
        return f'{_decimal(x, prec)}{suffix}'
    exp = int(math.floor(math.log10(abs_x)))
    mant = x / (10**exp)
    return rf'${mant:.{prec}f} \times 10^{{{exp}}}${suffix}'


_UNICODE_SUP = str.maketrans('0123456789-+', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺')


def sci_fmt_plain(x: float, *, prec: int = 2, unit: str = '') -> str:
    """Like :func:`sci_fmt` but emits Unicode superscripts instead of
    matplotlib mathtext.

    Use this inside monospace text blocks (e.g. summary tables) where
    mathtext would break column alignment. The output stays a plain
    ASCII / Unicode string with ``×`` and superscript digits, e.g.
    ``4.39 × 10¹⁹ kg``.
    """
    suffix = f' {unit}' if unit else ''
    if x == 0 or not math.isfinite(x):
        return f'0{suffix}'
    abs_x = abs(x)
    if 1e-2 <= abs_x < 1e4:
        return f'{_decimal(x, prec)}{suffix}'
    exp = int(math.floor(math.log10(abs_x)))
    mant = x / (10**exp)
    return f'{mant:.{prec}f} × 10{str(exp).translate(_UNICODE_SUP)}{suffix}'


FIGURE_DIR = (
    Path(__file__).resolve().parent.parent.parent / 'docs' / 'assets' / 'figures' / 'tutorials'
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

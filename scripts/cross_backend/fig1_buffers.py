"""Figure 1: The IW buffer divergence.

Top panel: log10 fO2 at the IW buffer vs temperature for the three
parameterisations used in the FWL ecosystem (O'Neill & Eggins 2002 the
CALLIOPE default; Fischer et al. 2011 the CALLIOPE alternative;
Hirschmann composite the atmodeller default).

Bottom panel: difference between Hirschmann and the two CALLIOPE
parameterisations, in dex.

Pure analytical evaluation: no chemistry solver, no random seed.
"""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np

from . import buffers
from .plot_style import COLOR_ATM, COLOR_CAL, COLOR_FIS, apply_style, panel_label, save

log = logging.getLogger('cross_backend.fig1')


def make_figure() -> dict:
    apply_style()

    T = np.linspace(800.0, 3500.0, 271)
    f_oneill = buffers.oneill(T)
    f_fischer = buffers.fischer(T)
    f_hirsch = buffers.hirschmann_composite(T)
    d_hirsch_oneill = f_hirsch - f_oneill
    d_hirsch_fischer = f_hirsch - f_fischer

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(6.4, 6.4), sharex=True,
        gridspec_kw={'height_ratios': [2.0, 1.0], 'hspace': 0.12},
    )

    ax_top.plot(T, f_fischer, color=COLOR_CAL, label='Fischer et al. 2011 (CALLIOPE default)')
    ax_top.plot(T, f_oneill, color=COLOR_FIS, linestyle='--', label="O'Neill & Eggins 2002 (CALLIOPE legacy)")
    ax_top.plot(T, f_hirsch, color=COLOR_ATM, label='Hirschmann composite (atmodeller default)')

    # Annotate the Hirschmann composite switchover. Anchor inside the
    # data range with a y just above the curves, not at the panel
    # edge, so the text never gets clipped against the frame.
    ax_top.axvline(1000.0, color='k', alpha=0.25, linestyle=':')
    y_lo, y_hi = ax_top.get_ylim()
    ax_top.text(1010, y_lo + 0.85 * (y_hi - y_lo),
                'H08 / H21\nswitchover',
                fontsize=8.5, va='top', ha='left', alpha=0.6)

    ax_top.set_ylabel(r'$\log_{10} f_{\mathrm{O}_2}$ at IW buffer')
    ax_top.set_title('Iron-wustite buffer parameterisations across magma-ocean temperatures')
    ax_top.legend(loc='lower right', framealpha=0.0)
    panel_label(ax_top, '(a)')

    ax_bot.axhline(0.0, color='k', alpha=0.4, linewidth=0.7)
    ax_bot.plot(T, d_hirsch_oneill, color=COLOR_ATM, label="Hirschmann − O'Neill")
    ax_bot.plot(T, d_hirsch_fischer, color=COLOR_FIS, linestyle='--', label='Hirschmann − Fischer')
    ax_bot.axvline(1000.0, color='k', alpha=0.25, linestyle=':')

    ax_bot.set_xlabel(r'$T$ [K]')
    ax_bot.set_ylabel(r'$\Delta\log_{10} f_{\mathrm{O}_2}$ [dex]')
    ax_bot.legend(loc='lower right')
    panel_label(ax_bot, '(b)')

    # Mark a few characteristic magma-ocean temperatures with vertical guides.
    for T_mark, label in [(1500, '1500 K'), (2000, '2000 K'), (3000, '3000 K')]:
        ax_bot.axvline(T_mark, color='k', alpha=0.1, linewidth=0.7)
        ax_top.axvline(T_mark, color='k', alpha=0.1, linewidth=0.7)

    paths = save(fig, 'fig1_buffer_divergence')
    plt.close(fig)
    return paths


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    out = make_figure()
    for ext, path in out.items():
        print(f'  {ext}: {path}')

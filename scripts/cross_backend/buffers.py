"""Analytical IW (iron-wustite) buffer formulae.

The buffer divergence is the single largest cross-backend systematic
and can be evaluated without any chemistry solver, so the harness
computes it from the published formulae directly.

- O'Neill & Eggins (2002): monolithic fit, used by CALLIOPE by default.
- Fischer et al. (2011): 1-bar reference, CALLIOPE alternative.
- Hirschmann composite: Hirschmann (2008) below 1000 K, Hirschmann
  (2021) above 1000 K. Used by atmodeller by default. The H21 branch is
  a multi-coefficient Saxena-style polynomial; we delegate to
  atmodeller's evaluator (`IronWustiteBufferHirschmann`) rather than
  reimplement it, so any future correction in atmodeller propagates
  here automatically.

All formulae return log10 fO2_IW at the buffer, in bar. Pressure is in
bar.
"""

from __future__ import annotations

import numpy as np


def oneill(T: np.ndarray) -> np.ndarray:
    """O'Neill & Eggins (2002) IW buffer, monolithic fit.

    Reproduces the CALLIOPE implementation in
    `calliope.oxygen_fugacity.oneill`. The 8.31441 J/mol/K gas-constant
    value reproduces O'Neill & Eggins Eq. 11 verbatim; do not substitute
    a more modern CODATA value.
    """
    T = np.asarray(T, dtype=float)
    R = 8.31441
    return 2.0 * (-244118.0 + 115.559 * T - 8.474 * T * np.log(T)) / (np.log(10.0) * R * T)


def fischer(T: np.ndarray) -> np.ndarray:
    """Fischer et al. (2011) IW buffer, 1-bar reference."""
    T = np.asarray(T, dtype=float)
    return 6.94059 - 28.1808e3 / T


def hirschmann_composite(T: np.ndarray, P_bar: float = 1.0) -> np.ndarray:
    """Hirschmann composite IW buffer: H08 below 1000 K, H21 above.

    Delegates to atmodeller's `IronWustiteBufferHirschmann` evaluator at
    `evaluation_pressure = P_bar`. Returns a numpy array.
    """
    from atmodeller.thermodata import IronWustiteBuffer

    buf = IronWustiteBuffer(log10_shift=0.0, evaluation_pressure=P_bar)
    T_arr = np.atleast_1d(np.asarray(T, dtype=float))
    out = np.array([float(buf.log10_fugacity_buffer(t, P_bar)) for t in T_arr])
    return out if T_arr.shape == np.asarray(T).shape else out


def hirschmann_minus_oneill_offset(T: np.ndarray, P_bar: float = 1.0) -> np.ndarray:
    """log10(fO2_Hirschmann) - log10(fO2_ONeill) at the buffer (no shift).

    This is the analytical correction that lets you subtract the buffer-
    convention contribution from a cross-backend Delta-IW comparison.
    """
    return hirschmann_composite(T, P_bar) - oneill(T)

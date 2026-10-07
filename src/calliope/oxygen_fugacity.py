# fO2 buffers
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger('fwl.' + __name__)

# Single source of truth for the default fO2 buffer. ``OxygenFugacity`` and
# ``chemistry.ModifiedKeq`` both default to this name so the two cannot drift
# out of step. Fischer is the more recent IW parameterisation and tracks the
# atmodeller Hirschmann composite to within ~0.2 dex across magma-ocean
# temperatures (see docs/Explanations/cross_backend_comparison.md).
DEFAULT_FO2_MODEL = 'fischer'

# Hirschmann (2021, GCA 313, 74) Table 1. Each parameter is
# m = m0 + m1*P + m2*P**2 + m3*P**3 + m4*P**0.5 with P in GPa, and
# log10 fO2 = a + b*T + c*T*ln(T) + d/T. Rows are (m0, m1, m2, m3, m4).
_H21_FCC_BCC = (
    (6.844864, 1.175691e-01, 1.143873e-03, 0.0, 0.0),  # a
    (5.791364e-04, -2.891434e-04, -2.737171e-07, 0.0, 0.0),  # b
    (-7.971469e-05, 3.198005e-05, 0.0, 1.059554e-10, 2.014461e-07),  # c
    (-2.769002e04, 5.285977e02, -2.919275e00, 0.0, 0.0),  # d
)
_H21_HCP = (
    (8.463095, -3.000307e-03, 7.213445e-05, 0.0, 0.0),  # e
    (1.148738e-03, -9.352312e-05, 5.161592e-07, 0.0, 0.0),  # f
    (-7.448624e-04, -6.329325e-06, 0.0, -1.407339e-10, 1.830014e-04),  # g
    (-2.782082e04, 5.285977e02, -8.473231e-01, 0.0, 0.0),  # h
)
# fcc-hcp iron boundary, Hirschmann (2021) Eq. 18: P_GPa = x0 + x1*T + x2*T**2
_H21_FCC_HCP = (-18.64, 0.04359, -5.069e-06)
_BAR_TO_GPA = 1e-4


class OxygenFugacity:
    """log10 oxygen fugacity as a function of temperature"""

    def __init__(self, model=DEFAULT_FO2_MODEL):
        self.callmodel = getattr(self, model)

    def __call__(self, T, fO2_shift=0):
        """Return log10 fO2"""
        if T <= 0:
            raise ValueError(
                f'Temperature must be positive (K), got T={T}. '
                'The IW formulae diverge at T<=0 via 1/T and T*log(T) terms.'
            )
        return self.callmodel(T) + fO2_shift

    def fischer(self, T):
        """Fischer et al. (2011) IW (FeO equation of state, EPSL 304, 496).

        The coefficients are the p = 1 bar isoline of their Figure 6,
        obtained by integrating their Equation 2.
        """
        return 6.94059 - 28.1808 * 1e3 / T

    def oneill(self, T):
        """O'Neill and Eggins (2002) IW"""
        # 8.31441 reproduces O'Neill & Eggins (2002) Eq. 11 verbatim;
        # do not replace with constants.R_gas (8.31446...).
        return 2 * (-244118 + 115.559 * T - 8.474 * T * np.log(T)) / (np.log(10) * 8.31441 * T)

    def hirschmann(self, T, P_bar=1.0):
        """Hirschmann (2021) IW (GCA 313, 74), empirical fit of Table 1.

        Accounts for variable wustite stoichiometry and for pressure.
        Calibrated over 1000-3000 K and 1 bar to 100 GPa; the paper
        advises against extrapolation below 1000 K or above 100 GPa.
        The dispatcher evaluates it at the 1 bar reference, the same
        reference as the other buffers, so ``fO2_shift`` keeps its
        meaning. Call this method directly with ``P_bar`` for the
        buffer at pressure. The fcc/bcc branch applies below the
        fcc-hcp iron boundary (Eq. 18), the hcp branch above it.
        """
        T = np.asarray(T, dtype=float)
        P = np.asarray(P_bar, dtype=float) * _BAR_TO_GPA
        if np.any(P < 0):
            raise ValueError(f'Pressure must be non-negative (bar), got P_bar={P_bar}.')

        def evaluate(coeffs):
            a, b, c, d = (
                m[0] + m[1] * P + m[2] * P**2 + m[3] * P**3 + m[4] * np.sqrt(P) for m in coeffs
            )
            return a + b * T + c * T * np.log(T) + d / T

        x0, x1, x2 = _H21_FCC_HCP
        hcp = P > x0 + x1 * T + x2 * T**2
        return np.where(hcp, evaluate(_H21_HCP), evaluate(_H21_FCC_BCC))[()]

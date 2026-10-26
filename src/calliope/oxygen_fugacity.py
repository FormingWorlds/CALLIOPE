# fO2 buffers
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger('fwl.' + __name__)

# Default fO2 buffer of both OxygenFugacity and chemistry.ModifiedKeq. Hirschmann (2021)
# accounts for non-stoichiometric wustite and is atmodeller's IW buffer above 1000 K
# (docs/Explanations/cross_backend_comparison.md).
DEFAULT_FO2_MODEL = 'hirschmann'

# Names accepted for the IW buffer, one per OxygenFugacity model method.
FO2_MODELS = ('hirschmann', 'fischer', 'oneill')

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
# Eq. 18 is evaluated at T held to [1000 K (calibration floor), vertex (~4300 K)]. Unheld,
# the parabola falls below 1 bar under ~451 K and above ~8148 K and would put the 1 bar
# buffer on the hcp branch, ~42 dex off at 8150 K.
_H21_FCC_HCP_T_RANGE = (1000.0, -_H21_FCC_HCP[1] / (2.0 * _H21_FCC_HCP[2]))
_BAR_TO_GPA = 1e-4
# The boundary stays above 19.8 GPa over that range, so at 1 bar the buffer is always
# the fcc/bcc branch; its (a, b, c, d) at 1 bar serve the dispatcher.
_H21_1BAR = tuple(
    m[0]
    + m[1] * _BAR_TO_GPA
    + m[2] * _BAR_TO_GPA**2
    + m[3] * _BAR_TO_GPA**3
    + m[4] * np.sqrt(_BAR_TO_GPA)
    for m in _H21_FCC_BCC
)


class OxygenFugacity:
    """log10 oxygen fugacity as a function of temperature"""

    def __init__(self, model=DEFAULT_FO2_MODEL):
        if model not in FO2_MODELS:
            raise ValueError(f'Unknown IW buffer {model!r}; expected one of {FO2_MODELS}.')
        self.callmodel = (
            self._hirschmann_1bar if model == 'hirschmann' else getattr(self, model)
        )

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

        Fischer et al. (2011) calculate the 1-bar Fe-FeO curve from Chase (1998)
        (Fig. 6) and add the pressure integral of their Eq. 2; their supplement
        (Table S2) gives a P-T polynomial. The source of this 2-term form is not
        stated in the paper body.
        """
        return 6.94059 - 28.1808 * 1e3 / T

    def oneill(self, T):
        """O'Neill and Eggins (2002) IW"""
        # DrG of Fe + 0.5 O2 = FeO(liq): O'Neill & Eggins (2002), Chem. Geol. 186,
        # Table 6, p. 163. Keep the coded gas constant 8.31441 J/mol/K; a change
        # moves every IW value.
        return 2 * (-244118 + 115.559 * T - 8.474 * T * np.log(T)) / (np.log(10) * 8.31441 * T)

    def _hirschmann_1bar(self, T):
        """``hirschmann(T)`` at 1 bar, without the pressure and branch evaluation."""
        a, b, c, d = _H21_1BAR
        return a + b * T + c * T * np.log(T) + d / T

    def hirschmann(self, T, P_bar=1.0):
        """Hirschmann (2021) IW (GCA 313, 74), empirical fit of Table 1.

        Accounts for variable wustite stoichiometry and for pressure.
        Calibrated over 1000-3000 K and 1 bar to 100 GPa; the paper
        advises against extrapolation below 1000 K or above 100 GPa.
        The dispatcher evaluates it at the 1 bar reference, the same
        reference as the other buffers, so ``fO2_shift`` keeps its
        meaning. Call this method directly with ``P_bar`` for the
        buffer at pressure. The fcc/bcc branch applies below the
        fcc-hcp iron boundary (Eq. 18), the hcp branch above it. The
        boundary is evaluated at T held to 1000 K - 4300 K, so it never
        drops to 1 bar; outside that range it stays at its end value.
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
        T_b = np.clip(T, *_H21_FCC_HCP_T_RANGE)
        hcp = P > x0 + x1 * T_b + x2 * T_b**2
        return np.where(hcp, evaluate(_H21_HCP), evaluate(_H21_FCC_BCC))[()]

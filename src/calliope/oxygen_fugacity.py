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


class OxygenFugacity:
    """log10 oxygen fugacity as a function of temperature"""

    def __init__(self, model=DEFAULT_FO2_MODEL):
        self.callmodel = getattr(self, model)

    def __call__(self, T, fO2_shift=0):
        """Return log10 fO2"""
        if T <= 0:
            raise ValueError(
                f'Temperature must be positive (K), got T={T}. '
                'Both IW formulae diverge at T<=0 via 1/T and T*log(T) terms.'
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

# Equilibrium chemistry
from __future__ import annotations

import logging

from .oxygen_fugacity import DEFAULT_FO2_MODEL, OxygenFugacity

log = logging.getLogger('fwl.' + __name__)


class ModifiedKeq:
    """Modified equilibrium constant (includes fO2)"""

    def __init__(self, Keq_model, fO2_model=DEFAULT_FO2_MODEL):
        self.fO2 = OxygenFugacity(fO2_model)
        self.callmodel = getattr(self, Keq_model)

    def __call__(self, T, fO2_shift):
        fO2 = self.fO2(T, fO2_shift)
        Keq, fO2_stoich = self.callmodel(T)
        Geq = 10 ** (Keq - fO2_stoich * fO2)
        return Geq

    def schaefer_CH4(self, T):
        """Schaefer log10Keq for CO2 + 2H2 = CH4 + fO2"""
        # second argument returns stoichiometry of O2
        return (-16276 / T - 5.4738, 1)

    def schaefer_C(self, T):
        """Schaefer log10Keq for CO2 = CO + 0.5 fO2"""
        return (-14787 / T + 4.5472, 0.5)

    def schaefer_H(self, T):
        """Schaefer log10Keq for H2O = H2 + 0.5 fO2"""
        return (-12794 / T + 2.7768, 0.5)

    def janaf_CO(self, T):
        """JANAF log10Keq for CO2 = CO + 0.5 fO2 (CALLIOPE use range 1500-3000 K)"""
        return (-14467.511400133637 / T + 4.348135473316284, 0.5)

    def janaf_H2(self, T):
        """JANAF log10Keq for H2O = H2 + 0.5 fO2 (CALLIOPE use range 1500-3000 K)"""
        return (-13152.477779978302 / T + 3.038586383273608, 0.5)

    def janaf_SO2(self, T):
        # JANAF log10Keq for S2 + 2 O2 = 2 SO2 (doubled form)
        # Coefficients = 2x the formation constant for 0.5 S2 + O2 = SO2
        # fO2_stoich=0: O2 dependence handled by p_O2^2 in the sqrt expression
        # https://doi.org/10.1016/j.gca.2022.08.032
        return (37774.0 / T - 7.6128, 0)

    def janaf_H2S(self, T):
        # JANAF log10Keq for S2 + 2 H2 = 2 H2S (doubled form)
        # Coefficients = 2x the formation constant for 0.5 S2 + H2 = H2S
        # See notebook in `tools/`
        return (13462.03094 / T - 7.24546062, 0)

    def janaf_NH3(self, T):
        # JANAF log10Keq for N2 + 3 H2 = 2 NH3 (doubled form)
        # Coefficients = 2x the formation constant
        # See notebook in `tools/`
        return (5328.031246 / T - 11.98476092, 0)

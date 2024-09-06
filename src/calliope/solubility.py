# Solubility laws
from __future__ import annotations

import logging

import numpy as np

from .oxygen_fugacity import OxygenFugacity

log = logging.getLogger("fwl."+__name__)

class Solubility:
    """Solubility base class.  All p in bar"""

    def __init__(self, composition):
        self.callmodel = getattr(self, composition)

    def power_law(self, p, const, exponent):
        return const*p**exponent

    def __call__(self, p, *args):
        '''Dissolved concentration in ppmw in the melt'''
        return self.callmodel(p, *args)

class SolubilityH2O(Solubility):
    """H2O solubility models"""

    # below default gives the default model used
    def __init__(self, composition='peridotite'):
        super().__init__(composition)

    def anorthite_diopside(self, p):
        '''Newcombe et al. (2017)'''
        return self.power_law(p, 727, 0.5)

    def peridotite(self, p):
        '''Sossi et al. (2022)'''
        return self.power_law(p, 524, 0.5)

    def basalt_dixon(self, p):
        '''Dixon et al. (1995) refit by Paolo Sossi'''
        return self.power_law(p, 965, 0.5)

    def basalt_wilson(self, p):
        '''Hamilton (1964) and Wilson and Head (1981)'''
        return self.power_law(p, 215, 0.7)

    def lunar_glass(self, p):
        '''Newcombe et al. (2017)'''
        return self.power_law(p, 683, 0.5)


class SolubilityS2(Solubility):
    """S2 solubility models"""

    # below default gives the default model used
    def __init__(self, composition='gaillard'):
        self.fO2_model = OxygenFugacity()
        super().__init__(composition)

    def gaillard(self, p, temp, fO2_shift):
        # Gaillard et al., 2022
        # https://doi.org/10.1016/j.epsl.2021.117255
        # https://ars.els-cdn.com/content/image/1-s2.0-S0012821X21005112-mmc1.pdf

        if p < 1.0e-20:
            return 0.0

        # melt composition [wt%]
        x_FeO  = 10.0

        # calculate fO2 [bar]
        fO2 = 10**self.fO2_model(temp, fO2_shift)

        # calculate log(Ss)
        out = 13.8426 - 26.476e3/temp + 0.124*x_FeO + 0.5*np.log(p/fO2)

        # convert to concentration ppmw
        out = np.exp(out) #* 10000.0

        return out


class SolubilityCO2(Solubility):
    """CO2 solubility models"""

    def __init__(self, composition='basalt_dixon'):
        super().__init__(composition)

    def basalt_dixon(self, p, temp):
        '''Dixon et al. (1995)'''
        ppmw = (3.8E-7)*p*np.exp(-23*(p-1)/(83.15*temp))
        ppmw = 1.0E4*(4400*ppmw) / (36.6-44*ppmw)
        return ppmw


class SolubilityN2(Solubility):
    """N2 solubility models"""

    def __init__(self, composition='libourel'):
        super().__init__(composition)

        # melt composition
        x_SiO2  = 0.56
        x_Al2O3 = 0.11
        x_TiO2  = 0.01
        self.dasfac_2 = np.exp(4.67 + 7.11*x_SiO2 - 13.06*x_Al2O3 - 120.67*x_TiO2)

    def libourel(self, p):
        '''Libourel et al. (2003)'''
        ppmw = self.power_law(p, 0.0611, 1.0)
        return ppmw

    def dasgupta(self, p, ptot, temp, fO2_shift):
        '''Dasgupta et al. (2022)'''

        # convert bar to GPa
        pb_N2  = p * 1.0e-4
        pb_tot = ptot * 1.0e-4

        pb_tot = max(pb_tot, 1e-15)

        # calculate N2 concentration in melt
        ppmw  = pb_N2**0.5 * np.exp(5908.0 * pb_tot**0.5/temp - 1.6*fO2_shift)
        ppmw += pb_N2 * self.dasfac_2

        return ppmw

class SolubilityCH4(Solubility):
    """CH4 solubility models"""

    def __init__(self, composition='basalt_ardia'):
        super().__init__(composition)

    def basalt_ardia(self, p, p_total):
        '''Ardia 2013'''
        p_total *= 1e-4  # Convert to GPa
        p *= 1e-4 # Convert to GPa
        ppmw = p*np.exp(4.93 - (0.000193 * p_total))
        return ppmw


class SolubilityCO(Solubility):
    """CO solubility models"""

    def __init__(self, composition='mafic_armstrong'):
        super().__init__(composition)

    def mafic_armstrong(self, p, p_total):
        '''Armstrong 2015'''
        ppmw = 10 ** (-0.738 + 0.876 * np.log10(p) - 5.44e-5 * p_total)
        return ppmw

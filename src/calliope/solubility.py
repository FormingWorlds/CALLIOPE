# Solubility laws
from __future__ import annotations

import logging

import numpy as np

from .constants import molar_mass, noble_gases
from .oxygen_fugacity import OxygenFugacity

log = logging.getLogger('fwl.' + __name__)

# Jambon, Weill & Braun (1986), doi:10.1016/0016-7037(86)90193-6. Henry's-law
# solubility constants for noble gases in tholeiitic basalt melt, in units of
# cm3 STP per gram of melt per bar of partial pressure. Measured at 1 bar and
# 1250-1600 C. These are the same primitive numbers atmodeller uses for its
# `<gas>_basalt_jambon86` models, so the two backends produce identical Henry
# constants and agree to solver tolerance at matched conditions.
JAMBON86_STP_HENRY = {  # cm3 STP / g / bar
    'He': 56e-5,
    'Ne': 25e-5,
    'Ar': 5.9e-5,
    'Kr': 3.0e-5,
    'Xe': 1.7e-5,
}

# Molar volume of an ideal gas at standard temperature and pressure
# [cm3 / mol]. Converts the STP-volume Henry constants above into a molar
# basis. The 2.24e4 literal matches atmodeller's conversion exactly.
STP_MOLAR_VOLUME_CM3 = 2.24e4


class Solubility:
    """Solubility base class.

    Pressures are in bar; subclasses return dissolved concentration in
    ppmw (parts-per-million by weight in the silicate melt).
    """

    def __init__(self, composition):
        self.callmodel = getattr(self, composition)

    def power_law(self, p, const, exponent):
        return const * p**exponent

    def __call__(self, p, *args):
        """Dissolved concentration in ppmw in the melt"""
        return self.callmodel(p, *args)


class SolubilityH2O(Solubility):
    """H2O solubility models"""

    # below default gives the default model used
    def __init__(self, composition='peridotite'):
        super().__init__(composition)

    def anorthite_diopside(self, p):
        """Newcombe et al. (2017)"""
        return self.power_law(p, 727, 0.5)

    def peridotite(self, p):
        """Sossi et al. (2023)"""
        return self.power_law(p, 524, 0.5)

    def basalt_dixon(self, p):
        """Dixon et al. (1995) refit by Paolo Sossi"""
        return self.power_law(p, 965, 0.5)

    def basalt_wilson(self, p):
        """Hamilton (1964) and Wilson and Head (1981)"""
        return self.power_law(p, 215, 0.7)

    def lunar_glass(self, p):
        """Newcombe et al. (2017)"""
        return self.power_law(p, 683, 0.5)


class SolubilityS2(Solubility):
    """S2 solubility models.

    Parameters
    ----------
    composition : str, default 'gaillard'
        Solubility-law name (currently only 'gaillard' is implemented).
    x_FeO : float, default 10.0
        Melt FeO content [wt%] used by the Gaillard et al. (2022) law.
        The default value matches the Earth-mantle reference adopted in
        prior CALLIOPE releases; override for non-Earth bulk
        compositions.
    """

    def __init__(self, composition='gaillard', x_FeO=10.0):
        self.fO2_model = OxygenFugacity()
        self.x_FeO = x_FeO
        super().__init__(composition)

    def gaillard(self, p, temp, fO2_shift):
        # Gaillard et al., 2022
        # https://doi.org/10.1016/j.epsl.2021.117255
        # https://ars.els-cdn.com/content/image/1-s2.0-S0012821X21005112-mmc1.pdf

        if p < 1.0e-20:
            return 0.0

        # calculate fO2 [bar]
        fO2 = 10 ** self.fO2_model(temp, fO2_shift)

        # calculate log(Ss); x_FeO [wt%] is set on the instance, default
        # 10.0 wt% (Earth-mantle reference).
        out = 13.8426 - 26.476e3 / temp + 0.124 * self.x_FeO + 0.5 * np.log(p / fO2)

        # convert to concentration ppmw
        out = np.exp(out)  # * 10000.0

        return out


class SolubilityCO2(Solubility):
    """CO2 solubility models"""

    def __init__(self, composition='basalt_dixon'):
        super().__init__(composition)

    def basalt_dixon(self, p, temp):
        """Dixon et al. (1995)"""
        ppmw = (3.8e-7) * p * np.exp(-23 * (p - 1) / (83.15 * temp))
        ppmw = 1.0e4 * (4400 * ppmw) / (36.6 - 44 * ppmw)
        return ppmw


class SolubilityN2(Solubility):
    """N2 solubility models.

    Parameters
    ----------
    composition : str, default 'libourel'
        Solubility-law name. 'libourel' selects the linear Henry's-law
        form of Libourel et al. (2003); 'dasgupta' selects the
        physical-state-dependent form of Dasgupta et al. (2022).
    x_SiO2, x_Al2O3, x_TiO2 : float
        Melt mole fractions used by the Dasgupta et al. (2022) law to
        compute the molecular-N2 prefactor `dasfac_2`. Defaults
        (0.56, 0.11, 0.01) match the Earth-mantle reference adopted in
        prior CALLIOPE releases; override for non-Earth compositions.
        These kwargs have no effect when `composition='libourel'`.
    """

    def __init__(self, composition='libourel', x_SiO2=0.56, x_Al2O3=0.11, x_TiO2=0.01):
        super().__init__(composition)

        # Stored on the instance so callers can introspect them; defaults
        # match the Earth-mantle reference used by Dasgupta et al. (2022).
        self.x_SiO2 = x_SiO2
        self.x_Al2O3 = x_Al2O3
        self.x_TiO2 = x_TiO2
        # dasfac_2 is only consumed by the dasgupta() path, so skip the
        # exp(...) precompute when libourel is selected. This avoids a
        # spurious RuntimeWarning at construction time when a libourel
        # caller passes extreme composition values that would overflow
        # the exponent (the dasgupta path is the only one that cares).
        if composition == 'dasgupta':
            self.dasfac_2 = np.exp(4.67 + 7.11 * x_SiO2 - 13.06 * x_Al2O3 - 120.67 * x_TiO2)
        else:
            self.dasfac_2 = None

    def libourel(self, p):
        """Libourel et al. (2003)"""
        ppmw = self.power_law(p, 0.0611, 1.0)
        return ppmw

    def dasgupta(self, p, ptot, temp, fO2_shift):
        """Dasgupta et al. (2022)"""

        # convert bar to GPa
        pb_N2 = p * 1.0e-4
        pb_tot = ptot * 1.0e-4

        pb_tot = max(pb_tot, 1e-15)

        # calculate N2 concentration in melt
        ppmw = pb_N2**0.5 * np.exp(5908.0 * pb_tot**0.5 / temp - 1.6 * fO2_shift)
        ppmw += pb_N2 * self.dasfac_2

        return ppmw


class SolubilityCH4(Solubility):
    """CH4 solubility models"""

    def __init__(self, composition='basalt_ardia'):
        super().__init__(composition)

    def basalt_ardia(self, p, p_total):
        """Ardia 2013"""
        p_total *= 1e-4  # Convert to GPa
        p *= 1e-4  # Convert to GPa
        ppmw = p * np.exp(4.93 - (1.93 * p_total))
        return ppmw


class SolubilityCO(Solubility):
    """CO solubility models"""

    def __init__(self, composition='mafic_armstrong'):
        super().__init__(composition)

    def mafic_armstrong(self, p, p_total):
        """Armstrong 2015"""
        ppmw = 10 ** (-0.738 + 0.876 * np.log10(p) - 5.44e-5 * p_total)
        return ppmw


def jambon86_ppmw_per_bar(gas: str) -> float:
    """Henry's-law solubility constant for a noble gas [ppmw / bar].

    Converts the Jambon et al. (1986) STP-volume Henry constant for `gas`
    into parts-per-million by weight of dissolved gas per bar of partial
    pressure, using the melt-independent chain

        const [ppmw/bar] = (k_STP / V_STP) [mol/g/bar]
                           * M [g/mol]
                           * 1e6 [ppmw per mass fraction]

    where `k_STP` is the tabulated constant in cm3 STP/g/bar, `V_STP` is the
    molar volume of an ideal gas at STP, and `M` is the molar mass. The
    result is the linear coefficient of the `ppmw = const * p` Henry law.

    Parameters
    ----------
    gas : str
        Noble gas symbol; one of `He`, `Ne`, `Ar`, `Kr`, `Xe`.

    Returns
    -------
    float
        Solubility constant in ppmw per bar.

    Raises
    ------
    KeyError
        If `gas` is not a Jambon et al. (1986) noble gas.
    """
    k_stp = JAMBON86_STP_HENRY[gas]
    molar_mass_g = molar_mass[gas] * 1.0e3  # kg/mol -> g/mol
    return (k_stp / STP_MOLAR_VOLUME_CM3) * molar_mass_g * 1.0e6


class SolubilityNobleGas(Solubility):
    """Noble gas solubility by Henry's law, Jambon et al. (1986).

    Each noble gas dissolves in silicate melt in proportion to its partial
    pressure, `ppmw = const * p`, with no melt-composition, temperature, or
    redox dependence in this parameterization. The linear (exponent 1) form
    is the defining property of Henry's law and distinguishes the noble
    gases from the square-root CHNOS laws in this module.

    Parameters
    ----------
    gas : str
        Noble gas symbol; one of `He`, `Ne`, `Ar`, `Kr`, `Xe`.

    Notes
    -----
    The calibration is tholeiitic basalt at 1 bar and 1250-1600 C. Applying
    it at the high surface pressures of a noble-gas-rich atmosphere is an
    extrapolation of a 1-bar, linear Henry law with no saturation term.
    """

    def __init__(self, gas: str):
        if gas not in noble_gases:
            raise ValueError(
                f"SolubilityNobleGas: '{gas}' is not a noble gas. "
                f'Expected one of {noble_gases}.'
            )
        self.gas = gas
        self.const = jambon86_ppmw_per_bar(gas)
        super().__init__('jambon86')

    def jambon86(self, p):
        """Jambon et al. (1986) linear Henry's law: `ppmw = const * p`."""
        return self.power_law(p, self.const, 1.0)

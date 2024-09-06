# Simple structure model
from __future__ import annotations

import logging

import numpy as np

from .constants import M_earth, R_earth

log = logging.getLogger("fwl."+__name__)


def calculate_mantle_mass(radius:float, mass:float, corefrac:float)->float:
    '''
    A very simple interior structure model.

    This calculates mantle mass given planetary mass, radius, and core fraction. This
    assumes a core density equal to that of Earth's, and that the planet mass is simply
    the sum of mantle and core.
    '''

    earth_fr = 0.55     # earth core radius fraction
    earth_fm = 0.325    # earth core mass fraction  (https://arxiv.org/pdf/1708.08718.pdf)

    core_rho = (3.0 * earth_fm * M_earth) / (4.0 * np.pi * ( earth_fr * R_earth )**3.0 )  # core density [kg m-3]
    log.debug("Core density = %.2f kg m-3" % core_rho)

    # Calculate mantle mass by subtracting core from total
    core_mass = core_rho * 4.0/3.0 * np.pi * (radius * corefrac )**3.0
    mantle_mass = mass - core_mass
    log.info("Total mantle mass = %.2e kg" % mantle_mass)
    if (mantle_mass <= 0.0):
        raise Exception("Something has gone wrong (mantle mass is negative)")

    return mantle_mass

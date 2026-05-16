# Simple structure model
from __future__ import annotations

import logging
import warnings

import numpy as np

from .constants import M_earth, R_earth

log = logging.getLogger('fwl.' + __name__)


def calculate_mantle_mass(
    radius: float,
    mass: float,
    core_frac: float | None = None,
    *,
    corefrac: float | None = None,
) -> float:
    """
    A very simple interior structure model.

    This calculates mantle mass given planetary mass, radius, and core fraction. This
    assumes a core density equal to that of Earth's, and that the planet mass is simply
    the sum of mantle and core.

    The legacy ``corefrac`` keyword is accepted as a deprecated alias for
    ``core_frac`` to keep older callers working; it will be removed in a
    future release.
    """
    # Backwards-compatibility shim for the corefrac -> core_frac rename.
    if corefrac is not None:
        if core_frac is not None:
            raise TypeError(
                "calculate_mantle_mass() received both 'core_frac' and "
                "'corefrac'; pass only 'core_frac'."
            )
        warnings.warn(
            "The 'corefrac' keyword is deprecated and will be removed; "
            "use 'core_frac' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        core_frac = corefrac
    if core_frac is None:
        raise TypeError("calculate_mantle_mass() missing required argument: 'core_frac'")

    earth_fr = 0.55  # earth core radius fraction
    # earth core mass fraction (32.5 +/- 0.3 wt%) from
    # Wang, Lineweaver & Ireland (2017), arxiv:1708.08718
    earth_fm = 0.325

    core_rho = (3.0 * earth_fm * M_earth) / (
        4.0 * np.pi * (earth_fr * R_earth) ** 3.0
    )  # core density [kg m-3]
    log.debug('Core density = %.2f kg m-3' % core_rho)

    # Calculate mantle mass by subtracting core from total
    core_mass = core_rho * 4.0 / 3.0 * np.pi * (radius * core_frac) ** 3.0
    mantle_mass = mass - core_mass
    log.info('Total mantle mass = %.2e kg' % mantle_mass)
    if mantle_mass <= 0.0:
        raise Exception('Something has gone wrong (mantle mass is negative)')

    return mantle_mass

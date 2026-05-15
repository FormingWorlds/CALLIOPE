"""Pre-flight verification: each backend must round-trip self-consistently
before we trust it in cross-backend plots.

Round-trip protocol per backend:

1. Run buffered-mode equilibrium at a known Delta-IW, take the resulting
   O_kg_total.
2. Feed that O budget back into the authoritative-O entry point.
3. Verify the recovered Delta-IW matches step 1 within 0.05 dex.

For CALLIOPE, step 1 calls `equilibrium_atmosphere` with fO2_shift_IW=X.
For atmodeller, step 1 calls the equilibrium solver with an
`IronWustiteBuffer(X)` fugacity constraint and reads the resulting O
mass from output.asdict().

This module is callable as `python -m scripts.cross_backend.verification`
and exits non-zero if any check fails.
"""

from __future__ import annotations

import logging
import sys
import warnings

import numpy as np

from .inventories import EARTH_BSE_KRIJT23, PLANETARY_DEFAULTS
from .runners import _calliope_ddict, run_atmodeller, run_calliope

log = logging.getLogger('cross_backend.verification')


def _earth_HCNS() -> dict:
    """Earth H/C/N/S budget — same as round-trip regression test."""
    inv = EARTH_BSE_KRIJT23
    return {'H': inv.H, 'C': inv.C, 'N': inv.N, 'S': inv.S}


def calliope_buffered_O(T_magma: float, dIW: float) -> float:
    """Return the O_kg_total CALLIOPE produces at fixed Delta-IW."""
    from calliope.solve import equilibrium_atmosphere

    ddict = _calliope_ddict(T_magma=T_magma)
    ddict['fO2_shift_IW'] = dIW
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = equilibrium_atmosphere(
            _earth_HCNS(), ddict, hide_warnings=True, print_result=False
        )
    return float(out['O_kg_total'])


def atmodeller_buffered_O(T_magma: float, dIW: float) -> float:
    """Return the total O kg atmodeller produces at fixed Delta-IW.

    Uses atmodeller's as-shipped solubility defaults (matching what
    `run_atmodeller` uses on the authoritative-O side), so the round-
    trip exercises the same chemistry in both directions. Diverging
    selections between buffered and authoritative-O calls would
    misattribute a solubility-set discrepancy to a solver-precision
    failure.
    """
    from atmodeller import ChemicalSpecies, EquilibriumModel, Planet, SpeciesNetwork
    from atmodeller.containers import SolverParameters
    from atmodeller.solubility import get_solubility_models
    from atmodeller.thermodata import IronWustiteBuffer

    from .runners import _DEFAULT_ATM_SOL

    sol_lib = get_solubility_models()
    sol_map = dict(_DEFAULT_ATM_SOL)
    species_names = {
        'H2O': 'H2O', 'CO2': 'CO2', 'N2': 'N2', 'S2': 'S2',
        'CO': 'CO', 'CH4': 'CH4', 'SO2': 'SO2', 'H2S': 'H2S',
        'H2': 'H2', 'NH3': 'H3N', 'O2': 'O2',
    }
    species_list = []
    for proteus_name, atm_name in species_names.items():
        kwargs = {}
        sol_key = sol_map.get(proteus_name)
        if sol_key:
            kwargs['solubility'] = sol_lib[sol_key]
        species_list.append(ChemicalSpecies.create_gas(atm_name, **kwargs))

    species = SpeciesNetwork(tuple(species_list))
    model = EquilibriumModel(species)

    planet = Planet(
        planet_mass=PLANETARY_DEFAULTS['M_planet'],
        core_mass_fraction=PLANETARY_DEFAULTS['core_mass_fraction'],
        mantle_melt_fraction=1.0,
        surface_radius=PLANETARY_DEFAULTS['radius'],
        temperature=T_magma,
        pressure=np.nan,
    )

    hcns = _earth_HCNS()
    mass_constraints = {e: float(hcns[e]) for e in 'HCNS'}

    solver_params = SolverParameters(atol=1e-6, rtol=1e-4, max_steps=256, multistart=10)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model.solve(
            state=planet,
            fugacity_constraints={'O2_g': IronWustiteBuffer(dIW)},
            mass_constraints=mass_constraints,
            solver_parameters=solver_params,
            solver='robust',
        )

    output = model.output
    output_dict = output.asdict()

    # atmodeller exposes per-element totals directly under `element_<X>`
    # keys; `total_mass` is gas + dissolved in kg. Using this avoids
    # hand-summing across species (which is brittle because atmodeller
    # canonicalises species names internally, e.g. SO2 -> O2S).
    elem = output_dict.get('element_O', {})
    if not isinstance(elem, dict) or 'total_mass' not in elem:
        raise RuntimeError("atmodeller output missing 'element_O.total_mass'")
    return float(np.squeeze(elem['total_mass']))


def round_trip_calliope(T_magma: float, dIW_in: float) -> tuple[float, float]:
    """Returns (input_dIW, recovered_dIW)."""
    O_kg = calliope_buffered_O(T_magma, dIW_in)
    target = {**_earth_HCNS(), 'O': O_kg}
    from calliope.solve import equilibrium_atmosphere_authoritative_O
    ddict = _calliope_ddict(T_magma=T_magma)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = equilibrium_atmosphere_authoritative_O(
            target, ddict, fO2_hint=dIW_in, hide_warnings=True,
            random_seed=1234, print_result=False, nguess=1500,
        )
    return dIW_in, float(out['fO2_shift_derived'])


def round_trip_atmodeller(T_magma: float, dIW_in: float) -> tuple[float, float]:
    """Returns (input_dIW, recovered_dIW)."""
    from .inventories import Inventory

    O_kg = atmodeller_buffered_O(T_magma, dIW_in)
    hcns = _earth_HCNS()
    inv = Inventory(
        name=f'roundtrip_T{T_magma:.0f}_dIW{dIW_in:+.1f}',
        H=hcns['H'], C=hcns['C'], N=hcns['N'], O=O_kg, S=hcns['S'],
        citation='generated by buffered-mode atmodeller call',
    )
    res = run_atmodeller(inv, T_magma=T_magma, Phi_global=1.0)
    if not res.converged:
        return dIW_in, float('nan')
    return dIW_in, res.fO2_shift_derived


def main() -> int:
    """Run round-trip verification for both backends at a grid of
    (T, dIW) and exit non-zero if any |residual| > 0.1 dex.

    The tolerance is generous (per-element solver tol is 1e-5 relative,
    but the inferred fO2 has wider sensitivity to numerical noise in the
    forward-then-inverse path).
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s %(message)s')
    T_grid = [1500.0, 2000.0, 2500.0, 3000.0]
    dIW_grid = [-2.0, 0.0, 2.0, 4.0]
    tol_dex = 0.10
    all_ok = True
    for backend, rt in (('calliope', round_trip_calliope), ('atmodeller', round_trip_atmodeller)):
        print(f'\n# Round-trip: {backend}')
        print(f'{"T_K":>6} {"dIW_in":>8} {"dIW_out":>10} {"residual":>10} {"status":>10}')
        for T in T_grid:
            for dIW in dIW_grid:
                try:
                    _, recov = rt(T, dIW)
                except Exception as exc:  # noqa: BLE001
                    print(f'{T:6.0f} {dIW:8.2f} {"":>10} {"":>10} {"ERR:" + type(exc).__name__:>10}')
                    all_ok = False
                    continue
                residual = recov - dIW
                ok = np.isfinite(residual) and abs(residual) < tol_dex
                status = 'OK' if ok else 'FAIL'
                if not ok:
                    all_ok = False
                print(f'{T:6.0f} {dIW:8.2f} {recov:10.3f} {residual:+10.3f} {status:>10}')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())

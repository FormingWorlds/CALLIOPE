"""Backend runners: CALLIOPE and atmodeller authoritative-O entry points.

Each runner accepts a normalised `Inventory`, planetary state, and
backend configuration knobs, runs the corresponding solver, and returns
a `BackendResult` with the converged Delta-IW, surface partial
pressures, dissolved masses, and convergence status. Failures are
caught and reported in the result (not raised), so a sweep over a
parameter grid can survive isolated non-convergence and report
coverage honestly.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np

from .inventories import PLANETARY_DEFAULTS, Inventory

log = logging.getLogger('cross_backend.runners')


@dataclass
class BackendResult:
    """Output of one backend call.

    Attributes
    ----------
    backend : str
        'calliope' or 'atmodeller'.
    converged : bool
        True if the solver returned a usable result. False if it raised
        or signalled non-convergence; the caller should mask this point
        from any quantitative comparison.
    fO2_shift_derived : float
        log10 IW-buffer offset the solver converged to. NaN if not
        converged.
    p_bar : dict
        Partial pressures keyed by species name (`H2O`, `CO2`, ...) in bar.
        Empty if not converged.
    dissolved_kg : dict
        Per-species dissolved mass in kg. Empty if not converged.
    total_P_bar : float
        Sum of partial pressures.
    error : str
        Empty on success, exception message on failure.
    inputs : dict
        Echo of the inputs that produced this result, for provenance.
    """

    backend: str
    converged: bool
    fO2_shift_derived: float = float('nan')
    p_bar: dict = field(default_factory=dict)
    dissolved_kg: dict = field(default_factory=dict)
    total_P_bar: float = float('nan')
    error: str = ''
    inputs: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# CALLIOPE
# ---------------------------------------------------------------------------


def _calliope_ddict(
    T_magma: float,
    Phi_global: float = 1.0,
    species_set: tuple = ('H2O', 'CO2', 'H2', 'CO', 'CH4', 'N2', 'NH3', 'S2', 'SO2', 'H2S'),
    **planet_overrides,
) -> dict:
    """Build the CALLIOPE coupler-options dict consumed by the solver."""
    p = dict(PLANETARY_DEFAULTS)
    p.update(planet_overrides)
    ddict = {
        'M_mantle': p['M_mantle'],
        'gravity': p['gravity'],
        'radius': p['radius'],
        'Phi_global': Phi_global,
        'T_magma': T_magma,
        'fO2_shift_IW': 0.0,
    }
    all_species = ('H2O', 'CO2', 'H2', 'CO', 'CH4', 'N2', 'NH3', 'S2', 'SO2', 'H2S')
    for sp in all_species:
        ddict[f'{sp}_included'] = 1 if sp in species_set else 0
        ddict[f'{sp}_initial_bar'] = 0.0
    return ddict


def _with_calliope_buffer(buffer: str):
    """Context manager that temporarily pins CALLIOPE's default IW
    buffer to ``buffer`` ('fischer' or 'oneill'). Restores the prior
    default on exit. Use this to compare both buffer choices from the
    same harness regardless of which one is the library default.
    """
    from contextlib import contextmanager

    from calliope.chemistry import ModifiedKeq
    from calliope.oxygen_fugacity import OxygenFugacity

    @contextmanager
    def _ctx():
        of_old = OxygenFugacity.__init__.__defaults__
        mk_old = ModifiedKeq.__init__.__defaults__
        OxygenFugacity.__init__.__defaults__ = (buffer,)
        ModifiedKeq.__init__.__defaults__ = (buffer,)
        try:
            yield
        finally:
            OxygenFugacity.__init__.__defaults__ = of_old
            ModifiedKeq.__init__.__defaults__ = mk_old

    return _ctx()


def run_calliope(
    inventory: Inventory,
    T_magma: float,
    fO2_hint: float = 2.0,
    Phi_global: float = 1.0,
    random_seed: int | None = 1234,
    nguess: int = 1500,
    print_result: bool = False,
    buffer: str = 'fischer',
    **planet_overrides,
) -> BackendResult:
    """Run CALLIOPE's authoritative-O solver.

    Suppresses fsolve convergence-warning chatter for cleaner sweep
    output. The `random_seed` is fixed by default so a re-run of the
    figure scripts produces identical numbers. ``buffer`` selects the
    IW parameterisation; default 'fischer' matches the CALLIOPE
    library default; 'oneill' reproduces the legacy CALLIOPE behaviour.
    """
    from calliope.solve import equilibrium_atmosphere_authoritative_O

    ddict = _calliope_ddict(T_magma=T_magma, Phi_global=Phi_global, **planet_overrides)
    target_d = inventory.asdict()
    inputs = dict(
        inventory=inventory.name,
        T_magma=T_magma,
        Phi_global=Phi_global,
        fO2_hint=fO2_hint,
        seed=random_seed,
        target_d=target_d,
        buffer=buffer,
    )
    try:
        with warnings.catch_warnings(), _with_calliope_buffer(buffer):
            warnings.simplefilter('ignore')
            out = equilibrium_atmosphere_authoritative_O(
                target_d,
                ddict,
                fO2_hint=fO2_hint,
                hide_warnings=True,
                random_seed=random_seed,
                nguess=nguess,
                print_result=print_result,
            )
    except Exception as exc:  # noqa: BLE001  (sweep must survive isolated failures)
        return BackendResult(
            backend='calliope',
            converged=False,
            error=f'{type(exc).__name__}: {exc}',
            inputs=inputs,
        )

    p_bar = {
        sp: float(out[f'{sp}_bar'])
        for sp in ('H2O', 'CO2', 'H2', 'CO', 'CH4', 'N2', 'NH3', 'S2', 'SO2', 'H2S')
        if f'{sp}_bar' in out
    }
    dissolved = {sp: float(out[f'{sp}_kg_liquid']) for sp in p_bar if f'{sp}_kg_liquid' in out}
    return BackendResult(
        backend='calliope',
        converged=True,
        fO2_shift_derived=float(out['fO2_shift_derived']),
        p_bar=p_bar,
        dissolved_kg=dissolved,
        total_P_bar=float(sum(p_bar.values())),
        inputs=inputs,
    )


# ---------------------------------------------------------------------------
# atmodeller
# ---------------------------------------------------------------------------


_ATM_SPECIES_MAP = {
    # PROTEUS-name -> atmodeller canonical-name (alphabetical-by-element).
    # SO2 is renamed to O2S internally; NH3 to H3N. Get the canonical
    # name wrong and the species silently drops out of output.asdict().
    'H2O': 'H2O',
    'H2': 'H2',
    'CO2': 'CO2',
    'CO': 'CO',
    'CH4': 'CH4',
    'N2': 'N2',
    'NH3': 'H3N',
    'S2': 'S2',
    'SO2': 'O2S',
    'H2S': 'H2S',
    'O2': 'O2',
}


_DEFAULT_ATM_SOL = {
    'H2O': 'H2O_peridotite_sossi23',
    'CO2': 'CO2_basalt_dixon95',
    'H2': 'H2_basalt_hirschmann12',
    'N2': 'N2_basalt_dasgupta22',
    'S2': 'S2_sulfide_basalt_boulliung23',
    'CO': 'CO_basalt_yoshioka19',
    'CH4': 'CH4_basalt_ardia13',
}


_CALLIOPE_ALIGNED_ATM_SOL = {
    # CALLIOPE has explicit solubility for H2O / CO2 / N2 / S2 only; the
    # other species inherit zero dissolved mass via Bower 2022 §2.2.3.
    # Matching that selection at the atmodeller side isolates the
    # remaining sources of disagreement.
    'H2O': 'H2O_peridotite_sossi23',
    'CO2': 'CO2_basalt_dixon95',
    'H2': None,
    'N2': 'N2_basalt_dasgupta22',
    'S2': 'S2_sulfide_basalt_boulliung23',  # Gaillard22 absent from atmodeller library
    'CO': None,
    'CH4': None,
}


def run_atmodeller(
    inventory: Inventory,
    T_magma: float,
    Phi_global: float = 1.0,
    solubility_map: dict | None = None,
    eos_map: dict | None = None,
    include_condensates: bool = False,
    solver_multistart: int = 10,
    solver_atol: float = 1e-6,
    solver_rtol: float = 1e-4,
    solver_max_steps: int = 256,
    **planet_overrides,
) -> BackendResult:
    """Run atmodeller's authoritative-O solver via its direct API.

    This bypasses the PROTEUS wrapper so the harness has no PROTEUS
    runtime dependency. The species list, solubility selection, EOS
    selection, and the no-fugacity-constraint authoritative-O path
    mirror what the PROTEUS wrapper would build at runtime.
    """
    from atmodeller import ChemicalSpecies, EquilibriumModel, Planet, SpeciesNetwork
    from atmodeller.containers import SolverParameters
    from atmodeller.solubility import get_solubility_models

    sol_map = solubility_map if solubility_map is not None else _DEFAULT_ATM_SOL
    eos_map = eos_map or {}
    sol_lib = get_solubility_models()

    p = dict(PLANETARY_DEFAULTS)
    p.update(planet_overrides)

    inputs = dict(
        inventory=inventory.name,
        T_magma=T_magma,
        Phi_global=Phi_global,
        solubility_map=dict(sol_map),
        eos_map=dict(eos_map),
        include_condensates=include_condensates,
        target_d=inventory.asdict(),
    )

    try:
        # Only include species whose constituent elements all have
        # non-zero budgets. Matches the wrapper's active-elements logic.
        species_elements = {
            'H2O': {'H'},
            'H2': {'H'},
            'CO2': {'C'},
            'CO': {'C'},
            'CH4': {'H', 'C'},
            'N2': {'N'},
            'NH3': {'H', 'N'},
            'S2': {'S'},
            'SO2': {'S'},
            'H2S': {'H', 'S'},
            'O2': set(),
        }
        budgets = inventory.asdict()
        active_elements = {e for e in 'HCNSO' if budgets.get(e, 0.0) > 0.0}
        active_species = {
            sp for sp, req in species_elements.items() if req.issubset(active_elements)
        }

        species_list = []
        for proteus_name, atm_name in _ATM_SPECIES_MAP.items():
            if proteus_name not in active_species:
                continue
            kwargs = {}
            sol_key = sol_map.get(proteus_name)
            if sol_key and sol_key in sol_lib:
                kwargs['solubility'] = sol_lib[sol_key]
            elif sol_key:
                raise ValueError(f'Unknown atmodeller solubility key: {sol_key!r}')
            eos_key = eos_map.get(proteus_name)
            if eos_key:
                from atmodeller.eos import get_eos_models

                eos_lib = get_eos_models()
                if eos_key in eos_lib:
                    kwargs['activity'] = eos_lib[eos_key]
                else:
                    raise ValueError(f'Unknown atmodeller EOS key: {eos_key!r}')
            species_list.append(ChemicalSpecies.create_gas(atm_name, **kwargs))

        if include_condensates and 'C' in active_elements:
            try:
                species_list.append(ChemicalSpecies.create_condensed('C'))
            except Exception:  # noqa: BLE001
                pass

        species = SpeciesNetwork(tuple(species_list))
        model = EquilibriumModel(species)

        planet = Planet(
            planet_mass=p['M_planet'],
            core_mass_fraction=p['core_mass_fraction'],
            mantle_melt_fraction=Phi_global,
            surface_radius=p['radius'],
            temperature=T_magma,
            pressure=np.nan,
        )

        mass_constraints = {e: float(budgets[e]) for e in 'HCNSO' if budgets.get(e, 0.0) > 0.0}

        solver_params = SolverParameters(
            atol=solver_atol,
            rtol=solver_rtol,
            max_steps=solver_max_steps,
            multistart=solver_multistart,
        )

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model.solve(
                state=planet,
                mass_constraints=mass_constraints,
                solver_parameters=solver_params,
                solver='robust',
            )

        output = model.output
        quick_look = output.quick_look()
        output_dict = output.asdict()

        reverse_map = {v: k for k, v in _ATM_SPECIES_MAP.items()}
        p_bar = {}
        for atm_name, p_val in quick_look.items():
            stripped = atm_name.replace('_g', '')
            proteus_name = reverse_map.get(stripped)
            if proteus_name is None:
                continue
            p_bar[proteus_name] = float(np.squeeze(p_val))

        dissolved = {}
        for proteus_name, atm_name in _ATM_SPECIES_MAP.items():
            key = f'{atm_name}_g'
            sd = output_dict.get(key, {})
            if isinstance(sd, dict):
                dm = sd.get('dissolved_mass')
                if dm is not None:
                    try:
                        dissolved[proteus_name] = max(0.0, float(np.squeeze(dm)))
                    except (TypeError, ValueError):
                        pass

        o2 = output_dict.get('O2_g', {})
        log10dIW = None
        if isinstance(o2, dict):
            v = o2.get('log10dIW_1_bar')
            if v is not None:
                log10dIW = float(np.squeeze(v))

        if log10dIW is None or not np.isfinite(log10dIW):
            return BackendResult(
                backend='atmodeller',
                converged=False,
                error='atmodeller returned no log10dIW_1_bar',
                inputs=inputs,
            )

        return BackendResult(
            backend='atmodeller',
            converged=True,
            fO2_shift_derived=log10dIW,
            p_bar=p_bar,
            dissolved_kg=dissolved,
            total_P_bar=float(sum(p_bar.values())),
            inputs=inputs,
        )
    except Exception as exc:  # noqa: BLE001
        return BackendResult(
            backend='atmodeller',
            converged=False,
            error=f'{type(exc).__name__}: {exc}',
            inputs=inputs,
        )

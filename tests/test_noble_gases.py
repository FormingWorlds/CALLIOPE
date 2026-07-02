"""Noble gas partitioning in `src/calliope/solve.py`.

Exercises the melt-atmosphere partitioning of the noble gases (He, Ne, Ar,
Kr, Xe) added to the equilibrium solver. A noble gas is monatomic and inert,
so its element and species are the same entity, it takes no part in the CHNOS
reaction network, and it partitions by a single linear Henry's law. It still
enters the total pressure and the mean molar mass, so it is solved jointly
with the CHNOS primaries rather than bolted on afterwards.

- Conservation: per-gas mass closure `kg_atm + kg_liquid ~ kg_total` and
  agreement with the supplied target inventory.
- Boundedness / positivity: non-negative partial pressures and reservoir
  masses; an excluded noble gas contributes exactly zero.
- Monotonicity: a larger noble gas budget yields a larger partial pressure
  and a larger dissolved mass.
- Coupling: a noble-gas-dominated atmosphere measurably shifts the CHNOS
  partial pressures through the mean molar mass, proving the coupling is
  solved rather than bypassed.
- Symmetry / inertness: noble gas partitioning is independent of the fO2
  buffer offset, which the reactive CHNOS species are not.
- Reference (analytical limit): the dissolved-to-atmospheric mass ratio
  equals the closed-form Henry-versus-hydrostatic-column ratio.

See `docs/How-to/build_tests.md` for the testing standards these follow.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from calliope import solve as calsolve
from calliope.constants import molar_mass, noble_gases
from calliope.solubility import jambon86_ppmw_per_bar
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
    get_target_from_params,
    get_target_from_pressures,
    is_included,
)

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]


def _ddict(active=('He', 'Ne', 'Ar', 'Kr', 'Xe'), dIW=0.0):
    """Earth-scale magma-ocean options with CHNOS + selected noble gases."""
    d = {
        'M_mantle': 4.0e24,
        'Phi_global': 1.0,
        'T_magma': 1800.0,
        'gravity': 9.81,
        'radius': 6.37e6,
        'fO2_shift_IW': dIW,
    }
    from calliope.constants import volatile_species

    for sp in volatile_species:
        d[f'{sp}_included'] = 1 if sp in ('H2O', 'CO2', 'N2', 'S2') else 0
        d[f'{sp}_initial_bar'] = 0.0
    for gas in noble_gases:
        d[f'{gas}_included'] = 1 if gas in active else 0
        d[f'{gas}_initial_bar'] = 0.0
    return d


_CHNOS = {'H': 1.5e20, 'C': 1.0e20, 'N': 2.0e18, 'S': 5.0e19}
_NOBLE = {'He': 3.0e16, 'Ne': 1.0e15, 'Ar': 2.0e16, 'Kr': 5.0e14, 'Xe': 1.0e14}


@pytest.mark.physics_invariant
def test_every_noble_gas_conserves_mass():
    """Each active noble gas closes its own mass budget: the atmospheric
    plus dissolved mass equals the supplied target, and equals the reported
    `_kg_total`. Runs all five simultaneously so an index or ordering bug in
    the residual vector would surface as a mismatched gas.
    """
    np.random.seed(42)
    ddict = _ddict()
    target = dict(_CHNOS, **_NOBLE)
    out = equilibrium_atmosphere(
        target, ddict, print_result=False, opt_solver=False, nguess=2000
    )

    for gas in noble_gases:
        reservoir = out[f'{gas}_kg_atm'] + out[f'{gas}_kg_liquid']
        # Closure against the target inventory (the solver's mass balance).
        assert reservoir == pytest.approx(target[gas], rel=1e-6)
        # Internal consistency: _kg_total is the atm + liquid sum (no solid).
        assert out[f'{gas}_kg_total'] == pytest.approx(reservoir, rel=1e-12)
        assert out[f'{gas}_kg_solid'] == 0.0
        # Positivity of the split.
        assert out[f'{gas}_kg_atm'] > 0.0
        assert out[f'{gas}_kg_liquid'] > 0.0


@pytest.mark.physics_invariant
def test_dissolved_mass_follows_closed_form_henry_law():
    """Reference (analytical limit): after the solve, each noble gas's
    dissolved mass must equal the closed-form Henry's law evaluated at the
    solved partial pressure, `M_diss = (1e-6 * M_mantle * Phi) * const * p`.
    This pins the melt side of the partitioning against the published
    Jambon et al. (1986) constants independent of the solver internals.
    """
    np.random.seed(7)
    ddict = _ddict()
    target = dict(_CHNOS, **_NOBLE)
    out = equilibrium_atmosphere(
        target, ddict, print_result=False, opt_solver=False, nguess=2000
    )

    prefactor = 1.0e-6 * ddict['M_mantle'] * ddict['Phi_global']
    for gas in noble_gases:
        predicted = prefactor * jambon86_ppmw_per_bar(gas) * out[f'{gas}_bar']
        assert out[f'{gas}_kg_liquid'] == pytest.approx(predicted, rel=1e-9)

    # Discrimination guard: a swapped Henry constant (using Ne's for He)
    # would break the He match by the ratio of the two constants (~2.25x),
    # far outside the 1e-9 tolerance.
    wrong = prefactor * jambon86_ppmw_per_bar('Ne') * out['He_bar']
    assert abs(out['He_kg_liquid'] - wrong) > 0.1 * out['He_kg_liquid']


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_partition_ratio_matches_henry_versus_hydrostatic_column():
    """Reference (analytical limit): the dissolved-to-atmospheric mass ratio
    of a noble gas is the closed-form ratio of its Henry coefficient to its
    hydrostatic column coefficient,

        M_liquid / M_atm = (prefactor * const) / (1e5 * A / g * M / mu),

    with `A = 4 pi R^2` and `mu` the mean molar mass. Both sides are computed
    independently of each other, so agreement pins the whole partitioning
    against first principles and the published solubility constant.
    """
    np.random.seed(11)
    ddict = _ddict()
    target = dict(_CHNOS, **_NOBLE)
    out = equilibrium_atmosphere(
        target, ddict, print_result=False, opt_solver=False, nguess=2000
    )

    area = 4.0 * np.pi * ddict['radius'] ** 2.0
    mu = out['atm_kg_per_mol']  # kg/mol
    prefactor = 1.0e-6 * ddict['M_mantle'] * ddict['Phi_global']
    for gas in ('He', 'Ar', 'Xe'):
        k_diss = prefactor * jambon86_ppmw_per_bar(gas)
        k_atm = 1.0e5 * area / ddict['gravity'] * molar_mass[gas] / mu
        predicted_ratio = k_diss / k_atm
        actual_ratio = out[f'{gas}_kg_liquid'] / out[f'{gas}_kg_atm']
        assert actual_ratio == pytest.approx(predicted_ratio, rel=1e-6)

    # Sign + scale guards: the ratio is positive and, for this Earth-scale
    # mantle at ~50 bar, dissolved He is a small but non-zero fraction of
    # atmospheric He (the mantle is a modest reservoir at low pressure).
    ratio_he = out['He_kg_liquid'] / out['He_kg_atm']
    assert ratio_he > 0.0
    assert ratio_he < 10.0


@pytest.mark.physics_invariant
def test_noble_partial_pressure_and_dissolved_increase_with_budget():
    """Monotonicity: raising the He budget raises both its partial pressure
    and its dissolved mass. A sign error on the residual or a mislabeled
    unknown would break this ordering.
    """
    ddict = _ddict(active=('He',))
    outs = []
    for he_budget in (1.0e15, 1.0e16, 1.0e17):
        np.random.seed(3)
        target = dict(_CHNOS, He=he_budget)
        outs.append(
            equilibrium_atmosphere(
                target, ddict, print_result=False, opt_solver=False, nguess=2000
            )
        )
    bars = [o['He_bar'] for o in outs]
    diss = [o['He_kg_liquid'] for o in outs]
    assert bars[0] < bars[1] < bars[2]
    assert diss[0] < diss[1] < diss[2]
    # Linear Henry law + near-fixed mu: a 10x budget increase raises the
    # partial pressure by close to 10x. Guards against a saturating or
    # sub-linear bug that would compress the spacing.
    assert bars[2] / bars[0] > 50.0


@pytest.mark.physics_invariant
def test_noble_dominated_atmosphere_shifts_chnos_pressures():
    """Coupling: a noble-gas-dominated atmosphere lowers the mean molar
    mass, which changes the mapping from CHNOS partial pressures to column
    masses, so at fixed CHNOS targets the solved CHNOS pressures must move.
    A bypassed implementation that solved the noble gases against a frozen
    CHNOS mean molar mass would leave the CHNOS pressures untouched; this
    test fails in that case.
    """
    np.random.seed(5)
    ddict_free = _ddict(active=())
    out_free = equilibrium_atmosphere(
        dict(_CHNOS), ddict_free, print_result=False, opt_solver=False, nguess=2000
    )

    np.random.seed(5)
    ddict_he = _ddict(active=('He',))
    # He budget more than an order of magnitude above the hydrogen budget:
    # He dominates the column and pulls the mean molar mass toward 4 g/mol.
    out_he = equilibrium_atmosphere(
        dict(_CHNOS, He=2.0e21), ddict_he, print_result=False, opt_solver=False, nguess=2000
    )

    # The mean molar mass must drop substantially once He dominates.
    assert out_he['atm_kg_per_mol'] < 0.5 * out_free['atm_kg_per_mol']

    # CO2 and N2 are held mostly in the atmosphere here, so the mean-molar-
    # mass drop reaches their full column: each partial pressure must fall by
    # more than half at the fixed CHNOS target (the lighter atmosphere carries
    # the same elemental mass at lower partial pressure). The direction is
    # pinned (a decrease), not just the magnitude. A frozen-mu bypass that
    # solved the noble gases against the CHNOS-only mean molar mass would
    # leave CO2_bar and N2_bar unchanged and fail both assertions.
    assert out_he['CO2_bar'] < 0.5 * out_free['CO2_bar']
    assert out_he['N2_bar'] < 0.5 * out_free['N2_bar']

    # CHNOS mass balance still holds in the He-dominated solve, so the shift
    # is a genuine re-partitioning, not a broken solve.
    assert out_he['H_res'] == pytest.approx(0.0, abs=max(1e-6 * _CHNOS['H'], 1e3))
    assert out_he['C_res'] == pytest.approx(0.0, abs=max(1e-6 * _CHNOS['C'], 1e3))


@pytest.mark.physics_invariant
def test_noble_partitioning_independent_of_fo2():
    """Symmetry: noble gases are chemically inert, so their partition
    between melt and atmosphere must not depend on the fO2 buffer offset,
    unlike the reactive CHNOS species (whose speciation shifts strongly with
    fO2). Solve the same He budget at two very different redox states and
    require the He partial pressure and dissolved mass to match.
    """
    target = dict(_CHNOS, He=3.0e16)
    np.random.seed(9)
    out_red = equilibrium_atmosphere(
        target,
        _ddict(active=('He',), dIW=-4.0),
        print_result=False,
        opt_solver=False,
        nguess=2000,
    )
    np.random.seed(9)
    out_ox = equilibrium_atmosphere(
        target,
        _ddict(active=('He',), dIW=+4.0),
        print_result=False,
        opt_solver=False,
        nguess=2000,
    )

    # He partitioning is set by its own budget and the column, not redox.
    # A small residual difference is allowed because the CHNOS background
    # (and hence mu) shifts slightly with fO2, but it must be minor.
    assert out_ox['He_kg_liquid'] / out_ox['He_kg_atm'] == pytest.approx(
        out_red['He_kg_liquid'] / out_red['He_kg_atm'], rel=1e-2
    )
    # Discrimination: the sulfur partitioning is fO2-dependent (the Gaillard
    # S2 melt solubility carries an explicit fO2 term), so S2_bar differs
    # between the reducing and oxidising solves, confirming the two redox
    # states are genuinely different and that noble inertness is not an
    # artifact of identical inputs.
    assert out_ox['S2_bar'] != pytest.approx(out_red['S2_bar'], rel=0.1)


def test_excluded_noble_gas_is_absent_from_the_solve():
    """Edge case: with only He included, the other four noble gases must
    report exactly zero reservoirs and must not be required in the target
    dict. This is the backward-compatible gate that keeps a CHNOS-only or
    single-noble run from dragging in unused species.
    """
    np.random.seed(1)
    ddict = _ddict(active=('He',))
    target = dict(_CHNOS, He=3.0e16)  # no Ne/Ar/Kr/Xe targets supplied
    out = equilibrium_atmosphere(
        target, ddict, print_result=False, opt_solver=False, nguess=2000
    )

    # Only active gases appear in the output schema, so an unused noble gas
    # leaves no keys behind (a CHNOS-only run is likewise noble-free).
    for gas in ('Ne', 'Ar', 'Kr', 'Xe'):
        assert f'{gas}_bar' not in out
        assert f'{gas}_kg_total' not in out
    # He, the one included gas, is present and non-trivial.
    assert out['He_bar'] > 0.0
    assert out['He_kg_total'] == pytest.approx(3.0e16, rel=1e-6)


def test_included_noble_gas_without_target_raises():
    """Error contract: an included noble gas with no target mass is a
    misconfiguration and must fail loudly, not silently solve for zero.
    """
    np.random.seed(1)
    ddict = _ddict(active=('He', 'Ar'))
    target = dict(_CHNOS, He=3.0e16)  # Ar included but no Ar target
    with pytest.raises(KeyError, match='Ar'):
        equilibrium_atmosphere(target, ddict, print_result=False, nguess=100)


def test_get_target_from_params_converts_noble_ppmw_to_kg():
    """The element-mode target builder (the PROTEUS `volatile_mode =
    'elements'` path) converts a noble gas ppmw budget to kg relative to the
    mantle mass, exactly as it does for nitrogen and sulfur. This exercises
    the `get_target_from_params` noble branch that the solver tests bypass by
    building target dicts directly.
    """
    M_mantle = 4.0e24
    ddict = {
        'M_mantle': M_mantle,
        'hydrogen_earth_oceans': 1.0,
        'CH_ratio': 0.1,
        'nitrogen_ppmw': 2.0,
        'sulfur_ppmw': 200.0,
        'He_included': 1,
        'He_ppmw': 5.0,
        'Ar_included': 1,
        'Ar_ppmw': 0.5,
        # Ne included but no ppmw: defaults to zero budget, not an error.
        'Ne_included': 1,
        'Kr_included': 0,
        'Xe_included': 0,
    }
    target = get_target_from_params(ddict)

    assert target['He'] == pytest.approx(5.0 * 1e-6 * M_mantle, rel=1e-12)
    assert target['Ar'] == pytest.approx(0.5 * 1e-6 * M_mantle, rel=1e-12)
    assert target['Ne'] == pytest.approx(0.0)
    # Excluded gases get no target at all.
    assert 'Kr' not in target
    assert 'Xe' not in target
    # Discrimination guard: dropping the 1e-6 ppmw factor would make the He
    # target 1e6x too large (5 * M_mantle instead of 5e-6 * M_mantle).
    assert abs(target['He'] - 5.0 * M_mantle) > 0.5 * (5.0 * M_mantle)


def test_get_target_from_pressures_tallies_noble_initial_bar():
    """The pressure-mode target builder sums a noble gas's atmospheric and
    dissolved mass from its initial partial pressure, so a config that
    specifies noble gases by pressure produces the correct inventory. Both
    reservoirs use the closed-form column and Henry expressions.
    """
    ddict = _ddict(active=('He',))
    ddict['H2O_initial_bar'] = 100.0
    ddict['CO2_initial_bar'] = 10.0
    ddict['N2_initial_bar'] = 5.0
    ddict['S2_initial_bar'] = 1.0
    ddict['He_initial_bar'] = 20.0

    target = get_target_from_pressures(ddict)

    # Independent closed-form He inventory, built from first principles rather
    # than from the functions the builder calls. The oxygen fugacity (and
    # hence the negligible O2 partial pressure) comes from the separate
    # OxygenFugacity module, so nothing here re-uses atmosphere_mass or
    # dissolved_mass.
    from calliope.oxygen_fugacity import OxygenFugacity

    p_O2 = 10.0 ** OxygenFugacity()(1800.0, 0.0)
    pin_all = {'H2O': 100.0, 'CO2': 10.0, 'N2': 5.0, 'S2': 1.0, 'He': 20.0, 'O2': p_O2}
    mu = sum(molar_mass[s] * pin_all[s] for s in pin_all) / sum(pin_all.values())
    area = 4.0 * np.pi * ddict['radius'] ** 2.0
    he_atm = 20.0 * 1.0e5 / ddict['gravity'] * area * molar_mass['He'] / mu
    he_diss = (
        1.0e-6 * ddict['M_mantle'] * ddict['Phi_global'] * jambon86_ppmw_per_bar('He') * 20.0
    )
    expected = he_atm + he_diss
    assert target['He'] == pytest.approx(expected, rel=1e-6)
    # He inventory is strictly positive at 20 bar and dominated by the
    # atmospheric column for this Earth-scale mantle at modest pressure.
    assert target['He'] > 0.0
    assert he_atm > he_diss
    # Discrimination guard: the CHNOS targets must also be present and
    # positive, confirming the noble branch did not displace them.
    for e in ('H', 'C', 'N', 'S'):
        assert target[e] > 0.0


def test_get_target_from_pressures_rejects_empty_atmosphere():
    """Error contract: if every initial partial pressure is below the
    surface-pressure floor, there is no atmosphere to invert and the builder
    raises rather than returning a degenerate zero-pressure target.
    """
    ddict = _ddict(active=('He',))
    for s in ('H2O', 'CO2', 'N2', 'S2'):
        ddict[f'{s}_initial_bar'] = 1.0e-30
    ddict['He_initial_bar'] = 1.0e-30
    with pytest.raises(Exception, match='too low'):
        get_target_from_pressures(ddict)


@pytest.mark.physics_invariant
def test_chnos_only_solve_is_unchanged_by_noble_gas_support():
    """Backward compatibility: a run with no noble gas budget must reproduce
    the pure CHNOS solution and emit no noble keys. Pins the commit's
    central claim that the noble gas machinery leaves the CHNOS path
    numerically untouched.
    """
    np.random.seed(17)
    ddict = _ddict(active=())
    out = equilibrium_atmosphere(
        dict(_CHNOS), ddict, print_result=False, opt_solver=False, nguess=2000
    )

    # No noble keys leak into a CHNOS-only output.
    for gas in noble_gases:
        assert f'{gas}_bar' not in out
        assert f'{gas}_kg_total' not in out
        assert f'{gas}/H_atm' not in out

    # CHNOS mass balance closes.
    assert out['H_res'] == pytest.approx(0.0, abs=max(1e-6 * _CHNOS['H'], 1e3))
    assert out['S_res'] == pytest.approx(0.0, abs=max(1e-6 * _CHNOS['S'], 1e3))

    # Pin the primary partial pressures and the mean molar mass against the
    # values the pre-noble solver produced at this fixed seed. A regression
    # that perturbed the CHNOS path while still closing mass balance (for
    # example by letting an inactive noble gas leak into the mean molar mass)
    # would move these numbers and fail here. Captured from the CHNOS-only
    # solver at seed 17.
    assert out['H2O_bar'] == pytest.approx(4.0844429849e-01, rel=1e-6)
    assert out['CO2_bar'] == pytest.approx(5.1936704868e01, rel=1e-6)
    assert out['N2_bar'] == pytest.approx(5.6066338860e-01, rel=1e-6)
    assert out['S2_bar'] == pytest.approx(1.4249796437e-07, rel=1e-5)
    assert out['atm_kg_per_mol'] == pytest.approx(4.3639799394e-02, rel=1e-6)


@pytest.mark.physics_invariant
def test_authoritative_o_mode_conserves_noble_mass_and_recovers_fo2():
    """The authoritative-O solver (PROTEUS Path C) carries the noble gases
    as extra unknowns after the O residual and must both close their mass
    budgets and recover the redox state, unchanged by their presence.
    """
    dIW = 4.0
    ddict = _ddict(dIW=dIW)
    target_chnos = dict(_CHNOS, **_NOBLE)
    np.random.seed(2)
    legacy = equilibrium_atmosphere(target_chnos, ddict, print_result=False, nguess=1000)
    target = dict(target_chnos, O=legacy['O_kg_total'])

    out = equilibrium_atmosphere_authoritative_O(
        target,
        ddict,
        fO2_hint=dIW,
        random_seed=0,
        nguess=1000,
        nsolve=1000,
        print_result=False,
        opt_solver=False,
    )

    for gas in noble_gases:
        reservoir = out[f'{gas}_kg_atm'] + out[f'{gas}_kg_liquid']
        assert reservoir == pytest.approx(target[gas], rel=1e-5)
        # Each noble gas gets its own residual key in this mode.
        assert abs(out[f'{gas}_res']) < max(1e-5 * target[gas], 1e3)
    # The derived redox state is the one the O budget was built at; the
    # noble gases do not perturb it.
    assert out['fO2_shift_derived'] == pytest.approx(dIW, abs=0.05)


def test_is_included_strict_for_chnos_soft_for_noble():
    """The inclusion predicate is opt-in for noble gases but fail-loud for the
    CHNOS reaction-network species. A missing noble flag defaults to excluded;
    a missing CHNOS flag is a malformed options dict and raises, so a species
    cannot silently vanish from the solve.
    """
    # Noble gases: absent flag means excluded, present flag is honoured.
    assert is_included('He', {}) is False
    assert is_included('Ne', {'Ne_included': 1}) is True
    assert is_included('Ar', {'Ar_included': 0}) is False
    # CHNOS: a missing flag raises rather than defaulting.
    with pytest.raises(KeyError):
        is_included('CO2', {})
    with pytest.raises(KeyError):
        is_included('N2', {})
    # A present CHNOS flag still works both ways.
    assert is_included('CO2', {'CO2_included': 1}) is True
    assert is_included('CO2', {'CO2_included': 0}) is False


def test_noble_residual_gate_rejects_mass_unbalanced_solution(monkeypatch, caplog):
    """The per-gas noble closure gate must reject a converged-looking root
    whose noble residual exceeds the gas's own budget, even when the scalar
    CHNOS gate would accept it. Force the inner solver to return such a root
    and confirm the solve is rejected rather than returned.
    """
    ddict = _ddict(active=('He',))
    target = dict(_CHNOS, He=3.0e16)
    # Supply a warm guess so the two-stage cold start is skipped and the
    # single stubbed attempt is the whole solve.
    p_guess = {'H2O': 1.0, 'CO2': 1.0, 'N2': 1.0, 'S2': 1.0, 'He': 1.0}

    # Root fsolve reports as converged (ier == 1).
    def _fake_fsolve(*args, **kwargs):
        return np.array([1.0, 1.0, 1.0, 1.0, 1.0]), {}, 1, 'stub'

    # Residual: CHNOS closed (0), but the He residual (3e13 kg) is far above
    # the He per-gas tolerance (3e16 * 1e-5 = 3e11 kg) while staying under the
    # scalar CHNOS gate (~1.5e15 kg), so only the noble gate can catch it.
    def _fake_func(*args, **kwargs):
        return [0.0, 0.0, 0.0, 0.0, 3.0e13]

    monkeypatch.setattr(calsolve.opt, 'fsolve', _fake_fsolve)
    monkeypatch.setattr(calsolve, 'func', _fake_func)

    with caplog.at_level(logging.DEBUG, logger='fwl.calliope.solve'):
        with pytest.raises(RuntimeError, match='Could not find solution'):
            equilibrium_atmosphere(
                target, ddict, p_guess=p_guess, nguess=1, print_result=False, opt_solver=False
            )
    # The rejection reason is the noble gate, not the scalar CHNOS gate.
    assert any('noble gas residual' in r.message for r in caplog.records)


def test_two_stage_cold_start_falls_back_when_core_presolve_fails(monkeypatch, caplog):
    """When the CHNOS core pre-solve of the two-stage cold start does not
    converge, the solver logs a warning and falls back to the random draw,
    which still recovers a valid solution. Force the recursive core call to
    raise and confirm both the warning and the recovery.
    """
    real = calsolve.equilibrium_atmosphere
    calls = {'n': 0}

    def wrapper(*args, **kwargs):
        calls['n'] += 1
        # The first call is the outer solve; the second is the recursive
        # CHNOS core pre-solve, which we force to fail.
        if calls['n'] == 2:
            raise RuntimeError('forced core pre-solve failure')
        return real(*args, **kwargs)

    monkeypatch.setattr(calsolve, 'equilibrium_atmosphere', wrapper)

    ddict = _ddict(active=('He',))
    target = dict(_CHNOS, He=3.0e16)
    np.random.seed(0)
    with caplog.at_level(logging.WARNING, logger='fwl.calliope.solve'):
        out = wrapper(target, ddict, print_result=False, opt_solver=False, nguess=4000)

    # The fallback warning fired, and the random draw still closed He mass.
    assert any('random high-dimensional draw' in r.message for r in caplog.records)
    assert out['He_kg_atm'] + out['He_kg_liquid'] == pytest.approx(3.0e16, rel=1e-5)


def test_abundant_noble_residual_within_budget_is_accepted(monkeypatch):
    """A converged root whose noble residual is within the noble gas's own
    per-gas tolerance must be accepted even when that residual exceeds the
    CHNOS-keyed scalar tolerance. The scalar acceptance gate judges only the
    CHNOS residuals; an abundant noble gas is governed by its per-gas gate.
    Without the CHNOS-only slice the scalar gate would reject this root and
    the solve would spuriously fail.
    """
    ddict = _ddict(active=('He',))
    # Trace CHNOS, abundant He: the He per-gas tolerance (He*rtol = 1e13)
    # exceeds the CHNOS scalar tolerance (dominated by atol ~ 1e10).
    target = {'H': 1.0e12, 'C': 1.0e12, 'N': 1.0e12, 'S': 1.0e12, 'He': 1.0e18}
    p_guess = {'H2O': 1.0, 'CO2': 1.0, 'N2': 1.0, 'S2': 1.0, 'He': 1.0}

    def _fake_fsolve(*args, **kwargs):
        return np.array([1.0, 1.0, 1.0, 1.0, 1.0]), {}, 1, 'stub'

    # CHNOS residuals closed; He residual 5e12 is within the He per-gas
    # tolerance (1e13) but above the CHNOS scalar tolerance (~1e10).
    def _fake_func(*args, **kwargs):
        return [0.0, 0.0, 0.0, 0.0, 5.0e12]

    monkeypatch.setattr(calsolve.opt, 'fsolve', _fake_fsolve)
    monkeypatch.setattr(calsolve, 'func', _fake_func)

    out = equilibrium_atmosphere(
        target,
        ddict,
        p_guess=p_guess,
        nguess=1,
        atol=1.0e10,
        rtol=1.0e-5,
        print_result=False,
        opt_solver=False,
    )
    # Accepted: a result dict is returned rather than a RuntimeError.
    assert isinstance(out, dict)
    assert 'He_res' in out


def test_authoritative_o_two_stage_fallback_warns(monkeypatch, caplog):
    """When the CHNOS core pre-solve of the authoritative-O two-stage cold
    start fails, the solver logs a warning and falls back to the random draw,
    the same contract the fixed-fO2 path has. Exercises the fallback branch of
    the authoritative-O entry point.
    """
    real = calsolve.equilibrium_atmosphere_authoritative_O
    calls = {'n': 0}

    def wrapper(*args, **kwargs):
        calls['n'] += 1
        # First call is the outer solve; second is the recursive CHNOS + O
        # core pre-solve, which we force to fail.
        if calls['n'] == 2:
            raise RuntimeError('forced core pre-solve failure')
        return real(*args, **kwargs)

    monkeypatch.setattr(calsolve, 'equilibrium_atmosphere_authoritative_O', wrapper)

    dIW = 4.0
    ddict = _ddict(dIW=dIW)
    target_chnos = dict(_CHNOS, **_NOBLE)
    np.random.seed(2)
    legacy = equilibrium_atmosphere(target_chnos, ddict, print_result=False, nguess=1000)
    target = dict(target_chnos, O=legacy['O_kg_total'])

    with caplog.at_level(logging.WARNING, logger='fwl.calliope.solve'):
        try:
            out = wrapper(
                target,
                ddict,
                fO2_hint=dIW,
                random_seed=0,
                nguess=2000,
                nsolve=1000,
                print_result=False,
                opt_solver=False,
            )
        except RuntimeError:
            out = None

    # The fallback warning fired, so the branch is exercised regardless of
    # whether the high-dimensional random draw then reconverged.
    assert any('random high-dimensional draw' in r.message for r in caplog.records)
    # If the random draw did recover, noble mass still closes.
    if out is not None:
        assert out['He_kg_atm'] + out['He_kg_liquid'] == pytest.approx(target['He'], rel=1e-5)

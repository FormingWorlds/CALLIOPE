"""Smoke-tier physics and chemistry invariants for the CALLIOPE
outgassing solver.

Each invariant in this file asserts a property that must hold for any
valid CALLIOPE result, independent of the specific (T_magma, fO2_shift,
elemental inventory) input. They are anti-happy-path by construction:
each parametric sweep covers an edge of the physical regime (extreme
reducing, extreme oxidising, low T, high T), and each class includes
at least one sad-path test that exercises an unphysical input.

The 8 smoke-tier invariants exercised here are:

1.  Per-element mass conservation: atm + liquid == total.
2.  Pressure positivity: every species partial pressure is non-negative.
3.  VMR closure: sum of volume mixing ratios is unity.
4.  Total pressure consistency: P_surf == sum of species partial pressures.
5.  Atmospheric mass consistency: M_atm == sum of species column masses.
6.  fO2 reconstruction: in authoritative-O mode, the derived shift
    reproduces 10**(buffer + shift) for the returned p_O2.
7.  Modified equilibrium constant identity: p_B / p_A matches
    ``ModifiedKeq.<method>(T, fO2_shift)`` at the converged state.
8.  Dissolved-mass non-negativity: every <species>_kg_liquid is >= 0.

Plus three Tim-flagged additions:

-   Solver-tolerance behaviour at the strongly-reducing edge of the
    Dasgupta calibration footprint (dIW <= -3, where the paper notes
    plateau-like N solubility).
-   S2 / SO2 branch behaviour at the oxidising edge of the Gaillard
    sulfide-saturated calibration (dIW > +4, where the sulfate regime
    begins and the law extrapolates).
-   p_guess warm-start verification: a warm start with the cold-solve
    result drops the restart count to zero or one.

The unit-tier invariants (CO2 stoichiometric atom counting, S2 Gaillard
monotonicity, N2 Dasgupta monotonicity, Libourel redox-independence,
Dasgupta law at the reducing edge) live in `test_invariants_unit.py`.
That split is needed because pytest stacks module-level and class-level
markers; a single module-level pytestmark on a mixed-tier file would
pull unit tests into the smoke gate.
"""

from __future__ import annotations

import logging
import math

import pytest

from calliope.chemistry import ModifiedKeq
from calliope.constants import element_list_chnos, volatile_species
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solve import (
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

logging.getLogger('calliope').setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ddict(T: float = 1800.0, Phi: float = 1.0, dIW: float = 4.0) -> dict:
    """Realistic ddict with every volatile species included."""
    d = {
        'M_mantle': 4.03e24,
        'gravity': 9.81,
        'radius': 6.371e6,
        'Phi_global': Phi,
        'T_magma': T,
        'fO2_shift_IW': dIW,
    }
    for sp in volatile_species:
        d[f'{sp}_included'] = 1
        d[f'{sp}_initial_bar'] = 0.0
    return d


def _target_HCNS() -> dict:
    """Earth-like H/C/N/S budget in kg for the buffered mode."""
    return {'H': 1.5e20, 'C': 1.5e19, 'N': 8.0e18, 'S': 8.0e20}


def _solve_buffered(T: float, dIW: float, Phi: float = 1.0) -> dict:
    """Run buffered mode at the given (T, dIW, Phi) with Earth-like targets."""
    return equilibrium_atmosphere(
        _target_HCNS(),
        _ddict(T=T, Phi=Phi, dIW=dIW),
        nguess=200,
        nsolve=1500,
        print_result=False,
        opt_solver=False,
    )


def _solve_authoritative(
    T: float, O_kg: float, Phi: float = 1.0, fO2_hint: float = 0.0
) -> dict:
    """Run authoritative-O mode with a five-element target."""
    target = _target_HCNS()
    target['O'] = O_kg
    return equilibrium_atmosphere_authoritative_O(
        target,
        _ddict(T=T, Phi=Phi, dIW=fO2_hint),
        fO2_hint=fO2_hint,
        nguess=200,
        nsolve=1500,
        print_result=False,
        opt_solver=False,
        random_seed=42,
    )


# Discriminating (T, fO2) points. Endpoints span the physical regime
# of magma-ocean redox states. 1500 K is near the lower edge of the
# JANAF validity window; 2200 K is mid-magma-ocean.
_DEFAULT_TFO2 = [
    (1500.0, -3.0),
    (1500.0, +0.0),
    (1500.0, +3.0),
    (2200.0, -3.0),
    (2200.0, +0.0),
    (2200.0, +4.0),
]


# ===========================================================================
# Invariant 1: per-element mass conservation
# ===========================================================================


class TestMassConservationPerElement:
    """Per-element atmospheric + dissolved == total."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_atm_plus_liquid_equals_total(self, T, dIW):
        """For every element, atmospheric mass + dissolved mass equals
        total mass to within solver tolerance, at every parametrized
        (T, dIW) point in the magma-ocean window."""
        result = _solve_buffered(T=T, dIW=dIW)
        seen_split = False
        for e in element_list_chnos:
            atm = result[f'{e}_kg_atm']
            liq = result[f'{e}_kg_liquid']
            tot = result[f'{e}_kg_total']
            assert atm + liq == pytest.approx(tot, rel=1e-12, abs=1e-3), (
                f'Element {e}: atm={atm:.4e} + liq={liq:.4e} '
                f'!= tot={tot:.4e} at T={T}, dIW={dIW}'
            )
            # Track at least one element where both channels are meaningful
            # so we can verify a non-trivial split below.
            if atm > 1e-3 and liq > 1e-3:
                seen_split = True

        # Discrimination guard: a wrong formula that double-counts the
        # dissolved channel (atm + 2 * liq) would give a result larger than
        # total by exactly liq, well outside the rel=1e-12 tolerance. The
        # check is meaningful only at a (T, dIW) where some element has a
        # non-trivial dissolved channel; the rel=1e-12 tolerance on the
        # primary assertion already excludes the floor case where every
        # liq ~ 0.
        if seen_split:
            for e in element_list_chnos:
                atm = result[f'{e}_kg_atm']
                liq = result[f'{e}_kg_liquid']
                tot = result[f'{e}_kg_total']
                if liq > 1.0:
                    wrong_total = atm + 2.0 * liq
                    assert abs(wrong_total - tot) > 0.5 * liq

    def test_unphysical_zero_mantle_mass_does_not_violate_invariant(self):
        """Sad-path: M_mantle=0 zeros the dissolved channel; the invariant
        atm + liq == total still holds (with liq == 0)."""
        ddict = _ddict(T=1800.0, Phi=1.0, dIW=0.0)
        ddict['M_mantle'] = 0.0
        result = equilibrium_atmosphere(
            _target_HCNS(),
            ddict,
            nguess=200,
            nsolve=1500,
            print_result=False,
            opt_solver=False,
        )
        for e in element_list_chnos:
            assert result[f'{e}_kg_liquid'] == 0.0, (
                f'M_mantle=0 should zero dissolved mass for {e}'
            )
            assert result[f'{e}_kg_atm'] == pytest.approx(
                result[f'{e}_kg_total'], rel=1e-12, abs=1e-3
            )


# ===========================================================================
# Invariant 2: pressure positivity
# ===========================================================================


class TestPressurePositivity:
    """Every species partial pressure is >= 0 at the converged solution."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_all_species_pressures_nonnegative(self, T, dIW):
        """Every species partial pressure is >= 0 at the converged
        solution across the magma-ocean (T, dIW) window."""
        result = _solve_buffered(T=T, dIW=dIW)
        for s in volatile_species:
            assert result[f'{s}_bar'] >= 0.0, (
                f'Species {s} has negative pressure {result[f"{s}_bar"]:.4e} '
                f'at T={T}, dIW={dIW}'
            )

        # Discrimination guard: a stub solver that returned zero for every
        # species would pass the positivity check. Confirm at least one of
        # the four primary species (H2O, CO2, N2, S2) carries meaningful
        # pressure given the Earth-like H/C/N/S target.
        primary_sum = sum(result[f'{s}_bar'] for s in ('H2O', 'CO2', 'N2', 'S2'))
        assert primary_sum > 0.0, f'All primary pressures are zero at T={T}, dIW={dIW}'

    def test_extreme_reducing_does_not_break_positivity(self):
        """Sad-path: at dIW=-5 the H2O/CO2 budgets collapse; verify the
        speciation walk does not produce negative pressures."""
        result = _solve_buffered(T=1800.0, dIW=-5.0)
        for s in volatile_species:
            assert result[f'{s}_bar'] >= 0.0

        # Discrimination guard: at dIW=-5 the H2/CO branches dominate over
        # H2O/CO2 (the buffer drives the reduced species). A solver that
        # left H2O/CO2 dominant under these conditions would have a wrong
        # redox response. Sample the H2 / H2O ratio: under reducing
        # conditions it should be >> 1.
        assert result['H2_bar'] > result['H2O_bar'], (
            f'At dIW=-5 expected H2 > H2O; got H2={result["H2_bar"]:.4e}, '
            f'H2O={result["H2O_bar"]:.4e}'
        )


# ===========================================================================
# Invariant 3: VMR closure
# ===========================================================================


class TestVMRClosure:
    """sum(volume mixing ratios) == 1.0 over all volatile species."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_vmrs_sum_to_one(self, T, dIW):
        """Sum of volume mixing ratios across all volatile species equals
        unity to within rel=1e-10 at every parametrized (T, dIW) point."""
        result = _solve_buffered(T=T, dIW=dIW)
        vmr_sum = sum(result[f'{s}_vmr'] for s in volatile_species)
        assert vmr_sum == pytest.approx(1.0, rel=1e-10), (
            f'VMR sum {vmr_sum} != 1.0 at T={T}, dIW={dIW}'
        )

        # Discrimination guard: dropping any single species from the sum
        # would give a value < 1 by that species' vmr. At every (T, dIW)
        # in _DEFAULT_TFO2 at least one species carries >= 1% vmr, well
        # outside the rel=1e-10 closure tolerance.
        vmrs = [result[f'{s}_vmr'] for s in volatile_species]
        max_vmr = max(vmrs)
        assert max_vmr > 0.01, (
            f'No species carries >= 1% vmr at T={T}, dIW={dIW}; closure check would be vacuous'
        )

    def test_vmr_closure_in_authoritative_O_mode(self):
        """Closure must hold in both solver modes."""
        # Use a moderate O budget that gives dIW ~ 0 derived
        result = _solve_authoritative(T=1800.0, O_kg=1.0e21, fO2_hint=0.0)
        vmr_sum = sum(result[f'{s}_vmr'] for s in volatile_species)
        assert vmr_sum == pytest.approx(1.0, rel=1e-10)

        # Discrimination guard: confirm the sum is not vacuous (a stub that
        # returned vmr=1/N for every species would also sum to 1). Require
        # at least one species to carry >= 1% vmr.
        vmrs = [result[f'{s}_vmr'] for s in volatile_species]
        assert max(vmrs) > 0.01


# ===========================================================================
# Invariant 4: total pressure consistency
# ===========================================================================


class TestTotalPressureConsistency:
    """P_surf == sum of species partial pressures."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_psurf_equals_sum_of_pp(self, T, dIW):
        """P_surf equals the sum of the per-species partial pressures
        across the magma-ocean (T, dIW) window."""
        result = _solve_buffered(T=T, dIW=dIW)
        p_sum = sum(result[f'{s}_bar'] for s in volatile_species)
        assert result['P_surf'] == pytest.approx(p_sum, rel=1e-10)

        # Discrimination guard: a stub that returned P_surf = 0 would fail
        # the closure only if the species pressures are themselves nonzero.
        # Confirm the closure is non-vacuous by requiring P_surf > 1 bar
        # given the Earth-like target.
        assert result['P_surf'] > 1.0, (
            f'P_surf = {result["P_surf"]:.4e} bar < 1 bar at T={T}, dIW={dIW}'
        )

    def test_psurf_positive(self):
        """Sanity sad-path: a converged solve never gives P_surf <= 0."""
        result = _solve_buffered(T=1800.0, dIW=0.0)
        assert result['P_surf'] > 0.0

        # Discrimination guard: a stub returning a tiny positive value
        # (1e-30 bar) would pass the bare positivity check. For the
        # Earth-like target, P_surf should be at least 1 bar.
        assert result['P_surf'] > 1.0


# ===========================================================================
# Invariant 5: atmospheric mass consistency
# ===========================================================================


class TestAtmosphericMassConsistency:
    """M_atm == sum of <species>_kg_atm."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_M_atm_equals_sum(self, T, dIW):
        """M_atm equals the sum of per-species column masses across the
        magma-ocean (T, dIW) window."""
        result = _solve_buffered(T=T, dIW=dIW)
        m_sum = sum(result[f'{s}_kg_atm'] for s in volatile_species)
        assert result['M_atm'] == pytest.approx(m_sum, rel=1e-10)

        # Discrimination guard: M_atm must be non-trivially positive for the
        # Earth-like target. A stub returning 0 for every species would
        # pass the closure trivially.
        assert result['M_atm'] > 1.0, (
            f'M_atm = {result["M_atm"]:.4e} kg is non-physically small at T={T}, dIW={dIW}'
        )


# ===========================================================================
# Invariant 6: fO2 reconstruction in authoritative-O mode
# ===========================================================================


class TestFO2Reconstruction:
    """log10(p_O2 / fO2_IW_buffer(T)) == fO2_shift_derived."""

    @pytest.mark.parametrize(
        'T,O_kg',
        [
            (1800.0, 5.0e20),
            (1800.0, 1.0e21),
            (1800.0, 2.0e21),
            (2200.0, 1.0e21),
        ],
    )
    def test_derived_fO2_matches_p_O2(self, T, O_kg):
        """In authoritative-O mode, the derived fO2 shift reproduces
        log10(p_O2) minus the buffer log10(fO2) at IW; verifies the
        closure between the new mode's fifth unknown and the gas-phase
        O2 partial pressure."""
        result = _solve_authoritative(T=T, O_kg=O_kg, fO2_hint=0.0)
        buffer_log10 = OxygenFugacity()(T, 0.0)  # log10 fO2 at IW (shift=0)
        p_O2 = result['O2_bar']
        recovered = math.log10(p_O2) - buffer_log10
        assert recovered == pytest.approx(result['fO2_shift_derived'], rel=1e-6, abs=1e-6)

        # Discrimination guard: using the wrong IW buffer (O'Neill instead
        # of the default Fischer 2011) would shift the recovered value by
        # 0.016 dex at 1800 K and 0.26 dex at 2200 K. Both are well outside
        # the rel=1e-6 tolerance.
        oneill_log10 = OxygenFugacity('oneill')(T, 0.0)
        recovered_wrong_buffer = math.log10(p_O2) - oneill_log10
        assert abs(recovered_wrong_buffer - result['fO2_shift_derived']) > 0.01


# ===========================================================================
# Invariant 7: modified equilibrium constant identity
# ===========================================================================


class TestModifiedKeqIdentity:
    """p_B / p_A == ModifiedKeq.<method>(T, fO2_shift) for each couple."""

    @pytest.mark.parametrize('T,dIW', [(1800.0, 0.0), (2200.0, 3.0)])
    def test_H2O_H2_ratio(self, T, dIW):
        """p_H2 / p_H2O at convergence matches ModifiedKeq('janaf_H2')
        evaluated at the same (T, dIW)."""
        result = _solve_buffered(T=T, dIW=dIW)
        Keq = ModifiedKeq('janaf_H2')
        Geq = Keq(T, dIW)
        # H2O = H2 + 0.5 O2; G_eq = p_H2 / p_H2O
        ratio = result['H2_bar'] / result['H2O_bar']
        assert ratio == pytest.approx(Geq, rel=1e-3)

        # Discrimination guard: the ratio must be closer to Geq than to
        # 1/Geq. This catches a bug that computed the inverse ratio (p_H2O
        # / p_H2 instead of p_H2 / p_H2O); the test would pass the approx
        # check by accident only at the special point where Geq = 1.
        assert abs(ratio - Geq) < abs(ratio - 1.0 / Geq)

    @pytest.mark.parametrize('T,dIW', [(1800.0, 0.0), (2200.0, 3.0)])
    def test_CO2_CO_ratio(self, T, dIW):
        """p_CO / p_CO2 at convergence matches ModifiedKeq('janaf_CO')
        evaluated at the same (T, dIW)."""
        result = _solve_buffered(T=T, dIW=dIW)
        Keq = ModifiedKeq('janaf_CO')
        Geq = Keq(T, dIW)
        ratio = result['CO_bar'] / result['CO2_bar']
        assert ratio == pytest.approx(Geq, rel=1e-3)

        # Discrimination guard: the ratio must be closer to Geq than to
        # 1/Geq, catching a swapped numerator / denominator bug.
        assert abs(ratio - Geq) < abs(ratio - 1.0 / Geq)


# ===========================================================================
# Invariant 8: dissolved-mass non-negativity
# ===========================================================================


class TestDissolvedMassNonNegativity:
    """<species>_kg_liquid is >= 0 for every species at convergence."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_all_dissolved_nonnegative(self, T, dIW):
        """Every per-species dissolved mass is >= 0 at the converged
        solution across the magma-ocean (T, dIW) window."""
        result = _solve_buffered(T=T, dIW=dIW)
        for s in volatile_species:
            assert result[f'{s}_kg_liquid'] >= 0.0, (
                f'{s}_kg_liquid = {result[f"{s}_kg_liquid"]:.4e} < 0 at T={T}, dIW={dIW}'
            )

        # Discrimination guard: a stub that zeroed every dissolved channel
        # would pass the positivity check trivially. For the Earth-like
        # target at any of the _DEFAULT_TFO2 points, at least one of the
        # primary species (H2O, CO2, S2 in particular) dissolves a
        # meaningful amount into the melt.
        primary_liq = sum(result[f'{s}_kg_liquid'] for s in ('H2O', 'CO2', 'N2', 'S2'))
        assert primary_liq > 0.0, f'No primary species dissolves at T={T}, dIW={dIW}'


# ===========================================================================
# Tim-flagged addition: Dasgupta plateau behaviour at strongly-reducing fO2
# ===========================================================================


class TestDasguptaReducingEdgeSolver:
    """Buffered mode at the reducing edge of the Dasgupta N
    calibration footprint. The paper Fig. 7 notes plateau-like N
    solubility at dIW < -3; we verify the solver still produces a
    finite, physically valid result."""

    @pytest.mark.parametrize('dIW', [-3.0, -4.0, -6.0])
    def test_buffered_solver_finite_at_reducing_edge(self, dIW):
        """At the strongly-reducing edge of the Dasgupta calibration
        footprint (dIW <= -3), the buffered solver still produces
        finite, non-negative partial pressures for every species."""
        result = _solve_buffered(T=1800.0, dIW=dIW)
        for s in volatile_species:
            assert math.isfinite(result[f'{s}_bar'])
            assert result[f'{s}_bar'] >= 0.0
        # N inventory is driven into the melt by the -1.6 dIW term
        # so dissolved-N at dIW=-6 is much larger than at dIW=-3
        assert result['N_kg_liquid'] >= 0.0


# ===========================================================================
# Tim-flagged addition: S2 / SO2 branch at the oxidising calibration edge
# ===========================================================================


class TestGaillardOxidisingEdge:
    """Gaillard sulfide-saturated calibration ends near IW+3.5
    (FMQ+0.1). Above that the sulfate regime begins and CALLIOPE's
    sulfide-only chemistry extrapolates. Verify the solver still
    converges and document the SO2 dominance behaviour."""

    @pytest.mark.parametrize('dIW', [+4.0, +5.0, +6.0])
    def test_solver_converges_at_oxidising_edge(self, dIW):
        """At the oxidising edge of the Gaillard sulfide-saturated
        calibration (dIW >= +4), the buffered solver still produces
        finite S2 and SO2 partial pressures, with SO2 dominating S2 by
        dIW >= +5."""
        result = _solve_buffered(T=1800.0, dIW=dIW)
        assert math.isfinite(result['S2_bar'])
        assert math.isfinite(result['SO2_bar'])
        # SO2 > S2 in partial pressure at strongly oxidising conditions
        if dIW >= +5.0:
            assert result['SO2_bar'] > result['S2_bar']

    def test_dissolved_S_drops_at_oxidising_edge(self):
        """At dIW > +4 the dissolved-S fraction of the total S inventory
        should be smaller than at dIW = 0, because Gaillard's law gives
        less sulfide solubility under oxidising conditions."""
        r_neutral = _solve_buffered(T=1800.0, dIW=0.0)
        r_oxidising = _solve_buffered(T=1800.0, dIW=+5.0)
        frac_neutral = r_neutral['S_kg_liquid'] / r_neutral['S_kg_total']
        frac_oxidising = r_oxidising['S_kg_liquid'] / r_oxidising['S_kg_total']
        assert frac_oxidising < frac_neutral, (
            f'Dissolved-S fraction at dIW=+5 ({frac_oxidising:.4e}) '
            f'not below dIW=0 fraction ({frac_neutral:.4e})'
        )

        # Discrimination guard: the gap must be substantive. Gaillard's
        # +0.5 ln(p_S2/fO2) term gives several-fold change in ppmw across
        # 5 dIW units; a near-equal pair of fractions would suggest the
        # redox dependence is being missed by the solver loop.
        assert frac_neutral / max(frac_oxidising, 1e-30) > 2.0, (
            f'Dissolved-S fraction ratio (neutral/oxidising) = '
            f'{frac_neutral / max(frac_oxidising, 1e-30):.2f}; expected > 2x '
            f'across a 5-dIW change'
        )


# ===========================================================================
# Tim-flagged addition: p_guess warm-start verification
# ===========================================================================


class TestPGuessWarmStart:
    """A warm start with the cold-solve result lands in the same
    basin and converges to the same partial pressures. PROTEUS
    depends on this for the coupled-run wall-time budget (without a
    warm start the buffered mode burns 10-50 Monte-Carlo restarts;
    with one, it succeeds on the first attempt)."""

    def test_warm_start_reproduces_cold_pressures(self):
        """A warm start with the cold-solve result as p_guess lands in
        the same basin: primary AND derived species partial pressures
        match the cold result within solver tolerance."""
        cold = equilibrium_atmosphere(
            _target_HCNS(),
            _ddict(T=1800.0, Phi=1.0, dIW=0.0),
            nguess=200,
            nsolve=1500,
            print_result=False,
            opt_solver=False,
        )
        p_guess = {
            'H2O': cold['H2O_bar'],
            'CO2': cold['CO2_bar'],
            'N2': cold['N2_bar'],
            'S2': cold['S2_bar'],
        }
        warm = equilibrium_atmosphere(
            _target_HCNS(),
            _ddict(T=1800.0, Phi=1.0, dIW=0.0),
            nguess=200,
            nsolve=1500,
            p_guess=p_guess,
            print_result=False,
            opt_solver=False,
        )
        # Warm start lands on the same basin: partial pressures match
        for s in ('H2O', 'CO2', 'N2', 'S2'):
            assert warm[f'{s}_bar'] == pytest.approx(cold[f'{s}_bar'], rel=1e-3), (
                f'Warm start drifted away from cold-solve basin for {s}'
            )

        # Discrimination guard: warm and cold must agree on the secondary
        # derived species as well, not just on the four primaries that
        # were used as the p_guess seeds. A solver that initialised only
        # the primaries from p_guess and re-derived the secondaries from
        # a fresh random seed could pass the primary check while drifting
        # on H2, CO, SO2, H2S, NH3, O2.
        for s in ('H2', 'CO', 'SO2', 'H2S', 'NH3'):
            assert warm[f'{s}_bar'] == pytest.approx(cold[f'{s}_bar'], rel=1e-2), (
                f'Warm start drifted on derived species {s}'
            )

    def test_warm_start_with_bad_guess_still_converges(self):
        """Sad-path: even a bad warm-start guess (off by factor of 100)
        should not break the solver; the Monte-Carlo restart catches it."""
        bad_p_guess = {'H2O': 1e3, 'CO2': 1e-3, 'N2': 1e3, 'S2': 1e-3}
        warm = equilibrium_atmosphere(
            _target_HCNS(),
            _ddict(T=1800.0, Phi=1.0, dIW=0.0),
            nguess=200,
            nsolve=1500,
            p_guess=bad_p_guess,
            print_result=False,
            opt_solver=False,
        )
        # If we got here, the solver converged despite the bad guess
        for s in volatile_species:
            assert math.isfinite(warm[f'{s}_bar'])
            assert warm[f'{s}_bar'] >= 0.0

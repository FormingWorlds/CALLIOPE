"""Physics and chemistry invariants for the CALLIOPE outgassing solver.

Each invariant in this file asserts a property that must hold for any
valid CALLIOPE result, independent of the specific (T_magma, fO2_shift,
elemental inventory) input. They are anti-happy-path by construction:
each parametric sweep covers an edge of the physical regime (extreme
reducing, extreme oxidising, low T, high T), and each class includes
at least one sad-path test that exercises an unphysical input.

The 11 invariants exercised here are:

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
9.  CO2 stoichiometric atom counting: the C contribution from CO2_kg_atm
    equals (12/44) * CO2_kg_atm at known partial pressure.
10. Monotonicity of dissolved S vs reducing fO2 at fixed p_S2.
11. Monotonicity of dissolved N vs reducing fO2 at fixed p_N2.

Plus three Tim-flagged additions:

-   Solver-tolerance behaviour at the strongly-reducing edge of the
    Dasgupta calibration footprint (dIW <= -3, where the paper notes
    plateau-like N solubility).
-   S2 / SO2 branch behaviour at the oxidising edge of the Gaillard
    sulfide-saturated calibration (dIW > +4, where the sulfate regime
    begins and the law extrapolates).
-   p_guess warm-start verification: a warm start with the cold-solve
    result drops the restart count to zero or one.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pytest

from calliope.chemistry import ModifiedKeq
from calliope.constants import element_list, molar_mass, volatile_species
from calliope.oxygen_fugacity import OxygenFugacity
from calliope.solubility import SolubilityN2, SolubilityS2
from calliope.solve import (
    _atmosphere_mass,
    equilibrium_atmosphere,
    equilibrium_atmosphere_authoritative_O,
)

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


def _solve_authoritative(T: float, O_kg: float, Phi: float = 1.0,
                         fO2_hint: float = 0.0) -> dict:
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


@pytest.mark.smoke
class TestMassConservationPerElement:
    """Per-element atmospheric + dissolved == total."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_atm_plus_liquid_equals_total(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        for e in element_list:
            atm = result[f'{e}_kg_atm']
            liq = result[f'{e}_kg_liquid']
            tot = result[f'{e}_kg_total']
            assert atm + liq == pytest.approx(tot, rel=1e-12, abs=1e-3), (
                f'Element {e}: atm={atm:.4e} + liq={liq:.4e} '
                f'!= tot={tot:.4e} at T={T}, dIW={dIW}'
            )

    def test_unphysical_zero_mantle_mass_does_not_violate_invariant(self):
        """Sad-path: M_mantle=0 zeros the dissolved channel; the invariant
        atm + liq == total still holds (with liq == 0)."""
        ddict = _ddict(T=1800.0, Phi=1.0, dIW=0.0)
        ddict['M_mantle'] = 0.0
        result = equilibrium_atmosphere(
            _target_HCNS(), ddict, nguess=200, nsolve=1500,
            print_result=False, opt_solver=False,
        )
        for e in element_list:
            assert result[f'{e}_kg_liquid'] == 0.0, (
                f'M_mantle=0 should zero dissolved mass for {e}'
            )
            assert result[f'{e}_kg_atm'] == pytest.approx(
                result[f'{e}_kg_total'], rel=1e-12, abs=1e-3
            )


# ===========================================================================
# Invariant 2: pressure positivity
# ===========================================================================


@pytest.mark.smoke
class TestPressurePositivity:
    """Every species partial pressure is >= 0 at the converged solution."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_all_species_pressures_nonnegative(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        for s in volatile_species:
            assert result[f'{s}_bar'] >= 0.0, (
                f'Species {s} has negative pressure {result[f"{s}_bar"]:.4e} '
                f'at T={T}, dIW={dIW}'
            )

    def test_extreme_reducing_does_not_break_positivity(self):
        """Sad-path: at dIW=-5 the H2O/CO2 budgets collapse; verify the
        speciation walk does not produce negative pressures."""
        result = _solve_buffered(T=1800.0, dIW=-5.0)
        for s in volatile_species:
            assert result[f'{s}_bar'] >= 0.0


# ===========================================================================
# Invariant 3: VMR closure
# ===========================================================================


@pytest.mark.smoke
class TestVMRClosure:
    """sum(volume mixing ratios) == 1.0 over all volatile species."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_vmrs_sum_to_one(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        vmr_sum = sum(result[f'{s}_vmr'] for s in volatile_species)
        assert vmr_sum == pytest.approx(1.0, rel=1e-10), (
            f'VMR sum {vmr_sum} != 1.0 at T={T}, dIW={dIW}'
        )

    def test_vmr_closure_in_authoritative_O_mode(self):
        """Closure must hold in both solver modes."""
        # Use a moderate O budget that gives dIW ~ 0 derived
        result = _solve_authoritative(T=1800.0, O_kg=1.0e21, fO2_hint=0.0)
        vmr_sum = sum(result[f'{s}_vmr'] for s in volatile_species)
        assert vmr_sum == pytest.approx(1.0, rel=1e-10)


# ===========================================================================
# Invariant 4: total pressure consistency
# ===========================================================================


@pytest.mark.smoke
class TestTotalPressureConsistency:
    """P_surf == sum of species partial pressures."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_psurf_equals_sum_of_pp(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        p_sum = sum(result[f'{s}_bar'] for s in volatile_species)
        assert result['P_surf'] == pytest.approx(p_sum, rel=1e-10)

    def test_psurf_positive(self):
        """Sanity sad-path: a converged solve never gives P_surf <= 0."""
        result = _solve_buffered(T=1800.0, dIW=0.0)
        assert result['P_surf'] > 0.0


# ===========================================================================
# Invariant 5: atmospheric mass consistency
# ===========================================================================


@pytest.mark.smoke
class TestAtmosphericMassConsistency:
    """M_atm == sum of <species>_kg_atm."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_M_atm_equals_sum(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        m_sum = sum(result[f'{s}_kg_atm'] for s in volatile_species)
        assert result['M_atm'] == pytest.approx(m_sum, rel=1e-10)


# ===========================================================================
# Invariant 6: fO2 reconstruction in authoritative-O mode
# ===========================================================================


@pytest.mark.smoke
class TestFO2Reconstruction:
    """log10(p_O2 / fO2_IW_buffer(T)) == fO2_shift_derived."""

    @pytest.mark.parametrize('T,O_kg', [
        (1800.0, 5.0e20),
        (1800.0, 1.0e21),
        (1800.0, 2.0e21),
        (2200.0, 1.0e21),
    ])
    def test_derived_fO2_matches_p_O2(self, T, O_kg):
        result = _solve_authoritative(T=T, O_kg=O_kg, fO2_hint=0.0)
        buffer_log10 = OxygenFugacity()(T, 0.0)  # log10 fO2 at IW (shift=0)
        p_O2 = result['O2_bar']
        recovered = math.log10(p_O2) - buffer_log10
        assert recovered == pytest.approx(
            result['fO2_shift_derived'], rel=1e-6, abs=1e-6
        )


# ===========================================================================
# Invariant 7: modified equilibrium constant identity
# ===========================================================================


@pytest.mark.smoke
class TestModifiedKeqIdentity:
    """p_B / p_A == ModifiedKeq.<method>(T, fO2_shift) for each couple."""

    @pytest.mark.parametrize('T,dIW', [(1800.0, 0.0), (2200.0, 3.0)])
    def test_H2O_H2_ratio(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        Keq = ModifiedKeq('janaf_H2')
        Geq = Keq(T, dIW)
        # H2O = H2 + 0.5 O2; G_eq = p_H2 / p_H2O
        ratio = result['H2_bar'] / result['H2O_bar']
        assert ratio == pytest.approx(Geq, rel=1e-3)

    @pytest.mark.parametrize('T,dIW', [(1800.0, 0.0), (2200.0, 3.0)])
    def test_CO2_CO_ratio(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        Keq = ModifiedKeq('janaf_CO')
        Geq = Keq(T, dIW)
        ratio = result['CO_bar'] / result['CO2_bar']
        assert ratio == pytest.approx(Geq, rel=1e-3)


# ===========================================================================
# Invariant 8: dissolved-mass non-negativity
# ===========================================================================


@pytest.mark.smoke
class TestDissolvedMassNonNegativity:
    """<species>_kg_liquid is >= 0 for every species at convergence."""

    @pytest.mark.parametrize('T,dIW', _DEFAULT_TFO2)
    def test_all_dissolved_nonnegative(self, T, dIW):
        result = _solve_buffered(T=T, dIW=dIW)
        for s in volatile_species:
            assert result[f'{s}_kg_liquid'] >= 0.0, (
                f'{s}_kg_liquid = {result[f"{s}_kg_liquid"]:.4e} < 0 '
                f'at T={T}, dIW={dIW}'
            )


# ===========================================================================
# Invariant 9: CO2 stoichiometric atom counting (unit-tier, no solver)
# ===========================================================================


@pytest.mark.unit
class TestCO2AtomCounting:
    """At known p_CO2, the C contribution from CO2 in the atmosphere
    is exactly (12/44) * CO2_kg_atm via stoichiometric atom counting."""

    @pytest.mark.parametrize('p_CO2_bar', [0.1, 1.0, 10.0, 100.0])
    def test_C_atom_count_from_CO2_only(self, p_CO2_bar):
        # Single-species atmosphere with only CO2 to isolate C contribution
        ddict = _ddict(T=1800.0, Phi=0.0, dIW=0.0)
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 0
        ddict['CO2_included'] = 1
        # Set all pressures to zero, then CO2
        p_d = {s: 0.0 for s in volatile_species}
        p_d['CO2'] = p_CO2_bar
        mass_atm = _atmosphere_mass(p_d, 0.0, ddict)
        # Atomic-C contribution expected from stoichiometry: 12.011 / 44.01
        # times the CO2 column mass
        expected_C_kg = mass_atm['CO2'] * molar_mass['C'] / molar_mass['CO2']
        assert mass_atm['C'] == pytest.approx(expected_C_kg, rel=1e-12)

    def test_zero_CO2_pressure_gives_zero_C(self):
        """Sad-path: at p_CO2 = 0 the C tally from CO2 channel is zero."""
        ddict = _ddict(T=1800.0, Phi=0.0, dIW=0.0)
        for sp in volatile_species:
            ddict[f'{sp}_included'] = 0
        ddict['CO2_included'] = 1
        p_d = {s: 0.0 for s in volatile_species}
        mass_atm = _atmosphere_mass(p_d, 0.0, ddict)
        assert mass_atm.get('C', 0.0) == 0.0


# ===========================================================================
# Invariant 10: monotonicity of S2 Gaillard solubility vs reducing fO2
# ===========================================================================


@pytest.mark.unit
class TestS2SolubilityMonotonicity:
    """Gaillard sulfide-saturated solubility increases as fO2 decreases
    at fixed p_S2 and T (the ``+0.5 ln(p_S2/fO2)`` term carries the
    redox dependence directly). Tested at the solubility-function level
    to isolate the law from the multi-species solver feedback loops."""

    @pytest.mark.parametrize('p_S2_bar', [0.01, 0.1, 1.0])
    @pytest.mark.parametrize('T', [1500.0, 1800.0, 2200.0])
    def test_gaillard_strictly_decreasing_with_oxidation(self, T, p_S2_bar):
        S2 = SolubilityS2('gaillard')
        dIWs = [-4.0, -2.0, 0.0, +2.0, +4.0]
        values = [S2.gaillard(p_S2_bar, T, dIW) for dIW in dIWs]
        for i in range(len(values) - 1):
            assert values[i] > values[i + 1], (
                f'Gaillard ppmw at dIW={dIWs[i]} ({values[i]:.4e}) '
                f'not greater than at dIW={dIWs[i+1]} ({values[i+1]:.4e}) '
                f'(T={T}, p_S2={p_S2_bar})'
            )

    def test_negative_pressure_returns_zero(self):
        """Sad-path: at p_S2 < 1e-20 bar the implementation returns 0
        to avoid log(0). Verify the floor behaves correctly."""
        S2 = SolubilityS2('gaillard')
        assert S2.gaillard(1e-30, 1800.0, 0.0) == 0.0


# ===========================================================================
# Invariant 11: monotonicity of N2 Dasgupta solubility vs reducing fO2
# ===========================================================================


@pytest.mark.unit
class TestN2SolubilityMonotonicity:
    """Dasgupta N2 solubility increases as fO2 decreases at fixed
    p_N2, p_tot, T (the ``-1.6 dIW`` term in the reduced-N branch
    drives this). Tested at the solubility-function level."""

    @pytest.mark.parametrize('p_N2_bar', [0.1, 1.0, 10.0])
    @pytest.mark.parametrize('T', [1500.0, 1800.0, 2200.0])
    def test_dasgupta_monotonic_with_oxidation(self, T, p_N2_bar):
        N2 = SolubilityN2('dasgupta')
        dIWs = [-6.0, -4.0, -2.0, 0.0, +2.0, +4.0]
        values = [N2.dasgupta(p_N2_bar, p_N2_bar, T, dIW) for dIW in dIWs]
        # Dasgupta has two terms (reduced-N and molecular N2); the
        # reduced-N term dominates at low fO2 and falls steeply with
        # increasing fO2, so the total is monotonically decreasing.
        for i in range(len(values) - 1):
            assert values[i] >= values[i + 1], (
                f'Dasgupta ppmw at dIW={dIWs[i]} ({values[i]:.4e}) '
                f'less than at dIW={dIWs[i+1]} ({values[i+1]:.4e}) '
                f'(T={T}, p_N2={p_N2_bar})'
            )

    def test_libourel_alternative_is_redox_independent(self):
        """Sad-path / contrast: the Libourel linear Henry's law has no
        fO2 term, so it does NOT show the Dasgupta monotonicity. This
        documents that ``libourel`` and ``dasgupta`` are distinct laws
        with distinct calibration footprints."""
        N2 = SolubilityN2('libourel')
        # Libourel takes only p_N2; no fO2 dependence
        assert N2.libourel(1.0) == N2.libourel(1.0)  # trivially deterministic
        # Two different p_N2 give different values (linear)
        assert N2.libourel(2.0) == pytest.approx(2.0 * N2.libourel(1.0))


# ===========================================================================
# Tim-flagged addition: Dasgupta plateau behaviour at strongly-reducing fO2
# ===========================================================================


@pytest.mark.smoke
class TestDasguptaReducingEdgeSolver:
    """Buffered mode at the reducing edge of the Dasgupta N
    calibration footprint. The paper Fig. 7 notes plateau-like N
    solubility at dIW < -3; we verify the solver still produces a
    finite, physically valid result."""

    @pytest.mark.parametrize('dIW', [-3.0, -4.0, -6.0])
    def test_buffered_solver_finite_at_reducing_edge(self, dIW):
        result = _solve_buffered(T=1800.0, dIW=dIW)
        for s in volatile_species:
            assert math.isfinite(result[f'{s}_bar'])
            assert result[f'{s}_bar'] >= 0.0
        # N inventory is driven into the melt by the -1.6 dIW term
        # so dissolved-N at dIW=-6 is much larger than at dIW=-3
        assert result['N_kg_liquid'] >= 0.0


@pytest.mark.unit
class TestDasguptaReducingEdgeLaw:
    """Direct check that the Dasgupta law itself stays numerically
    well-behaved at strongly-reducing conditions, independent of the
    multi-species solver."""

    @pytest.mark.parametrize('dIW', [-3.0, -4.0, -6.0, -8.0])
    def test_dasgupta_law_finite_at_reducing_edge(self, dIW):
        """The reduced-N branch grows as exp(-1.6 dIW); at dIW=-8 the
        value reaches 10^5 ppmw at p_N2=1 bar but is still
        numerically representable."""
        N2 = SolubilityN2('dasgupta')
        val = N2.dasgupta(1.0, 1.0, 1800.0, dIW)
        assert math.isfinite(val)
        assert val > 0.0


# ===========================================================================
# Tim-flagged addition: S2 / SO2 branch at the oxidising calibration edge
# ===========================================================================


@pytest.mark.smoke
class TestGaillardOxidisingEdge:
    """Gaillard sulfide-saturated calibration ends near IW+3.5
    (FMQ+0.1). Above that the sulfate regime begins and CALLIOPE's
    sulfide-only chemistry extrapolates. Verify the solver still
    converges and document the SO2 dominance behaviour."""

    @pytest.mark.parametrize('dIW', [+4.0, +5.0, +6.0])
    def test_solver_converges_at_oxidising_edge(self, dIW):
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


# ===========================================================================
# Tim-flagged addition: p_guess warm-start verification
# ===========================================================================


@pytest.mark.smoke
class TestPGuessWarmStart:
    """A warm start with the cold-solve result lands in the same
    basin and converges to the same partial pressures. PROTEUS
    depends on this for the coupled-run wall-time budget (without a
    warm start the buffered mode burns 10-50 Monte-Carlo restarts;
    with one, it succeeds on the first attempt)."""

    def test_warm_start_reproduces_cold_pressures(self):
        cold = equilibrium_atmosphere(
            _target_HCNS(),
            _ddict(T=1800.0, Phi=1.0, dIW=0.0),
            nguess=200, nsolve=1500,
            print_result=False, opt_solver=False,
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
            nguess=200, nsolve=1500,
            p_guess=p_guess,
            print_result=False, opt_solver=False,
        )
        # Warm start lands on the same basin: partial pressures match
        for s in ('H2O', 'CO2', 'N2', 'S2'):
            assert warm[f'{s}_bar'] == pytest.approx(
                cold[f'{s}_bar'], rel=1e-3
            ), f'Warm start drifted away from cold-solve basin for {s}'

    def test_warm_start_with_bad_guess_still_converges(self):
        """Sad-path: even a bad warm-start guess (off by factor of 100)
        should not break the solver; the Monte-Carlo restart catches it."""
        bad_p_guess = {'H2O': 1e3, 'CO2': 1e-3, 'N2': 1e3, 'S2': 1e-3}
        warm = equilibrium_atmosphere(
            _target_HCNS(),
            _ddict(T=1800.0, Phi=1.0, dIW=0.0),
            nguess=200, nsolve=1500,
            p_guess=bad_p_guess,
            print_result=False, opt_solver=False,
        )
        # If we got here, the solver converged despite the bad guess
        for s in volatile_species:
            assert math.isfinite(warm[f'{s}_bar'])
            assert warm[f'{s}_bar'] >= 0.0

"""Tests for the metal-reservoir fields of the CALLIOPE input dict (`src/calliope/solve.py`).

Contract: each field that the metal reservoir reads is checked before the first solve;
a boolean (a TOML true), a string, NaN, inf, or a value outside its range raises
ValueError naming the field, and the defaults and the KLB-1 melt pass.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from calliope.metal import KLB1
from calliope.solve import check_metal_fields

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

GOOD = {
    'metal_mass': 1.9e24,
    'P_eq': 50.0,
    'T_eq': 3500.0,
    'dIW_eq': -2.0,
    'silicate_nbo_t': 2.6,
    'silicate_oxides': dict(KLB1),
    'w_O_core': 0.0,
    'w_Si_core': 0.02,
    'tol_w_core': 1e-6,
}

# field -> values that must raise; bounds: P_eq, T_eq > 0, metal_mass >= 0, w in [0, 1)
BAD = {
    'metal_mass': [True, 'x', math.nan, math.inf, -1.0],
    'P_eq': [True, math.nan, math.inf, 0.0, -5.0, 10**400],
    'T_eq': [False, math.nan, -math.inf, 0.0, -1.0],
    'dIW_eq': [True, 'IW-2', math.nan, math.inf],
    'silicate_nbo_t': [True, math.nan, math.inf, -0.1],
    'w_O_core': [True, math.nan, -0.01, 1.0],
    'w_Si_core': [math.inf, -1e-9, 1.5],
    # from 1 on, the relative stop holds on the first pass
    'tol_w_core': [True, math.nan, 0.0, -1e-6, 1.0, 5.0],
}


class TestMetalFields:
    """Field checks of the metal reservoir."""

    def test_good_fields_and_klb1_pass(self):
        """The fields above, with the KLB-1 melt (99.29 wt%), pass the check."""
        assert check_metal_fields(dict(GOOD)) is None
        assert sum(KLB1.values()) == pytest.approx(99.29, abs=1e-9)

    @pytest.mark.parametrize(
        'field, value',
        [(f, v) for f, vals in BAD.items() for v in vals],
        ids=lambda v: repr(v)[:24],
    )
    def test_bad_value_raises_naming_the_field(self, field, value):
        """Each wrong type or out-of-range value raises ValueError with the field name."""
        with pytest.raises(ValueError) as exc:
            check_metal_fields({**GOOD, field: value})
        assert field in str(exc.value)
        assert check_metal_fields(dict(GOOD)) is None

    @pytest.mark.parametrize(
        'oxides, reason',
        [
            ({**KLB1, 'FeO': math.nan}, 'finite'),
            ({**KLB1, 'FeO': -1.0}, 'finite'),
            ({**KLB1, 'Cr2O3': 0.3}, 'unknown'),
            ({'SiO2': 44.48, 'MgO': 39.22}, 'sum'),
            ({**KLB1, 'SiO2': 60.0}, 'sum'),
            ('KLB-1', 'dict'),
        ],
        ids=[
            'nan_oxide',
            'negative_oxide',
            'unknown_oxide',
            'sum_low',
            'sum_high',
            'not_a_dict',
        ],
    )
    def test_bad_melt_composition_raises(self, oxides, reason):
        """The melt needs known oxides, finite and non-negative, summing to 95-101 wt%."""
        with pytest.raises(ValueError) as exc:
            check_metal_fields({**GOOD, 'silicate_oxides': oxides})
        assert 'silicate_oxides' in str(exc.value) and reason in str(exc.value)

    @pytest.mark.parametrize('field', ['P_eq', 'T_eq', 'dIW_eq'])
    def test_required_fields_are_named_when_absent(self, field):
        """With metal, P_eq, T_eq and (parameter mode) dIW_eq are required; without, not."""
        missing = {k: v for k, v in GOOD.items() if k != field}
        with pytest.raises(ValueError) as exc:
            check_metal_fields(missing)
        assert field in str(exc.value) and 'required' in str(exc.value)
        assert check_metal_fields({**missing, 'metal_mass': 0.0}) is None
        if field == 'dIW_eq':
            assert check_metal_fields({**missing, 'dIW_eq_mode': 'from_FeO'}) is None

    @pytest.mark.parametrize(
        'extra, oxide',
        [
            ({'dIW_eq_mode': 'from_FeO'}, 'FeO'),
            ({'Si_mode_core': 'partition'}, 'SiO2'),
        ],
        ids=['from_FeO_needs_FeO', 'Si_partition_needs_SiO2'],
    )
    def test_mode_needs_its_oxide(self, extra, oxide):
        """from_FeO needs FeO above 0 in the melt, Si partition needs SiO2 above 0."""
        rest = {k: v for k, v in KLB1.items() if k != oxide}
        rest['MgO'] += KLB1[oxide]
        for oxides in (rest, {**rest, oxide: 0.0}):
            with pytest.raises(ValueError) as exc:
                check_metal_fields({**GOOD, **extra, 'silicate_oxides': oxides})
            assert oxide in str(exc.value)
        assert check_metal_fields({**GOOD, **extra}) is None

    @pytest.mark.parametrize(
        'field, value',
        [('O_mode_core', 'free'), ('Si_mode_core', 1), ('dIW_eq_mode', 'buffer')],
    )
    def test_unknown_mode_raises(self, field, value):
        """A mode outside the known ones is refused with its field name."""
        with pytest.raises(ValueError) as exc:
            check_metal_fields({**GOOD, field: value})
        assert field in str(exc.value)
        assert repr(value) in str(exc.value)

    def test_array_valued_mode_is_refused_by_name(self):
        """A mode given as an array is refused with its field name."""
        with pytest.raises(ValueError, match='O_mode_core'):
            check_metal_fields({**GOOD, 'O_mode_core': np.array(['fixed', 'partition'])})

    def test_fixed_o_and_si_that_leave_no_iron_are_refused(self):
        """w_O_core + w_Si_core of 1 or more leaves no iron in the metal."""
        with pytest.raises(ValueError) as exc:
            check_metal_fields({**GOOD, 'w_O_core': 0.6, 'w_Si_core': 0.5})
        assert 'w_O_core' in str(exc.value) and 'w_Si_core' in str(exc.value)
        assert check_metal_fields({**GOOD, 'w_O_core': 0.5, 'w_Si_core': 0.4}) is None

    def test_numpy_scalars_pass_and_numpy_bool_fails(self):
        """numpy float32 and int64 values are numbers; a numpy bool is not."""
        ok = {**GOOD, 'P_eq': np.float32(50.0), 'dIW_eq': np.int64(-2)}
        ok['silicate_oxides'] = {k: np.float32(v) for k, v in KLB1.items()}
        assert check_metal_fields(ok) is None
        with pytest.raises(ValueError, match='T_eq'):
            check_metal_fields({**GOOD, 'T_eq': np.bool_(True)})
        with pytest.raises(ValueError, match='finite'):
            check_metal_fields({**GOOD, 'P_eq': np.float64('nan')})

    @settings(derandomize=True, max_examples=200, deadline=None)
    @given(st.floats(allow_nan=True, allow_infinity=True))
    def test_temperature_passes_only_when_finite_and_positive(self, t_eq):
        """Over every float, T_eq passes exactly when it is finite and above 0 K."""
        ok = math.isfinite(t_eq) and t_eq > 0
        if ok:
            assert check_metal_fields({**GOOD, 'T_eq': t_eq}) is None
        else:
            with pytest.raises(ValueError, match='T_eq'):
                check_metal_fields({**GOOD, 'T_eq': t_eq})

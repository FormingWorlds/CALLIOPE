"""Tests for the metal-reservoir fields of the CALLIOPE input dict (`src/calliope/solve.py`).

Contract: each field that the metal reservoir reads is checked before the first solve;
a boolean (a TOML true), a string, NaN, inf, or a value outside its range raises
ValueError naming the field, and the defaults and the KLB-1 melt pass.
"""

from __future__ import annotations

import math

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
    'P_eq': [True, math.nan, math.inf, 0.0, -5.0],
    'T_eq': [False, math.nan, -math.inf, 0.0, -1.0],
    'dIW_eq': [True, 'IW-2', math.nan, math.inf],
    'silicate_nbo_t': [True, math.nan, math.inf, -0.1],
    'w_O_core': [True, math.nan, -0.01, 1.0],
    'w_Si_core': [math.inf, -1e-9, 1.5],
    'tol_w_core': [True, math.nan, 0.0, -1e-6],
}


class TestMetalFields:
    """Field checks of the metal reservoir."""

    def test_good_fields_and_klb1_pass(self):
        """The fields above, with the KLB-1 melt (99.29 wt%), pass the check."""
        assert check_metal_fields(dict(GOOD)) is None
        assert sum(KLB1.values()) == pytest.approx(99.29, abs=1e-9)

    @pytest.mark.parametrize(
        'field, value', [(f, v) for f, vals in BAD.items() for v in vals], ids=repr
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

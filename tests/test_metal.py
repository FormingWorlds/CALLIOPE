"""Tests for `src/calliope/metal.py`: the metal-silicate partition coefficient file.

Contract clauses: a file with schema_version 1 and one section per element (H, C, N,
S, O, Si) loads with every field typed and checked; each malformed field (a TOML
boolean or string where a number belongs, NaN, inf, an unknown form, quantity or
basis, a mole basis without its silicate basis, a reversed or non-finite range, a
pressure unit other than GPa, an interaction term for an unknown element, a
covariance that is not square, not symmetric, not positive semidefinite or of the
wrong size) raises PartitionTableError naming the element and field; the packaged
placeholder file loads and says that it is a placeholder.
"""

from __future__ import annotations

import math

import pytest

from calliope.metal import (
    ELEMENTS,
    PartitionTableError,
    load_partition_table,
    log10_D,
    out_of_range,
)

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Tagawa et al. 2021 Eq. 2 for H, the shape every entry of the file has
_H = {
    'form': '"log10_linear"',
    'quantity': '"K_D"',
    'basis': '"mole_xprime_noHC"',
    'silicate_basis': '"single_cation"',
    'a': '2.42',
    'b': '-2892.0',
    'c': '-32.0',
    'd': '-0.25',
    'e': '0.0',
    'dIW_factor': '1.0',
    'eps': '{}',
    'covariance': '[]',
    'P_range_GPa': '[30.0, 60.0]',
    'T_range_K': '[3100.0, 4600.0]',
    'dIW_range': '[-1.33, -0.88]',
    'P_unit': '"GPa"',
    'source': '"Tagawa et al. 2021, Nat. Commun. 12, 2588, Eq. 2"',
    'doi': '"10.1038/s41467-021-22035-0"',
}


def _text(entries=None, header='schema_version = 1'):
    """TOML text with a meta table and one section per element of `entries`."""
    entries = {'H': _H} if entries is None else entries
    out = [header, '[meta]', 'name = "test"', 'placeholder = true']
    for el, fields in entries.items():
        out.append(f'[{el}]')
        out += [f'{k} = {v}' for k, v in fields.items()]
    return '\n'.join(out) + '\n'


def _write(tmp_path, text):
    path = tmp_path / 'partition.toml'
    path.write_text(text)
    return path


def _with(**changes):
    """The H entry with fields replaced (value None removes the field)."""
    fields = dict(_H)
    for k, v in changes.items():
        if v is None:
            fields.pop(k)
        else:
            fields[k] = v
    return {'H': fields}


class TestValidFile:
    """A well-formed file loads with typed fields."""

    def test_tagawa_entry_loads_with_numbers_and_ranges(self, tmp_path):
        """Every coefficient is a float, the ranges are pairs and the meta table is kept."""
        table = load_partition_table(_write(tmp_path, _text()))
        h = table['H']
        assert h['a'] == pytest.approx(2.42)
        assert h['c'] == pytest.approx(-32.0)
        assert h['P_range_GPa'] == (pytest.approx(30.0), pytest.approx(60.0))
        assert h['eps'] == {}
        assert table['meta']['placeholder'] is True
        assert set(table) == {'schema_version', 'meta', 'H'}

    def test_interaction_terms_and_matching_covariance_load(self, tmp_path):
        """Two eps terms need a 7 x 7 covariance in the order a, b, c, d, e, eps by name."""
        cov = [[1.0 if i == j else 0.0 for j in range(7)] for i in range(7)]
        entry = _with(eps='{ S = 9.6, O = -19.5 }', covariance=str(cov))
        table = load_partition_table(_write(tmp_path, _text(entry)))
        assert table['H']['eps'] == {'O': pytest.approx(-19.5), 'S': pytest.approx(9.6)}
        assert len(table['H']['covariance']) == 7

    def test_empty_range_means_an_unknown_range(self, tmp_path):
        """An empty list is a range the source does not give: it loads as None."""
        table = load_partition_table(
            _write(tmp_path, _text(_with(T_range_K='[]', dIW_range='[]')))
        )
        assert table['H']['T_range_K'] is None
        assert table['H']['dIW_range'] is None
        assert table['H']['P_range_GPa'] == (pytest.approx(30.0), pytest.approx(60.0))

    def test_packaged_placeholder_covers_every_element(self):
        """The default file is the placeholder set, flagged as such, with all 6 elements."""
        table = load_partition_table()
        assert table['meta']['placeholder'] is True
        assert set(ELEMENTS) <= set(table)
        assert ELEMENTS == ('H', 'C', 'N', 'S', 'O', 'Si')
        for el in ELEMENTS:
            assert table[el]['source'] and math.isfinite(table[el]['a'])


class TestMalformedFile:
    """Each malformed field raises PartitionTableError with the element and the field."""

    @pytest.mark.parametrize(
        'entry, field',
        [
            (_with(a='true'), 'a'),
            (_with(b='"-2892"'), 'b'),
            (_with(c='nan'), 'c'),
            (_with(d='inf'), 'd'),
            (_with(e=None), 'e'),
            (_with(dIW_factor='-inf'), 'dIW_factor'),
            (_with(form='"power_law"'), 'form'),
            (_with(quantity='"Kd"'), 'quantity'),
            (_with(basis='"molar"'), 'basis'),
            (_with(silicate_basis=None), 'silicate_basis'),
            (_with(silicate_basis='"oxide"'), 'silicate_basis'),
            (_with(P_range_GPa='[60.0, 30.0]'), 'P_range_GPa'),
            (_with(T_range_K='[3100.0, nan]'), 'T_range_K'),
            (_with(dIW_range='[-1.0]'), 'dIW_range'),
            (_with(P_unit='"Pa"'), 'P_unit'),
            (_with(eps='{ Xe = 1.0 }'), 'eps'),
            (_with(eps='{ S = true }'), 'eps'),
            (_with(covariance='[[1.0, 0.0], [0.0, 1.0]]'), 'covariance'),
            (_with(covariance='[[1.0, 2.0, 0.0, 0.0, 0.0]]'), 'covariance'),
            (_with(source='""'), 'source'),
        ],
        ids=[
            'bool_coefficient',
            'string_coefficient',
            'nan_coefficient',
            'inf_coefficient',
            'missing_coefficient',
            'infinite_diw_factor',
            'unknown_form',
            'unknown_quantity',
            'unknown_basis',
            'mole_basis_without_silicate_basis',
            'unknown_silicate_basis',
            'reversed_pressure_range',
            'nan_temperature_range',
            'short_diw_range',
            'pressure_in_pa',
            'eps_unknown_element',
            'eps_bool',
            'covariance_wrong_size',
            'covariance_not_square',
            'empty_source',
        ],
    )
    def test_malformed_field_raises_with_its_name(self, tmp_path, entry, field):
        """The error names the element and the field, and is a ValueError."""
        with pytest.raises(PartitionTableError) as exc:
            load_partition_table(_write(tmp_path, _text(entry)))
        assert "'H'" in str(exc.value) and f"'{field}'" in str(exc.value)
        assert isinstance(exc.value, ValueError)

    @pytest.mark.parametrize(
        'cov, reason',
        [
            (
                [
                    [1.0, 0.5, 0, 0, 0],
                    [0.4, 1.0, 0, 0, 0],
                    [0, 0, 1.0, 0, 0],
                    [0, 0, 0, 1.0, 0],
                    [0, 0, 0, 0, 1.0],
                ],
                'symmetric',
            ),
            (
                [
                    [1.0, 2.0, 0, 0, 0],
                    [2.0, 1.0, 0, 0, 0],
                    [0, 0, 1.0, 0, 0],
                    [0, 0, 0, 1.0, 0],
                    [0, 0, 0, 0, 1.0],
                ],
                'semidefinite',
            ),
        ],
        ids=['asymmetric', 'negative_eigenvalue'],
    )
    def test_covariance_must_be_a_covariance(self, tmp_path, cov, reason):
        """An asymmetric matrix or one with a negative eigenvalue (here -1) is refused."""
        with pytest.raises(PartitionTableError) as exc:
            load_partition_table(_write(tmp_path, _text(_with(covariance=str(cov)))))
        assert reason in str(exc.value)
        assert "'covariance'" in str(exc.value)

    @pytest.mark.parametrize(
        'header, message',
        [
            ('', 'schema_version'),
            ('schema_version = 2', 'schema_version'),
            ('schema_version = true', 'schema_version'),
        ],
        ids=['missing_version', 'future_version', 'bool_version'],
    )
    def test_schema_version_must_be_one(self, tmp_path, header, message):
        """Only schema version 1 (an integer, not a boolean) is read."""
        with pytest.raises(PartitionTableError) as exc:
            load_partition_table(_write(tmp_path, _text(header=header)))
        assert message in str(exc.value)
        assert load_partition_table(_write(tmp_path, _text()))['schema_version'] == 1

    def test_unknown_element_section_and_missing_meta_raise(self, tmp_path):
        """A section that is no element is refused, as is a file without its meta table."""
        with pytest.raises(PartitionTableError, match="'Xe'"):
            load_partition_table(_write(tmp_path, _text({'Xe': _H})))
        with pytest.raises(PartitionTableError, match='meta'):
            load_partition_table(_write(tmp_path, 'schema_version = 1\n[H]\na = 1.0\n'))

    def test_missing_file_raises_file_not_found(self, tmp_path):
        """A path that does not exist is an OSError, not a malformed table."""
        with pytest.raises(FileNotFoundError):
            load_partition_table(tmp_path / 'none.toml')
        assert not (tmp_path / 'none.toml').exists()


@settings(derandomize=True, max_examples=100, deadline=None)
@given(st.floats(allow_nan=True, allow_infinity=True))
def test_any_float_coefficient_loads_only_when_finite(tmp_path_factory, value):
    """Over every float, a coefficient loads exactly when it is finite."""
    path = tmp_path_factory.mktemp('h') / 'partition.toml'
    path.write_text(_text(_with(a=repr(value))))
    if math.isfinite(value):
        assert load_partition_table(path)['H']['a'] == value
    else:
        with pytest.raises(PartitionTableError, match="'a'"):
            load_partition_table(path)


class TestLaw:
    """The log10 D law of an entry and its range flag."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    @pytest.mark.parametrize(
        'p_gpa, t_k, expected',
        [(20.0, 3000.0, 1.2427), (135.0, 4000.0, 0.6170), (300.0, 5000.0, -0.0784)],
    )
    def test_tagawa_eq2_values(self, p_gpa, t_k, expected):
        """Tagawa et al. 2021 Eq. 2: log10 K_D = 2.42 - 2892/T - 32.0 P/T (dIW 0), 4 decimals."""
        h = load_partition_table()['H']
        val = log10_D(h, p_gpa, t_k, 0.0)
        assert val == pytest.approx(expected, abs=5e-5)
        # the dIW term is -dIW/4: one log unit more reduced raises log10 D by 0.25
        assert log10_D(h, p_gpa, t_k, -1.0) - val == pytest.approx(0.25, abs=1e-12)
        # Tagawa Eq. 10 (0.692 - 4590/T + 141 P/T) differs by more than 0.3 at each point
        assert abs(val - (0.692 - 4590 / t_k + 141 * p_gpa / t_k)) > 0.3

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_fischer_eq1_single_stage_value(self):
        """Fischer et al. 2020 p3: Eq. 1 gives D_C = 0.550 at 54 GPa, 4100 K, NBO/T 2.6, IW-2.2, X_O 0.06, X_S 0.03."""
        c = load_partition_table()['C']
        d_c = 10 ** log10_D(c, 54.0, 4100.0, -2.2, x_metal={'O': 0.06, 'S': 0.03}, nbo_t=2.6)
        assert d_c == pytest.approx(0.550, abs=5e-4)
        # Suer SA prints the dIW coefficient as 0.328: that law lands far outside
        wrong = 10 ** (log10_D(c, 54.0, 4100.0, -2.2, {'O': 0.06, 'S': 0.03}, 2.6) + 0.09 * 2.2)
        assert abs(wrong - d_c) > 0.1
        # without the O and S terms D_C is 0.447 (derivation section B)
        assert 10 ** log10_D(c, 54.0, 4100.0, -2.2, nbo_t=0.0) == pytest.approx(0.447, abs=2e-3)

    def test_out_of_range_flags_each_axis_and_unknown_ranges(self):
        """A point outside any range, or on an axis whose range is unknown, is flagged."""
        table = load_partition_table()
        h, o = table['H'], table['O']
        assert not out_of_range(h, 45.0, 4000.0, -1.0)
        assert out_of_range(h, 135.0, 4000.0, -1.0)
        assert out_of_range(h, 45.0, 2500.0, -1.0)
        assert out_of_range(h, 45.0, 4000.0, -2.3)
        assert out_of_range(o, 50.0, 4000.0, -2.0)

    def test_non_finite_condition_raises(self):
        """A NaN pressure, a non-positive temperature or a metal fraction of 1 is refused."""
        c = load_partition_table()['C']
        for args in ((math.nan, 4000.0, -2.0), (50.0, 0.0, -2.0), (50.0, 4000.0, math.inf)):
            with pytest.raises(ValueError):
                log10_D(c, *args)
        with pytest.raises(ValueError, match='X_O'):
            log10_D(c, 50.0, 4000.0, -2.0, x_metal={'O': 1.0})

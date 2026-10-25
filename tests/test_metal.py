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
    core_light_elements,
    load_partition_table,
    log10_D,
    mass_D,
    metal_masses,
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


class TestBasis:
    """Molar D on the H- and C-free metal basis converted to core mass fractions."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    @pytest.mark.parametrize(
        'd_h, d_c, w_o, w_si, w_s, core_h, core_c',
        [
            (92, 67, 0.00, 3.43, 2.00, 0.68, 0.75),
            (56, 96, 2.52, 4.42, 2.00, 0.45, 1.15),
            (44, 107, 1.49, 7.78, 2.00, 0.35, 1.29),
        ],
        ids=['wade_wood_40gpa', 'siebert_50gpa', 'fischer_54gpa'],
    )
    def test_tsutsumi_single_stage_core(self, d_h, d_c, w_o, w_si, w_s, core_h, core_c):
        """Tsutsumi et al. 2025 SI Table 1: BSE 710 ppm H2O, 120 ppm C, KLB-1 melt (Takahashi 1986)."""
        w_melt = {'H': 710e-6 * 2 * 1.008 / 18.01528, 'C': 120e-6}
        fixed = {'O': w_o / 100, 'Si': w_si / 100, 'S': w_s / 100}
        w = core_light_elements({'H': d_h, 'C': d_c}, w_melt, fixed)
        assert round(100 * w['H'], 2) == pytest.approx(core_h, abs=1e-9)
        assert round(100 * w['C'], 2) == pytest.approx(core_c, abs=1e-9)
        assert 0 < w['H'] < 0.01 and 0 < w['C'] < 0.02

    def test_plain_mass_d_misses_the_wade_wood_case(self):
        """Without the basis conversion, D times the melt H gives 0.73 wt%, not 0.68 (40 GPa case)."""
        w_melt = {'H': 710e-6 * 2 * 1.008 / 18.01528, 'C': 120e-6}
        plain = 92 * w_melt['H'] * 100
        w = core_light_elements({'H': 92, 'C': 67}, w_melt, {'Si': 0.0343, 'S': 0.02})
        assert round(plain, 2) != pytest.approx(0.68, abs=1e-9)
        assert abs(plain / (100 * w['H']) - 1) > 0.06

    def test_lighter_metal_holds_more_h_and_c_at_fixed_mole_ratio(self):
        """O and Si are lighter than Fe: a kg of metal holds more H- and C-free moles, so more H."""
        w_melt = {'H': 7.94e-5, 'C': 1.2e-4}
        clean = core_light_elements({'H': 50, 'C': 50}, w_melt, {})
        loaded = core_light_elements({'H': 50, 'C': 50}, w_melt, {'O': 0.05, 'Si': 0.05})
        assert loaded['H'] > clean['H'] and loaded['C'] > clean['C']
        # the ratio is the moles per kg of the H- and C-free metal, about 1.17 here
        assert loaded['H'] / clean['H'] == pytest.approx(loaded['C'] / clean['C'], rel=1e-3)
        with pytest.raises(ValueError, match='no iron'):
            core_light_elements({'H': 50, 'C': 50}, w_melt, {'O': 0.6, 'Si': 0.5})


class TestMetalReservoir:
    """The metal term of the element balance: core_i = D_i M_metal w_i(melt)."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    @pytest.mark.parametrize('d', [0.1, 1.0, 10.0, 30.0, 100.0])
    def test_gas_free_core_fraction_is_dx_over_one_plus_dx(self, d):
        """Melt and metal alone: the core holds f = Dx/(1 + Dx) of an element, x = 0.325/0.675."""
        m_planet, x = 6e24, 0.325 / 0.675
        ddict = {
            'M_mantle': 0.675 * m_planet,
            'Phi_global': 1.0,
            'metal_mass': 0.325 * m_planet,
        }
        total = 1e20
        melt = total / (1 + d * x)  # the gas-free balance melt + core = total
        core = metal_masses({'H': melt}, ddict, {'H': d})['H']
        f = core / (melt + core)
        assert f == pytest.approx(d * x / (1 + d * x), rel=1e-12)
        assert melt + core == pytest.approx(total, rel=1e-12)
        if (
            d >= 1.0
        ):  # f = Dx (no melt back-reaction) would exceed the right value by a factor 1 + Dx
            assert abs(f - d * x) > 0.3 * f

    def test_no_metal_gives_no_core_term_and_no_melt_raises(self):
        """metal_mass 0 gives zero core masses; metal with Phi_global 0 has no melt to equilibrate."""
        ddict = {'M_mantle': 4e24, 'Phi_global': 1.0, 'metal_mass': 0.0}
        assert metal_masses({'H': 1e20, 'C': 1e19}, ddict, {'H': 50.0, 'C': 5.0}) == {
            'H': 0.0,
            'C': 0.0,
        }
        ddict.update(metal_mass=2e24, Phi_global=0.0)
        with pytest.raises(ValueError, match='melt'):
            metal_masses({'H': 1e20}, ddict, {'H': 50.0})
        ddict.update(Phi_global=0.5)
        with pytest.raises(ValueError, match='D'):
            metal_masses({'H': 1e20}, ddict, {'H': -1.0})


class TestMassD:
    """The mass-basis D that the element balance uses."""

    @staticmethod
    def _const(log10_d, basis):
        entry = dict(load_partition_table()['H'])
        entry.update(a=log10_d, b=0.0, c=0.0, d=0.0, e=0.0, eps={}, basis=basis)
        return entry

    @pytest.mark.physics_invariant
    def test_xprime_basis_reproduces_the_closed_form_core(self):
        """On the self-consistent metal, D_w times the melt H equals the closed-form core H."""
        w_melt = {'H': 710e-6 * 2 * 1.008 / 18.01528, 'C': 120e-6}
        fixed = {'Si': 0.0343, 'S': 0.02}
        core = core_light_elements({'H': 92, 'C': 67}, w_melt, fixed)
        w_metal = {**fixed, **core}
        d_w = mass_D(
            self._const(math.log10(92), 'mole_xprime_noHC'), 45.0, 3750.0, -2.3, w_metal, w_melt
        )
        assert d_w * w_melt['H'] == pytest.approx(core['H'], rel=1e-12)
        # the conversion factor is not 1: plain D would put 7 % more H in the core
        assert abs(92 / d_w - 1) > 0.06

    def test_mass_basis_is_the_law_itself(self):
        """A mass-basis entry gives 10**log10_D whatever the metal and melt."""
        entry = self._const(1.5, 'mass')
        d_w = mass_D(entry, 45.0, 3750.0, -2.3, {'O': 0.05}, {'H': 1e-4})
        assert d_w == pytest.approx(10**1.5, rel=1e-15)
        assert mass_D(entry, 45.0, 3750.0, -2.3, {}, {}) == pytest.approx(d_w, rel=1e-15)

    def test_metal_with_no_iron_raises(self):
        """Light elements that sum to 1 leave no iron for the metal basis."""
        entry = self._const(1.0, 'mole_xprime_noHC')
        with pytest.raises(ValueError, match='no iron'):
            mass_D(entry, 45.0, 3750.0, -2.3, {'O': 0.6, 'Si': 0.4}, {'H': 1e-4})
        assert mass_D(entry, 45.0, 3750.0, -2.3, {'O': 0.1}, {'H': 1e-4}) > 0

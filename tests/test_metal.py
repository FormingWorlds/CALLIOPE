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
    exchange_metal,
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
            (_with(covariance='[1.0, 2.0]'), 'covariance'),
            (_with(covariance='[[1.0, 0.0], [0.0]]'), 'covariance'),
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
            'covariance_flat_list',
            'covariance_ragged',
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
        'el, changes, field',
        [
            ('H', {'form': '"exchange_KD"'}, 'form'),
            ('O', {'form': '"log10_linear"', 'quantity': '"D"', 'basis': '"mass"'}, 'form'),
            ('O', {'quantity': '"D"'}, 'quantity'),
            ('Si', {'basis': '"mass"'}, 'basis'),
            ('Si', {'basis': '"mole_xprime_noHC"'}, 'basis'),
            ('O', {'d': '0.5'}, 'd'),
            ('O', {'e': '-0.1'}, 'e'),
            ('Si', {'eps': '{ O = -5.0 }'}, 'eps'),
            ('O', {'dIW_factor': '2.0'}, 'dIW_factor'),
        ],
        ids=[
            'H_exchange_form',
            'O_linear_form',
            'exchange_quantity_D',
            'exchange_basis_mass',
            'exchange_basis_xprime',
            'exchange_d_term',
            'exchange_e_term',
            'exchange_eps_term',
            'exchange_diw_factor',
        ],
    )
    def test_unread_field_values_are_refused(self, tmp_path, el, changes, field):
        """A value the evaluator of the element does not read is refused, naming the field."""
        exchange = {
            **_H,
            'form': '"exchange_KD"',
            'quantity': '"K_D"',
            'basis': '"mole"',
            'a': '0.986',
            'b': '-3275.0',
            'c': '0.0',
            'd': '0.0',
        }
        entry = dict(_H if el == 'H' else exchange)
        entry.update(changes)
        with pytest.raises(PartitionTableError) as exc:
            load_partition_table(_write(tmp_path, _text({el: entry})))
        assert f"'{el}'" in str(exc.value) and f"'{field}'" in str(exc.value)
        # the unchanged entry of the same element loads
        base = _H if el == 'H' else exchange
        assert load_partition_table(_write(tmp_path, _text({el: base})))[el]['form']

    def test_evaluators_refuse_an_entry_of_the_other_form(self):
        """A table built in code bypasses the loader: mass_D and exchange_metal check the form."""
        table = load_partition_table()
        exchange_h = dict(table['H'], form='exchange_KD')
        with pytest.raises(PartitionTableError, match="'form'"):
            mass_D(exchange_h, 50.0, 3500.0, -2.0, {}, {})
        linear_o = dict(table, O=dict(table['O'], form='log10_linear'))
        with pytest.raises(PartitionTableError, match="'form'"):
            exchange_metal(linear_o, 50.0, 3500.0, -2.0)
        assert mass_D(table['H'], 50.0, 3500.0, -2.0, {}, {}) > 0
        assert 0 < exchange_metal(table, 50.0, 3500.0, -2.0, elements=('Si',))['Si'] < 1

    def test_element_entry_must_be_a_table(self, tmp_path):
        """An element given as a value instead of a table is refused with its name."""
        text = 'schema_version = 1\nH = 3.0\n[meta]\nname = "x"\nplaceholder = true\n'
        with pytest.raises(PartitionTableError) as exc:
            load_partition_table(_write(tmp_path, text))
        assert "'H'" in str(exc.value) and 'table' in str(exc.value)
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

    def test_elements_outside_their_roles_are_refused(self):
        """D_x is for H and C, the fixed fractions for O, Si, S and N; anything else raises."""
        w_melt = {'H': 7.94e-5, 'C': 1.2e-4}
        for d_x, fixed in (
            ({'H': 50, 'N': 5}, {}),
            ({'H': 50}, {'C': 0.01}),
            ({}, {'H': 0.01}),
        ):
            with pytest.raises(ValueError, match='d_x takes H and C'):
                core_light_elements(d_x, w_melt, fixed)
        # the edge: every allowed role filled is accepted and leaves iron
        w = core_light_elements(
            {'H': 50, 'C': 50}, w_melt, {'O': 0.01, 'Si': 0.01, 'S': 0.01, 'N': 0.01}
        )
        assert 0 < w['H'] < 0.01 and 0 < w['C'] < 0.01

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
        # f = Dx (no melt back-reaction) would exceed the right value by a factor 1 + Dx
        if d >= 1.0:
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


class TestExchange:
    """O and Si in the metal from the exchange K_D of FeO = Fe + O and SiO2 + 2 Fe = Si + 2 FeO."""

    @pytest.mark.physics_invariant
    @pytest.mark.reference_pinned
    def test_ideal_exchange_mole_fractions(self):
        """x_O = K_D 10^(dIW/2) and x_Si = K_D x_SiO2 10^(-dIW), Tagawa SI Table 3 (refs. 2, 57)."""
        table = load_partition_table()
        t_k, d_iw = 4000.0, -2.0
        x = exchange_metal(table, 50.0, t_k, d_iw)
        k_o = 10 ** (0.986 - 3275.0 / t_k)
        k_si = 10 ** (0.549 - 12324.0 / t_k)
        x_sio2 = (44.48 / 60.083) / 1.96965  # KLB-1 single-cation basis, mol per 100 g
        assert x['O'] == pytest.approx(k_o * 10 ** (d_iw / 2), rel=1e-12)
        assert x['Si'] == pytest.approx(k_si * x_sio2 * 10 ** (-d_iw), rel=1e-3)
        # read as the K of the half reaction, K_D is squared: x_Si about 3e-4, not 0.11
        assert x['Si'] > 100 * k_si**2 * x_sio2 * 10 ** (-d_iw)

    def test_more_oxidised_metal_takes_more_o_and_less_si(self):
        """Raising dIW by 1 multiplies x_O by 10^0.5 and x_Si by 10^-1; too much raises."""
        table = load_partition_table()
        lo = exchange_metal(table, 50.0, 4000.0, -2.0)
        hi = exchange_metal(table, 50.0, 4000.0, -1.0)
        assert hi['O'] / lo['O'] == pytest.approx(10**0.5, rel=1e-12)
        assert hi['Si'] / lo['Si'] == pytest.approx(0.1, rel=1e-12)
        with pytest.raises(ValueError, match='no iron'):
            exchange_metal(table, 50.0, 4000.0, -6.0)

    def test_one_partitioned_element_ignores_the_other_law(self):
        """With only O partitioned the Si law is not evaluated, so a reduced dIW passes."""
        table = load_partition_table()
        both = exchange_metal(table, 50.0, 4000.0, -2.0)
        only_o = exchange_metal(table, 50.0, 4000.0, -2.0, elements=('O',))
        assert set(only_o) == {'O'}
        assert only_o['O'] == pytest.approx(both['O'], rel=1e-15)
        # at dIW -8 the Si law alone exceeds 1, the O law gives a small x_O
        assert exchange_metal(table, 50.0, 4000.0, -8.0, elements=('O',))['O'] < 1e-3
        with pytest.raises(ValueError, match='no iron'):
            exchange_metal(table, 50.0, 4000.0, -8.0, elements=('Si',))


# Tagawa et al. 2021 SI Table 1: P [GPa], T [K], dIW, measured log10 K_D and its 1 sigma
TAGAWA_RUNS = [
    (46.0, 3920.0, -0.88, 1.32, 0.06),
    (48.0, 3450.0, -0.99, 1.10, 0.06),
    (57.0, 3860.0, -1.08, 1.22, 0.05),
    (60.0, 4560.0, -1.17, 1.37, 0.05),
    (47.0, 4230.0, -1.27, 1.35, 0.06),
    (30.0, 3080.0, -1.33, 1.18, 0.06),
]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_tagawa_runs_lie_within_one_sigma_of_eq2(validation_plot):
    """Tagawa et al. 2021 SI Table 1: each of the 6 runs within its 1 sigma of the H law."""
    h = load_partition_table()['H']
    fit = [log10_D(h, p, t, 0.0) for p, t, _, _, _ in TAGAWA_RUNS]
    ratio = [(m - f) / s for (_, _, _, m, s), f in zip(TAGAWA_RUNS, fit)]
    assert max(abs(r) for r in ratio) <= 1.0
    # Tagawa Eq. 10 (fit with the multi-anvil data) misses every run by more than 2 sigma
    eq10 = [0.692 - 4590 / t + 141 * p / t for p, t, _, _, _ in TAGAWA_RUNS]
    assert min(abs(m - e) / s for (_, _, _, m, s), e in zip(TAGAWA_RUNS, eq10)) > 2.0
    # the ranges of the H entry are those of the 6 runs, so every run is in range
    assert not any(out_of_range(h, p, t, d) for p, t, d, _, _ in TAGAWA_RUNS)
    assert out_of_range(h, 30.0, 3079.0, -1.33) and out_of_range(h, 60.0, 4561.0, -1.17)

    def figure():
        import matplotlib.pyplot as plt

        meas = [r[3] for r in TAGAWA_RUNS]
        fig, ax = plt.subplots()
        ax.errorbar(fit, meas, yerr=[r[4] for r in TAGAWA_RUNS], fmt='o', label='runs #1-#6')
        ax.plot([1.0, 1.5], [1.0, 1.5], ls='--', color=plt.rcParams['text.color'], label='1:1')
        ax.set(xlabel=r'Eq. 2 $\log_{10} K_D$', ylabel=r'measured $\log_{10} K_D$')
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'metal_tagawa2021_h_runs')


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_gas_free_core_fraction_curve(validation_plot):
    """The core term on a D grid lies on Dx/(1 + Dx); f = Dx lies more than 10 % above it."""
    x = 0.325 / 0.675
    ddict = {'M_mantle': 0.675, 'Phi_global': 1.0, 'metal_mass': 0.325}
    ds = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
    f = []
    for d in ds:
        melt = 1.0 / (1 + d * x)
        core = metal_masses({'H': melt}, ddict, {'H': d})['H']
        f.append(core / (melt + core))
    assert f == pytest.approx([d * x / (1 + d * x) for d in ds], rel=1e-12)
    # f = Dx overshoots by (Dx)^2/(1 + Dx): more than 10 % of f from D = 0.3 on
    assert all(d * x - fi > 0.1 * fi for fi, d in zip(f[1:], ds[1:]))

    def figure():
        import matplotlib.pyplot as plt

        grid = [10 ** (k / 20) for k in range(-25, 46)]
        fig, ax = plt.subplots()
        ax.plot(grid, [d * x / (1 + d * x) for d in grid], label=r'$Dx/(1+Dx)$')
        ax.plot(
            grid, [min(d * x, 1.0) for d in grid], ls='--', label=r'$Dx$ (no back-reaction)'
        )
        ax.plot(ds, f, 'o', label='metal_masses')
        ax.set(xscale='log', xlabel='D (metal / melt, mass)', ylabel='core fraction')
        ax.minorticks_off()
        ax.legend()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'metal_gas_free_core_fraction')


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_tsutsumi_core_values_plot(validation_plot):
    """Core H and C of the 3 Tsutsumi et al. 2025 single-stage cases against SI Table 1."""
    cases = [
        ((92, 67), (0.00, 3.43), (0.68, 0.75)),
        ((56, 96), (2.52, 4.42), (0.45, 1.15)),
        ((44, 107), (1.49, 7.78), (0.35, 1.29)),
    ]
    w_melt = {'H': 710e-6 * 2 * 1.008 / 18.01528, 'C': 120e-6}
    got = []
    for (d_h, d_c), (o, si), _ in cases:
        w = core_light_elements(
            {'H': d_h, 'C': d_c}, w_melt, {'O': o / 100, 'Si': si / 100, 'S': 0.02}
        )
        got.append((100 * w['H'], 100 * w['C']))
    for (_, _, paper), (h, c) in zip(cases, got):
        assert abs(h - paper[0]) <= 0.005 + 1e-9
        assert abs(c - paper[1]) <= 0.005 + 1e-9
    # D times the melt H without the basis conversion misses the 40 GPa case (0.73 wt%)
    assert abs(100 * 92 * w_melt['H'] - cases[0][2][0]) > 0.005

    def figure():
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        labels = ['40 GPa', '50 GPa', '54 GPa']
        for k, (name, marker) in enumerate((('H', 'o'), ('C', 's'))):
            ax.plot(
                labels,
                [c[2][k] for c in cases],
                marker,
                mfc='none',
                ms=10,
                label=f'{name}, SI Table 1',
            )
            ax.plot(labels, [g[k] for g in got], marker, ms=4, label=f'{name}, CALLIOPE')
        ax.set(ylabel='core mass fraction [wt%]')
        ax.legend(fontsize='small')
        fig.tight_layout()
        return fig

    validation_plot(figure, 'metal_tsutsumi2025_core_h_c')

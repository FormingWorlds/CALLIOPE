"""Metal-silicate partitioning of volatiles into core-forming metal.

The partition laws live in a TOML file with ``schema_version = 1``, a ``[meta]`` table
(``name``, ``placeholder``) and one table per element of :data:`ELEMENTS`. Each entry
gives ``log10 D = a + b/T + c P/T + d dIW_factor dIW + e NBO/T + sum_j eps_j
log10(1 - X_j)`` with T in K, P in GPa and X_j the mole fraction of element j in the
metal, together with the quantity the source fits (``D`` or ``K_D``), the metal
``basis`` of D (``mass``, ``mole`` or ``mole_xprime_noHC``, the mole fraction over the
metal without H and C), the ``silicate_basis`` of a mole basis (``single_cation``),
the ranges of the source data (an empty T or dIW range means unknown), an optional covariance of the coefficients in the
order a, b, c, d, e and the eps terms by element name, and the source and DOI. The
packaged file ``data/partition_coefficients_placeholder.toml`` holds placeholder laws.
"""

from __future__ import annotations

import math
import tomllib
from importlib import resources
from pathlib import Path

import numpy as np

from .constants import molar_mass, molar_mass_rock

ELEMENTS = ('H', 'C', 'N', 'S', 'O', 'Si')
EXCHANGE_ELEMENTS = ('O', 'Si')
SCHEMA_VERSION = 1
FORMS = ('log10_linear', 'exchange_KD')
QUANTITIES = ('D', 'K_D')
BASES = ('mass', 'mole', 'mole_xprime_noHC')
SILICATE_BASES = ('single_cation',)
COEFFICIENTS = ('a', 'b', 'c', 'd', 'e', 'dIW_factor')
RANGES = ('P_range_GPa', 'T_range_K', 'dIW_range')
KEYS = frozenset(('form', 'quantity', 'basis', 'silicate_basis', 'P_unit', 'eps', 'covariance'))
KEYS = KEYS | set(COEFFICIENTS) | set(RANGES) | {'source', 'doi'}
PLACEHOLDER_FILE = 'partition_coefficients_placeholder.toml'


class PartitionTableError(ValueError):
    """A partition coefficient file that breaks the schema."""


def _number(value, el, field):
    """Return a finite float; a boolean, string or non-finite value raises."""
    if (
        isinstance(value, bool)
        or not isinstance(value, int | float)
        or not math.isfinite(value)
    ):
        raise PartitionTableError(f'{el!r}: {field!r} must be a finite number, got {value!r}')
    return float(value)


def _choice(entry, el, field, options):
    value = entry.get(field)
    if value not in options:
        raise PartitionTableError(f'{el!r}: {field!r} must be one of {options}, got {value!r}')
    return value


def _entry(el, raw):
    """Check and type one element entry."""
    if not isinstance(raw, dict):
        raise PartitionTableError(f'{el!r}: the entry must be a table')
    unknown = sorted(set(raw) - KEYS)
    if unknown:
        raise PartitionTableError(f'{el!r}: unknown keys {unknown}')
    out = {
        'element': el,
        'form': _choice(raw, el, 'form', FORMS),
        'quantity': _choice(raw, el, 'quantity', QUANTITIES),
        'basis': _choice(raw, el, 'basis', BASES),
        'P_unit': _choice(raw, el, 'P_unit', ('GPa',)),
    }
    out['silicate_basis'] = (
        _choice(raw, el, 'silicate_basis', SILICATE_BASES) if out['basis'] != 'mass' else None
    )
    for field in COEFFICIENTS:
        out[field] = _number(raw.get(field), el, field)
    for field in RANGES:
        pair = raw.get(field)
        if pair == [] and field != 'P_range_GPa':
            out[field] = None  # not given by the source: every use is out of range
            continue
        if not isinstance(pair, list) or len(pair) != 2:
            raise PartitionTableError(f'{el!r}: {field!r} must be [low, high], got {pair!r}')
        lo, hi = (_number(v, el, field) for v in pair)
        if lo > hi:
            raise PartitionTableError(f'{el!r}: {field!r} low {lo} above high {hi}')
        if lo < 0 and field != 'dIW_range':
            raise PartitionTableError(f'{el!r}: {field!r} must not be negative, got {lo}')
        out[field] = (lo, hi)
    eps = raw.get('eps', {})
    if not isinstance(eps, dict) or any(k not in ELEMENTS for k in eps):
        raise PartitionTableError(f'{el!r}: {"eps"!r} keys must be elements of {ELEMENTS}')
    out['eps'] = {k: _number(v, el, 'eps') for k, v in sorted(eps.items())}
    out['covariance'] = _covariance(raw.get('covariance', []), 5 + len(eps), el)
    for field in ('source', 'doi'):
        if not isinstance(raw.get(field), str) or not raw[field].strip():
            raise PartitionTableError(f'{el!r}: {field!r} must be a non-empty string')
        out[field] = raw[field]
    _check_read_values(el, out)
    return out


def _need_form(entry, form):
    """Refuse an entry whose form is not the one the evaluator implements."""
    if entry['form'] != form:
        raise PartitionTableError(f'{"form"!r} must be {form!r} here, got {entry["form"]!r}')


def _check_read_values(el, entry):
    """Refuse a value that the evaluator of the element does not read.

    H, C, N and S use :func:`mass_D` (``log10_linear``, every basis); O and Si use
    :func:`exchange_metal` (``exchange_KD``, ``K_D`` on the mole basis, with only the
    a, b and c terms).
    """
    if el in EXCHANGE_ELEMENTS:
        need = {'form': 'exchange_KD', 'quantity': 'K_D', 'basis': 'mole', 'dIW_factor': 1.0}
        need |= {'d': 0.0, 'e': 0.0}
    else:
        need = {'form': 'log10_linear'}
    for field, value in need.items():
        if entry[field] != value:
            raise PartitionTableError(
                f'{el!r}: {field!r} must be {value!r} for this element, got {entry[field]!r}'
            )
    if el in EXCHANGE_ELEMENTS and entry['eps']:
        raise PartitionTableError(f'{el!r}: {"eps"!r} must be empty for an exchange law')


def _covariance(cov, size, el):
    """An empty list, or a symmetric positive semidefinite size x size matrix."""
    if isinstance(cov, list) and not cov:
        return []
    try:
        mat = np.array(
            [[_number(v, el, 'covariance') for v in row] for row in cov], dtype=float
        )
    except (TypeError, ValueError) as err:
        raise PartitionTableError(f'{el!r}: {"covariance"!r} must be a matrix') from err
    if mat.shape != (size, size):
        raise PartitionTableError(
            f'{el!r}: {"covariance"!r} must be {size} x {size} (a, b, c, d, e, eps), '
            f'got shape {mat.shape}'
        )
    if not np.allclose(mat, mat.T, rtol=1e-12, atol=0.0):
        raise PartitionTableError(f'{el!r}: {"covariance"!r} must be symmetric')
    if np.linalg.eigvalsh(mat).min() < -1e-12 * max(1.0, np.abs(mat).max()):
        raise PartitionTableError(f'{el!r}: {"covariance"!r} must be positive semidefinite')
    return mat.tolist()


def load_partition_table(path=None):
    """Read and check a partition coefficient file.

    Parameters
    ----------
    path : str or Path or None
        TOML file; None reads the packaged placeholder file.

    Returns
    -------
    dict
        ``schema_version``, ``meta`` and one checked entry per element in the file.

    Raises
    ------
    PartitionTableError
        A field that breaks the schema, named with its element.
    FileNotFoundError
        No file at ``path``.
    tomllib.TOMLDecodeError
        A file that is not valid TOML (a ``ValueError``).
    """
    if path is None:
        text = resources.files('calliope.data').joinpath(PLACEHOLDER_FILE).read_text()
    else:
        text = Path(path).read_text()
    raw = tomllib.loads(text)
    version = raw.pop('schema_version', None)
    if type(version) is not int or version != SCHEMA_VERSION:
        raise PartitionTableError(f'schema_version must be {SCHEMA_VERSION}, got {version!r}')
    meta = raw.pop('meta', None)
    name = meta.get('name') if isinstance(meta, dict) else None
    if not (
        isinstance(name, str) and name.strip() and isinstance(meta.get('placeholder'), bool)
    ):
        raise PartitionTableError(
            'the file needs a [meta] table with a non-empty name and a boolean placeholder'
        )
    unknown = [k for k in raw if k not in ELEMENTS]
    if unknown:
        raise PartitionTableError(f'unknown element sections {unknown!r}; known: {ELEMENTS}')
    table = {'schema_version': version, 'meta': dict(meta)}
    table.update({el: _entry(el, entry) for el, entry in raw.items()})
    return table


def _finite(name, value, positive=False):
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(
            f'{name} must be finite{" and positive" if positive else ""}, got {value!r}'
        )


def log10_D(entry, p_gpa, t_k, d_iw, x_metal=None, nbo_t=2.6):
    """log10 of the partition coefficient of one entry, in the entry's basis.

    Parameters
    ----------
    entry : dict
        One element entry of :func:`load_partition_table`.
    p_gpa : float
        Equilibration pressure in GPa.
    t_k : float
        Equilibration temperature in K.
    d_iw : float
        Oxygen fugacity relative to the iron-wustite buffer, log10 units.
    x_metal : dict or None
        Mole fractions of the elements of the eps terms in the metal; absent is 0.
    nbo_t : float
        NBO/T of the silicate melt.

    Returns
    -------
    float
        ``a + b/T + c P/T + d dIW_factor dIW + e NBO/T + sum_j eps_j log10(1 - X_j)``.

    Raises
    ------
    ValueError
        P, dIW or NBO/T not finite, T not finite and positive, or an X_j of an eps term
        outside [0, 1).
    """
    _finite('P [GPa]', p_gpa)
    _finite('T [K]', t_k, positive=True)
    _finite('dIW', d_iw)
    _finite('NBO/T', nbo_t)
    x_metal = x_metal or {}
    val = entry['a'] + entry['b'] / t_k + entry['c'] * p_gpa / t_k
    val += entry['d'] * entry['dIW_factor'] * d_iw + entry['e'] * nbo_t
    for el, eps in entry['eps'].items():
        x = x_metal.get(el, 0.0)
        if not 0.0 <= x < 1.0:
            raise ValueError(f'X_{el} in the metal must be in [0, 1), got {x!r}')
        val += eps * math.log10(1.0 - x)
    return val


def _pow10(x, label, p_gpa, t_k, d_iw):
    """10**x, or ValueError naming ``label`` and the conditions beyond the float range."""
    if not x < 308.0:
        raise ValueError(
            f'{label}: 10**{x:.4g} at P {p_gpa:g} GPa, T {t_k:g} K, dIW {d_iw:g} overflows'
        )
    return 10.0**x


def out_of_range(entry, p_gpa, t_k, d_iw):
    """True when (P, T, dIW) lies outside a range of the entry or on an unknown range.

    Parameters
    ----------
    entry : dict
        One element entry of :func:`load_partition_table`.
    p_gpa : float
        Equilibration pressure in GPa, against ``P_range_GPa``.
    t_k : float
        Equilibration temperature in K, against ``T_range_K``.
    d_iw : float
        Oxygen fugacity relative to the iron-wustite buffer, log10 units, against
        ``dIW_range``.

    Returns
    -------
    bool
        True when a value lies outside the closed range ``[low, high]`` of the entry, or
        when the T or dIW range is unknown (an empty range in the file).
    """
    ranges = (entry['P_range_GPa'], entry['T_range_K'], entry['dIW_range'])
    return any(r is None or not r[0] <= v <= r[1] for v, r in zip((p_gpa, t_k, d_iw), ranges))


# Takahashi 1986, JGR 91, 9367, Table 1, KLB-1 bulk: the 7 major oxides [wt%]
KLB1 = {
    'SiO2': 44.48,
    'TiO2': 0.16,
    'Al2O3': 3.59,
    'FeO': 8.10,
    'MgO': 39.22,
    'CaO': 3.44,
    'Na2O': 0.30,
}
# cation, cations and O atoms per formula unit
_OXIDES = {
    'SiO2': ('Si', 1, 2),
    'TiO2': ('Ti', 1, 2),
    'Al2O3': ('Al', 2, 3),
    'FeO': ('Fe', 1, 1),
    'MgO': ('Mg', 1, 1),
    'CaO': ('Ca', 1, 1),
    'Na2O': ('Na', 2, 1),
}
_OXIDE_MASS = {
    ox: n_cat * molar_mass_rock[cat] + n_o * molar_mass['O']
    for ox, (cat, n_cat, n_o) in _OXIDES.items()
}


def silicate_units_per_kg(oxides=KLB1):
    """Moles of single-cation units (SiO2, TiO2, AlO1.5, FeO, MgO, CaO, NaO0.5) per kg.

    Parameters
    ----------
    oxides : dict
        Melt composition in oxide wt%, keyed by the oxides of :data:`KLB1`; used as
        given, not normalised.

    Returns
    -------
    float
        Moles of single-cation units per kg of melt [mol/kg].
    """
    return sum(wt / 100.0 * _OXIDES[ox][1] / _OXIDE_MASS[ox] for ox, wt in oxides.items())


def core_light_elements(d_x, w_melt, w_fixed, oxides=KLB1):
    """Core H and C mass fractions from molar D on the H- and C-free metal basis.

    ``x'_i = D_i x_i`` with the melt mole fraction ``x_i`` on the single-cation basis,
    H counted as HO0.5 and C as one unit each, and ``x'_i = N_i / sum_{k != H, C} N_k``
    in the metal, whose other light elements are fixed and whose rest is Fe. The metal
    composition depends on H and C only through Fe, so the solve is closed form.

    Parameters
    ----------
    d_x : dict
        Molar D for 'H' and/or 'C'.
    w_melt : dict
        Mass fraction of H and C in the melt.
    w_fixed : dict
        Mass fractions of O, Si, S and N in the metal.
    oxides : dict
        Melt composition in oxide wt%.

    Returns
    -------
    dict
        Mass fraction of each element of ``d_x`` in the metal.

    Raises
    ------
    ValueError
        An element outside H and C in ``d_x``, outside O, Si, S, N in ``w_fixed``, or a
        metal with no iron left.
    """
    if not set(d_x) <= {'H', 'C'} or not set(w_fixed) <= {'O', 'Si', 'S', 'N'}:
        raise ValueError(f'd_x takes H and C, w_fixed O, Si, S and N; got {d_x}, {w_fixed}')
    m_fe = molar_mass_rock['Fe']
    rest = 1.0 - sum(w_fixed.values())
    n0 = rest / m_fe + sum(w / _mass(el) for el, w in w_fixed.items())
    n_sil = silicate_units_per_kg(oxides) + sum(w / molar_mass[el] for el, w in w_melt.items())
    xp = {el: d * (w_melt[el] / molar_mass[el]) / n_sil for el, d in d_x.items()}
    n_metal = n0 / (1.0 + sum(molar_mass[el] * x for el, x in xp.items()) / m_fe)
    out = {el: molar_mass[el] * x * n_metal for el, x in xp.items()}
    if rest - sum(out.values()) <= 0:
        raise ValueError(f'the metal has no iron left: fixed {w_fixed}, H and C {out}')
    return out


def _mass(el):
    return molar_mass_rock[el] if el in molar_mass_rock else molar_mass[el]


def metal_masses(mass_int_d, ddict, d_mass):
    """Element masses in the equilibrated metal [kg]: D_i M_metal w_i(melt).

    Parameters
    ----------
    mass_int_d : dict
        Element masses dissolved in the melt [kg].
    ddict : dict
        Reads ``metal_mass`` [kg] (0 or absent: no metal), ``M_mantle`` [kg] and
        ``Phi_global``; the melt mass is ``M_mantle * Phi_global``.
    d_mass : dict
        Mass-basis D per element (metal over melt mass fraction).

    Returns
    -------
    dict
        Metal mass of each element of ``d_mass``.

    Raises
    ------
    ValueError
        ``metal_mass`` negative or not finite, metal with no melt, or a D that is
        negative or not finite.
    """
    m_metal = ddict.get('metal_mass', 0.0)
    if not (math.isfinite(m_metal) and m_metal >= 0):
        raise ValueError(f'metal_mass must be finite and >= 0, got {m_metal!r}')
    if m_metal == 0.0:
        return {el: 0.0 for el in d_mass}
    m_melt = ddict['M_mantle'] * ddict['Phi_global']
    if not m_melt > 0:
        raise ValueError(
            f'metal_mass {m_metal!r} kg needs melt, but M_mantle * Phi_global is {m_melt!r}'
        )
    out = {}
    for el, d in d_mass.items():
        if not (math.isfinite(d) and d >= 0):
            raise ValueError(f'D of {el} must be finite and >= 0, got {d!r}')
        out[el] = d * m_metal * mass_int_d[el] / m_melt
    return out


def _metal_moles(w_metal, exclude=()):
    """Moles per kg of metal of the elements outside ``exclude``, the rest being Fe."""
    rest = 1.0 - sum(w_metal.values())
    if rest <= 0:
        raise ValueError(f'the metal has no iron left: {w_metal}')
    n = {el: w / _mass(el) for el, w in w_metal.items() if el not in exclude}
    return rest / molar_mass_rock['Fe'] + sum(n.values()), n


def mass_D(entry, p_gpa, t_k, d_iw, w_metal, w_melt, oxides=KLB1, nbo_t=2.6):
    """Mass-basis D (metal over melt mass fraction) of one entry.

    Parameters
    ----------
    entry : dict
        One element entry of :func:`load_partition_table`, form ``log10_linear``.
    p_gpa, t_k, d_iw : float
        Equilibration pressure [GPa], temperature [K] and dIW [log10 units].
    w_metal : dict
        Light-element mass fractions of the metal (the rest is Fe).
    w_melt : dict
        Mass fractions of H and C in the melt, counted as HO0.5 and C units.
    oxides : dict
        Melt composition in oxide wt%.
    nbo_t : float
        NBO/T of the melt.

    Returns
    -------
    float
        ``D_w = D N_metal / N_melt`` for a mole basis, with N the moles per kg of the
        metal (without H and C for ``mole_xprime_noHC``) and of the melt units; ``D``
        for the mass basis.

    Raises
    ------
    PartitionTableError
        An entry whose form is not ``log10_linear``.
    ValueError
        A metal with no iron left, a D beyond the float range (naming the law and its P,
        T and dIW), or an invalid input of :func:`log10_D`.

    Notes
    -----
    The X_j of the eps terms are ``x'_j = N_j / sum_{k != H, C} N_k``, the metal basis of
    Tagawa et al. (2021, Methods Eq. 4) and of the H law. Fischer et al. (2020, p. 8744)
    fit X_O and X_S as mole fractions in metal without H and do not name the elements of
    the denominator; x' is the choice of this module.
    """
    _need_form(entry, 'log10_linear')
    total, n = _metal_moles(w_metal)
    free = total - n.get('H', 0.0) - n.get('C', 0.0)
    x_metal = {el: v / free for el, v in n.items()}
    label = f'{entry["element"]} law'
    d = _pow10(log10_D(entry, p_gpa, t_k, d_iw, x_metal, nbo_t), label, p_gpa, t_k, d_iw)
    if entry['basis'] == 'mass':
        return d
    n_metal = free if entry['basis'] == 'mole_xprime_noHC' else total
    n_melt = silicate_units_per_kg(oxides) + sum(w / molar_mass[el] for el, w in w_melt.items())
    return d * n_metal / n_melt


def cation_fraction(oxides, name):
    """Mole fraction of oxide ``name`` on the single-cation basis of the melt.

    Parameters
    ----------
    oxides : dict
        Melt composition in oxide wt%, keyed by the oxides of :data:`KLB1`.
    name : str
        One oxide of ``oxides``, for example ``'FeO'`` or ``'SiO2'``.

    Returns
    -------
    float
        Cation units of ``name`` over all cation units of the melt, from 0 to 1 (an
        oxide with 2 cations, such as Al2O3, counts 2 units of AlO1.5).
    """
    return (
        oxides[name]
        / 100.0
        / _OXIDE_MASS[name]
        * _OXIDES[name][1]
        / silicate_units_per_kg(oxides)
    )


def exchange_metal(table, p_gpa, t_k, d_iw, oxides=KLB1, elements=('O', 'Si'), check=True):
    """Metal mole fractions of O and Si from their exchange K_D, ideal activities.

    ``FeO = Fe + O``: ``x_O = K_D x_FeO / x_Fe = K_D 10^(dIW/2)``; ``SiO2 + 2 Fe = Si +
    2 FeO``: ``x_Si = K_D x_SiO2 (x_Fe / x_FeO)^2 = K_D x_SiO2 10^(-dIW)``, with
    ``dIW = 2 log10(x_FeO / x_Fe)``. Melt fractions are on the single-cation basis; metal
    fractions are x' over the metal without H and C (Tagawa et al. 2021, Methods Eq. 4).
    Activity coefficients are 1; the fits of Siebert et al. come with activity terms
    that this release does not use.

    Only the laws of ``elements`` are evaluated. ``check=False`` returns finite fractions
    that leave no iron as they are, for a fixed-point pass that shortens its step.

    Parameters
    ----------
    table : dict
        Table of :func:`load_partition_table` with an ``exchange_KD`` law (``K_D`` on the
        mole basis) for each element of ``elements``.
    p_gpa : float
        Equilibration pressure in GPa.
    t_k : float
        Equilibration temperature in K.
    d_iw : float
        Oxygen fugacity of the metal-melt equilibrium relative to the iron-wustite
        buffer, log10 units.
    oxides : dict
        Melt composition in oxide wt%; its SiO2 sets ``x_SiO2`` of the Si exchange.
    elements : tuple of str
        The exchanged elements to evaluate, a subset of ('O', 'Si').
    check : bool
        Raise when the O and Si fractions sum to 1 or more.

    Returns
    -------
    dict
        Mole fraction of each element of ``elements`` in the metal, on the basis of the
        metal without H and C (x').

    Raises
    ------
    PartitionTableError
        A law of ``elements`` whose form is not ``exchange_KD``.
    ValueError
        P not finite, T not finite and positive, dIW not finite, a K_D or exchange term
        beyond the float range (naming the law and its P, T and dIW), or O and Si that
        leave no iron in the metal (with ``check``) or are not finite.
    """
    for el in elements:
        _need_form(table[el], 'exchange_KD')
    x = {}
    if 'O' in elements:
        k_d = _pow10(log10_D(table['O'], p_gpa, t_k, 0.0), 'O law', p_gpa, t_k, d_iw)
        x['O'] = k_d * _pow10(d_iw / 2.0, 'O exchange', p_gpa, t_k, d_iw)
    if 'Si' in elements:
        x_sio2 = cation_fraction(oxides, 'SiO2')
        k_d = _pow10(log10_D(table['Si'], p_gpa, t_k, 0.0), 'Si law', p_gpa, t_k, d_iw)
        x['Si'] = k_d * x_sio2 * _pow10(-d_iw, 'Si exchange', p_gpa, t_k, d_iw)
    total = sum(x.values())
    if not math.isfinite(total):
        raise ValueError(
            f'the O and Si exchange is not finite ({x}) at P {p_gpa:g} GPa, T {t_k:g} K, '
            f'dIW {d_iw:g}'
        )
    if check and total >= 1.0:
        raise ValueError(f'O and Si leave the metal with no iron: {x}')
    return x

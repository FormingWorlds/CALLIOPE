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

ELEMENTS = ('H', 'C', 'N', 'S', 'O', 'Si')
SCHEMA_VERSION = 1
FORMS = ('log10_linear', 'exchange_KD')
QUANTITIES = ('D', 'K_D')
BASES = ('mass', 'mole', 'mole_xprime_noHC')
SILICATE_BASES = ('single_cation',)
COEFFICIENTS = ('a', 'b', 'c', 'd', 'e', 'dIW_factor')
RANGES = ('P_range_GPa', 'T_range_K', 'dIW_range')
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
    out = {
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
    return out


def _covariance(cov, size, el):
    """An empty list, or a symmetric positive semidefinite size x size matrix."""
    if isinstance(cov, list) and not cov:
        return []
    try:
        mat = np.array(
            [[_number(v, el, 'covariance') for v in row] for row in cov], dtype=float
        )
    except TypeError as err:
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
    """
    if path is None:
        text = resources.files('calliope.data').joinpath(PLACEHOLDER_FILE).read_text()
    else:
        text = Path(path).read_text()
    raw = tomllib.loads(text)
    version = raw.pop('schema_version', None)
    if isinstance(version, bool) or version != SCHEMA_VERSION:
        raise PartitionTableError(f'schema_version must be {SCHEMA_VERSION}, got {version!r}')
    meta = raw.pop('meta', None)
    if not isinstance(meta, dict) or not isinstance(meta.get('placeholder'), bool):
        raise PartitionTableError('the file needs a [meta] table with a boolean placeholder')
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


def log10_D(entry, p_gpa, t_k, d_iw, x_metal=None, nbo_t=0.0):
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


def out_of_range(entry, p_gpa, t_k, d_iw):
    """True when (P, T, dIW) lies outside a range of the entry or on an unknown range."""
    for value, rng in zip(
        (p_gpa, t_k, d_iw), (entry['P_range_GPa'], entry['T_range_K'], entry['dIW_range'])
    ):
        if rng is None or not rng[0] <= value <= rng[1]:
            return True
    return False

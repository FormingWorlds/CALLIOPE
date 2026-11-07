"""Core and mantle split of a bulk element inventory at one equilibration point.

The metal of mass ``M_metal`` equilibrates with a silicate of mass ``M_silicate`` at
``(P_eq, T_eq, dIW)``. For a mass-basis D_i the silicate mass fraction is
``w_i = M_i / (M_silicate + D_i M_metal)``, the core holds ``D_i M_metal w_i`` and the
mantle ``M_silicate w_i``. D_i depends on the metal and melt composition (interaction
terms, mole bases), so the split is a fixed point on the metal composition. The laws come
from :func:`calliope.metal.get_law` by name (:data:`DEFAULT_LAWS`) or as entries. ``dIW``
is the one redox input: a law with an FeO term (the S law) takes the melt FeO mole fraction
``X_FeO = x'_Fe 10^(dIW/2)`` of ideal mixing from the metal of the pass, and the melt
``oxides`` set only the silicate units of the mole bases.

A law outside its (P, T, dIW) range is evaluated at the given point, not clamped, and
flagged with its distance per axis. O and Si are not split: their bulk is the O and Si of
the rock, which CALLIOPE does not hold; their metal mass fractions enter as ``w_fixed``.
The exchange-K_D split of O and Si is not part of this module.
"""

from __future__ import annotations

import logging
import math

from .metal import (
    ELEMENTS,
    EXCHANGE_ELEMENTS,
    KLB1,
    PartitionTableError,
    _metal_moles,
    _real,
    check_law,
    get_law,
    mass_D,
    range_distance,
)

log = logging.getLogger('fwl.' + __name__)

DEFAULT_LAWS = {
    'H': 'tagawa2021_eq2',
    'C': 'fischer2020_eq1',
    'N': 'suer2023_line',
    'S': 'suer2017_eq7',
}
P_EQ_GPA = 60.0  # Suer et al. 2023 Table 1, deep magma ocean column
T_EQ_K = 4000.0
MAX_PASSES = 50
TOL_W = 1e-12  # relative step of the metal mass fractions that ends the fixed point
TOL_W_FLOOR = 1e-15  # absolute floor of that step, for an element near 0
LEDGER_RTOL = 1e-12
LEDGER_ATOL = 2 * math.ulp(0.0)  # round-off of a subnormal bulk [kg]

# (element, law) pairs whose out-of-range warning was given in this process
_WARNED = set()


def _check(name, value, positive=False, signed=False):
    """A finite real number (not a boolean), >= 0 (or > 0, or of any sign)."""
    ok = _real(value) and (signed or value > 0 or (value == 0 and not positive))
    if not ok:
        rule = 'finite' if signed else ('finite and > 0' if positive else 'finite and >= 0')
        raise ValueError(f'{name} must be {rule}, got {value!r}')
    return float(value)


def _masses(bulk, d, m_metal, m_sil):
    """Core and mantle masses [kg] of each element at fixed mass-basis D."""
    total = {el: m_sil + d[el] * m_metal for el in bulk}
    for el in (el for el, t in total.items() if not math.isfinite(t)):
        raise ValueError(f'{el!r}: M_silicate + D M_metal is not finite (D {d[el]!r})')
    core = {el: d[el] * m_metal / total[el] * m for el, m in bulk.items()}
    return core, {el: m_sil / total[el] * m for el, m in bulk.items()}


def split_core_mantle(
    bulk_kg,
    M_metal,
    M_silicate,
    *,
    dIW,
    P_eq=P_EQ_GPA,
    T_eq=T_EQ_K,
    laws=None,
    w_fixed=None,
    oxides=KLB1,
    nbo_t=2.6,
):
    """Split the bulk mass of each element between core and mantle.

    Parameters
    ----------
    bulk_kg : dict
        Bulk mass [kg] of each element to split, H, C, N or S (or any element of
        :data:`calliope.metal.ELEMENTS` with a D law); O and Si are refused.
    M_metal : float
        Mass of the equilibrating metal [kg], >= 0.
    M_silicate : float
        Mass of the equilibrating silicate [kg], > 0.
    dIW : float
        Oxygen fugacity of the equilibrium relative to the iron-wustite buffer [log10].
    P_eq, T_eq : float
        Equilibration pressure [GPa] and temperature [K], > 0.
    laws : dict or None
        Law per element, a name for :func:`calliope.metal.get_law` or an entry; elements
        not given take :data:`DEFAULT_LAWS`.
    w_fixed : dict or None
        Metal mass fractions of elements that are not split (O, Si), each in [0, 1) with a
        sum below 1.
    oxides : dict
        Silicate composition in oxide wt%: the silicate units of the mole bases and the
        melt a callable law gets (the FeO of a law with an FeO term comes from ``dIW``).
    nbo_t : float
        NBO/T of the silicate.

    Returns
    -------
    core_kg, mantle_kg : dict
        Core and mantle mass [kg] of each element of ``bulk_kg``; their sum is the bulk to
        1e-12 relative plus 2 ulp of 0 (asserted).
    flags : dict
        Per element: ``out_of_range`` (1 when a distance is not 0) and the distances
        ``P_GPa``, ``T_K``, ``dIW`` of :func:`calliope.metal.range_distance`. Under
        ``'metal'``: the mole fractions ``x_light`` (all light elements) and ``x_S`` over
        the whole metal, ``S_over_Fe`` (moles), and ``beyond_FeS`` (1 when the metal holds
        more S than Fe, past FeS; no clamp).

    Raises
    ------
    ValueError
        An input outside its domain (named), O or Si in ``bulk_kg``, a key of ``laws``
        outside ``bulk_kg``, a D that overflows with ``M_metal``, or a metal with no iron
        left (naming the split).
    calliope.metal.PartitionTableError
        A law name that is not held, or a law in ``laws`` that is not a D law of its
        element.
    RuntimeError
        A metal composition that does not settle in :data:`MAX_PASSES` passes, or a split
        that breaks the element ledger.
    """
    bulk = {}
    for el, m in bulk_kg.items():
        if el in EXCHANGE_ELEMENTS:
            raise ValueError(
                f'{el} is not split: its bulk is the rock {el}, which CALLIOPE does not hold; '
                f'give the metal {el} mass fraction in w_fixed'
            )
        if el not in ELEMENTS:
            raise ValueError(f'bulk element {el!r} has no partition law; known: {ELEMENTS}')
        bulk[el] = _check(f'bulk_kg[{el!r}]', m)
    m_metal = _check('M_metal', M_metal)
    m_sil = _check('M_silicate', M_silicate, positive=True)
    p_eq = _check('P_eq', P_eq, positive=True)
    t_eq = _check('T_eq', T_eq, positive=True)
    d_iw = _check('dIW', dIW, signed=True)
    fixed = {}
    for el, w in (w_fixed or {}).items():
        if el not in ELEMENTS or el in bulk:
            raise ValueError(f'w_fixed: {el!r} must be an element of {ELEMENTS} not in bulk_kg')
        fixed[el] = _check(f'w_fixed[{el!r}]', w)
        if fixed[el] >= 1.0:
            raise ValueError(f'w_fixed[{el!r}] = {w!r} must be below 1')
    if sum(fixed.values()) >= 1.0:
        raise ValueError(f'w_fixed {fixed} leaves no iron in the metal')
    if not isinstance(laws or {}, dict):
        raise ValueError(f'laws must be a dict of element to law, got {laws!r}')
    extra = sorted(set(laws or {}) - set(bulk))
    if extra:
        raise ValueError(f'laws names elements outside bulk_kg: {extra}')
    names, entries = {**DEFAULT_LAWS, **(laws or {})}, {}
    for el in bulk:
        spec = names[el]
        try:
            entries[el] = check_law(get_law(el, spec) if isinstance(spec, str) else spec, el)
        except PartitionTableError as err:
            raise PartitionTableError(f'laws[{el!r}] = {spec!r:.60}: {err}') from err
    w, w_melt = {**fixed, **dict.fromkeys(bulk, 0.0)}, {}
    for _ in range(MAX_PASSES):
        d = {}
        for el in bulk:
            try:
                d[el] = mass_D(entries[el], p_eq, t_eq, d_iw, w, w_melt, oxides, nbo_t)
            except ValueError as err:
                raise ValueError(f'split of {el}: {err}') from err
        core, mantle = _masses(bulk, d, m_metal, m_sil)
        new = {**fixed, **{el: core[el] / m_metal if m_metal else 0.0 for el in bulk}}
        w_melt = {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}
        settled = all(abs(new[el] - w[el]) <= TOL_W * abs(new[el]) + TOL_W_FLOOR for el in bulk)
        w = new
        if settled:
            break
    else:
        raise RuntimeError(f'the metal of the split did not settle in {MAX_PASSES} passes: {w}')
    for el, m in bulk.items():
        if not (
            core[el] >= 0.0
            and mantle[el] >= 0.0
            and abs(core[el] + mantle[el] - m) <= LEDGER_RTOL * m + LEDGER_ATOL
        ):
            raise RuntimeError(
                f'{el} ledger of the split fails: core {core[el]!r} + mantle {mantle[el]!r} kg '
                f'against bulk {m!r} kg'
            )
    total, n = _metal_moles(w)
    n_s, n_fe = n.get('S', 0.0), total - sum(n.values())
    flags = {
        'metal': {
            'x_light': sum(n.values()) / total,
            'x_S': n_s / total,
            'S_over_Fe': n_s / n_fe,
            'beyond_FeS': int(n_s > n_fe),
        }
    }
    for el, entry in entries.items():
        dist = range_distance(entry, p_eq, t_eq, d_iw)
        flags[el] = {'out_of_range': int(any(v != 0.0 for v in dist.values())), **dist}
        label = names[el] if isinstance(names[el], str) else f'entry ({entry["source"][:40]})'
        key = (el, label) + tuple(entry[r] for r in ('P_range_GPa', 'T_range_K', 'dIW_range'))
        if flags[el]['out_of_range'] and key not in _WARNED:
            _WARNED.add(key)
            shown = {k: round(v, 3) for k, v in dist.items()}
            log.warning(
                f'{el} partition law {label} evaluated outside its range: P_eq {p_eq:g} GPa, '
                f'T_eq {t_eq:g} K, dIW {d_iw:g} (distance {shown})'
            )
    return core, mantle, flags

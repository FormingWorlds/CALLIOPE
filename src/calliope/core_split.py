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
    _mass,
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
MAX_RELAXED_PASSES = 1000
DAMPED_LAM = 0.02  # step factor of the passes that look for a root with more iron
MAX_DAMPED_PASSES = 20000
IRON_FLOOR = 1e-8  # iron mass fraction below which a relaxed pass raises ValueError
GRID_RTOL = 1e-9  # step of a steep g that the float grid of w cannot resolve further
TOL_W = 1e-12  # relative step of the metal mass fractions that ends the fixed point
TOL_W_FLOOR = 1e-15  # absolute floor of that step, for an element near 0
LEDGER_RTOL = 1e-12
LEDGER_ATOL = 2 * math.ulp(0.0)  # round-off of a subnormal bulk [kg]

# (element, law, source, ranges) of the out-of-range warnings given in this process
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


class _Unsettled(Exception):
    """Plain substitution left the metal with no iron or did not settle."""


class _IronOut(ValueError):
    """A relaxed pass left the metal with iron below IRON_FLOOR."""


def _settled(g, w):
    return all(abs(g[el] - w[el]) <= TOL_W * abs(g[el]) + TOL_W_FLOOR for el in g)


def _plain_fixed_point(evaluate, start, fixed):
    """Plain substitution w <- g(w) until the metal and the melt H and C settle.

    An error after the first pass or no settling in MAX_PASSES raises _Unsettled.
    """
    w, w_melt = start, {}
    for k in range(MAX_PASSES):
        try:
            core, mantle, g, new_melt = evaluate(w, w_melt)
        except ValueError:
            if k == 0:
                raise
            raise _Unsettled from None
        melt_settled = _settled(new_melt, w_melt) if w_melt else not new_melt
        settled, w, w_melt = _settled(g, w) and melt_settled, {**fixed, **g}, new_melt
        if settled:
            return core, mantle, w
    raise _Unsettled


def _relaxed_fixed_point(evaluate, start, fixed):
    """Relaxed passes, first with a step factor per element, then with one for all elements.

    Each element moves lam = 1/(1 + |s|) of its step g - w, s being the secant slope of g against
    w of the element over its last 2 passes (as in the metal reservoir of calliope.solve for a
    falling g; a rising slope from the coupling of elements is damped the same way). If those
    passes do not settle or leave no iron, they start again from ``start`` with the smallest lam
    for every element, which keeps coupled moves in step. Raises ValueError when the common
    passes leave no iron, or do not settle after the first passes left no iron, and RuntimeError
    when neither settles. A return at w with the metal beyond FeS goes to _more_iron; the 4th
    item is 1 when its damped root replaced the first root, else 0.
    """
    iron_out = None
    try:
        out = _relaxed_passes(evaluate, start, fixed, False)
    except _IronOut as err:
        out, iron_out = None, err
    if out is None:
        try:
            out = _relaxed_passes(evaluate, start, fixed, True)
        except _IronOut as err:
            raise ValueError(str(err)) from None
    if out is not None:
        return _more_iron(evaluate, start, fixed, out)
    if iron_out is not None:
        raise ValueError(
            f'{iron_out}; the common passes then did not settle in {MAX_RELAXED_PASSES} passes'
        ) from None
    raise RuntimeError(
        f'the metal of the split did not settle in 2 x {MAX_RELAXED_PASSES} relaxed passes'
    )


def _iron(w):
    return 1.0 - sum(w.values())


def _more_iron(evaluate, start, fixed, out):
    """The root of damped passes from ``start`` when it holds more iron than a return at w.

    A return at w (``out[3]``) whose metal holds more S than Fe in moles can sit on a root
    with little iron where an iron-rich root exists too; passes at lam DAMPED_LAM for every
    element (melt undamped) reach the iron-rich root. Their root, settled or returned at w,
    replaces ``out`` when it holds more iron. This prefers that root; it does not prove the
    root unique. Damped passes that do not settle or raise ValueError (no iron, a law) keep out.
    """
    core, mantle, w, at_w = out
    if not at_w or w.get('S', 0.0) / _mass('S') <= _iron(w) / _mass('Fe'):
        return core, mantle, w, 0
    try:
        damped = _relaxed_passes(evaluate, start, fixed, True, DAMPED_LAM, MAX_DAMPED_PASSES)
    except ValueError:
        damped = None
    if damped is not None and _iron(damped[2]) > _iron(w):
        return (*damped[:3], 1)
    return core, mantle, w, 0


def _bracketed(evaluate, w, w_melt, el, ulps=4):
    """g - w of element el changes sign from + to - across ``ulps`` floats each side of w[el]."""

    def excess(x):
        try:
            return evaluate({**w, el: x}, w_melt)[2][el] - x
        except ValueError:
            return math.nan

    step = ulps * math.ulp(w[el])
    return excess(w[el] - step) > 0.0 > excess(w[el] + step)


def _relaxed_passes(evaluate, start, fixed, common, lam_all=None, passes=None):
    """Relaxed passes from ``start``: (core, mantle, w, at_w), or None when they do not settle.

    ``common`` moves every element at one lam, ``lam_all`` or else the smallest own lam;
    ``passes`` caps the passes (default MAX_RELAXED_PASSES).
    The melt H and C of the pass must have settled too, since they set the silicate units of
    the mole bases. A pass adds at most half of the iron left. An element has settled when its
    step meets the settle test, or when g falls with w, its move at its own lam is within 2 ulp
    of w and its step within GRID_RTOL of g (the float grid of a steep g), and g - w changes
    sign within 4 ulp each side of the evaluated w; such a pass is returned at w (core mass
    fractions w, mantle the rest of the bulk) when no mantle is negative. A pass that leaves
    iron below IRON_FLOOR raises _IronOut.
    """
    w, w_melt, last = start, {}, {}
    for _ in range(MAX_RELAXED_PASSES if passes is None else passes):
        core, mantle, g, new_melt = evaluate(w, w_melt)
        melt_settled = _settled(new_melt, w_melt) if w_melt else not new_melt
        w_melt_prev, w_melt = w_melt, new_melt
        slope = {}
        for el, gv in g.items():
            slope[el] = 0.0
            if el in last and w[el] != last[el][0]:
                slope[el] = (gv - last[el][1]) / (w[el] - last[el][0])
            last[el] = (w[el], gv)
        own = {el: 1.0 / (1.0 + abs(s)) for el, s in slope.items()}
        lam = own
        if common:
            lam = dict.fromkeys(own, min(own.values()) if lam_all is None else lam_all)
        move = {el: lam[el] * (gv - w[el]) for el, gv in g.items()}
        if melt_settled and all(
            _settled({el: g[el]}, {el: w[el]})
            or (
                slope[el] < 0.0
                and own[el] * abs(g[el] - w[el]) <= 2.0 * math.ulp(w[el])
                and abs(g[el] - w[el]) <= GRID_RTOL * abs(g[el])
            )
            for el in g
        ):
            if _settled(g, w):
                return core, mantle, {**fixed, **g}, False
            if all(
                _settled({el: g[el]}, {el: w[el]}) or _bracketed(evaluate, w, w_melt_prev, el)
                for el in g
            ):
                at_w = {el: w[el] * core[el] / g[el] if g[el] else 0.0 for el in g}
                rest = {el: core[el] + mantle[el] - at_w[el] for el in g}
                if min(rest.values()) >= 0.0:
                    return at_w, rest, {**fixed, **{el: w[el] for el in g}}, True
        while sum(move.values()) > 0.5 * (1.0 - sum(w.values())):
            move = {el: 0.5 * v for el, v in move.items()}
        w = {**w, **{el: w[el] + v for el, v in move.items()}}
        if 1.0 - sum(w.values()) < IRON_FLOOR:
            raise _IronOut(
                f'the metal has no iron left: a pass leaves iron below {IRON_FLOOR:g} by mass '
                f'for this bulk and metal mass ({w})'
            )
    return None


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
        Mass of the equilibrating metal [kg], >= 0; 0 evaluates no law and returns the bulk
        as mantle.
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
        more S than Fe, past FeS; no clamp). A relaxed split past FeS prefers the root with
        more iron (see the metal page); ``root_rule`` is 1 when that rule replaced the first
        root, else 0. Per element also ``D_residual``: |w - D w_sil| / w
        at the returned split, with w = core / M_metal and w_sil = mantle / M_silicate (0 for
        no core mass).

    Raises
    ------
    ValueError
        An input outside its domain (named), O or Si in ``bulk_kg``, a key of ``laws``
        outside ``bulk_kg``, a D that overflows with ``M_metal``, or a metal with no iron
        left (naming the split), also when the common relaxed passes do not settle after
        the per-element passes ran out of iron.
    calliope.metal.PartitionTableError
        A law name that is not held, or a law in ``laws`` that is not a D law of its
        element.
    RuntimeError
        A metal composition that settles neither in :data:`MAX_PASSES` plain passes nor in
        2 x :data:`MAX_RELAXED_PASSES` relaxed passes, or a split that breaks the element
        ledger.
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
    if laws is not None and not isinstance(laws, dict):
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
    start = {**fixed, **dict.fromkeys(bulk, 0.0)}

    def evaluate(w, w_melt):
        d = {}
        for el in bulk:
            try:
                d[el] = mass_D(entries[el], p_eq, t_eq, d_iw, w, w_melt, oxides, nbo_t)
            except ValueError as err:
                raise ValueError(f'split of {el}: {err}') from err
        core, mantle = _masses(bulk, d, m_metal, m_sil)
        g = {el: core[el] / m_metal for el in bulk}
        return core, mantle, g, {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}

    root_rule = 0
    if m_metal == 0.0:
        core, mantle, w = dict.fromkeys(bulk, 0.0), dict(bulk), start
    else:
        try:
            core, mantle, w = _plain_fixed_point(evaluate, start, fixed)
        except _Unsettled:
            core, mantle, w, root_rule = _relaxed_fixed_point(evaluate, start, fixed)
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
    w_ret = {**fixed, **{el: core[el] / m_metal for el in bulk}} if m_metal else w
    total, n = _metal_moles(w_ret)
    n_s, n_fe = n.get('S', 0.0), (1.0 - sum(w_ret.values())) / _mass('Fe')
    flags = {
        'metal': {
            'x_light': sum(n.values()) / total,
            'x_S': n_s / total,
            'S_over_Fe': n_s / n_fe,
            'beyond_FeS': int(n_s > n_fe),
            'root_rule': root_rule,
        }
    }
    w_melt = {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}
    for el, entry in entries.items():
        dist = range_distance(entry, p_eq, t_eq, d_iw)
        flags[el] = {'out_of_range': int(any(v != 0.0 for v in dist.values())), **dist}
        flags[el]['D_residual'] = 0.0
        if core[el] > 0.0:
            d = mass_D(entry, p_eq, t_eq, d_iw, w_ret, w_melt, oxides, nbo_t)
            flags[el]['D_residual'] = (
                abs(core[el] - d * mantle[el] * m_metal / m_sil) / core[el]
            )
        label = names[el] if isinstance(names[el], str) else f'entry ({entry["source"][:40]})'
        ranges = tuple(entry[r] for r in ('P_range_GPa', 'T_range_K', 'dIW_range'))
        key = (el, label, entry['source']) + ranges
        if flags[el]['out_of_range'] and key not in _WARNED:
            _WARNED.add(key)
            shown = {k: round(v, 3) for k, v in dist.items()}
            log.warning(
                f'{el} partition law {label} evaluated outside its range: P_eq {p_eq:g} GPa, '
                f'T_eq {t_eq:g} K, dIW {d_iw:g} (distance {shown})'
            )
    return core, mantle, flags

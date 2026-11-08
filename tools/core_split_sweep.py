"""Random splits of ``split_core_mantle``: the sampler behind the sweep counts in the metal docs.

Usage, with the calliope under test first on ``sys.path`` (for example ``PYTHONPATH=src``)::

    python tools/core_split_sweep.py run broad out.json
    python tools/core_split_sweep.py run small-silicate out.json
    python tools/core_split_sweep.py compare old.json new.json

``run`` writes one JSON row per split (input, outcome, core, mantle, path, D residual) and prints
the outcome counts. The D residual of an element is |w - D(w, w_melt) w_sil| / w at the returned
split, with w = core / M_metal and w_sil = mantle / M_silicate. ``compare`` prints, for 2 runs of
the same stratum, the splits that change, that stop converging and that start to converge.
"""

from __future__ import annotations

import collections
import json
import sys

import numpy as np

from calliope import core_split
from calliope.core_split import DEFAULT_LAWS, P_EQ_GPA, T_EQ_K, split_core_mantle
from calliope.metal import get_law, mass_D

# name: (seed, splits, log10 range of M_silicate [kg])
STRATA = {'broad': (11, 3000, (16.0, 24.6)), 'small-silicate': (12, 2500, (18.0, 20.0))}


def sample(stratum):
    """Yield (bulk, M_metal, M_silicate, dIW, w_fixed) of each split of a stratum."""
    seed, n, (lo, hi) = STRATA[stratum]
    rng = np.random.default_rng(seed)
    for _ in range(n):
        m_metal = 10 ** rng.uniform(20, 24.5)
        m_sil = 10 ** rng.uniform(lo, hi)
        els = [el for el in ('H', 'C', 'N', 'S') if rng.random() < 0.6] or ['S']
        bulk = {el: 10 ** rng.uniform(-6, 0.4) * m_metal for el in els}
        diw = rng.uniform(-5, -0.2)
        w_fixed = {'O': float(rng.choice([0.0, 0.02, 0.1]))}
        yield bulk, m_metal, m_sil, diw, w_fixed


def run_one(bulk, m_metal, m_sil, diw, w_fixed):
    """Outcome row of one split."""
    relaxed = []
    real = core_split._relaxed_fixed_point
    core_split._relaxed_fixed_point = lambda *a: (relaxed.append(1), real(*a))[1]
    try:
        core, mantle, _ = split_core_mantle(bulk, m_metal, m_sil, dIW=diw, w_fixed=w_fixed)
    except Exception as err:  # older versions of the split raise other types too
        return {'outcome': type(err).__name__, 'message': str(err), 'relaxed': bool(relaxed)}
    finally:
        core_split._relaxed_fixed_point = real
    w = {**w_fixed, **{el: core[el] / m_metal for el in bulk}}
    w_melt = {el: mantle[el] / m_sil for el in ('H', 'C') if el in bulk}
    res = {}
    for el in bulk:
        d = mass_D(get_law(el, DEFAULT_LAWS[el]), P_EQ_GPA, T_EQ_K, diw, w, w_melt)
        res[el] = abs(w[el] - d * mantle[el] / m_sil) / w[el] if w[el] else 0.0
    return {
        'outcome': 'ok',
        'core': core,
        'mantle': mantle,
        'iron': 1.0 - sum(w.values()),
        'relaxed': bool(relaxed),
        'residual': res,
    }


def run(stratum, out):
    rows = [{'input': list(case), **run_one(*case)} for case in sample(stratum)]
    with open(out, 'w') as fh:
        fh.write(
            '[\n' + ',\n'.join(json.dumps(r, separators=(',', ':')) for r in rows) + '\n]\n'
        )
    print(stratum, dict(collections.Counter(r['outcome'] for r in rows)))


def compare(old_path, new_path):
    old, new = json.load(open(old_path)), json.load(open(new_path))
    assert [r['input'] for r in old] == [r['input'] for r in new], 'not the same stratum'
    changed, lost, gained, other = [], [], [], []
    for k, (a, b) in enumerate(zip(old, new)):
        if a['outcome'] != b['outcome'] and 'ok' not in (a['outcome'], b['outcome']):
            other.append((k, a['outcome'], b['outcome']))
        if a['outcome'] == 'ok' and b['outcome'] == 'ok':
            rel = max(
                abs(b[p][el] - a[p][el]) / a[p][el]
                for p in ('core', 'mantle')
                for el in a[p]
                if a[p][el]
            )
            if rel > 0.0:
                changed.append((k, rel))
        elif a['outcome'] == 'ok':
            lost.append((k, b['outcome']))
        elif b['outcome'] == 'ok':
            gained.append((k, a['outcome']))
    print(
        f'both converge {sum(a["outcome"] == b["outcome"] == "ok" for a, b in zip(old, new))}; '
        f'changed {len(changed)} (above 1e-9: {sum(r > 1e-9 for _, r in changed)}, '
        f'above 1e-3: {sum(r > 1e-3 for _, r in changed)}); lost {lost}; gained {gained}; '
        f'error type changed {other}'
    )


if __name__ == '__main__':
    {'run': run, 'compare': compare}[sys.argv[1]](*sys.argv[2:])

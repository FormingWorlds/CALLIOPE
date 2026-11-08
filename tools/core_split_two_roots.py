"""Round-off starts of the 2-root S-C-H-N split: the root of each, with and without the rule.

Usage, with the calliope under test first on ``sys.path`` (for example ``PYTHONPATH=src``)::

    python tools/core_split_two_roots.py

The split (silicate 1.8e-5 of the metal mass, metal O 0.02) has roots at iron 5.0e-2 and
6.2e-7 by mass. Each of 30 starts (seed 3) scales every input by 1 + e, e uniform in
[-4e-16, 4e-16]; start 0 is the unperturbed input. The run without the rule replaces
``core_split._more_iron`` by the first relaxed return.
"""

from __future__ import annotations

import numpy as np

from calliope import core_split
from calliope.core_split import split_core_mantle

BULK = {
    'S': 1.6351025151692482e24,
    'H': 4.026924495995128e22,
    'C': 2.0258800008695458e24,
    'N': 3.113694670982211e23,
}
M_METAL, M_SIL, DIW, W_FIXED = (
    2.6777700956688316e24,
    4.869342993603356e19,
    -0.9403058010523511,
    {'O': 0.02},
)


def irons(starts=30, seed=3):
    """Iron mass fraction of the metal returned by each start."""
    rng, out = np.random.default_rng(seed), []
    for k in range(starts):
        eps = rng.uniform(-1, 1, 7) * 4e-16 if k else np.zeros(7)
        bulk = {el: m * (1 + e) for (el, m), e in zip(BULK.items(), eps)}
        m_metal = M_METAL * (1 + eps[4])
        core, _, _ = split_core_mantle(
            bulk, m_metal, M_SIL * (1 + eps[5]), dIW=DIW * (1 + eps[6]), w_fixed=W_FIXED
        )
        out.append(1.0 - sum(W_FIXED.values()) - sum(core.values()) / m_metal)
    return out


if __name__ == '__main__':
    with_rule = irons()
    real = core_split._more_iron
    core_split._more_iron = lambda evaluate, start, fixed, out: (*out, 0)
    try:
        without = irons()
    finally:
        core_split._more_iron = real
    for name, fe in (('with the rule', with_rule), ('without the rule', without)):
        rich, poor = sum(f > 1e-2 for f in fe), sum(f < 1e-4 for f in fe)
        print(f'{name}: {rich} of {len(fe)} starts return iron > 1e-2, {poor} iron < 1e-4')

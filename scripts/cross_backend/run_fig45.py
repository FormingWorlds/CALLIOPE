"""Driver that runs Fig 4 (attribution) and Fig 5 (Earth anchor) in
the same Python process so they share the JAX warmup cost.

Calling each fig script independently means a cold JAX compile per
process (~10 min on a Mac). Both figs only need 2-5 atmodeller calls,
so amortising the warmup across both is a 2x speedup.
"""

from __future__ import annotations

import logging

from . import fig4_attribution, fig5_earth_anchor


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')
    print('=== Fig 5 (Earth anchor) ===')
    out5 = fig5_earth_anchor.make_figure()
    for ext, path in out5.items():
        print(f'  {ext}: {path}')

    print('=== Fig 4 (attribution) ===')
    out4 = fig4_attribution.make_figure()
    for ext, path in out4.items():
        print(f'  {ext}: {path}')


if __name__ == '__main__':
    main()

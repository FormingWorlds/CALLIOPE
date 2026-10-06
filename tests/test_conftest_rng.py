"""The autouse fixture seeds the global NumPy RNG with 0 before every test."""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.parametrize('run', [1, 2])
def test_each_test_starts_from_the_seeded_state(run):
    """Every test, also the second of two, starts from the seed-0 state and its draws."""
    seed0 = np.random.RandomState(0)
    np.testing.assert_array_equal(np.random.get_state()[1], seed0.get_state()[1])
    np.testing.assert_array_equal(np.random.random_sample(5), seed0.random_sample(5))

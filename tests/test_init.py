from __future__ import annotations

import pytest

from calliope import __version__

pytestmark = pytest.mark.unit


def test_version():
    assert __version__

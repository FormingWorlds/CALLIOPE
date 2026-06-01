from __future__ import annotations

import pytest

from calliope import __version__

pytestmark = pytest.mark.unit


def test_version():
    """`calliope.__version__` must be a non-empty PEP 440 version string."""
    assert __version__
    # Discrimination guard: a stub that returned a non-string truthy value
    # (e.g. True or an int) would pass the bare assert. The version must
    # be a string with at least one dot (so 'x.y' or 'x.y.z' shape).
    assert isinstance(__version__, str)
    assert '.' in __version__

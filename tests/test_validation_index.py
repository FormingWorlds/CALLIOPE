"""Tests for the case table of `docs/Validation/index.md`.

Contract: every test that carries `reference_pinned` (on the function, its class or its
module) appears in the table, and every test id in the table names such a test. A
shorthand id `::name` replaces the last part of the id before it in the same cell.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / 'docs' / 'Validation' / 'index.md'

# Carries the marker of its module but checks only an environment setting.
NOT_A_CASE = {'tests/test_cross_code_consistency.py::test_required_env_guards'}


def _pinned(nodes) -> bool:
    """True when any expression in ``nodes`` names ``reference_pinned``."""
    return any(
        isinstance(n, ast.Attribute) and n.attr == 'reference_pinned'
        for node in nodes
        for n in ast.walk(node)
    )


def _pytestmark(body) -> list[ast.expr]:
    """Values assigned to ``pytestmark`` in a module or class body."""
    return [
        s.value
        for s in body
        if isinstance(s, ast.Assign)
        and any(getattr(t, 'id', '') == 'pytestmark' for t in s.targets)
    ]


def pinned_tests(tests_dir: Path) -> set[str]:
    """Return ``file::[Class::]function`` of every reference-pinned test under ``tests_dir``."""
    found = set()
    for path in sorted(tests_dir.glob('test_*.py')):
        tree = ast.parse(path.read_text())
        rel = path.relative_to(tests_dir.parent).as_posix()
        mod = _pinned(_pytestmark(tree.body))
        funcs = [(f'{rel}::', mod, n) for n in tree.body]
        for c in tree.body:
            if isinstance(c, ast.ClassDef) and c.name.startswith('Test'):
                cls = mod or _pinned(c.decorator_list) or _pinned(_pytestmark(c.body))
                funcs += [(f'{rel}::{c.name}::', cls, n) for n in c.body]
        for prefix, inherited, n in funcs:
            if isinstance(n, ast.FunctionDef | ast.AsyncFunctionDef) and n.name.startswith(
                'test'
            ):
                if inherited or _pinned(n.decorator_list):
                    found.add(prefix + n.name)
    return found


def index_ids(text: str) -> set[str]:
    """Return the test ids of the case table, with each ``::name`` shorthand expanded.

    Raises
    ------
    ValueError
        A shorthand id has no full id before it in its cell.
    """
    ids = set()
    for line in text.splitlines():
        cells = line.split('|')
        if not line.startswith('| [') or len(cells) < 5:
            continue
        last = None
        for token in re.findall(r'`([^`]+)`', cells[3]):
            if token.startswith('::'):
                if last is None:
                    raise ValueError(f'shorthand {token} without a full id: {line}')
                token = last.rsplit('::', 1)[0] + token
            last = token
            ids.add(token)
    return ids


def test_index_lists_every_reference_pinned_test():
    """The table names each reference-pinned test once in some row, and nothing else."""
    pinned = pinned_tests(ROOT / 'tests')
    listed = index_ids(INDEX.read_text())
    # 28 pinned test functions exist; an empty scan would pass the set checks vacuously
    assert len(pinned) >= 28
    assert NOT_A_CASE <= pinned, 'stale NOT_A_CASE entry'
    assert sorted(pinned - NOT_A_CASE - listed) == [], 'pinned tests missing from the index'
    assert sorted(listed - pinned) == [], 'index ids that name no reference-pinned test'


def test_scan_follows_function_class_and_module_marks(tmp_path):
    """A mark on the function, on the class or in pytestmark counts; an unmarked test does not."""
    tests = tmp_path / 'tests'
    tests.mkdir()
    (tests / 'test_a.py').write_text(
        'import pytest\n'
        '@pytest.mark.reference_pinned\n'
        'def test_f(): pass\n'
        'def test_plain(): pass\n'
        '@pytest.mark.reference_pinned\n'
        'class TestC:\n'
        '    def test_m(self): pass\n'
        'class TestD:\n'
        '    pytestmark = pytest.mark.reference_pinned\n'
        '    def test_n(self): pass\n'
        '    def helper(self): pass\n'
    )
    (tests / 'test_b.py').write_text(
        'import pytest\npytestmark = [pytest.mark.unit, pytest.mark.reference_pinned]\n'
        'def test_g(): pass\n'
    )
    (tests / 'test_c.py').write_text('import pytest\npytestmark = [pytest.mark.unit]\n')
    found = pinned_tests(tests)
    assert found == {
        'tests/test_a.py::test_f',
        'tests/test_a.py::TestC::test_m',
        'tests/test_a.py::TestD::test_n',
        'tests/test_b.py::test_g',
    }
    assert 'tests/test_a.py::test_plain' not in found


def test_index_parser_expands_shorthand_and_refuses_a_bare_one():
    """`::name` takes the class of the id before it; a cell that starts with one raises."""
    row = '| [P](p.md) | case | `tests/t.py::TestA::test_x`, `::test_y`, `tests/u.py::test_z` | ref | unit |'
    assert index_ids(row) == {
        'tests/t.py::TestA::test_x',
        'tests/t.py::TestA::test_y',
        'tests/u.py::test_z',
    }
    # header, separator and prose lines carry no ids
    assert (
        index_ids('| Process | Case | Test id |\n|---|---|---|\nSee `tests/t.py::test_q`.')
        == set()
    )
    with pytest.raises(ValueError, match='without a full id'):
        index_ids('| [P](p.md) | case | `::test_y` | ref | unit |')

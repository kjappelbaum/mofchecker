# -*- coding: utf-8 -*-
"""Test composiion checks"""
from mofchecker.checks.global_structure import HasCarbon, HasHydrogen


def test_no_c(get_no_c):
    """Check we can correctly flag a case without C"""
    for structure in get_no_c:
        has_c = HasCarbon(structure)
        assert not has_c.is_ok


def test_no_h(get_no_h):
    """Check we correctly flag a case without H"""
    for structure in get_no_h:
        has_h = HasHydrogen(structure)
        assert not has_h.is_ok

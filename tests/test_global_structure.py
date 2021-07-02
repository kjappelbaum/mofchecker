# -*- coding: utf-8 -*-
from mofchecker.checks.global_structure import HasCarbon, HasHydrogen


def test_no_c(get_no_c):
    for structure in get_no_c:
        has_c = HasCarbon(structure)
        assert not has_c.is_ok


def test_no_h(get_no_h):
    for structure in get_no_h:
        has_h = HasHydrogen(structure)
        assert not has_h.is_ok

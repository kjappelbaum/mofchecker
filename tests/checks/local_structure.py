# -*- coding: utf-8 -*-
from mofchecker.checks.local_structure.overlapping_atoms import AtomicOverlapCheck


def test_clashing(get_clashing_structures):
    for structure in get_clashing_structures:
        overlap_check = AtomicOverlapCheck(structure)
        assert not overlap_check.is_ok

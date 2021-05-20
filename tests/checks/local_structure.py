# -*- coding: utf-8 -*-
from pymatgen.core import Structure

from mofchecker.checks.local_structure.overlapping_atoms import AtomicOverlapCheck
from mofchecker.checks.local_structure.undercoordinated_carbon import (
    UnderCoordinatedCarbonCheck,
)


def test_clashing(get_clashing_structures):
    for structure in get_clashing_structures:
        overlap_check = AtomicOverlapCheck(structure)
        assert not overlap_check.is_ok


def check_undercoordinated_carbon():
    ...

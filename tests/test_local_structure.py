# -*- coding: utf-8 -*-
"""Test the checks of the local coordination environment"""
import os

from pymatgen.core import Structure

from mofchecker.checks.local_structure.false_oxo import FalseOxoCheck
from mofchecker.checks.local_structure.overlapping_atoms import AtomicOverlapCheck
from mofchecker.checks.local_structure.undercoordinated_rare_earth import (
    UnderCoordinatedRareEarthCheck,
)
from mofchecker.graph import get_structure_graph

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def test_clashing(get_clashing_structures):
    """Test the check for atomic overlaps"""
    for structure in get_clashing_structures:
        overlap_check = AtomicOverlapCheck(structure)
        assert not overlap_check.is_ok


def test_false_oxo():
    """Test the check Andrew recommended for "false" terminal oxo"""
    structure = Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "false_terminal_oxo_ca.cif")
    )
    structure_graph = get_structure_graph(structure, "vesta")

    checker = FalseOxoCheck(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 1


def test_undercoordinated_rare_metal():
    """Test our check for rare earth metals that likely miss some
    neighbors"""
    structure = Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "GADRAH_Ce_clean.cif")
    )
    structure_graph = get_structure_graph(structure, "vesta")

    checker = UnderCoordinatedRareEarthCheck(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 4

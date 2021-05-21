# -*- coding: utf-8 -*-
import os

from pymatgen.core import Structure

from mofchecker.checks.local_structure.false_oxo import FalseOxoCheck
from mofchecker.checks.local_structure.overlapping_atoms import AtomicOverlapCheck
from mofchecker.graph import get_structure_graph

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def test_clashing(get_clashing_structures):
    for structure in get_clashing_structures:
        overlap_check = AtomicOverlapCheck(structure)
        assert not overlap_check.is_ok


def test_false_oxo():
    structure = Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "false_terminal_oxo_ca.cif")
    )
    structure_graph = get_structure_graph(structure, "vesta")

    checker = FalseOxoCheck(structure, structure_graph)
    print(checker.is_ok)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 1

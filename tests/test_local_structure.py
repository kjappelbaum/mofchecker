# -*- coding: utf-8 -*-
"""Testing the checks on the local chemical environment."""
import os

from pymatgen.core import Structure
from structuregraph_helpers.create import get_structure_graph

from mofchecker.checks.local_structure.false_oxo import FalseOxoCheck
from mofchecker.checks.local_structure.geometrically_exposed_metal import GeometricallyExposedMetal
from mofchecker.checks.local_structure.overlapping_atoms import AtomicOverlapCheck
from mofchecker.checks.local_structure.undercoordinated_alkaline import (
    UnderCoordinatedAlkaliAlkaline,
)
from mofchecker.checks.local_structure.undercoordinated_rare_earth import (
    UnderCoordinatedRareEarthCheck,
)

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def test_clashing(get_clashing_structures):
    """Testing the check for overlapping atoms."""
    for structure in get_clashing_structures:
        overlap_check = AtomicOverlapCheck(structure)
        assert not overlap_check.is_ok


def test_false_oxo():
    """Testing the check for suspicious oxo groups."""
    structure = Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "false_terminal_oxo_ca.cif")
    )
    structure_graph = get_structure_graph(structure, "vesta")

    checker = FalseOxoCheck(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 1


def test_undercoordinated_rare_earth_check():
    """Testing the check for undercoordinated rare earth metals."""
    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "GADRAH_Ce_clean.cif"))
    structure_graph = get_structure_graph(structure, "vesta")

    checker = UnderCoordinatedRareEarthCheck(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 4


def test_undercoordinated_alkali_alkaline_check():
    """Testing the check for undercoordinated alkali/alkaline earth metals."""
    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "MOTMAK_clean.cif"))
    structure_graph = get_structure_graph(structure, "vesta")

    checker = UnderCoordinatedAlkaliAlkaline(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 2


def test_geometrically_exposed_metal():
    """Testing the check for geometrically exposed metals."""
    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    structure_graph = get_structure_graph(structure, "vesta")

    checker = GeometricallyExposedMetal(structure, structure_graph)
    assert checker.is_ok
    assert len(checker.flagged_indices) == 0

    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "MOTMAK_clean.cif"))
    structure_graph = get_structure_graph(structure, "vesta")

    checker = GeometricallyExposedMetal(structure, structure_graph)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 2

    structure = Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "046_flu+N131+126_charge.cif")
    )
    structure_graph = get_structure_graph(structure, "vesta")

    checker = GeometricallyExposedMetal(structure, structure_graph, tight=True)
    assert not checker.is_ok
    assert len(checker.flagged_indices) == 8

    checker = GeometricallyExposedMetal(structure, structure_graph, tight=False)
    assert checker.is_ok

# -*- coding: utf-8 -*-
"""Test the symmetry subpackage"""
import os

from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation

from mofchecker import MOFChecker
from mofchecker.symmetry import get_symmetry_hash

from .conftest import THIS_DIR


def test_get_symmetry_hash():
    """Test the symmetry hashes"""
    zif3 = MOFChecker.from_cif(os.path.join(THIS_DIR, "test_files", "ZIF-3.cif"))
    print(type(zif3.structure))
    hash_a = get_symmetry_hash(zif3.structure, tight=True)
    assert len(hash_a) == 275
    hash_b = get_symmetry_hash(zif3.structure)
    assert len(hash_b) == 47
    assert hash_a[-3:] == hash_b[-3:]

    zif4 = MOFChecker.from_cif(os.path.join(THIS_DIR, "test_files", "ZIF-4.cif"))
    hash_zif4_a = get_symmetry_hash(zif4.structure, tight=True)
    hash_zif4_b = get_symmetry_hash(zif4.structure)

    assert hash_zif4_a != hash_a
    assert hash_zif4_b != hash_b

    assert zif4.symmetry_hash == hash_zif4_b

    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    original_hash = get_symmetry_hash(MOFChecker(structure).structure)

    # rotate structure
    rotation_transformer = RotationTransformation([1, 0, 0], 10)
    rotated_structure = rotation_transformer.apply_transformation(structure)
    assert get_symmetry_hash(MOFChecker(rotated_structure).structure) == original_hash

    # create supercell
    structure.make_supercell([1, 2, 1])
    print(type(structure))
    assert get_symmetry_hash(MOFChecker(structure).structure) == original_hash

# -*- coding: utf-8 -*-
"""Testing the hash functions."""
import os

import pytest

from mofchecker import MOFChecker
from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation

from .conftest import THIS_DIR


def test_graph_hash():
    """Basic check that the function call works"""
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    )
    assert isinstance(mofchecker.graph_hash, str)


def test_graph_hash_robustness():  # pylint: disable=too-many-locals
    """Check that duplicating or rotating the structure produces the same hash."""
    structure = Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    original_hash = MOFChecker(structure).graph_hash

    # rotate structure
    rotation_transformer = RotationTransformation([1, 0, 0], 10)
    rotated_structure = rotation_transformer.apply_transformation(structure)
    assert MOFChecker(rotated_structure).graph_hash == original_hash

    # create supercell
    structure.make_supercell([1, 2, 1])
    assert MOFChecker(structure).graph_hash == original_hash

    # check the MOF-74 structures
    mohgoi_checker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "MOHGOI.cif"))
    )

    todyuj_checker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "TODYUJ.cif"))
    )

    vogtiv_checker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "VOGTIV.cif"))
    )

    # There water on TODYUJ
    assert mohgoi_checker.graph_hash != todyuj_checker.graph_hash
    # one is the supercell of the other
    assert mohgoi_checker.graph_hash == vogtiv_checker.graph_hash

    # MOF-74-Zr.cif
    mof_74_zr = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "MOF-74-Zr.cif"))
    )
    assert mof_74_zr.graph_hash != todyuj_checker.graph_hash
    assert mof_74_zr.graph_hash != vogtiv_checker.graph_hash

    # # MOF-74-Zn.cif
    mof_74_zn = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "MOF-74-Zn.cif"))
    )
    assert mof_74_zr.scaffold_hash == mof_74_zn.scaffold_hash
    assert mof_74_zr.graph_hash != mof_74_zn.graph_hash

    # # MOF-5 is not ZIF-8
    mof_5 = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "mof-5_cellopt.cif"))
    )
    zif_8 = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ZIF-8-RASPA.cif"))
    )
    assert mof_5.graph_hash != zif_8.graph_hash

    # # Mn-MOF-74 and UiO-67
    coknun = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "coknun01.cif"))
    )
    wizmac = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "WIZMAV02_auto.cif"))
    )
    assert coknun.graph_hash != wizmac.graph_hash


@pytest.mark.past_issue
def test_graph_hash_false_positives():
    """Test the hash on past issues of false positives.

    Cases where different structures gave the same hash.
    """
    # issue 130
    cof_18141N2 = MOFChecker.from_cif(  # pylint: disable=invalid-name
        os.path.join(THIS_DIR, "test_files", "18141N2.cif")
    )
    cof_20211N2 = MOFChecker.from_cif(  # pylint: disable=invalid-name
        os.path.join(THIS_DIR, "test_files", "20211N2.cif")
    )
    assert cof_18141N2.graph_hash != cof_20211N2.graph_hash

    # # Daniele's report
    mmpf7 = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "943643.cif"))
    )
    mmpf8 = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "943644.cif"))
    )

    assert mmpf7.graph_hash != mmpf8.graph_hash
    assert mmpf7.scaffold_hash != mmpf8.scaffold_hash

    # ZIF3/4
    # zif3 = MOFChecker.from_cif(os.path.join(THIS_DIR, "test_files", "ZIF-3.cif"))
    # zif4 = MOFChecker.from_cif(os.path.join(THIS_DIR, "test_files", "ZIF-4.cif"))
    # assert zif3.graph_hash != zif4.graph_hash


@pytest.mark.past_issue
def test_graph_hash_false_negatives():
    """Test the hash on past issues of false negatives.

    Cases of structures that should match but the graph hash did not.
    """
    # issue 107: Mn-MOF-74: ASR ORIWET and COKNUN… give different hash
    oriwet = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ORIWET.cif"))
    )

    coknun = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "COKNUN.cif"))
    )

    assert oriwet.graph_hash == coknun.graph_hash
    assert oriwet.scaffold_hash == coknun.scaffold_hash

# -*- coding: utf-8 -*-
"""Tests of mofchecker"""
# pylint: disable=missing-function-docstring,singleton-comparison,invalid-name
import os

from pymatgen import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation

from mofchecker import MOFChecker

from .conftest import THIS_DIR


def test_has_oms(get_cn4_structre, get_cn5_paddlewheel_structure):

    omsdetector = MOFChecker(get_cn4_structre)
    assert omsdetector.has_oms == True

    omsdetector = MOFChecker(get_cn5_paddlewheel_structure)
    assert omsdetector.has_oms == False


def test_has_oms_multiple(get_testdict):
    for k, v in get_testdict.items():
        print(k)
        omsdetector = MOFChecker.from_cif(k)
        assert omsdetector.has_oms == v


def test_clashing(get_clashing_structures):
    for structure in get_clashing_structures:
        mofchecker = MOFChecker(structure)
        assert mofchecker.has_atomic_overlaps == True


def test_name():
    s = os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif")
    mofchecker = MOFChecker.from_cif(s)
    assert mofchecker.name == "ABAVIJ_clean"


def test_no_h(get_no_h):
    for structure in get_no_h:
        mofchecker = MOFChecker(structure)
        assert mofchecker.has_hydrogen == False


def test_no_c(get_no_c):
    for structure in get_no_c:
        mofchecker = MOFChecker(structure)
        assert mofchecker.has_carbon == False


def test_overvalent_c(get_overvalent_c_structures):
    for structure in get_overvalent_c_structures:
        mofchecker = MOFChecker(structure)
        assert mofchecker.has_overvalent_c == True

    mofchecker = MOFChecker.from_cif(
        os.path.join(THIS_DIR, "test_files", "XIGFOJ_manual.cif")
    )
    assert mofchecker.has_overvalent_c == False

    # alkine ligand
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "RUDQUD_clean.cif"))
    )
    assert mofchecker.has_overvalent_c == False

    # based on issue https://github.com/kjappelbaum/mofchecker/issues/63
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AGARUW_clean.cif"))
    )
    assert mofchecker.has_overvalent_c == False


def test_lone_atom():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif"))
    )
    assert mofchecker.has_lone_atom == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "HKUST_floating.cif"))
    )
    assert mofchecker.has_lone_atom == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "FEZTIP_clean.cif"))
    )
    assert mofchecker.has_lone_atom == False


def test_lone_molecule():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif"))
    )
    assert mofchecker.has_lone_molecule == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "HKUST_floating.cif"))
    )
    assert mofchecker.has_lone_molecule == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "UiO_66_water.cif"))
    )
    assert mofchecker.has_lone_molecule == True


def test_undercoordinated_c():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_c == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AHOKIR_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_c == True

    # alkine ligand
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "RUDQUD_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_c == False


def test_undercoordinated_n():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "VUGYED_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "mil-53-al-nh2.cif"))
    )
    assert mofchecker.has_undercoordinated_n == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "1246903.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False

    mofchecker = MOFChecker(
        Structure.from_file(
            os.path.join(THIS_DIR, "test_files", "1246903_missing_H.cif")
        )
    )
    assert mofchecker.has_undercoordinated_n == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABOVOF_FSR.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False


def test_chargecheck():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AMUFIZ_clean.cif"))
    )
    assert mofchecker.has_high_charges == False

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
    assert mohgoi_checker.graph_hash == vogtiv_checker.graph_hash


def test_graph_hash():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAXUZ.cif"))
    )
    assert isinstance(mofchecker.graph_hash, str)


def test_graph_hash_robustness():
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


def test_dicts():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "VUGYED_clean.cif"))
    )
    assert isinstance(mofchecker.check_descriptions, dict)
    assert isinstance(mofchecker.check_expected_values, dict)

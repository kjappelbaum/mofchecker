# -*- coding: utf-8 -*-
"""Tests of mofchecker"""
# pylint: disable=missing-function-docstring,singleton-comparison,invalid-name
import os

import pytest
from pymatgen import Structure

from mofchecker import MOFChecker

from .conftest import THIS_DIR


def test_partial_occupancy():
    with pytest.raises(NotImplementedError):
        MOFChecker.from_cif(os.path.join(THIS_DIR, "test_files", "ABUBIK.cif"))


def test_name():
    s = os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif")
    mofchecker = MOFChecker.from_cif(s)
    assert mofchecker.name == "ABAVIJ_clean"


def test_unknown_elements():
    """Parsing structure with unknown element raises warning for covalent radius."""
    with pytest.warns(UserWarning) as record:
        mofchecker = MOFChecker.from_cif(
            os.path.join(THIS_DIR, "test_files", "GUPQOA.cif")
        )
        mofchecker.get_mof_descriptors()
    assert len(record) >= 1


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

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "N_MOF_ASR.cif"))
    )
    assert mofchecker.has_undercoordinated_n == True


def test_chargecheck():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AMUFIZ_clean.cif"))
    )
    assert mofchecker.has_high_charges == False


def test_dicts():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "VUGYED_clean.cif"))
    )
    assert isinstance(mofchecker.check_descriptions, dict)
    assert isinstance(mofchecker.check_expected_values, dict)


def test_is_porous(get_cn5_paddlewheel_structure):
    mc = MOFChecker(get_cn5_paddlewheel_structure)
    assert mc.is_porous == True

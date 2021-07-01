# -*- coding: utf-8 -*-
"""Tests of mofchecker"""
# pylint: disable=missing-function-docstring,singleton-comparison,invalid-name
import os

import pytest
from ase.io import read
from pymatgen.core import Structure

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
    # with pytest.warns(UserWarning) as record:
    #     mofchecker = MOFCheƒ√cker.from_cif(
    #         os.path.join(THIS_DIR, "test_files", "GUPQOA.cif")
    #     )
    #     mofchecker.get_mof_descriptors()
    # assert len(record) >= 1

    with pytest.raises(NotImplementedError):
        mofchecker = MOFChecker.from_cif(
            os.path.join(THIS_DIR, "test_files", "GUPQOA.cif")
        )
        mofchecker.get_mof_descriptors()


def test_overvalent_h():
    mofchecker = MOFChecker.from_cif(
        os.path.join(THIS_DIR, "test_files", "overvalent_h.cif")
    )
    assert mofchecker.has_overvalent_h
    assert len(mofchecker.overvalent_h_indices) == 3

    mofchecker = MOFChecker.from_cif(
        os.path.join(THIS_DIR, "test_files", "XIGFOJ_manual.cif")
    )
    assert not mofchecker.has_overvalent_h


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


@pytest.mark.past_issue
def test_overvalent_c_past_issue():
    # based on issue https://github.com/kjappelbaum/mofchecker/issues/63
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AGARUW_clean.cif"))
    )
    assert mofchecker.has_overvalent_c == False


def test_lone_molecule():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif"))
    )
    assert mofchecker.has_lone_molecule == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "HKUST_floating.cif"))
    )
    assert mofchecker.has_lone_molecule == True

    assert mofchecker.lone_molecule_indices == [[144]]

    atoms = read(os.path.join(THIS_DIR, "test_files", "overvalent_h.cif"))
    mofchecker = MOFChecker.from_ase(atoms, primitive=False)
    assert len(mofchecker.lone_molecule_indices) == 3


@pytest.mark.past_issue
def test_lone_molecule_past_issue():
    atoms = read(os.path.join(THIS_DIR, "test_files", "floating_check.cif"))
    mofchecker = MOFChecker.from_ase(atoms)
    assert len(mofchecker.lone_molecule_indices) == 1
    assert len(mofchecker.lone_molecule_indices[0]) == 5
    species = []
    for ind in mofchecker.lone_molecule_indices[0]:
        species.append(str(mofchecker.structure[ind].specie))
    assert set(species) == set(["H", "H", "H", "H", "C"])
    assert mofchecker.has_lone_molecule == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "OTOXIF_clean.cif"))
    )
    assert mofchecker.has_lone_molecule == False


def test_undercoordinated_c():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_c == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AHOKIR_clean.cif"))
    )
    assert mofchecker.has_undercoordinated_c == True

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "missing_h_on_c.cif"))
    )
    assert mofchecker.has_undercoordinated_c == True
    assert len(mofchecker.undercoordinated_c_indices) == 2
    assert len(mofchecker.undercoordinated_c_candidate_positions) == 2


@pytest.mark.past_issue
def slow_check_undercoordinated_c_past_issue():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "DALVIY_clean.cif"))
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


@pytest.mark.past_issue
def test_undercoordinated_n_past_issue():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "axipee.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "ABOVOF_FSR.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False

    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "abiqae.cif"))
    )
    assert mofchecker.has_undercoordinated_n == False


def test_chargecheck():
    mofchecker = MOFChecker(
        Structure.from_file(os.path.join(THIS_DIR, "test_files", "AMUFIZ_clean.cif"))
    )
    assert mofchecker.has_high_charges == False


def test_is_porous(get_cn5_paddlewheel_structure):
    mc = MOFChecker(get_cn5_paddlewheel_structure)
    assert mc.is_porous == True

# -*- coding: utf-8 -*-
"""pytest fixtures"""
# pylint: disable=missing-function-docstring
import os

import pytest
from pymatgen.core import Structure

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture(scope="module")
def get_cn5_paddlewheel_structure():
    return Structure.from_file(
        os.path.join(THIS_DIR, "test_files", "paddlewheel_cn5.cif")
    )


@pytest.fixture(scope="module")
def get_cn4_structure():
    return Structure.from_file(os.path.join(THIS_DIR, "test_files", "cn4.cif"))


@pytest.fixture(scope="module")
def get_testdict():

    dct = {
        str(os.path.join(THIS_DIR, "test_files", "ABEXEM_clean.cif")): True,
        str(os.path.join(THIS_DIR, "test_files", "ABEXIQ_clean.cif")): True,
        str(os.path.join(THIS_DIR, "test_files", "ZUSNOS_clean.cif")): True,
        # ToDo: this is a trickier one (?)
        # str(os.path.join(THIS_DIR, 'test_files', 'ac403674p_si_001_clean.cif')):
        # False,
        str(os.path.join(THIS_DIR, "test_files", "ZUQBUK_clean.cif")): False,
        str(os.path.join(THIS_DIR, "test_files", "ABAVIJ_clean.cif")): False,
        str(os.path.join(THIS_DIR, "test_files", "FAPXIG_clean.cif")): False,
        # str(os.path.join(THIS_DIR, 'test_files', 'ADABIS_clean.cif')):
        # False,  # https://onlinelibrary.wiley.com/doi/epdf/10.1002/anie.201202992
        # -> I feel one there is actually True
        str(os.path.join(THIS_DIR, "test_files", "ALUJOH_clean.cif")): True,
        str(
            os.path.join(THIS_DIR, "test_files", "AMUFIZ_clean.cif")
        ): True,  # CN=2 error is caught here
        str(os.path.join(THIS_DIR, "test_files", "DEJCIF_clean.cif")): True,
        str(
            os.path.join(THIS_DIR, "test_files", "DAWWEF_clean.cif")
        ): True,  # CN=4 see-saw
        str(os.path.join(THIS_DIR, "test_files", "UFEXOT_clean.cif")): True,
        str(os.path.join(THIS_DIR, "test_files", "REHHIX_clean.cif")): True,
        str(os.path.join(THIS_DIR, "test_files", "ZADDAJ_clean.cif")): True,
        # str(os.path.join(THIS_DIR, 'test_files', 'ELIYUU_clean.cif')): False
        str(os.path.join(THIS_DIR, "test_files", "ZUWXUM_clean.cif")): False,
        str(
            os.path.join(THIS_DIR, "test_files", "VUGYED_clean.cif")
        ): True,  # wrong in CoRE
        str(
            os.path.join(THIS_DIR, "test_files", "TONTIB_clean.cif")
        ): False,  # wrong in CoRE
        str(
            os.path.join(THIS_DIR, "test_files", "ELUQIM13_clean.cif")
        ): False,  # wrong in CoRE
    }

    return dct


@pytest.fixture(scope="module")
def get_no_h():
    structures = []

    names = ["ZADDAJ_clean.cif"]

    for name in names:
        structures.append(
            Structure.from_file(os.path.join(THIS_DIR, "test_files", name))
        )

    return structures


@pytest.fixture(scope="module")
def get_no_c():
    structures = []

    names = ["VOXVEL_clean.cif"]

    for name in names:
        structures.append(
            Structure.from_file(os.path.join(THIS_DIR, "test_files", name))
        )

    return structures


@pytest.fixture(scope="module")
def get_clashing_structures():
    structures = []

    names = ["RUYGEZ_clean.cif", "PADHIM_clean.cif"]

    for name in names:
        structures.append(
            Structure.from_file(os.path.join(THIS_DIR, "test_files", name))
        )

    return structures


@pytest.fixture(scope="module")
def get_overvalent_c_structures():
    structures = []
    names = ["ZUQBUK_clean.cif"]

    for name in names:
        structures.append(
            Structure.from_file(os.path.join(THIS_DIR, "test_files", name))
        )

    return structures

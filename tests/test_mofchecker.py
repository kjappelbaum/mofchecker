# -*- coding: utf-8 -*-
import os

import pytest

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
    s = os.path.join(THIS_DIR, 'test_files', 'ABAVIJ_clean.cif')
    mofchecker = MOFChecker.from_cif(s)
    assert mofchecker.name == 'ABAVIJ_clean'


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
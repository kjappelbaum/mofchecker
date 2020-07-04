# -*- coding: utf-8 -*-
import pytest

from mofchecker import MOFChecker


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

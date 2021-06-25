# -*- coding: utf-8 -*-
"""Testing the code for the OMS detection"""
import pytest

from mofchecker import MOFChecker


@pytest.mark.skip
def test_has_oms(get_cn4_structure, get_cn5_paddlewheel_structure):
    """Test two specific structures"""
    omsdetector = MOFChecker(get_cn4_structure)
    assert omsdetector.has_oms

    omsdetector = MOFChecker(get_cn5_paddlewheel_structure)
    assert not omsdetector.has_oms


@pytest.mark.skip
def test_has_oms_multiple(get_testdict):
    """Test on a bunch of different structures"""
    for key, value in get_testdict.items():
        omsdetector = MOFChecker.from_cif(key)
        assert omsdetector.has_oms == value

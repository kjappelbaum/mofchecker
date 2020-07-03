# -*- coding: utf-8 -*-
import pytest

from omsdetector import OMSDetector


def test_has_oms(get_cn4_structre, get_cn5_paddlewheel_structure):

    omsdetector = OMSDetector(get_cn4_structre)
    assert omsdetector.has_oms == True

    omsdetector = OMSDetector(get_cn5_paddlewheel_structure)
    assert omsdetector.has_oms == False


def test_has_oms_multiple(get_testdict):
    for k, v in get_testdict.items():
        print(k)
        omsdetector = OMSDetector.from_cif(k)
        assert omsdetector.has_oms == v

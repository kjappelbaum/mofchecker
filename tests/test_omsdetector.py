# -*- coding: utf-8 -*-
import pytest

from omsdetector import OMSDetector


def test_has_oms(get_cn4_structre, get_cn5_paddlewheel_structure):

    omsdetector = OMSDetector(get_cn4_structre)
    assert omsdetector.has_oms == True

    omsdetector = OMSDetector(get_cn5_paddlewheel_structure)
    assert omsdetector.has_oms == False

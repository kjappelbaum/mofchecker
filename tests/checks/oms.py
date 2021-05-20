# -*- coding: utf-8 -*-
def test_has_oms(get_cn4_structure, get_cn5_paddlewheel_structure):

    omsdetector = MOFChecker(get_cn4_structure)
    assert omsdetector.has_oms == True

    omsdetector = MOFChecker(get_cn5_paddlewheel_structure)
    assert omsdetector.has_oms == False


def test_has_oms_multiple(get_testdict):
    for k, v in get_testdict.items():
        omsdetector = MOFChecker.from_cif(k)
        assert omsdetector.has_oms == v

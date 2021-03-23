# -*- coding: utf-8 -*-
"""Testing the zeo++ module"""
import pytest

from mofchecker.zeopp import check_if_porous, parse_zeopp, run_zeopp

TEST_LINE = "output_file.res    1.70107 0.95106  1.64805"


def test_parse_zeopp():
    """Simple parsing of the pore output"""
    res = parse_zeopp(TEST_LINE)
    assert res == {
        "lis": 1.70107,  # largest included sphere
        "lifs": 0.95106,  # largest free sphere
        "lifsp": 1.64805,  # largest included sphere along free sphere path
    }


def test_run_zeopp(get_cn5_paddlewheel_structure):
    """Running the full analysis starting from a pmg structure"""
    structure = get_cn5_paddlewheel_structure
    res = run_zeopp(structure)
    assert res["lis"] == pytest.approx(7.65505)
    assert res["lifs"] == pytest.approx(5.81104)
    assert res["lifsp"] == pytest.approx(7.63730)


def test_check_if_porous(get_cn5_paddlewheel_structure):
    """Run the porosity check"""
    structure = get_cn5_paddlewheel_structure
    assert check_if_porous(structure)

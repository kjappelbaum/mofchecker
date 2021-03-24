# -*- coding: utf-8 -*-
"""Testing the zeo++ module"""
import pytest

from mofchecker.checks.zeopp import _parse_zeopp, check_if_porous, run_zeopp

TEST_LINE = "output_file.res    1.70107 0.95106  1.64805"


def test__parse_zeopp():
    """Simple parsing of the pore output"""
    res = _parse_zeopp(TEST_LINE)
    assert res == {
        "lis": 1.70107,  # largest included sphere
        "lifs": 0.95106,  # largest free sphere
        "lifsp": 1.64805,  # largest included sphere along free sphere path
    }


def test_run_zeopp(get_cn5_paddlewheel_structure):
    """Running the full analysis starting from a pmg structure"""
    structure = get_cn5_paddlewheel_structure
    res = run_zeopp(structure)
    assert res["lis"] == pytest.approx(7.65505, abs=0.02)
    assert res["lifs"] == pytest.approx(5.81104, abs=0.02)
    assert res["lifsp"] == pytest.approx(7.63730, abs=0.02)


def test_check_if_porous(get_cn5_paddlewheel_structure):
    """Run the porosity check"""
    structure = get_cn5_paddlewheel_structure
    assert check_if_porous(structure)

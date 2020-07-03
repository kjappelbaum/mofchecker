# -*- coding: utf-8 -*-
import os

import pytest
from pymatgen import Structure

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture(scope='module')
def get_cn5_paddlewheel_structure():
    return Structure.from_file(
        os.path.join(THIS_DIR, 'test_files', 'paddlewheel_cn5.cif'))


@pytest.fixture(scope='module')
def get_cn4_structre():
    return Structure.from_file(os.path.join(THIS_DIR, 'test_files', 'cn4.cif'))

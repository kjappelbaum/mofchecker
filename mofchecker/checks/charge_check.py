# -*- coding: utf-8 -*-
from tempfile import NamedTemporaryFile

import numpy as np
from pyeqeq.main import run_on_cif

from .check_base import AbstractCheck


class ChargeCheck(AbstractCheck):
    def __init__(self, structure):
        self.structure = structure
        self.threshold = 3

    @property
    def name(self):
        return "High charges"

    @property
    def description(self):
        return f"Check that the charges of the structure are reasonable (abs(charge) not higher than {self.threshold})."

    def _run_check(self):
        with NamedTemporaryFile("w", suffix=".cif") as file:
            self.structure.to("cif", file.name)
            charges = run_on_cif(file.name)
            has_high_charges = np.sum(np.abs(charges) > self.threshold)

        return not has_high_charges

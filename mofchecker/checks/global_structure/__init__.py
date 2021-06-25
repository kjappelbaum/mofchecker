# -*- coding: utf-8 -*-
"""Checks on the composition"""
from ..check_base import AbstractCheck
from ..utils.get_indices import (
    get_c_indices,
    get_h_indices,
    get_metal_indices,
    get_n_indices,
)


class HasCarbon(AbstractCheck):
    "Checks if the structure has any carbon atom."

    def __init__(self, structure):
        self.structure = structure

    def _run_check(self):
        c_indices = get_c_indices(self.structure)
        return len(c_indices) > 0

    @property
    def description(self):
        return "Checks if the structure has any carbon atom."


class HasNitrogen(AbstractCheck):
    "Checks if the structure has any nitrogen atom."

    def __init__(self, structure):
        self.structure = structure

    def _run_check(self):
        n_indices = get_n_indices(self.structure)
        return len(n_indices) > 0

    @property
    def description(self):
        return "Checks if the structure has any nitrogen atom."


class HasHydrogen(AbstractCheck):
    "Checks if the structure has any hydrogen atom."

    def __init__(self, structure):
        self.structure = structure

    def _run_check(self):
        h_indices = get_h_indices(self.structure)
        return len(h_indices) > 0

    @property
    def description(self):
        return "Checks if the structure has any hydrogen atom."


class HasMetal(AbstractCheck):
    "Checks if the structure has any metal atom."

    def __init__(self, structure):
        self.structure = structure

    def _run_check(self):
        metal_indices = get_metal_indices(self.structure)
        return len(metal_indices) > 0

    @property
    def description(self):
        return "Checks if the structure has any metal atom."

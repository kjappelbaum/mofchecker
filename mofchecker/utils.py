# -*- coding: utf-8 -*-
"""Helper functions for the MOFChecker"""
import json
import pickle
from functools import cached_property

import pymatgen
from pymatgen.core.structure import Structure


def read_pickle(file):
    """Read a pickle file"""
    with open(file, "rb") as handle:
        return pickle.load(handle)


def _check_metal_coordination(site, coordination_number: int) -> bool:
    # Lanthanides like to have many neighbors
    # Low coordinatio number is usually only
    # possible with really bulky ligands
    # for this reason, a Lanthanide with low
    # coordination number, e.g., <= 4 can be considered "interesting"
    if (
        (site.specie.is_lanthanoid)
        or (site.specie.is_actinoid)
        or (site.specie.symbol in ("Mo", "Cr", "Hf", "Mb"))
    ):
        if coordination_number <= 4:
            return True

    # Also for the alkaline/alkaline earth metals,
    # I would find a low coordination number surprising
    # elif site.specie.is_alkali or site.specie.is_alkaline:
    #     if coordination_number <= 4:
    #         return True

    return False


def print_dict(dictionary):
    """Print a dictionary to stdout line by line."""
    for k, v in dictionary.items():  # pylint: disable=invalid-name
        print(k, v)


def _check_if_ordered(structure):
    if not structure.is_ordered:
        raise NotImplementedError(
            "Support of unordered structures with partial occupancies \
                is not implemented (yet)."
        )
    for site in structure:
        if site.specie.atomic_radius is None:
            raise NotImplementedError(
                f"Pymatgen currently does not support this {str(site.specie)} element"
            )


class IStructure(pymatgen.core.structure.IStructure):
    """pymatgen IStructure with faster equality comparison.

    This dramatically speeds up lookups in the LRU cache when an object
    with the same __hash__ is already in the cache.
    """

    __hash__ = pymatgen.core.structure.IStructure.__hash__

    def __eq__(self, other):
        """Use specific, yet performant hash for equality comparison."""
        return self._dict_hash == other._dict_hash

    @cached_property
    def _dict_hash(self):
        """Specific, yet performant hash."""
        return hash(json.dumps(self.as_dict(), sort_keys=True))

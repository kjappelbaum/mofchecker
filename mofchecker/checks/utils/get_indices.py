# -*- coding: utf-8 -*-
import functools
from typing import Union

import pymatgen
from pymatgen.core.structure import IStructure, Structure

from ...definitions import METALS


def is_metal(site: pymatgen.core.Site) -> bool:
    """according to conquest help:
    transition metal, lanthanide, actinide,
    or Al, Ga, In, Tl, Ge, Sn, Pb, Sb, Bi, Po"""
    if str(site.specie) in METALS:
        return True
    return False


@functools.lru_cache(maxsize=2, typed=False)
def _get_indices(immutable_structure: IStructure) -> dict:
    return {
        "c": _get_c_indices(immutable_structure),
        "h": _get_h_indices(immutable_structure),
        "n": _get_n_indices(immutable_structure),
        "metal": _get_metal_indices(immutable_structure),
    }


def _get_c_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "C"]


def _get_h_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "H"]


def _get_n_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "N"]


def _get_metal_indices(structure):
    return [i for i, site in enumerate(structure) if is_metal(site)]


def get_h_indices(structure):
    return get_indices(structure)["h"]


def get_c_indices(structure):
    return get_indices(structure)["c"]


def get_n_indices(structure):
    return get_indices(structure)["n"]


def get_metal_indices(structure):
    return get_indices(structure)["metal"]


def get_indices(structure: Union[Structure, IStructure]):
    if isinstance(structure, Structure):
        structure = IStructure.from_sites(structure)
    return _get_indices(structure)

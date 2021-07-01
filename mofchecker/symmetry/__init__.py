# -*- coding: utf-8 -*-
"""Analyze symmetrized structures and hash them"""
from typing import Union

from pymatgen.core import IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.structure import SymmetrizedStructure

from .hash import hash_symmetrized_structure


def get_symmetrized_structure(
    structure: Union[Structure, IStructure]
) -> SymmetrizedStructure:
    """Constructs a SymmetrizedStructure, i.e. a structure
    where the spacegroup and symmetry operations are defined."""
    return SpacegroupAnalyzer(structure).get_symmetrized_structure()


def get_symmetry_hash(
    structure: Union[Structure, IStructure, SymmetrizedStructure], tight: bool = False
) -> str:
    """Hashes the symmetrical positions of the SymmetrizedStructure.
    The tight setting also considers the ordering,
    otherwise only the number and identity
    of the elements is considered.

    Args:
        structure (Union[Structure, IStructure, SymmetrizedStructure]):
            A structure for which the symmetry hash is calculated
        tight (bool, optional): If True, also consider the ordering. Defaults to False.

    Returns:
        str: hash
    """
    if not isinstance(structure, SymmetrizedStructure):
        structure = get_symmetrized_structure(structure)
    return hash_symmetrized_structure(structure, tight)

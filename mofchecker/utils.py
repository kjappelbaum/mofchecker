# -*- coding: utf-8 -*-
import numpy as np
from pymatgen import Structure
from scipy import sparse

from .definitions import COVALENT_RADII


class LowCoordinationNumber(KeyError):
    pass


class HighCoordinationNumber(KeyError):
    pass


class NoOpenDefined(KeyError):
    pass


class NoMetal(KeyError):
    pass


def compute_overlap_matrix(distance_matrix: np.array,
                           allatomtypes: list,
                           tolerance: float = 1.0):
    """
    Find atomic overlap based on pairwise distance and Covalent radii.
    Criterion: if dist < min (CovR_1,CovR_2) -> overlap (this function is used in molsimplify)
    """
    overlap_matrix = np.zeros(distance_matrix.shape)
    for i, e1 in enumerate(allatomtypes[:-1]):
        for j, e2 in enumerate(allatomtypes[i + 1:]):
            dist = distance_matrix[i, i + j + 1]
            # check for atomic overlap:
            if dist < tolerance * min(COVALENT_RADII[e1], COVALENT_RADII[e2]):
                overlap_matrix[i, i + j + 1] = 1
                overlap_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(overlap_matrix)


def get_overlaps(s: Structure) -> list:
    distance_matrix = s.distance_matrix
    atomtypes = [str(species) for species in s.species]
    overlap_matrix = compute_overlap_matrix(distance_matrix, atomtypes)
    overlap_atoms = []
    for at in set(sparse.find(overlap_matrix)[0]):
        overlap_atoms.append(at.item())
    return overlap_atoms

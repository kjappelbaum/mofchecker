"""Helpers for analysing the geometry of a structure."""
import numpy as np
from element_coder.encode import encode_many
from libconeangle import cone_angle
from numpy.linalg import matrix_rank
from pymatgen.analysis.graphs import StructureGraph


def are_coplanar(coords) -> bool:
    """Check if the given coordinates are coplanar."""
    coords = np.unique(np.array(coords), axis=0)
    coords -= coords.mean(axis=0)
    if matrix_rank(coords, tol=0.1) <= 2:
        return True
    return False


def _get_coords_and_elements_of_neighbors(graph, index):
    neighbors = graph.get_connected_sites(index)
    coords = []
    species = []
    coords.append(graph.structure.cart_coords[index])
    species.append(str(graph.structure[index].specie))
    for neighbor in neighbors:
        coords.append(neighbor.site.coords)
        species.append(str(neighbor.site.specie))
    return coords, species


def get_open_angle(graph: StructureGraph, index: int) -> float:
    """Get 360 - cone angle of the site with the given index."""
    coords, species = _get_coords_and_elements_of_neighbors(graph, index)
    encodings = encode_many(species, "van_der_waals_radius")
    print(are_coplanar(coords))
    try:
        angle, _, _ = cone_angle(coords, encodings, 0)
        return 360 - angle
    except ValueError:
        coords = np.unique(np.array(coords), axis=0)
        coords -= coords.mean(axis=0)
        if matrix_rank(coords) <= 2:
            return 180
        return np.nan


def has_open_angle(graph: StructureGraph, index: int, threshold: float = 80) -> bool:
    """Check if the site with the given index has an open angle.
    Nans are treated as False.
    """
    angle = get_open_angle(graph, index)
    if angle > threshold:
        return True
    return False

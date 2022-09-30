import numpy as np
from element_coder.encoders import encode_many
from libconeangle import cone_angle
from numpy.linalg import matrix_rank


def are_coplanar(coords):
    coords = np.unique(np.array(coords), axis=0)
    coords -= coords.mean(axis=0)
    if matrix_rank(coords, tol=.1) <= 2:
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

def get_open_angle(graph, index):
    coords, species = _get_coords_and_elements_of_neighbors(graph, index)
    encodings = encode_many(species, 'van_der_waals_radius')
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


def has_open_angle(graph, index, threshold=80):
    angle = get_open_angle(graph, index)
    if angle > threshold:
        return True
    return False
# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
from pymatgen import Structure
from pymatgen.core import Molecule
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


def print_dict(dictionary):
    for k, v in sorted(dictionary.items()):
        print(k, v)






def get_subgraphs_as_molecules_all(sg, use_weights=False):
    """Copied from http://pymatgen.org/_modules/pymatgen/analysis/graphs.html#StructureGraph.get_subgraphs_as_molecules and removed the duplicate check

    Args:
        sg ([type]): [description]
        use_weights (bool, optional): [description]. Defaults to False.

    Returns:
        [type]: [description]
    """

    # creating a supercell is an easy way to extract
    # molecules (and not, e.g., layers of a 2D crystal)
    # without adding extra logic

    supercell_sg = sg * (3, 3, 3)

    # make undirected to find connected subgraphs
    supercell_sg.graph = nx.Graph(supercell_sg.graph)

    # find subgraphs
    all_subgraphs = list(nx.connected_component_subgraphs(supercell_sg.graph))

    # discount subgraphs that lie across *supercell* boundaries
    # these will subgraphs representing crystals
    molecule_subgraphs = []
    for subgraph in all_subgraphs:
        intersects_boundary = any(
            [d['to_jimage'] != (0, 0, 0) for u, v, d in subgraph.edges(data=True)]
        )
        if not intersects_boundary:
            molecule_subgraphs.append(subgraph)

    # add specie names to graph to be able to test for isomorphism
    for subgraph in molecule_subgraphs:
        for n in subgraph:
            subgraph.add_node(n, specie=str(supercell_sg.structure[n].specie))

    # get Molecule objects for each subgraph
    molecules = []
    for subgraph in molecule_subgraphs:

        coords = [supercell_sg.structure[n].coords for n in subgraph.nodes()]
        species = [supercell_sg.structure[n].specie for n in subgraph.nodes()]

        molecule = Molecule(species, coords)

        molecules.append(molecule)

    return molecules

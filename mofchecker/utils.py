# -*- coding: utf-8 -*-
"""Helper functions for the MOFChecker"""
import networkx as nx
import numpy as np
import pymatgen
from pymatgen import Structure
from pymatgen.core import Molecule
from scipy import sparse

from .definitions import COVALENT_RADII


class LowCoordinationNumber(KeyError):
    """Error for low coordination number"""


class HighCoordinationNumber(KeyError):
    """Error for high coordination number"""


class NoOpenDefined(KeyError):
    """Error in case the open check is not defined
    for this coordination numberF"""


class NoMetal(KeyError):
    """Error in case there is no metal in structure"""


def compute_overlap_matrix(
    distance_matrix: np.array, allatomtypes: list, tolerance: float = 1.0
):
    """
    Find atomic overlap based on pairwise distance and Covalent radii.
    Criterion: if dist < min (CovR_1,CovR_2) -> overlap (this function is used in molsimplify)
    """
    overlap_matrix = np.zeros(distance_matrix.shape)
    for i, elem_1 in enumerate(allatomtypes[:-1]):
        for j, elem_2 in enumerate(allatomtypes[i + 1 :]):
            dist = distance_matrix[i, i + j + 1]
            # check for atomic overlap:
            if dist < tolerance * min(COVALENT_RADII[elem_1], COVALENT_RADII[elem_2]):
                overlap_matrix[i, i + j + 1] = 1
                overlap_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(overlap_matrix)


def get_overlaps(s: Structure) -> list:
    """Find overlapping atoms in a structure."""
    distance_matrix = s.distance_matrix
    atomtypes = [str(species) for species in s.species]
    overlap_matrix = compute_overlap_matrix(distance_matrix, atomtypes)
    overlap_atoms = []
    for atom in set(sparse.find(overlap_matrix)[0]):
        overlap_atoms.append(atom.item())
    return overlap_atoms


def print_dict(dictionary):
    """Print a dictionary to stdout line by line."""
    for k, v in sorted(dictionary.items()):
        print(k, v)


def get_subgraphs_as_molecules_all(
    structure_graph: pymatgen.analysis.graphs.StructureGraph,
):
    """Copied from
    http://pymatgen.org/_modules/pymatgen/analysis/graphs.html#StructureGraph.get_subgraphs_as_molecules
    and removed the duplicate check

    Args:
        structure_graph ( pymatgen.analysis.graphs.StructureGraph): Structuregraph

    Returns:
        List: list of molecules
    """

    # creating a supercell is an easy way to extract
    # molecules (and not, e.g., layers of a 2D crystal)
    # without adding extra logic
    supercell_sg = structure_graph * (3, 3, 3)

    # make undirected to find connected subgraphs
    supercell_sg.graph = nx.Graph(supercell_sg.graph)

    # find subgraphs
    all_subgraphs = [
        supercell_sg.graph.subgraph(c)
        for c in nx.connected_components(supercell_sg.graph)
    ]

    # discount subgraphs that lie across *supercell* boundaries
    # these will subgraphs representing crystals
    molecule_subgraphs = []
    for subgraph in all_subgraphs:
        intersects_boundary = any(
            [d["to_jimage"] != (0, 0, 0) for u, v, d in subgraph.edges(data=True)]
        )
        if not intersects_boundary:
            molecule_subgraphs.append(nx.MultiDiGraph(subgraph))

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

        # shift so origin is at center of mass
        molecule = molecule.get_centered_molecule()

        molecules.append(molecule)

    return molecules

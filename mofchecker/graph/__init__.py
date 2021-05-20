# -*- coding: utf-8 -*-
"""Helper functions for the graph hash calculation"""
import networkx as nx
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import (
    BrunnerNN_relative,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    JmolNN,
    MinimumDistanceNN,
    VoronoiNN,
)
from pymatgen.core import Structure, Molecule

VESTA_NN = CutOffDictNN.from_preset("vesta_2019")


def construct_clean_graph(
    structure: Structure, structure_graph: StructureGraph
) -> nx.Graph:
    """Creates a networkx graph with atom numbers as node labels"""
    edges = {(u, v) for u, v, d in structure_graph.graph.edges(keys=False, data=True)}
    graph = nx.Graph()
    graph.add_edges_from(edges)

    for node in graph.nodes:

        graph.nodes[node]["specie"] = str(structure[node].specie)

    return graph


def _get_cn(structure_graph, site_index):
    return len(structure_graph.get_connected_sites(site_index))


def get_local_env_method(method):
    method = method.lower()

    if method.lower() == "crystalnn":
        return CrystalNN()
    elif method.lower() == "econnn":
        return EconNN()
    elif method.lower() == "brunnernn":
        return BrunnerNN_relative()
    elif method.lower() == "minimumdistance":
        return MinimumDistanceNN()
    elif method.lower() == "vesta":
        return VESTA_NN
    elif method.lower() == "voronoinn":
        return VoronoiNN()
    else:
        return JmolNN()


def get_structure_graph(structure, method):
    return StructureGraph.with_local_env_strategy(
        structure, get_local_env_method(method)
    )



def get_subgraphs_as_molecules_all(
    structure_graph: StructureGraph,
):
    """Copied from
    http://pymatgen.org/_modules/pymatgen/analysis/graphs.html#StructureGraph.get_subgraphs_as_molecules
    and removed the duplicate check

    Args:
        structure_graph (StructureGraph): Structuregraph

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
            (d["to_jimage"] != (0, 0, 0) for u, v, d in subgraph.edges(data=True))
        )
        if not intersects_boundary:
            molecule_subgraphs.append(nx.MultiDiGraph(subgraph))

    # add specie names to graph to be able to test for isomorphism
    for subgraph in molecule_subgraphs:
        for node in subgraph:
            subgraph.add_node(node, specie=str(supercell_sg.structure[node].specie))

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
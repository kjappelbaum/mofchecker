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
from pymatgen.core import Structure

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

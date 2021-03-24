# -*- coding: utf-8 -*-
"""Helper functions for the graph hash calculation"""
import networkx as nx
from pymatgen import Structure
from pymatgen.analysis.graphs import StructureGraph


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

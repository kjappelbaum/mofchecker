# -*- coding: utf-8 -*-
"""Check that attempts to find connected components in the cell that do not cross PBC"""
import networkx as nx

from ..graph import get_subgraphs_as_molecules
from .check_base import AbstractIndexCheck


class FloatingSolventCheck(AbstractIndexCheck):
    """Checks if there is any non-periodic connected component in the cell,
    which can be a solvent, charge compensating counter ion, etc."""

    def __init__(self, structure_graph):
        self.structure_graph = structure_graph
        nx.set_node_attributes(
            self.structure_graph.graph,
            name="idx",
            values=dict(zip(range(len(structure_graph)), range(len(structure_graph)))),
        )

    @property
    def name(self):
        return "Floating atom or molecule"

    @classmethod
    def from_mofchecker(cls, mofchecker):
        """Initialize a checker from a mofchecker instance"""
        checker = cls(mofchecker.graph)
        return checker

    def _run_check(self):
        _, _, idx, _, _ = get_subgraphs_as_molecules(
            self.structure_graph, return_unique=False
        )
        return len(idx) == 0, idx

    @property
    def description(self):
        return "Checks if there is any non-periodic connected component in the cell,\
             which can be a solvent, charge compensating counter ion, etc."

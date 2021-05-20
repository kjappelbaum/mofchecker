from .check_base import AbstractIndexCheck
from ..graph import get_subgraphs_as_molecules

class FloatingSolventCheck(AbstractIndexCheck):
    def __init__(self, structure_graph):
        self.structure_graph = structure_graph

    @classmethod
    def from_mofchecker(cls, mofchecker):
        checker = cls(mofchecker.graph)
        return checker

    def _run_check(self):
        _, _, idx, _, _ = get_subgraphs_as_molecules(self.structure_graph)
        return len(idx) == 0, idx
        
    @property
    def description(self):
        return "Checks if there is any non-periodic connected component in the cell, which can be a solvent, charge compensating counter ion, etc."
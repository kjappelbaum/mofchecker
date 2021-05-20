# -*- coding: utf-8 -*-
from ..utils.get_indices import _is_any_neighbor_metal, get_h_indices
from .base_coordination_check import BaseCoordinationCheck


class OverCoordinatedHydrogenCheck(BaseCoordinationCheck):
    def __init__(self, structure, structure_graph):
        self.structure = structure
        self.h_indices = get_h_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def description(self):
        return "Checks, using geometric heuristics, if there are any hydrogen that are likely overcoordinated (i.e., CN>1)."

    def _run_check(self):
        overcoordinated_hydrogens = self._get_overcoordinated_hydrogens()
        return len(overcoordinated_hydrogens) == 0, overcoordinated_hydrogens

    def _get_overcoordinated_hydrogens(self):
        """Check for all H if CN>1, ignore metal bonds."""
        overcoordinated_hydrogens = []

        for site_index in self.h_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 1:
                if not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                    overcoordinated_hydrogens.append(site_index)
        return overcoordinated_hydrogens

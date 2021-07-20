# -*- coding: utf-8 -*-
"""Checks, using geometric heuristics if there are any carbons
that are likely overcoordinated (i.e., CN>4)"""
from ..utils.get_indices import _is_any_neighbor_metal, get_n_indices
from .base_coordination_check import BaseCoordinationCheck


class OverCoordinatedNitrogenCheck(BaseCoordinationCheck):
    """Checks, using geometric heuristics if there are any carbons
    that are likely overcoordinated (i.e., CN>4)"""

    def __init__(
        self, structure, structure_graph
    ):  # pylint: disable=super-init-not-called
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        return "Overcoordinated nitrogen"

    @property
    def description(self):
        return "Checks, using geometric heuristics,\
             if there are any carbons that are likely overcoordinated (i.e., CN>4)."

    def _run_check(self):
        overcoordinated_nitrogen = self._get_overcoordinated_nitrogen()
        return len(overcoordinated_nitrogen) == 0, overcoordinated_nitrogen

    def _get_overcoordinated_nitrogen(self):
        """Check for all N if CN>4, ignore metal bonds."""
        overcoordinated_nitrogen = []

        for site_index in self.n_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 4:
                if not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                    overcoordinated_nitrogen.append(site_index)

        return overcoordinated_nitrogen

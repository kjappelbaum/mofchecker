# -*- coding: utf-8 -*-
import numpy as np

from ...graph import _get_cn
from ..utils.get_indices import _vdw_radius_neighbors, get_c_indices, is_metal
from .base_coordination_check import BaseCoordinationCheck
from .geometry import _maximum_angle


class UnderCoordinatedCarbonCheck(BaseCoordinationCheck):
    def __init__(self, structure, structure_graph):
        self.structure = structure
        self.c_indices = get_c_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def description(self):
        return "Checks, using geometric heuristics, if there are any carbons that are likely undercoordinated."

    def _run_check(self):
        undercoordinated_carbons = self._get_undercoordinated_carbons()
        return len(undercoordinated_carbons) == 0, undercoordinated_carbons

    def _get_undercoordinated_carbons(self, tolerance: int = 10):
        """Idea is that carbon should at least have three neighbors if it is not sp1.
        In sp1 case it is linear. So we can just check if there are carbons with
        non-linear coordination with less than three neighbors. An example in CoRE
        MOF would be AHOKIR. In principle this should also flag the quite common
        case of benzene rings with missing hydrogens.
        """
        undercoordinated_carbons = []

        for site_index in self.c_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn == 2:
                neighbors = self.get_connected_sites(site_index)
                angle = _maximum_angle(
                    self.structure.get_angle(
                        site_index, neighbors[0].index, neighbors[1].index
                    )
                )
                if (np.abs(180 - angle) > tolerance) or (np.abs(180 - 0) > tolerance):
                    if (not is_metal(neighbors[0].site)) or (
                        not is_metal(neighbors[1].site)
                    ):
                        if len(_vdw_radius_neighbors(self.structure, site_index)) <= 2:
                            undercoordinated_carbons.append(site_index)

        return undercoordinated_carbons

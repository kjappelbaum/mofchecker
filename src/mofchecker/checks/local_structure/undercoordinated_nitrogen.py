# -*- coding: utf-8 -*-
"""Check for undercoordinated nitrogens"""
from ..utils.get_indices import get_n_indices
from .base_missing_check import BaseMissingCheck
from .geometry import (
    _guess_underbound_nitrogen_cn2,
    _guess_underbound_nitrogen_cn3,
    add_sp2_hydrogen,
    add_sp3_hydrogen,
    add_sp_hydrogen,
)


class UnderCoordinatedNitrogenCheck(BaseMissingCheck):
    """Check for undercoordinated nitrogens"""

    def __init__(
        self, structure, structure_graph
    ):  # pylint: disable=super-init-not-called
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        return "Undercoordinated nitrogen"

    @property
    def description(self):
        return "Checks, using geometric heuristics,\
             if there are any nitrogens that are likely undercoordinated."

    def _run_check(self):
        undercoordinated_nitrogens, positions = self._get_undercoordinated_nitrogens()
        return (
            len(undercoordinated_nitrogens) == 0,
            undercoordinated_nitrogens,
            positions,
        )

    def _get_undercoordinated_nitrogens(self, tolerance: int = 25):
        """Attempts to captures missing hydrogens on nitrogen groups
        using heuristics
        """
        undercoordinated_nitrogens = []
        h_positions = []
        for site_index in self.n_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            neighbors = self.get_connected_sites(site_index)
            if cn == 1:
                # this is suspicous, but it also might a CN which is perfectly fine.
                # to check this, we first see if the neighbor is carbon
                # and then what its coordination number is. If it is greater than 2
                # then we likely do not have a CN for which the carbon should be a
                # linear sp one
                if (self.get_cn(neighbors[0].index) > 2) and not neighbors[
                    0
                ].site.specie.is_metal:
                    undercoordinated_nitrogens.append(site_index)
                    h_positions.append(
                        add_sp_hydrogen(self.structure[site_index], neighbors)
                    )
            elif cn == 2:
                undercoordinated_nitrogen = _guess_underbound_nitrogen_cn2(
                    self.structure,
                    site_index,
                    neighbors,
                    self.get_connected_sites(neighbors[0].index),
                    self.get_connected_sites(neighbors[1].index),
                    tolerance,
                )
                if undercoordinated_nitrogen:
                    undercoordinated_nitrogens.append(site_index)
                    h_positions.append(
                        add_sp2_hydrogen(self.structure[site_index], neighbors)
                    )
            elif cn == 3:
                undercoordinated_nitrogen = _guess_underbound_nitrogen_cn3(
                    self.structure, site_index, neighbors, tolerance
                )
                if undercoordinated_nitrogen:
                    undercoordinated_nitrogens.append(site_index)
                    h_positions.append(
                        add_sp3_hydrogen(self.structure[site_index], neighbors)
                    )
        return undercoordinated_nitrogens, h_positions

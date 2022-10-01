# -*- coding: utf-8 -*-
""""Check if there are any metals
that are sterically exposed"""
from .base_coordination_check import BaseCoordinationCheck
from ..utils.geometry import has_open_angle
from ..utils.get_indices import get_alkali_alkaline_indices, get_rare_earth_indices


class GeometricallyExposedMetal(BaseCoordinationCheck):
    """Check if there are any metals (alkali/alkaline earth or rare earth)
    that are likely geometrically exposed, i.e. which form
    a small cone angle with their binding partners"""

    def __init__(self, structure, structure_graph):  # pylint: disable=super-init-not-called
        self.structure = structure
        self.relevant_metals = get_alkali_alkaline_indices(structure) + get_rare_earth_indices(
            structure
        )
        self.structure_graph = structure_graph
        self.threshold = 80

    @property
    def name(self):
        return "Geometrically exposed metal."

    @property
    def description(self):
        return "Check if there are any metals (alkali/alkaline earth or rare earth) that are likely \
            geometrically exposed, i.e. which form a small cone angle \
                with their binding partners"

    def _run_check(self):
        exposed_metals = self._get_exposed_metals()
        return (
            len(exposed_metals) == 0,
            exposed_metals,
        )

    def _get_exposed_metals(self):
        """Check for all geometrically exposed metals."""
        geometrically_exposed_metals = []

        for site_index in self.relevant_metals:
            if has_open_angle(self.structure_graph, site_index, self.threshold):
                geometrically_exposed_metals.append(site_index)

        return geometrically_exposed_metals

# -*- coding: utf-8 -*-
import numpy as np
from pymatgen.core import Structure
from ..utils.get_indices import is_metal

def _maximum_angle(angle):
    diff_to_180 = np.abs(180 - angle)
    return max([angle, diff_to_180])




def _guess_underbound_nitrogen_cn3(
    structure: Structure, site_index: int, neighbors: list, tolerance: int = 10
) -> bool:
    """Check if there is a nitrogen with three neighbors
    that likely misses some coordination.

    Args:
        structure (Structure): pymatgen Structure object
        site_index (int): index of the central site that is check
        neighbors (list): list of neighboring sites
        tolerance (int, optional): Tolerance for angle checks in degree.
            Defaults to 10.

    Returns:
        bool: True if the nitrogen is likely missing some coordination partner
    """
    angle_a = structure.get_angle(site_index, neighbors[0].index, neighbors[1].index)
    angle_b = structure.get_angle(site_index, neighbors[0].index, neighbors[2].index)
    angle_c = structure.get_angle(site_index, neighbors[1].index, neighbors[2].index)
    min_angle = np.min([angle_a, angle_b, angle_c])

    any_metal = False
    for neighbor in neighbors:
        if is_metal(neighbor.site):
            any_metal = True

    num_h = 0
    for neighbor in neighbors:
        if str(neighbor.site.specie) == "H":
            num_h += 1

    if min_angle + tolerance < 110:
        # let's only do this if one of the neighbors is a metal.
        # sometimes the M-N bond is so long that it isn't correctly recognized
        # obviously, this now won't detect missing H on a floating NH3
        # but this is probably a rare situation
        if any_metal and (num_h == 2):
            return True

    return False


def _guess_underbound_nitrogen_cn2(  # pylint:disable=too-many-arguments
    structure: Structure,
    site_index: int,
    neighbors: list,
    connected_sites_a: list,
    connected_sites_b: list,
    tolerance: int = 10,
) -> bool:
    """Check if there is a nitrogen with CN 2 that probably misses
    some coordination.

    Args:
        structure (Structure): pymatgen Structure object
        site_index (int): Index of the site on which the check is performed
        neighbors (list): List of neighboring sites
        connected_sites_a (list): List of neighbor sites for first neighbor
        connected_sites_b (list): List of neighbor sites for second neighbor
        tolerance (int, optional): Tolerance for angle checks in degree.
             Defaults to 10.

    Returns:
        bool: True if there is a nitrogen that likely misses some coordination.
    """
    angle = structure.get_angle(site_index, neighbors[0].index, neighbors[1].index)
    neighbor_species = set(
        [str(neighbors[0].site.specie), str(neighbors[1].site.specie)]
    )

    if (np.abs(180 - angle) < tolerance) or (np.abs(0 - angle) < tolerance):
        # sp hybridization if the nitrogen is linear
        # this could be a nitride or a nitrosyl
        # usually, there is nothing to worry about if this is the case
        return False
    if angle < 115:
        # typically angle around 109.5 degree for sp3 hybridization
        # if we only have two neighbors but the nitrogen is likely
        # sp3 this is suspicious
        # to be sure we will check if it is planar (pyridine) or
        # not (piperazine) in the case the two neighbors are carbon
        # if neighbor_species == set(["C", "C"]):

        dihedral_a = structure.get_dihedral(
            site_index,
            neighbors[0].index,
            neighbors[1].index,
            connected_sites_a[0].index,
        )
        dihedral_b = structure.get_dihedral(
            site_index,
            neighbors[0].index,
            neighbors[1].index,
            connected_sites_b[0].index,
        )

        mean_dihedral = np.mean([dihedral_a, dihedral_b])
        if (np.abs(mean_dihedral - 180) < tolerance) or (
            np.abs(mean_dihedral - 0) < tolerance
        ):
            return False
        return True

    # larger angles should indicate sp2 hybridization
    # one case where MOFs might have an issue with sp2
    # is an NH2 group planar to the ring where one H is missing
    # the heuristic we use to catch this is if one of the neighbors
    # is H
    if "H" in neighbor_species:
        return True
    return False
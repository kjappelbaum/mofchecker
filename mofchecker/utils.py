# -*- coding: utf-8 -*-
"""Helper functions for the MOFChecker"""
import warnings

import networkx as nx
import numpy as np
import pymatgen
from pymatgen import Structure
from pymatgen.core import Molecule
from pymatgen.io.cif import CifWriter
from scipy import sparse

from .definitions import COVALENT_RADII


def _vdw_radius_neighbors(structure, site_index, tolerance: float = 1.5):
    radius = structure[site_index].specie.van_der_waals_radius
    return structure.get_neighbors(structure[site_index], tolerance * radius)


def _is_any_neighbor_metal(neighbors):
    for neighbor in neighbors:
        if neighbor.site.specie.is_metal:
            return True

    return False


def _check_metal_coordination(site, coordination_number: int) -> bool:
    # Lanthanides like to have many neighbors
    # Low coordinatio number is usually only
    # possible with really bulky ligands
    # for this reason, a Lanthanide with low
    # coordination number, e.g., <= 4 can be considered "interesting"
    if (
        (site.specie.is_lanthanoid)
        or (site.specie.is_actinoid)
        or (site.specie.symbol in ("Mo", "Cr", "Hf", "Mb"))
    ):
        if coordination_number <= 4:
            return True

    # Also for the alkaline/alkaline earth metals,
    # I would find a low coordination number surprising
    # elif site.specie.is_alkali or site.specie.is_alkaline:
    #     if coordination_number <= 4:
    #         return True

    return False


def _maximum_angle(angle):
    diff_to_180 = np.abs(180 - angle)
    return max([angle, diff_to_180])


def get_charges(structure: Structure):
    """Compute EqEq charges for a pymatgen structure"""
    try:
        from openbabel import pybel  # pylint:disable=import-outside-toplevel

        cif_structure = str(CifWriter(structure))
        mol = pybel.readstring("cif", cif_structure)
        mol.calccharges("eqeq")
        charges = [a.partialcharge for a in mol]
        return charges
    except ImportError:
        return None
    except Exception as execp:  # pylint:disable=broad-except
        warnings.warn(f"Exception occured during the charge calculation {execp}")
        return None


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
        if neighbor.site.specie.is_metal:
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
    print(angle)
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


class LowCoordinationNumber(KeyError):
    """Error for low coordination number"""


class HighCoordinationNumber(KeyError):
    """Error for high coordination number"""


class NoOpenDefined(KeyError):
    """Error in case the open check is not defined
    for this coordination numberF"""


class NoMetal(KeyError):
    """Error in case there is no metal in structure"""


def compute_overlap_matrix(
    distance_matrix: np.array, allatomtypes: list, tolerance: float = 1.0
):
    """
    Find atomic overlap based on pairwise distance and Covalent radii.

    Criterion: if dist < min (CovR_1,CovR_2) -> overlap
        (this function is used in molsimplify)
    """
    overlap_matrix = np.zeros(distance_matrix.shape)
    for i, elem_1 in enumerate(allatomtypes[:-1]):
        for j, elem_2 in enumerate(allatomtypes[i + 1 :]):
            dist = distance_matrix[i, i + j + 1]
            # check for atomic overlap:
            if dist < tolerance * min(COVALENT_RADII[elem_1], COVALENT_RADII[elem_2]):
                overlap_matrix[i, i + j + 1] = 1
                overlap_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(overlap_matrix)


def get_overlaps(s: Structure) -> list:  # pylint: disable=invalid-name
    """Find overlapping atoms in a structure."""
    distance_matrix = s.distance_matrix
    atomtypes = [str(species) for species in s.species]
    overlap_matrix = compute_overlap_matrix(distance_matrix, atomtypes)
    overlap_atoms = []
    for atom in set(sparse.find(overlap_matrix)[0]):
        overlap_atoms.append(atom.item())
    return overlap_atoms


def print_dict(dictionary):
    """Print a dictionary to stdout line by line."""
    for k, v in dictionary.items():  # pylint: disable=invalid-name
        print(k, v)


def get_subgraphs_as_molecules_all(
    structure_graph: pymatgen.analysis.graphs.StructureGraph,
):
    """Copied from
    http://pymatgen.org/_modules/pymatgen/analysis/graphs.html#StructureGraph.get_subgraphs_as_molecules
    and removed the duplicate check

    Args:
        structure_graph ( pymatgen.analysis.graphs.StructureGraph): Structuregraph

    Returns:
        List: list of molecules
    """

    # creating a supercell is an easy way to extract
    # molecules (and not, e.g., layers of a 2D crystal)
    # without adding extra logic
    supercell_sg = structure_graph * (3, 3, 3)

    # make undirected to find connected subgraphs
    supercell_sg.graph = nx.Graph(supercell_sg.graph)

    # find subgraphs
    all_subgraphs = [
        supercell_sg.graph.subgraph(c)
        for c in nx.connected_components(supercell_sg.graph)
    ]

    # discount subgraphs that lie across *supercell* boundaries
    # these will subgraphs representing crystals
    molecule_subgraphs = []
    for subgraph in all_subgraphs:
        intersects_boundary = any(
            (d["to_jimage"] != (0, 0, 0) for u, v, d in subgraph.edges(data=True))
        )
        if not intersects_boundary:
            molecule_subgraphs.append(nx.MultiDiGraph(subgraph))

    # add specie names to graph to be able to test for isomorphism
    for subgraph in molecule_subgraphs:
        for node in subgraph:
            subgraph.add_node(node, specie=str(supercell_sg.structure[node].specie))

    # get Molecule objects for each subgraph
    molecules = []
    for subgraph in molecule_subgraphs:
        coords = [supercell_sg.structure[n].coords for n in subgraph.nodes()]
        species = [supercell_sg.structure[n].specie for n in subgraph.nodes()]

        molecule = Molecule(species, coords)

        # shift so origin is at center of mass
        molecule = molecule.get_centered_molecule()

        molecules.append(molecule)

    return molecules

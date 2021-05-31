# -*- coding: utf-8 -*-
"""Helper functions for the graph hash calculation"""
import os
from collections import defaultdict
from typing import List, Tuple

import networkx as nx
import numpy as np
import yaml
from pymatgen.analysis.graphs import MoleculeGraph, StructureGraph
from pymatgen.analysis.local_env import (
    BrunnerNN_relative,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    JmolNN,
    MinimumDistanceNN,
    VoronoiNN,
)
from pymatgen.core import Molecule, Structure

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(THIS_DIR, "atom_typing_radii.yml"), "r") as handle:
    ATOM_TYPING_CUTOFFS = yaml.load(handle)


with open(os.path.join(THIS_DIR, "li_radii.yml"), "r") as handle:
    LI_TYPING_CUTOFFS = yaml.load(handle)

VESTA_NN = CutOffDictNN.from_preset("vesta_2019")
ATR_NN = CutOffDictNN(cut_off_dict=ATOM_TYPING_CUTOFFS)
LI_NN = CutOffDictNN(cut_off_dict=LI_TYPING_CUTOFFS)


def construct_clean_graph(
    structure: Structure, structure_graph: StructureGraph
) -> nx.Graph:
    """Creates a networkx graph with atom numbers as node labels"""
    edges = {(u, v) for u, v, d in structure_graph.graph.edges(keys=False, data=True)}
    graph = nx.Graph()
    graph.add_edges_from(edges)
    graph_descriptors = get_cycle_descriptors(graph)
    sg_centrality = nx.subgraph_centrality(graph)
    for node in graph.nodes:

        graph.nodes[node]["specie"] = str(structure[node].specie)
        graph.nodes[node]["specie-cn"] = (
            str(structure[node].specie)
            + "-"
            + str(structure_graph.get_coordination_of_site(node))
        )
        graph.nodes[node]["species-cycle"] = get_node_string(
            node, structure_graph, graph_descriptors["memberships"]
        )
        graph.nodes[node]["memership-hash"] = get_node_string(
            node,
            structure_graph,
            graph_descriptors["memberships"],
            membership_only=True,
        )
        graph.nodes[node]["species-sg-central"] = str(structure[node].specie) + str(
            sg_centrality[node]
        )
    return graph


def _get_cn(structure_graph, site_index):
    return len(structure_graph.get_connected_sites(site_index))


def get_local_env_method(method):
    method = method.lower()

    if method.lower() == "crystalnn":
        return CrystalNN(porous_adjustment=True, x_diff_weight=2)
    elif method.lower() == "econnn":
        return EconNN()
    elif method.lower() == "brunnernn":
        return BrunnerNN_relative()
    elif method.lower() == "minimumdistance":
        return MinimumDistanceNN()
    elif method.lower() == "vesta":
        return VESTA_NN
    elif method.lower() == "voronoinn":
        return VoronoiNN()
    elif method.lower() == "atr":
        return ATR_NN
    elif method.lower() == "li":
        return LI_NN
    else:
        return JmolNN()


def _is_in_cell(frac_coords):
    return all(frac_coords <= 1)


def _is_any_atom_in_cell(frac_coords):
    for row in frac_coords:
        if _is_in_cell(row):
            return True
    return False


def get_structure_graph(structure, method: str = "vesta"):
    return StructureGraph.with_local_env_strategy(
        structure, get_local_env_method(method)
    )


def _select_parts_in_cell(
    molecules: List[Molecule],
    graphs: List[MoleculeGraph],
    indices: List[List[int]],
    indices_here: List[List[int]],
    centers: List[np.ndarray],
    fractional_coordinates: np.ndarray,
    coordinates: np.ndarray,
) -> Tuple[List[Molecule], List[MoleculeGraph], List[List[int]]]:
    valid_indices = defaultdict(list)
    for i, ind in enumerate(indices_here):
        # change this check to having an atom in the cell
        frac_coords = fractional_coordinates[ind]

        if _is_any_atom_in_cell(frac_coords):
            sorted_idx = sorted(indices[i])
            valid_indices[str(sorted_idx)].append(i)

    molecules_ = []
    selected_indices = []
    graphs_ = []
    centers_ = []
    coordinates_ = []

    for _, v in valid_indices.items():
        for index in v:
            selected_indices.append(indices[index])
            molecules_.append(molecules[index])
            graphs_.append(graphs[index])
            centers_.append(centers[index])
            coordinates_.append(coordinates[index])

    return molecules_, graphs_, selected_indices, centers_, coordinates_


def get_subgraphs_as_molecules(
    structure_graph: StructureGraph,
    use_weights: bool = False,
    return_unique: bool = True,
    disable_boundary_crossing_check: bool = False,
    filter_in_cell: bool = True,
) -> Tuple[List[Molecule], List[MoleculeGraph], List[List[int]], List[np.ndarray]]:
    """Copied from
    http://pymatgen.org/_modules/pymatgen/analysis/graphs.html#StructureGraph.get_subgraphs_as_molecules
    and removed the duplicate check
    Args:
        structure_graph ( pymatgen.analysis.graphs.StructureGraph): Structuregraph
        use_weights (bool): If True, use weights for the edge matching
        return_unique (bool): If true, it only returns the unique molecules.
            If False, it will return all molecules that are completely included in the unit cell
            and fragments of the ones that are only partly in the cell
        filter_in_cell (bool): If True, it will only return molecules that have at least one atom
            in the cell
    Returns:
        Tuple[List[Molecule], List[MoleculeGraph], List[List[int]], List[np.ndarray]]
    """
    # creating a supercell is an easy way to extract
    # molecules (and not, e.g., layers of a 2D crystal)
    # without adding extra logic
    supercell_sg = structure_graph * (3, 3, 3)

    # make undirected to find connected subgraphs
    supercell_sg.graph = nx.Graph(supercell_sg.graph)

    # find subgraphs
    all_subgraphs = [
        supercell_sg.graph.subgraph(c).copy()
        for c in nx.connected_components(supercell_sg.graph)
    ]

    # discount subgraphs that lie across *supercell* boundaries
    # these will subgraphs representing crystals
    molecule_subgraphs = []

    for subgraph in all_subgraphs:
        if disable_boundary_crossing_check:
            molecule_subgraphs.append(nx.MultiDiGraph(subgraph))
        else:
            intersects_boundary = any(
                [d["to_jimage"] != (0, 0, 0) for u, v, d in subgraph.edges(data=True)]
            )
            if not intersects_boundary:
                molecule_subgraphs.append(nx.MultiDiGraph(subgraph))

    # add specie names to graph to be able to test for isomorphism
    for subgraph in molecule_subgraphs:
        for node in subgraph:
            subgraph.add_node(node, specie=str(supercell_sg.structure[node].specie))

    unique_subgraphs = []

    def node_match(n1, n2):
        return n1["specie"] == n2["specie"]

    def edge_match(e1, e2):
        if use_weights:
            return e1["weight"] == e2["weight"]
        else:
            return True

    if return_unique:
        for subgraph in molecule_subgraphs:
            already_present = [
                nx.is_isomorphic(
                    subgraph, g, node_match=node_match, edge_match=edge_match
                )
                for g in unique_subgraphs
            ]

            if not any(already_present):
                unique_subgraphs.append(subgraph)

    def make_mols(molecule_subgraphs=molecule_subgraphs, center=False):
        molecules = []
        indices = []
        indices_here = []
        mol_centers = []
        coordinates = []
        for subgraph in molecule_subgraphs:
            coords = [supercell_sg.structure[n].coords for n in subgraph.nodes()]
            species = [supercell_sg.structure[n].specie for n in subgraph.nodes()]
            # binding = [
            #     supercell_sg.structure[n].properties["binding"]
            #     for n in subgraph.nodes()
            # ]
            idx = [subgraph.nodes[n]["idx"] for n in subgraph.nodes()]
            idx_here = [n for n in subgraph.nodes()]
            molecule = Molecule(
                species, coords
            )  #  site_properties={"binding": binding}
            mol_centers.append(
                np.mean(supercell_sg.structure.cart_coords[idx_here], axis=0)
            )
            # shift so origin is at center of mass
            if center:
                molecule = molecule.get_centered_molecule()
            indices.append(idx)
            molecules.append(molecule)
            indices_here.append(idx_here)
            coordinates.append(coords)
        return molecules, indices, indices_here, mol_centers, coordinates

    def relabel_graph(multigraph):
        mapping = dict(zip(multigraph, range(0, len(multigraph.nodes()))))
        return nx.readwrite.json_graph.adjacency_data(
            nx.relabel_nodes(multigraph, mapping)
        )

    if return_unique:
        mol, idx, indices_here, centers, coordinates = make_mols(
            unique_subgraphs, center=True
        )
        return_subgraphs = unique_subgraphs
        return (
            mol,
            [
                MoleculeGraph(mol, relabel_graph(graph))
                for mol, graph in zip(mol, return_subgraphs)
            ],
            idx,
            centers,
            coordinates,
        )

    mol, idx, indices_here, centers, coordinates = make_mols(molecule_subgraphs)

    return_subgraphs = [
        MoleculeGraph(mol, relabel_graph(graph))
        for mol, graph in zip(mol, molecule_subgraphs)
    ]

    if filter_in_cell:
        mol, return_subgraphs, idx, centers, coordinates = _select_parts_in_cell(
            mol,
            return_subgraphs,
            idx,
            indices_here,
            centers,
            structure_graph.structure.lattice.get_fractional_coords(
                supercell_sg.structure.cart_coords
            ),
            coordinates,
        )

    return mol, return_subgraphs, idx, centers, coordinates


def get_cycle_lengths(graph):
    cycles = nx.minimum_cycle_basis(graph)
    lengths = count_sublist_lengths(cycles)
    return lengths


def count_sublist_lengths(list_of_lists):
    lengths = []

    for l in list_of_lists:
        lengths.append(len(l))

    return lengths


def get_cycle_memberships(graph, cycles, lengths):

    unique_lengths = np.unique(lengths)

    cycle_memberships = {}

    for node in graph.nodes():

        length_hash = dict(zip(unique_lengths, [0] * len(unique_lengths)))
        for cycle, length in zip(cycles, lengths):
            if node in cycle:
                length_hash[length] += 1

        cycle_memberships[node] = length_hash

    return cycle_memberships


def get_cycle_descriptors(graph):
    cycles = nx.minimum_cycle_basis(graph)

    lengths = count_sublist_lengths(cycles)

    memberships = get_cycle_memberships(graph, cycles, lengths)

    return {
        "cycles": cycles,
        "lengths": lengths,
        "memberships": memberships,
    }


def get_node_string(index, structure_graph, memberships, membership_only: bool = False):
    membership_hash = "".join([str(e) for e in memberships[index].values()])
    atom_label = str(structure_graph.structure[index].specie)
    if membership_only:
        return membership_hash
    return atom_label + membership_hash

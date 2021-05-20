# -*- coding: utf-8 -*-
"""MOFChecker: Basic sanity checks for MOFs"""
import logging
import os
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import List, Union

import networkx as nx
import numpy as np
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from pymatgen.analysis.graphs import ConnectedSite, StructureGraph
from pymatgen.analysis.local_env import (
    BrunnerNN_relative,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    JmolNN,
    MinimumDistanceNN,
    VoronoiNN,
)
from pymatgen.core import Structure
from pymatgen.core.structure import IStructure, Structure
from pymatgen.io.cif import CifParser

from ._version import get_versions
from .checks.global_structure import HasCarbon, HasHydrogen, HasMetal, HasNitrogen
from .checks.local_structure import AtomicOverlapCheck
from .checks.utils.get_indices import (
    get_c_indices,
    get_h_indices,
    get_metal_indices,
    get_n_indices,
)
from .checks.zeopp import check_if_porous
from .definitions import CHECK_DESCRIPTIONS, EXPECTED_CHECK_VALUES
from .graph import construct_clean_graph
from .utils import (
    _check_if_ordered,
    _guess_underbound_nitrogen_cn2,
    _guess_underbound_nitrogen_cn3,
    _is_any_neighbor_metal,
    _maximum_angle,
    _vdw_radius_neighbors,
    get_charges,
    get_subgraphs_as_molecules_all,
    is_metal,
)

__version__ = get_versions()["version"]
del get_versions

__all__ = ["__version__", "MOFChecker"]

MOFCheckLogger = logging.getLogger(__name__)
MOFCheckLogger.setLevel(logging.DEBUG)
VESTA_NN = CutOffDictNN.from_preset("vesta_2019")
try:
    from openbabel import pybel  # pylint:disable=import-outside-toplevel, unused-import

    HAS_OPENBABEL = True
except ImportError:
    warnings.warn(
        "For the charge check openbabel needs to be installed. \
    This can be done, for example using conda install openbabel"
    )
    HAS_OPENBABEL = False


class MOFChecker:  # pylint:disable=too-many-instance-attributes, too-many-public-methods
    """MOFChecker performs basic sanity checks for MOFs"""

    def __init__(self, structure: Structure, primitive: bool = True):
        """Class that can perform basic sanity checks for MOF structures

        Args:
            structure (Structure): pymatgen Structure object
            primitive (bool): If True, it will perform the analysis
                on the primitive structure

        Raises:
            NotImplementedError in the case of partial occupancies
        """
        self.structure = structure
        _check_if_ordered(structure)
        if primitive:
            self.structure = self.structure.get_primitive_structure()
        self.metal_indices = get_metal_indices(self.structure)

        self.porous_adjustment = False
        self.charges = None
        self._porous = ""
        self.metal_features = None
        self._cnn = None
        self._cnn_method = None
        self._filename = None
        self._name = None
        self.c_indices = get_c_indices(self.structure)
        self.h_indices = get_h_indices(self.structure)
        self.n_indices = get_n_indices(self.structure)
        self._overvalent_c = None
        self._overvalent_n = None
        self._overvalent_h = None
        self._undercoordinated_carbon = None
        self._undercoordinated_nitrogen = None
        self.check_expected_values = EXPECTED_CHECK_VALUES
        self.check_descriptions = CHECK_DESCRIPTIONS

        self._graph = None
        self._nx_graph = None
        self._connected_sites = {}
        self._cns = {}
        self._set_cnn()
        self._checks = {
            "has_c": HasCarbon(self.structure),
            "has_h": HasHydrogen(self.structure),
            "has_metal": HasMetal(self.structure),
            "has_nitrogen": HasNitrogen(self.structure),
            "has_atomic_overlaps": AtomicOverlapCheck(self.structure),
        }

    def _set_filename(self, path):
        self._filename = os.path.abspath(path)
        self._name = Path(path).stem

    def get_overlapping_indices(self):
        """Return the indices of overlapping atoms"""
        return self._checks["has_atomic_overlaps"].flagged_indices

    @property
    def graph_hash(self):
        """Return the Weisfeiler-Lehman graph hash.
        Hashes are identical for isomorphic graphs
        (taking the atomic kinds into account)
        and there are guarantees that non-isomorphic graphs will get different hashes.
        """
        return weisfeiler_lehman_graph_hash(self.nx_graph, node_attr="specie")

    @property
    def scaffold_hash(self):
        """Return the Weisfeiler-Lehman graph hash.
        Hashes are identical for isomorphic graphs and there are
        guarantees that non-isomorphic graphs will get different hashes.
        """
        return weisfeiler_lehman_graph_hash(self.nx_graph)

    @property
    def has_atomic_overlaps(self):
        """Check if there are any overlaps in the structure"""
        return self._checks["has_atomic_overlaps"].is_ok

    @property
    def name(self):
        """Return filename if the MOFChecker instance was created based on
        a histogram."""
        return self._name

    @property
    def has_carbon(self):
        """Check if there is any carbon atom in the structure"""
        return self._checks["has_c"].is_ok

    @property
    def has_hydrogen(self):
        """Check if there is any hydrogen atom in the structure"""
        return self._checks["has_h"].is_ok

    @property
    def density(self):
        """Density of structure"""
        return self.structure.density

    @property
    def volume(self):
        """Volume of structure in A^3"""
        return self.structure.volume

    @property
    def formula(self):
        """Return the chemical formula of the structure"""
        return self.structure.formula

    @property
    def has_overvalent_c(self) -> bool:
        """Returns true if there is some carbon in the structure
        that has more than 4 neighbors.

        Returns:
            [bool]: True if carbon with CN > 4 in structure.
        """
        if self._overvalent_c is not None:
            return self._overvalent_c

        self._has_overvalent_c()
        return self._overvalent_c

    @property
    def has_overvalent_h(self) -> bool:
        """Returns true if there is some hydrogen in the structure
        that has more than 1 neighbor.

        Returns:
            [bool]: True if hydrogen with CN > 1 in structure.
        """
        if self._overvalent_h is not None:
            return self._overvalent_h

        self._has_overvalent_h()
        return self._overvalent_h

    @property
    def has_undercoordinated_c(self) -> bool:
        """Check if there is a carbon that likely misses
        hydrogen"""
        if self._undercoordinated_carbon is not None:
            return self._undercoordinated_carbon

        self._has_undercoordinated_carbon()
        return self._undercoordinated_carbon

    @property
    def has_undercoordinated_n(self) -> bool:
        """Check if there is a nitrogen that likely misses
        hydrogen"""
        if self._undercoordinated_nitrogen is not None:
            return self._undercoordinated_nitrogen

        self._has_undercoordinated_nitrogen()
        return self._undercoordinated_nitrogen

    def _has_overvalent_c(self):
        overvalent_c = False
        for site_index in self.c_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 4:
                if not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                    overvalent_c = True
                    break
        self._overvalent_c = overvalent_c

    def _has_overvalent_h(self):
        overvalent_h = False
        for site_index in self.h_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 1:
                if not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                    overvalent_h = True
                    break
        self._overvalent_h = overvalent_h

    def _has_overvalent_n(self):
        overvalent_n = False
        for site_index in self.n_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn > 4:
                if not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                    overvalent_n = True
                    break
        self._overvalent_n = overvalent_n

    def _has_undercoordinated_carbon(self, tolerance: int = 10):
        """Idea is that carbon should at least have three neighbors if it is not sp1.
        In sp1 case it is linear. So we can just check if there are carbons with
        non-linear coordination with less than three neighbors. An example in CoRE
        MOF would be AHOKIR. In principle this should also flag the quite common
        case of benzene rings with missing hydrogens.
        """
        undercoordinated_carbon = False

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
                            undercoordinated_carbon = True
                            break
        self._undercoordinated_carbon = undercoordinated_carbon

    def _has_undercoordinated_nitrogen(self, tolerance: int = 15):
        """
        Attempts to captures missing hydrogens on nitrogen groups
        using heuristics
        """
        undercoordinated_nitrogen = False
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
                    undercoordinated_nitrogen = True
                    break
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
                    break
            elif cn == 3:
                undercoordinated_nitrogen = _guess_underbound_nitrogen_cn3(
                    self.structure, site_index, neighbors, tolerance
                )
                if undercoordinated_nitrogen:
                    break

        self._undercoordinated_nitrogen = undercoordinated_nitrogen

    def _has_high_charges(self, threshold=3) -> Union[bool, None]:
        if (self.charges is None) and HAS_OPENBABEL:
            self.charges = get_charges(self.structure)

        if isinstance(self.charges, list):
            if np.sum(np.abs(self.charges) > threshold):
                return True
        else:
            return None

        return False

    def _is_porous(self) -> Union[bool, None]:
        if self._porous == "":
            self._porous = check_if_porous(self.structure)
        return self._porous

    @property
    def is_porous(self) -> Union[bool, None]:
        """Returns True if the MOF is porous according to the CoRE-MOF definition.
        Returns None if the check could not be run successfully."""
        return self._is_porous()

    @property
    def has_high_charges(self) -> Union[bool, None]:
        """Check if the structure has unreasonably high EqEq charges.
        Returns None if the check could not be run successfully."""
        return self._has_high_charges()

    @property
    def nx_graph(self) -> nx.Graph:
        """Returns a networkx graph with atom numbers as node labels"""
        if self._nx_graph is None:
            _ = self.graph
        return self._nx_graph

    @property
    def graph(self) -> StructureGraph:
        """pymatgen structure graph."""
        if self._graph is None:
            self._graph = StructureGraph.with_local_env_strategy(
                self.structure, self._cnn
            )
            self._nx_graph = construct_clean_graph(self.structure, self._graph)
        return self._graph

    def get_connected_sites(self, site_index) -> List[ConnectedSite]:
        """Get connected sites for given index.

        Uses internal cache for speedup.
        """
        if site_index not in self._connected_sites:
            self._connected_sites[site_index] = self.graph.get_connected_sites(
                site_index
            )
        return self._connected_sites[site_index]

    def get_cn(self, site_index) -> int:
        """Get coordination number for site with CrystalNN method

        Uses internal cache for speedup.

        Args:
            site_index (int): index of site in pymatgen Structure

        Returns:
            int: Coordination number
        """
        if site_index not in self._cns:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._cns[site_index] = self._cnn.get_cn(self.structure, site_index)
        return self._cns[site_index]

    @property
    def has_overvalent_n(self) -> bool:
        """Returns true if there is some nitrogen in the structure
        that has more than 4 neighbors.

        Returns:
            [bool]: True if nitrogen with CN > 4 in structure.
        """
        if self._overvalent_n is not None:
            return self._overvalent_n

        self._has_overvalent_n()
        return self._overvalent_n

    @property
    def has_lone_atom(self) -> bool:
        """Returns True if there is a isolated floating atom"""
        return self._has_lone_atom()

    @property
    def has_lone_molecule(self) -> bool:
        """Returns true if there is a isolated floating atom or molecule"""
        return self._has_stray_molecules()

    def _has_lone_atom(self, safe: bool = True) -> bool:
        for site_index in range(len(self.structure)):
            nbr = self.get_connected_sites(site_index)
            if not nbr:
                lone = True
                # safe option checks if there is really nothing around 1.5 * VdW radius
                # this is useful as sometimes the NN method is off, one example for this
                # is FEZTIP
                if safe:
                    if _vdw_radius_neighbors(self.structure, site_index):
                        lone = False
                if lone:
                    return True
        return False

    @classmethod
    def _from_file(cls, path: str):
        structure = Structure.from_file(path)
        omscls = cls(structure)
        omscls._set_filename(path)  # pylint:disable=protected-access
        return omscls

    @classmethod
    def from_cif(cls, path: Union[str, Path], primitive: bool = True):
        """Create a MOFChecker instance from a CIF file

        Args:
            path (Union[str, Path]): Path to string file
            primitive (bool): If True, it will perform the analysis
                on the primitive structure

        Returns:
            MOFChecker: Instance of MOFChecker
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cifparser = CifParser(path)
            structure = cifparser.get_structures()[0]
            omscls = cls(structure, primitive=primitive)
            omscls._set_filename(path)  # pylint:disable=protected-access
            return omscls

    def _set_cnn(self, method="vesta", porous_adjustment: bool = False):
        if self._cnn_method == method.lower():
            return
        self._cnn_method = method.lower()

        self.porous_adjustment = porous_adjustment
        if method.lower() == "crystalnn":
            self._cnn = CrystalNN(porous_adjustment=self.porous_adjustment)
        elif method.lower() == "econnn":
            self._cnn = EconNN()
        elif method.lower() == "brunnernn":
            self._cnn = BrunnerNN_relative()
        elif method.lower() == "minimumdistance":
            self._cnn = MinimumDistanceNN()
        elif method.lower() == "vesta":
            self._cnn = VESTA_NN
        elif method.lower() == "voronoinn":
            self._cnn = VoronoiNN()
        else:
            self._cnn = JmolNN()

    def _has_stray_molecules(self) -> bool:
        molecules = get_subgraphs_as_molecules_all(self.graph)
        if len(molecules) > 0:
            return True
        return False

    def get_mof_descriptors(self) -> OrderedDict:
        """Run most of the sanity checks
        and get a dictionary with the result

        Returns:
            OrderedDict: result of overall checks
        """
        result_dict = OrderedDict(
            (
                ("name", self.name),
                ("graph_hash", self.graph_hash),
                ("formula", self.formula),
                ("path", self._filename),
                ("density", self.density),
                ("has_carbon", self.has_carbon),
                ("has_hydrogen", self.has_hydrogen),
                ("has_atomic_overlaps", self.has_atomic_overlaps),
                ("has_overcoordinated_c", self.has_overvalent_c),
                ("has_overcoordinated_n", self.has_overvalent_n),
                ("has_overcoordinated_h", self.has_overvalent_h),
                ("has_undercoordinated_c", self.has_undercoordinated_c),
                ("has_undercoordinated_n", self.has_undercoordinated_n),
                ("has_metal", self.has_metal),
                ("has_lone_atom", self.has_lone_atom),
                ("has_lone_molecule", self.has_lone_molecule),
                ("has_high_charges", self.has_high_charges),
                # ("has_undercoordinated_metal", self.has_undercoordinated_metal),
                ("is_porous", self.is_porous),
            )
        )
        return result_dict

    # def _has_low_metal_coordination(self):
    #     for site_index in self.metal_indices:
    #         if _check_metal_coordination(self.structure[site_index],
    #                                      self.get_cn(site_index)):
    #             if self.is_site_open(site_index):
    #                 return True
    #     return False

    # @property
    # def has_undercoordinated_metal(self):
    #     """Check if a metal has unusually low coordination"""
    #     return self._has_low_metal_coordination()

    @property
    def has_metal(self):
        """Checks if there is at least one metal in the structure"""
        self._checks["has_metal"].is_ok

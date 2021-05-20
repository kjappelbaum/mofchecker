# -*- coding: utf-8 -*-
"""MOFChecker: Basic sanity checks for MOFs"""
import os
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import List, Union

import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from pymatgen.analysis.graphs import ConnectedSite, StructureGraph
from pymatgen.core.structure import IStructure, Structure
from pymatgen.io.cif import CifParser

from ._version import get_versions
from .checks.charge_check import ChargeCheck
from .checks.floating_solvent import FloatingSolventCheck
from .checks.global_structure import HasCarbon, HasHydrogen, HasMetal, HasNitrogen
from .checks.local_structure import (
    AtomicOverlapCheck,
    OverCoordinatedCarbonCheck,
    OverCoordinatedHydrogenCheck,
    OverCoordinatedNitrogenCheck,
    UnderCoordinatedCarbonCheck,
    UnderCoordinatedNitrogenCheck,
)
from .checks.oms import MOFOMS
from .checks.utils.get_indices import (
    get_c_indices,
    get_h_indices,
    get_metal_indices,
    get_n_indices,
)
from .checks.zeopp import PorosityCheck
from .definitions import CHECK_DESCRIPTIONS, EXPECTED_CHECK_VALUES
from .graph import _get_cn, construct_clean_graph, get_structure_graph
from .utils import _check_if_ordered

__version__ = get_versions()["version"]
del get_versions

__all__ = ["__version__", "MOFChecker"]


class MOFChecker:  # pylint:disable=too-many-instance-attributes, too-many-public-methods
    """MOFChecker performs basic sanity checks for MOFs"""
    def __init__(self,
                 structure: Union[Structure, IStructure],
                 primitive: bool = True):
        """Class that can perform basic sanity checks for MOF structures

        Args:
            structure (Structure): pymatgen Structure object
            primitive (bool): If True, it will perform the analysis
                on the primitive structure

        Raises:
            NotImplementedError in the case of partial occupancies
        """
        if isinstance(structure, Structure):
            self.structure = IStructure.from_sites(structure)
        else:
            self.structure = structure
        _check_if_ordered(structure)
        if primitive:
            self.structure = self.structure.get_primitive_structure()
        self.metal_indices = get_metal_indices(self.structure)

        self.charges = None
        self._porous = ""
        self.metal_features = None
        self._cnn_method = "vesta"
        self._filename = None
        self._name = None
        self.c_indices = get_c_indices(self.structure)
        self.h_indices = get_h_indices(self.structure)
        self.n_indices = get_n_indices(self.structure)
        self._overvalent_c = None
        self._overvalent_n = None
        self._overvalent_h = None
        self.check_expected_values = EXPECTED_CHECK_VALUES
        self.check_descriptions = CHECK_DESCRIPTIONS

        self._graph = None
        self._nx_graph = None
        self._connected_sites = {}
        self._cns = {}
        self._checks = {
            "has_c":
            HasCarbon(self.structure),
            "has_h":
            HasHydrogen(self.structure),
            "has_metal":
            HasMetal(self.structure),
            "has_nitrogen":
            HasNitrogen(self.structure),
            "no_atomic_overlaps":
            AtomicOverlapCheck(self.structure),
            "no_undercoordinated_carbon":
            UnderCoordinatedCarbonCheck.from_mofchecker(self),
            "no_overcoordinated_carbon":
            OverCoordinatedCarbonCheck.from_mofchecker(self),
            "no_overcoordinated_hydrogen":
            OverCoordinatedHydrogenCheck.from_mofchecker(self),
            "no_overcoordinated_nitrogen":
            OverCoordinatedNitrogenCheck.from_mofchecker(self),
            "no_undercoordinated_nitrogen":
            UnderCoordinatedNitrogenCheck.from_mofchecker(self),
            "no_floating_molecule":
            FloatingSolventCheck.from_mofchecker(self),
            "no_high_charges":
            ChargeCheck(self.structure),
            "is_porous":
            PorosityCheck(self.structure),
            "no_oms":
            MOFOMS.from_mofchecker(self),
        }

    def _set_filename(self, path):
        self._filename = os.path.abspath(path)
        self._name = Path(path).stem

    def get_overlapping_indices(self):
        """Return the indices of overlapping atoms"""
        return self._checks["no_atomic_overlaps"].flagged_indices

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
        return self._checks["no_atomic_overlaps"].is_ok

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
        return not self._checks["no_overcoordinated_carbon"].is_ok

    @property
    def overvalent_c_indices(self) -> bool:
        """Returns indices of carbon in the structure
        that has more than 4 neighbors.

        Returns:
            [list]:
        """
        return self._checks["no_overcoordinated_carbon"].flagged_indices

    @property
    def has_overvalent_h(self) -> bool:
        """Returns true if there is some hydrogen in the structure
        that has more than 1 neighbor.

        Returns:
            [bool]: True if hydrogen with CN > 1 in structure.
        """
        return not self._checks["no_overcoordinated_hydrogen"].is_ok

    @property
    def overvalent_h_indices(self) -> bool:
        """Returns indices of hydrogen in the structure
        that has more than 1 neighbors.

        Returns:
            [list]:
        """
        return self._checks["no_overcoordinated_hydrogen"].flagged_indices

    @property
    def has_undercoordinated_c(self) -> bool:
        """Check if there is a carbon that likely misses
        hydrogen"""
        return not self._checks["no_undercoordinated_carbon"].is_ok

    @property
    def undercoordinated_c_indices(self) -> bool:
        """Returns indices of carbon in the structure
        that likely miss some neighbors.

        Returns:
            [list]:
        """
        return self._checks["no_undercoordinated_carbon"].flagged_indices

    @property
    def has_undercoordinated_n(self) -> bool:
        """Check if there is a nitrogen that likely misses
        hydrogen"""
        return not self._checks["no_undercoordinated_nitrogen"].is_ok

    @property
    def undercoordinated_n_indices(self) -> bool:
        """Returns indices of nitrogen in the structure
        that likely miss some neighbors.

        Returns:
            [list]:
        """
        return self._checks["no_undercoordinated_nitrogen"].flagged_indices

    @property
    def is_porous(self) -> Union[bool, None]:
        """Returns True if the MOF is porous according to the CoRE-MOF definition.
        Returns None if the check could not be run successfully."""
        return self._checks["is_porous"].is_ok

    @property
    def has_high_charges(self) -> Union[bool, None]:
        """Check if the structure has unreasonably high EqEq charges.
        Returns None if the check could not be run successfully."""
        return not self._checks["no_high_charges"].is_ok

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
            self._graph = get_structure_graph(self.structure, self._cnn_method)
            self._nx_graph = construct_clean_graph(self.structure, self._graph)
        return self._graph

    def get_connected_sites(self, site_index) -> List[ConnectedSite]:
        """Get connected sites for given index.

        Uses internal cache for speedup.
        """
        if site_index not in self._connected_sites:
            self._connected_sites[site_index] = self.graph.get_connected_sites(
                site_index)
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
                self._cns[site_index] = _get_cn(self.graph, site_index)
        return self._cns[site_index]

    @property
    def has_overvalent_n(self) -> bool:
        """Returns true if there is some nitrogen in the structure
        that has more than 4 neighbors.

        Returns:
            [bool]: True if nitrogen with CN > 4 in structure.
        """
        return not self._checks["no_overcoordinated_nitrogen"].is_ok

    @property
    def has_lone_molecule(self) -> bool:
        """Returns true if there is a isolated floating atom or molecule"""
        return not self._checks["no_floating_molecule"].is_ok

    @property
    def lone_molecule_indices(self):
        """Returns indices of non-periodic connected component in the structure"""
        return self._checks["no_floating_molecule"].flagged_indices

    @classmethod
    def _from_file(cls, path: str):
        structure = Structure.from_file(path)
        mofchecker = cls(structure)
        mofchecker._set_filename(path)  # pylint:disable=protected-access
        return mofchecker

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

    @property
    def has_oms(self):
        return not self._checks["no_oms"].is_ok

    def _set_cnn(self, method="vesta"):
        if self._cnn_method == method.lower():
            return
        self._cnn_method = method.lower()

    def get_mof_descriptors(self) -> OrderedDict:
        """Run most of the sanity checks
        and get a dictionary with the result

        Returns:
            OrderedDict: result of overall checks
        """
        result_dict = OrderedDict((
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
            ("has_lone_molecule", self.has_lone_molecule),
            ("has_high_charges", self.has_high_charges),
            # ("has_undercoordinated_metal", self.has_undercoordinated_metal),
            ("is_porous", self.is_porous),
        ))
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
        return self._checks["has_metal"].is_ok

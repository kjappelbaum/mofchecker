# -*- coding: utf-8 -*-
"""MOFChecker: Basic sanity checks for MOFs"""
import os
import warnings
from pathlib import Path
from typing import Union

import numpy as np
from pymatgen import Structure
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN, LocalStructOrderParams
from pymatgen.io.cif import CifParser

from ._version import get_versions
from .definitions import OP_DEF
from .utils import (
    HighCoordinationNumber,
    LowCoordinationNumber,
    NoMetal,
    NoOpenDefined,
    get_overlaps,
    get_subgraphs_as_molecules_all,
)

__version__ = get_versions()["version"]
del get_versions

__all__ = ["__version__", "MOFChecker"]


class MOFChecker:  # pylint:disable=too-many-instance-attributes
    """MOFChecker performs basic sanity checks for MOFs"""

    def __init__(self, structure: Structure, porous_adjustment: bool = True):
        """Class that can perform basic sanity checks for MOF structures

        Args:
            structure (Structure): pymatgen Structure object
            porous_adjustment (bool, optional): If true, porous adjustment
                is used for CrystalNN to find the coordination number. Defaults to True.
        """
        self.structure = structure
        self.metal_indices = [
            i for i, species in enumerate(self.structure.species) if species.is_metal
        ]
        self.porous_adjustment = porous_adjustment
        self.metal_features = None
        self._open_indices: set = set()
        self._has_oms = None
        self._cnn = None
        self._filename = None
        self._atomic_overlaps = None
        self._name = None
        self.c_indices = [
            i for i, species in enumerate(self.structure.species) if str(species) == "C"
        ]
        self.h_indices = [
            i for i, species in enumerate(self.structure.species) if str(species) == "H"
        ]
        self.n_indices = [
            i for i, species in enumerate(self.structure.species) if str(species) == "N"
        ]
        self._overvalent_c = None
        self._overvalent_n = None

    def _set_filename(self, path):
        self._filename = os.path.abspath(path)
        self._name = Path(path).stem

    def _get_atomic_overlaps(self):
        if self._atomic_overlaps is not None:
            return self._atomic_overlaps

        self._atomic_overlaps = get_overlaps(self.structure)
        return self._atomic_overlaps

    def get_overlapping_indices(self):
        return self._get_atomic_overlaps()

    @property
    def has_atomic_overlaps(self):
        atomic_overlaps = self._get_atomic_overlaps()
        return len(atomic_overlaps) > 0

    @property
    def name(self):
        return self._name

    @property
    def has_carbon(self):
        return len(self.c_indices) > 0

    @property
    def has_hydrogen(self):
        return len(self.h_indices) > 0

    @property
    def density(self):
        return self.structure.density

    @property
    def volume(self):
        return self.structure.volume

    @property
    def formula(self):
        return self.structure.formula

    @property
    def has_overvalent_c(self) -> bool:
        """Returns true if there is some carbon in the structure that has more than 4 neighbors.

        Returns:
            [bool]: True if carbon with CN > 4 in structure.
        """
        if self._overvalent_c is not None:
            return self._overvalent_c

        self._has_overvalent_c()
        return self._overvalent_c

    def _has_overvalent_c(self):
        overvalent_c = False
        for site_index in self.c_indices:
            cn = self.get_cn(site_index)
            if cn > 4:
                overvalent_c = True
                break
        self._overvalent_c = overvalent_c

    @property
    def has_overvalent_n(self) -> bool:
        """Returns true if there is some nitrogen in the structure that has more than 4 neighbors.

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

    def _has_lone_atom(self):
        self._set_cnn()
        graph = StructureGraph.with_local_env_strategy(self.structure, self._cnn)
        for site in range(len(self.structure)):
            nbr = graph.get_connected_sites(site)
            if not nbr:
                return True
        return False

    def _has_overvalent_n(self):
        overvalent_n = False
        for site_index in self.n_indices:
            cn = self.get_cn(site_index)
            if cn > 4:
                overvalent_n = True
                break
        self._overvalent_n = overvalent_n

    @classmethod
    def _from_file(cls, path: str, porous_adjustment: bool = True):
        s = Structure.from_file(path)
        omscls = cls(s, porous_adjustment)
        omscls._set_filename(path)  # pylint:disable=protected-access
        return omscls

    @classmethod
    def from_cif(cls, path: Union[str, Path], porous_adjustment: bool = True):
        """Create a MOFChecker instance from a CIF file

        Args:
            path (Union[str, Path]): Path to string file
            porous_adjustment (bool, optional):  If true, porous adjustment
                is used for CrystalNN to find the coordination number. Defaults to True.

        Returns:
            MOFChecker: Instance of MOFChecker
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cifparser = CifParser(path)
            s = cifparser.get_structures()[0]
            omscls = cls(s, porous_adjustment)
            omscls._set_filename(path)  # pylint:disable=protected-access
            return omscls

    def _set_cnn(self):
        if self._cnn is None:
            self._cnn = CrystalNN(porous_adjustment=self.porous_adjustment)

    def get_cn(self, site_index: int) -> int:
        """Compute coordination number (CN) for site with CrystalNN method

        Args:
            site_index (int): index of site in pymatgen Structure

        Returns:
            int: Coordination number
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._set_cnn()
            return self._cnn.get_cn(self.structure, site_index)

    def _get_ops_for_site(self, site_index):
        cn = self.get_cn(site_index)
        try:
            names = OP_DEF[cn]["names"]
            is_open = OP_DEF[cn]["open"]
            weights = OP_DEF[cn]["weights"]
            lsop = LocalStructOrderParams(names)
            return (
                cn,
                names,
                lsop.get_order_parameters(self.structure, site_index),
                is_open,
                weights,
            )
        except KeyError:
            # For a bit more fine grained error messages
            if cn <= 3:  # pylint:disable=no-else-raise
                raise LowCoordinationNumber(
                    "Coordination number {} is low and order parameters undefined".format(
                        cn
                    )
                )
            elif cn > 8:
                raise HighCoordinationNumber(
                    "Coordination number {} is high and order parameters undefined".format(
                        cn
                    )
                )

            return cn, None, None, None, None

    def is_site_open(self, site_index: int) -> bool:
        """Check for a site if is open (based on the values of
        some coordination geomeetry fingerprints)

        Args:
            site_index (int): Index of the site in the structure

        Returns:
            bool: True if site is open
        """
        if site_index not in self._open_indices:
            try:
                _, _, lsop, is_open, weights = self._get_ops_for_site(site_index)
                site_open = MOFChecker._check_if_open(lsop, is_open, weights)
                if site_open:
                    self._open_indices.add(site_index)
                return site_open
            except LowCoordinationNumber:
                return True
            except HighCoordinationNumber:
                return None
        return True

    @staticmethod
    def _check_if_open(lsop, is_open, weights, threshold: float = 0.5):
        if lsop is not None:
            if is_open is None:
                return False
            lsop = np.array(lsop) * np.array(weights)
            open_contributions = lsop[is_open].sum()
            close_contributions = lsop.sum() - open_contributions
            return (
                open_contributions / (open_contributions + close_contributions)
                > threshold
            )
        return None

    def _get_metal_descriptors_for_site(self, site_index: int):
        metal = str(self.structure[site_index].species)
        try:
            cn, names, lsop, is_open, weights = self._get_ops_for_site(site_index)
            site_open = MOFChecker._check_if_open(lsop, is_open, weights)
            if site_open:
                self._open_indices.add(site_index)
            descriptors = {
                "metal": metal,
                "lsop": dict(zip(names, lsop)),
                "open": site_open,
                "cn": cn,
            }
        except LowCoordinationNumber:
            descriptors = {"metal": metal, "lsop": None, "open": True, "cn": None}
        except HighCoordinationNumber:
            descriptors = {"metal": metal, "lsop": None, "open": None, "cn": None}
        return descriptors

    def _has_stray_molecules(self) -> bool:
        self._set_cnn()
        sgraph = StructureGraph.with_local_env_strategy(self.structure, self._cnn)
        molecules = get_subgraphs_as_molecules_all(sgraph)
        if len(molecules) > 0:
            return True
        return False

    def get_mof_descriptors(self) -> dict:
        """Run most of the sanity checks
        and get a dictionary with the result

        Returns:
            dict: result of overall checks
        """
        d = {
            "name": self.name,
            "path": self._filename,
            "has_oms": self.has_oms,
            "has_carbon": self.has_carbon,
            "has_hydrogen": self.has_hydrogen,
            "has_atomic_overlaps": self.has_atomic_overlaps,
            "has_overcoordinated_c": self.has_overvalent_c,
            "has_overcoordinated_n": self.has_overvalent_n,
            "has_metal": self.has_metal,
            "has_lone_atom": self.has_lone_atom,
            "has_lone_molecule": self.has_lone_molecule,
            "density": self.density,
        }
        return d

    def get_metal_descriptors_for_site(self, site_index: int) -> dict:
        if not self.has_metal:
            raise NoMetal
        return self._get_metal_descriptors_for_site(site_index)

    def _get_metal_descriptors(self):
        descriptordict = {}
        for site_index in self.metal_indices:
            descriptordict[site_index] = self._get_metal_descriptors_for_site(
                site_index
            )

        self.metal_features = descriptordict

        return descriptordict

    def get_metal_descriptors(self) -> dict:
        """Return local structure order parameters for coordination number (CN),
        element string and wheter site is open or not. Key is the site index.

        Raises:
            NoMetal: If no metal can be found in the structure

        Returns:
            dict: Key is the site index.
        """
        if not self.has_metal:
            raise NoMetal
        return self._get_metal_descriptors()

    @property
    def has_metal(self):
        if self.metal_indices:
            return True
        return False

    @property
    def has_oms(self) -> bool:
        """True if the structure contains open metal sites (OMS).
        Also returns True in case of low coordination numbers (CN <=3)
        which typically also means open coordination for MOFs.
        For high coordination numbers, for which we do not have a good order
        parameter for open structures. For this reason we return None even though
        this might change in a future release.

        Raises:
            NoMetal: Raised if the structure contains no metal

        Returns:
            [bool]: True if the structure contains OMS
        """
        if not self.has_metal:
            raise NoMetal("This structure does not contain a metal")
        if self._has_oms is not None:  # pylint:disable=no-else-return
            return self._has_oms
        else:
            for site_index in self.metal_indices:
                if self.is_site_open(site_index):
                    self._has_oms = True
                    return True
            self._has_oms = False
            return False

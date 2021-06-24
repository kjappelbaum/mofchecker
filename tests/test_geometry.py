# -*- coding: utf-8 -*-
import numpy as np
from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import PeriodicNeighbor, PeriodicSite

from mofchecker.checks.local_structure.geometry import (
    add_sp3_hydrogens_on_cn1,
    get_some_orthorgonal_vector,
    add_methylene_hydrogens,
    add_sp2_hydrogen
)


def test_get_some_orthorgonal_vector():
    some_vec = np.array([1, 1, 1])
    new_vec = get_some_orthorgonal_vector(some_vec)
    assert np.abs(np.dot(some_vec, new_vec)) < 0.01
    assert np.abs(np.linalg.norm(new_vec) - 1) < 0.01


def test_add_sp3_hydrogens_on_cn1():
    site_b_coord = [0, 0, 0]
    site_a_coord = [0.5, 0, 0]

    lattice = Lattice.from_parameters(1, 1, 1, 90, 90, 90)
    site_a = PeriodicSite(Composition("H"), site_a_coord, lattice)
    site_b = PeriodicNeighbor(Composition("H"), site_b_coord, lattice)

    vectors = add_sp3_hydrogens_on_cn1(site_a, [site_b])

    for vector in vectors:
        assert (np.linalg.norm(vector - np.array(site_a_coord)) - 1) < 0.01
        assert vector[0] > 0.5
        assert vector[1] != 0
        assert vector[2] != 0

def test_add_methylene_hydrogens():
    # simplified  test case of existing CN2 coordination 
    site_b_coord = [0.1, 0.1, 0]
    site_a_coord = [0, 0.5, 0]
    site_c_coord = [0.1,0.9,0]

    lattice = Lattice.from_parameters(1, 1, 1, 90, 90, 90)
    site_a = PeriodicSite(Composition("H"), site_a_coord, lattice)
    site_b = PeriodicNeighbor(Composition("H"), site_b_coord, lattice)
    site_c = PeriodicNeighbor(Composition("H"), site_c_coord,lattice)

    hydrogens = add_methylene_hydrogens(site_a, [site_b, site_c])

    assert len(hydrogens) == 2
    for vector in hydrogens:
        assert np.abs(np.linalg.norm(vector-np.array(site_a_coord)) - 1) < 0.01
        assert np.abs(vector[1] -0.5) < 0.01
        assert vector[0] < 0.5
        assert vector[2] <= 0.5 

def test_add_sp2_hydrogen():
    # simplified  test case of existing CN2 coordination 
    site_b_coord = [0.1, 0.1, 0]
    site_a_coord = [0, 0.5, 0]
    site_c_coord = [0.1,0.9,0]

    lattice = Lattice.from_parameters(1, 1, 1, 90, 90, 90)
    site_a = PeriodicSite(Composition("H"), site_a_coord, lattice)
    site_b = PeriodicNeighbor(Composition("H"), site_b_coord, lattice)
    site_c = PeriodicNeighbor(Composition("H"), site_c_coord,lattice)

    hydrogen = add_sp2_hydrogen(site_a, [site_b, site_c])
    assert len(hydrogen) == 3
    print(hydrogen)
    assert np.abs(np.linalg.norm(hydrogen-np.array(site_a_coord)) - 1) < 0.01
    assert np.abs(hydrogen[1]-0.5) < 0.001
    assert hydrogen[0] < 0
    assert np.abs(hydrogen[2]) < 0.001 

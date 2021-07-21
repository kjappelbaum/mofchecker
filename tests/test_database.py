# -*- coding: utf-8 -*-
"""Testing the database module"""
from mofchecker.database import MOFCheckerDB


def test_mofcheckerdb():
    """Test the database lookup"""
    database = MOFCheckerDB()
    assert len(database.lookup_composition("H54 C90 N18 O12")) == 1
    assert len(database.lookup_graph_hash("f75ab1f90320513bad85f2242e7c07ee")) == 1
    assert len(database.lookup_graph_hash("fdblkblabla")) == 0

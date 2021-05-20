# -*- coding: utf-8 -*-

import abc

from ...graph import _get_cn
from ..check_base import AbstractIndexCheck


class BaseCoordinationCheck(AbstractIndexCheck):
    @abc.abstractmethod
    def __init__(self, structure, structure_graph):
        self.structure = structure
        self.structure_graph = structure_graph

    def get_cn(self, index):
        return _get_cn(self.structure_graph, index)

    def get_connected_sites(self, index):
        return self.structure_graph.get_connected_sites(index)

    @classmethod
    def from_mofchecker(cls, mofchecker):
        checker = cls(mofchecker.structure, mofchecker.graph)
        checker.get_cn = mofchecker.get_cn
        checker.get_connected_sites = mofchecker.get_connected_sites
        return checker

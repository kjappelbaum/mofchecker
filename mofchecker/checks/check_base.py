# -*- coding: utf-8 -*-
# pylint: disable=missing-module-docstring,missing-function-docstring,missing-class-docstring
import abc

from backports.cached_property import cached_property


class AbstractCheck(abc.ABC):
    @property
    @abc.abstractmethod
    def description(self):
        pass

    @cached_property
    def is_ok(self):
        return self._run_check()

    @abc.abstractmethod
    def _run_check(self):
        pass


class AbstractIndexCheck(abc.ABC):
    @property
    @abc.abstractmethod
    def description(self):
        pass

    @cached_property
    def is_ok_and_indices(self):
        result, indices = self._run_check()
        return result, indices

    @cached_property
    def is_ok(self):
        result, _ = self.is_ok_and_indices
        return result

    @cached_property
    def flagged_indices(self):
        _, indices = self.is_ok_and_indices
        return indices

    @abc.abstractmethod
    def _run_check(self):
        pass


class AbstractMissingCheck(abc.ABC):
    @property
    @abc.abstractmethod
    def description(self):
        pass

    @property
    @abc.abstractmethod
    def name(self):
        pass

    @cached_property
    def is_ok_indices_positions(self):
        result, indices, positions = self._run_check()
        return result, indices, positions

    @property
    def is_ok(self):
        result, _, _ = self.is_ok_indices_positions
        return result

    @property
    def flagged_indices(self):
        _, indices, _ = self.is_ok_indices_positions
        return indices

    @property
    def candidate_positions(self):
        _, _, positions = self.is_ok_indices_positions
        return positions

    @abc.abstractmethod
    def _run_check(self):
        pass

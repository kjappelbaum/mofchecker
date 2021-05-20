# -*- coding: utf-8 -*-
import abc
from functools import cached_property


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

    @property
    def is_ok(self):
        result, _ = self.is_ok_and_indices
        return result

    @property
    def flagged_indices(self):
        _, indices = self.is_ok_and_indices
        return indices

    @abc.abstractmethod
    def _run_check(self):
        pass

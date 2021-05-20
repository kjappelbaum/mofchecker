# -*- coding: utf-8 -*-
import numpy as np


def _maximum_angle(angle):
    diff_to_180 = np.abs(180 - angle)
    return max([angle, diff_to_180])

# -*- coding: utf-8 -*-
"""Functions for running basic pore analysis with zeo++"""
import os
import subprocess
from tempfile import TemporaryDirectory

from pymatgen import Structure

# from .utils import is_tool

ZEOPP_BASE_COMMAND = ["network", "-ha", "-res"]


def run_zeopp(structure: Structure) -> dict:
    """Run zeopp with network -ha -res (http://www.zeoplusplus.org/examples.html)
    to find the pore diameters

    Args:
        structure (Structure): pymatgen Structure object

    Returns:
        dict: pore analysis results
    """
    with TemporaryDirectory() as tempdir:
        structure_path = os.path.join(tempdir, "structure.cif")
        result_path = os.path.join(tempdir, "result.res")
        structure.to("cif", structure_path)
        cmd = ZEOPP_BASE_COMMAND + [str(result_path), str(structure_path)]
        _ = subprocess.run(
            cmd,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )

        with open(result_path, "r") as handle:
            results = handle.read()

        zeopp_results = parse_zeopp(results)

        return zeopp_results


def parse_zeopp(filecontent: str) -> dict:
    """Parse the results line of a network call to zeopp

    Args:
        filecontent (str): results file

    Returns:
        dict: largest included sphere, largest free sphere,
            largest included sphera along free sphere path
    """
    first_line = filecontent.split("\n")[0]
    parts = first_line.split()

    results = {
        "lis": float(parts[1]),  # largest included sphere
        "lifs": float(parts[2]),  # largest free sphere
        "lifsp": float(parts[3]),  # largest included sphere along free sphere path
    }

    return results

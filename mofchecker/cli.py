# -*- coding: utf-8 -*-
"""Command line interface"""

import json

import click

from mofchecker import MOFChecker


@click.command()
@click.option(
    "--primitive/--no-primitive",
    default=True,
    help="Perform the analysis on the primitive structure",
    show_default=True,
)
@click.argument(
    "CIF_FILES", type=click.Path(exists=True, dir_okay=False), nargs=-1
)
def run(primitive, cif_files):
    """
    Check provided structures and print list of JSON objects with descriptors.
    """
    # Note: we want to see output as things progress,
    # thus this clumsy way of creating a JSON list
    print("[")
    for index, structure_file in enumerate(cif_files):
        mofchecker = MOFChecker.from_cif(structure_file, primitive=primitive)
        descriptors = mofchecker.get_mof_descriptors()

        string = json.dumps(descriptors, indent=2)
        if index != len(structure_files) - 1:
            string += ","
        print(string)
    print("]")

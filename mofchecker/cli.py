# -*- coding: utf-8 -*-
"""Command line interface tool for MOFChecker"""
import click

from . import MOFChecker
from .utils import print_dict


@click.command("cli")
@click.argument("structure", type=click.Path(exists=True))
@click.option("-c", "--check-mof", is_flag=True, help="run all sanity checks")
@click.option("-oms", "--has-oms", is_flag=True, help="check for open metal sites")
@click.option(
    "-hc",
    "--has-carbon",
    is_flag=True,
    help="check if there is carbon in the structure",
)
@click.option(
    "-ovc",
    "--has-overvalent-c",
    is_flag=True,
    help="check if the structure has overcoordinated carbon",
)
@click.option(
    "-ovn",
    "--has-overvalent-n",
    is_flag=True,
    help="check if the structure has overcoordinated nitrogen",
)
@click.option(
    "-clash",
    "--has-clashing",
    is_flag=True,
    help="check if the structure has clashing atoms",
)
@click.option(
    "-la",
    "--has-lone-atom",
    is_flag=True,
    help="check if the structure has floating atoms",
)
def main(  # pylint:disable=too-many-arguments
    structure,
    check_mof,
    has_oms,
    has_carbon,
    has_overvalent_c,
    has_overvalent_n,
    has_clashing,
    has_lone_atom,
):
    """Click CLI"""
    mofchecker = MOFChecker.from_cif(structure)

    if check_mof:
        mof_desc = mofchecker.get_mof_descriptors()
        print_dict(mof_desc)

    elif has_oms:
        print(mofchecker.has_oms)

    elif has_carbon:
        print(mofchecker.has_carbon)

    elif has_overvalent_c:
        print(mofchecker.has_overvalent_c)

    elif has_overvalent_n:
        print(mofchecker.has_overvalent_n)

    elif has_clashing:
        print(mofchecker.has_atomic_overlaps)

    elif has_lone_atom:
        print(mofchecker.has_lone_atom)


if __name__ == "__main__":
    main()  # pylint:disable=no-value-for-parameter

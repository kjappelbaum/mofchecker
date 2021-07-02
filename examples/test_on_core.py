# -*- coding: utf-8 -*-
"""Example that runs the checks on all structures in the CoRE MOF"""
import concurrent.futures
from glob import glob

import pandas as pd
from tqdm import tqdm  # pylint:disable=import-error

from mofchecker import MOFChecker
from mofchecker.errors import NoMetal

all_structures = glob("2019-11-01-ASR-public_12020/structure_10143/*.cif")


def get_feat_one_structure(cif):
    """Run the MOFCheck on one structure"""
    try:
        mofchecker = MOFChecker.from_cif(cif)
        descriptors = mofchecker.get_mof_descriptors()
        return descriptors
    except NoMetal:
        print("{} has no metal".format(cif))
        return None


def main():
    """Loops over all structures"""
    mof_features = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in tqdm(
            executor.map(get_feat_one_structure, all_structures),
            total=len(all_structures),
        ):
            if result is not None:
                mof_features.append(result)

    df = pd.DataFrame(mof_features)  # pylint:disable=invalid-name
    df.to_csv("mof_feat.csv", index=False)


if __name__ == "__main__":
    main()

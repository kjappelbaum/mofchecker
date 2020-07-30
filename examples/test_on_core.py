from glob import glob
from mofchecker import MOFChecker, NoMetal
import pandas as pd
import concurrent.futures
from tqdm import tqdm

all_structures = glob(
    '/Users/kevinmaikjablonka/Downloads/2019-11-01-ASR-public_12020/structure_10143/*.cif'
)


def get_feat_one_structure(cif):
    try:
        mofchecker = MOFChecker.from_cif(cif)
        descriptors = mofchecker.get_mof_descriptors()
        return descriptors
    except NoMetal:
        print('{} has no metal'.format(cif))
        return None


def main():
    mof_features = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in tqdm(executor.map(get_feat_one_structure,
                                        all_structures),
                           total=len(all_structures)):
            if result is not None:
                mof_features.append(result)

    df = pd.DataFrame(mof_features)
    df.to_csv('mof_feat.csv', index=False)


if __name__ == "__main__":
    main()
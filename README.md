# mofchecker

[//]: # "Badges"

![Python package](https://github.com/kjappelbaum/omsdetector/workflows/Python%20package/badge.svg)
[![codecov](https://codecov.io/gh/kjappelbaum/mofchecker/branch/master/graph/badge.svg?token=TQ82D3PFIU)](https://codecov.io/gh/kjappelbaum/mofchecker)
[![Documentation Status](https://readthedocs.org/projects/mofchecker/badge/?version=latest)](https://mofchecker.readthedocs.io/en/latest/?badge=latest)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mofchecker)

## What does it do?

`mofchecker` performs quick sanity checks on crystal structures of metal-organic frameworks (MOFs).

Sanity checks:

- Presence of at least one metal, carbon and hydrogen atom
- Overlapping atoms (distance between atoms above covalent *radius* of the smaller atom)
- Overvalent carbons (coordination number above 4), nitrogens (heuristics), or hydrogens (CN > 1)
- Missing hydrogen on common coordination geometries of C and N (heuristics)
- Atoms with excessive [EQeq partial charge](https://pubs.acs.org/doi/10.1021/jz3008485)

Basic analysis:
- Presence of floating atoms or molecules
- Hash of the atomic structure graph (useful to identify duplicates)

The sanity checks can be used to weed out really unreasonable structures (nothing too fancy).
The code is a rewrite of similar tools in [structure_comp](https://github.com/kjappelbaum/structure_comp).

## Installation

Development version:

```bash
pip install git+https://github.com/kjappelbaum/mofchecker.git
```

Latest stable release

```bash
pip install mofchecker
```

Note that you need to install [zeopp](https://anaconda.org/conda-forge/zeopp-lsmo) if you want to use the porosity features.

```bash
conda install -c conda-forge zeopp-lsmo
```



A web app is currently being developed [in another repository](https://github.com/kjappelbaum/webmofchecker) and deployed on [MatCloud](http://mofchecker.matcloud.xyz/).

## Usage

### Command line interface

```bash
mofchecker --help # list options
mofchecker structure1.cif structure2.cif  # prints JSON output
mofchecker -d has_metal -d has_atomic_overlaps *.cif  # compute only selected descriptors
```

### In Python

```python
from mofchecker import MOFChecker
mofchecker = MOFChecker.from_cif(<path_to_cif>)
# or: MOFChecker(structure=my_pymatgen_structure)

# Test for OMS
mofchecker.has_oms

# Test for clashing atoms
mofchecker.has_overlapping_atoms

# Run basic checks on a list of cif paths (sample_structures)
results = []

for structure in sample_structures:
    mofchecker = MOFChecker.from_cif(structure)
    results.append(mofchecker.get_mof_descriptors())
```

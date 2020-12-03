# mofchecker

[//]: # "Badges"

![Python package](https://github.com/kjappelbaum/omsdetector/workflows/Python%20package/badge.svg)
[![codecov](https://codecov.io/gh/kjappelbaum/mofchecker/branch/master/graph/badge.svg?token=TQ82D3PFIU)](https://codecov.io/gh/kjappelbaum/mofchecker)
[![Documentation Status](https://readthedocs.org/projects/mofchecker/badge/?version=latest)](https://mofchecker.readthedocs.io/en/latest/?badge=latest)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mofchecker)

## What does it do?

Perform quick sanity checks on your MOF:

- Find open metal sites (OMS) in metal-organic frameworks (MOFs).
- Find atomic overlaps.
- Find overvalent (CN>4) carbons, nitrogens, or hydrogen.
- Check if there is metal, carbon or hydrogen.
- Check if there is floating atoms or molecules.
- Check if there is missing hydrogen on common coordination geometries of C and N.

The idea is to have nothing to fancy but a fast tool that we can run to eliminate the really unreasonable structures. The code is basically a rewrite of the checking tools that we implemented in [structure_comp](https://github.com/kjappelbaum/structure_comp).

## Installation

Development version:

```bash
pip install git+https://github.com/kjappelbaum/mofchecker.git
```

Latest stable release

```bash
pip install mofchecker
```

## Usage

### In Python

```python
from mofchecker import MOFChecker
mofchecker = MOFChecker.from_cif(<path_to_cif>)

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

### CLI

For example, you can use

```bash
mofchecker <cif> --has-oms
```

You can get an overview over all options with

```bash
mofchecker --help
```

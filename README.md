# mofchecker

[//]: # "Badges"

![Python package](https://github.com/kjappelbaum/omsdetector/workflows/Python%20package/badge.svg)
[![codecov](https://codecov.io/gh/kjablonk/omsdetector/branch/master/graph/badge.svg)](https://codecov.io/gh/kjablonk/omsdetector/branch/master)

## What does it do?

Perform quick sanity checks on your MOF:

- Find open metal sites (OMS) in metal-organic frameworks (MOFs).
- Find atomic overlaps.
- Find overvalent (CN>4) carbons.
- Check if there is metal, carbon or hydrogen.

The idea is to have nothing to fancy but a fast tool that we can run to eliminate the really unreasonable structures.

A basic CLI is in development.

### OMS detection

If works similar to the [OMS detector used for the CoRE-MOF database](https://github.com/emmhald/open_metal_detector) but it attempts to be a bit more lightweight and correct in some cases (like paddlewheels).

We use [pymatgen](https://pymatgen.org) to calculate the [order-parameters proposed by Zimmermann and Jain](https://pubs.rsc.org/en/content/articlelanding/2020/RA/C9RA07755C#!divAbstract). Note that we introduce some weighting factors.

Identifying OMS is important as our [force-fields fail to describe them correctly](https://pubs.acs.org/doi/10.1021/acs.jpcc.7b02302).

### Overlapping atoms

Overlapping atoms are detected based on the pairwise distance matrix and the covalent radii.

### Why check for carbons?

I think it is a pretty fast check if there is something organic ...

## Installation

Development version:

```(bash)
pip install git+https://github.com/kjappelbaum/mofchecker.git
```

## Usage

### In Python

```(python)
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

```(bash)
mofchecker <cif> --has-oms
```

You can get an overview over all options with

```(bash)
mofchecker --help
```

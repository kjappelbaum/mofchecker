<!--
<p align="center">
  <img src="https://github.com/kjappelbaum/mofchecker/raw/main/docs/source/figures/logo.png" height="300">
</p> -->
<h1 align="center">
    mofchecker
</h1>
<p align="center">
    <a href="https://github.com/kjappelbaum/mofchecker/actions?query=workflow%3Apython_package">
        <img alt="Tests" src="https://github.com/kjappelbaum/mofchecker/actions/workflows/python-package.yml/badge.svg" />
    </a>
    <a href="https://pypi.org/project/mofchecker">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/mofchecker" />
    </a>
    <a href="https://pypi.org/project/mofchecker">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/mofchecker" />
    </a>
    <a href="https://github.com/kjappelbaum/mofchecker/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/mofchecker" />
    </a>
    <a href='https://mofchecker.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/mofchecker/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' />
    </a>
</p>


## What does it do?

`mofchecker` performs quick sanity checks on crystal structures of metal-organic frameworks (MOFs).

Try the live web app at https://github.com/kjappelbaum/webmofchecker !

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

## 🚀 Installation

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

## 💪 Getting Started

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
mofchecker.has_atomic_overlaps

# Run basic checks on a list of cif paths (sample_structures)
results = []

for structure in sample_structures:
    mofchecker = MOFChecker.from_cif(structure)
    results.append(mofchecker.get_mof_descriptors())
```


## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.rst](https://github.com/kjappelbaum/mofchecker/blob/master/CONTRIBUTING.rst) for more information on getting involved.


### ⚖️ License

The code in this package is licensed under the MIT License.


### 💰 Funding

The research was supported by the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme ([grant agreement 666983, MaGic](https://cordis.europa.eu/project/id/666983)), by the [NCCR-MARVEL](https://www.nccr-marvel.ch/), funded by the Swiss National Science Foundation, and by the Swiss National Science Foundation (SNSF) under Grant 200021_172759.


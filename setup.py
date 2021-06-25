# -*- coding: utf-8 -*-
"""Perform sanity checks on MOF structures"""
import sys

from setuptools import find_packages, setup

import versioneer

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

with open("README.md", "r") as handle:
    LONG_DESCRIPTION = handle.read()

setup(
    # Self-descriptive entries which should always be present
    name="mofchecker",
    author="Kevin M. Jablonka",
    author_email="kevin.jablonka@epfl.ch",
    description="Perform sanity checks on MOF structures",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://github.com/kjappelbaum/mofchecker",
    license="GPL",
    packages=find_packages(),
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=pytest_runner,
    install_requires=[
        "pymatgen>=2021.1,<2022",
        "click==8.*",
        "networkx>=2.5",
        "pyeqeq",
        "backports.cached-property",
        "ase",
        "ruamel.yaml",
    ],
    extras_require={
        "testing": ["pytest", "pytest-cov<2.12"],
        "docs": [
            "sphinx",
            "sphinx-book-theme",
            "sphinx-autodoc-typehints",
            "sphinx-copybutton",
        ],
        "dev": ["versioneer"],
        "pre-commit": ["pylint", "pre-commit"],
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)

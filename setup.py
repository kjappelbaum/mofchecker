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

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    # Self-descriptive entries which should always be present
    name="mofchecker",
    author="Kevin M. Jablonka",
    author_email="kevin.jablonka@epfl.ch",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=pytest_runner,
    install_requires=requirements,
    extras_require={
        "testing": ["pytest", "pytest-cov<2.11"],
        "docs": [
            "sphinx",
            "sphinx-book-theme",
            "sphinx-autodoc-typehints",
            "sphinx-copybutton",
        ],
        "dev": ["versioneer"],
    },
    entry_points={
        "console_scripts": [
            "mofchecker = mofchecker.cli:main",
        ],
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)

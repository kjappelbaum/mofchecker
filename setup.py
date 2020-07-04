# -*- coding: utf-8 -*-
import sys

from setuptools import find_packages, setup

import versioneer

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open('README.md', 'r') as handle:
        long_description = handle.read()
except:
    long_description = '\n'.join(short_description[2:])

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    # Self-descriptive entries which should always be present
    name='mofchecker',
    author='Kevin M. Jablonka',
    author_email='kevin.jablonka@epfl.ch',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=pytest_runner,
    install_requires=requirements,
    extras_require={
        'docs': ['guzzle_sphinx_theme'],
        'pre-commit': [
            'pre-commit',
            'yapf',
            'prospector',
            'pylint',
            'versioneer',
            'isort',
        ],
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)

Getting started with mofchecker
===================================

Installation
--------------

We recommend installing mofchecker in a clean virtual environment environment (e.g., a `conda environment <https://docs.conda.io/projects/conda/en/latest/index.html>`_)


You can install the latest stable release from PyPi using

.. code-block:: bash

    pip install mofchecker


or the latest development version using

.. code-block:: bash

    pip install git+https://github.com/kjappelbaum/mofchecker.git


If you want to use the charge check (based on the EqEq method), you'll need to install openbabel, for example using 

.. code-block:: bash

    conda install openbabel

Running checks
----------------

To run all checks you only need the following code snippet

.. code-block:: python

    from mofchecker import MOFChecker

    # starting from a pymatgen structure object
    checker = MOFChecker(<pymatgen_structure_object>)

    # alternatively, starting from a file
    checker = MOFChecker.from_file(<path_to_file>)

    check_result = checker.get_mof_descriptors()

:code:`check_result` will be a :code:`OrderedDict` in which the keys are the names of the checks and the values are the check results.

The "ideal"/"expected" values for the checks are defined in :py:attr:`~mofchecker.check_expected_values`.

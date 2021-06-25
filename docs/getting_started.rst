Getting started with mofchecker
===================================

Installation
--------------

We recommend installing mofchecker in a clean virtual environment (e.g., a `conda environment <https://docs.conda.io/projects/conda/en/latest/index.html>`_)


You can install the latest stable release from PyPi using

.. code-block:: bash

    pip install mofchecker


or the latest development version using

.. code-block:: bash

    pip install git+https://github.com/kjappelbaum/mofchecker.git


if you want to use the porosity check, you'll need to have the `network` binary of the zeo++ library in your PATH. You can install it using

.. code-block:: bash

  conda install -c conda-forge zeopp-lsmo

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

The “ideal”/”expected” values for the checks are defined in :py:attr:`~mofchecker.check_expected_values`.

Adding missing hydrogens
--------------------------

We currently do not provide a fully-automated routine as the algorithm for identification of candidate sites for hydrogen is still quite primitive and not tested in production.

.. code-block:: python

    from mofchecker import MOFChecker
    checker = MOFChecker.from_file(<path_to_file>)

    # mofchecker uses an immutable object internally to avoid side effects
    # and to allow hashing
    mutable_structure = Structure.from_sites(checker.structure.sites)

    # here we flatten a list that can potentially be a list of lists
    # of Cartesian coordinates
    h_positions = sum([], undercoordinated_c_candidate_positions)

    # now, we can append the new sites to the structure
    for h in h_positions:
        h = mutable_structure.lattice.get_fractional_coords(h)
        mutable_structure.append('H',  h)

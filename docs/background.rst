Background information about the checks
=========================================


Porosity check
----------------
`Potential porosity <https://blogs.rsc.org/ce/2013/01/08/iupac-provisional-recommendations-on-metal-organic-framework-and-coordination-polymer-terminology/?doing_wp_cron=1616568093.4138350486755371093750>`_ is thought to be a key defining criterion for MOFs.
For assembling the widely used CoRE-MOF database (see `first <https://pubs.acs.org/doi/10.1021/cm502594j>`_ and `second <https://pubs.acs.org/doi/10.1021/acs.jced.9b00835>`_ paper) Yongchul G. Chung et al. used a pore limiting diameter (PLD) of more than 2.4 Å (approximately the van der Waals diameter of a hydrogen molecule) as threshold for the distinction of porous and non-porous (or nanoporous) MOFs.

In mofchecker, we follow this definition and use `zeopp <http://www.zeoplusplus.org/>`_ with the high-accuracy flag and the default atom radii to compute the largest free sphere. If it is above or equal to 2.4 Å the :py:attr`~mofchecker.MOFChecker.is_porous` will return :code:`True`.

Graph hash
-----------
The structure graph hashes that mofchecker returns are calculated using the Weisfeiler Lehman (WL) algorithm.
It is an algorithm that iteratively performs neighborhood aggregation and is guaranteed to give different hashes for different structure and hash strong guarantees that it will produce different hashes for different structure. Unfortunately, there might be cases where different structures are mapped to the same hash.

What variants does mofchecker offer?
..........................................

Currently, mofchecker implements :py:attr:`~mofchecker.MOFChecker.graph_hash` and :py:attr:`~mofchecker.MOFChecker.scaffold_hash`. The :py:attr:`~mofchecker.MOFChecker.graph_hash` takes the atom labels into account and will, for example, return different hashes for Ni-MOF-74 and Mg-MOF-74. The :py:attr:`~mofchecker.MOFChecker.scaffold_hash` does not take atom labels into account and will return the same hash for structures with the same connectivity. That is, it will return the same hash for Ni-MOF-74 and Mg-MOF-74.

.. image:: _static/hash_comparison_mof_74.jpg
  :width: 600
  :alt: Comparison of graph and scaffold hash.


What can it be used for?
............................

The most important use case for the hashes is to identify duplicates in databases for which it is unfeasible to perform :math:`N^2` graph isomorphirsm checks.
Note that the definition of duplicate, under the :py:attr:`~mofchecker.MOFChecker.graph_hash`, is similar to the one proposed by `Barthel et al. <https://pubs.acs.org/doi/pdf/10.1021/acs.cgd.7b01663>`_ :

    However, from a MOF point of view two structures are considered identical if they share the same bond network, with respect to the atom types and their embedding:
    i.e., if two structures can in principle be deformed into each other without breaking and forming bonds.

The scaffold hash can be useful to find families of related MOFs. For example, all members of the unfunctionalized MOF-74 familiy would group under the same hash. Similarly, all functionalized UiO-66 structures would group under the same hash as long as the position and the size (e.g., one atom) of the functionalization is the same.

What can go wrong?
.......................

There are multiple ways in which the hash can yield results different from what one would expect.

1. Incorrect bonding network
2. Did not reduce to primitive cell
3. Unlucky hash clash (Weisfeiler Lehman has some `edge cases <https://informaconnect.com/beyond-weisfeiler-lehman-using-substructures-for-provably-expressive-graph-neural-networks/>`_)


How does it work?
....................

Under the hood, mofchecker uses `pymatgen <http://pymatgen.org/>`_ to analyze the bonding network to create a structure graph.
We then apply the Weisfeiler Lehman algorithm (as implemented in `networkx <https://networkx.org/>`_) to compute the hash.

The Weisfeiler Lehman algorithm is explained in the `English translation of the original paper <https://www.iti.zcu.cz/wl2018/pdf/wl_paper_translation.pdf>`_
and a `popular blog post <https://davidbieber.com/post/2019-05-10-weisfeiler-lehman-isomorphism-test/#:~:text=The%20core%20idea%20of%20the,used%20to%20check%20for%20isomorphism>`_. Briefly, the algorithm uses iterative recoloring to compute a hash.
In the first iteration on might label the structures by the number of neighbors (as for the scaffold hash) or the atom type. Then one builds new labels based on the past labels of the intermediate neighbors. The figure below (adopted from `Michael Bronstein's blog <https://towardsdatascience.com/expressive-power-of-graph-neural-networks-and-the-weisefeiler-lehman-test-b883db3c7c49>`_) illustrates this

.. image:: _static/wl_hash.png
  :width: 600
  :alt: Illustration of the WL hashing algorithm, based on https://towardsdatascience.com/expressive-power-of-graph-neural-networks-and-the-weisefeiler-lehman-test-b883db3c7c49.

The initial labeling here is based on the connectivity and in the next step we extend the labels with the labels of the nearest neighbors and we continue this process until self-consistency (or a maximum number of iterations) are reached. A histogram of the node colors can then be converted into a hash string.

Background information about the checks
=========================================


Porosity check
----------------
`Potential porosity <https://blogs.rsc.org/ce/2013/01/08/iupac-provisional-recommendations-on-metal-organic-framework-and-coordination-polymer-terminology/?doing_wp_cron=1616568093.4138350486755371093750>`_ is thought to be a key defining criterion for MOFs.
For assembling the widely used CoRE-MOF database (see `first <https://pubs.acs.org/doi/10.1021/cm502594j>`_ and `second <https://pubs.acs.org/doi/10.1021/acs.jced.9b00835>`_ paper) Yongchul G. Chung et al. used a pore limiting diameter (PLD) of more than 2.4 Å (approximately the van der Waals diameter of a hydrogen molecule) as threshold for the distinction of porous and non-porous (or nanoporous) MOFs.

In mofchecker, we follow this definition and use `zeopp <http://www.zeoplusplus.org/>`_ with the high-accuracy flag and the default atom radii to compute the largest free sphere. If it is above or equal to 2.4 Å the :py:attr`mofchecker.MOFChecker.is_porous` will return `True`.

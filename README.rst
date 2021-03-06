.. -*- mode: rst; coding: utf-8 -*-

==========================
 Path Similarity Analysis
==========================

|zenodo|

:Author:    Sean Seyler
:Year:      2015
:License:   GNU Public Licence, version 3 (or higher)
:Copyright: © 2015 Sean Seyler
:Citation:  SL Seyler, A Kumar, MF Thorpe, O Beckstein. 
            *Path Similarity Analysis: a Method for Quantifying Macromolecular Pathways.* 
            ArXiv e-prints (2015). `arXiv:1505.04807`_ [q-bio.QM]
            
.. |zenodo| image:: https://zenodo.org/badge/13219/Becksteinlab/PSAnalysisTutorial.svg
    :alt: doi: 10.5281/zenodo.17902
    :target: http://dx.doi.org/10.5281/zenodo.17902

Summary
=======

Path Similarity Analysis (PSA) comprises a computational framework designed to
enhance the quantitative comparison of macromolecular transition paths [Seyler2015]_. 
This tutorial provides two examples to demonstrate a simple comparison, using PSA, of
closed to open adenylate kinase (AdK) transition paths generated by a selection
of various algorithms [Seyler2014]_. Hierarchical clustering is used as a
simple, but powerful approach to exploratory data analysis by construction of a
heat map-dendrogram representation of the quantitative comparison.


Background
==========

PSA is based on measuring the geometric similarity of transition paths in
configuration space using the Hausdorff and Fréchet path metrics. PSA takes
advantage of MDAnalysis_ [Michaud-Agrawal2011]_ to provide a seamless interface
to the Python and NumPy arrays, and a mechanism for performing path comparisons
using arbitrary atom selections. MDAnalysis also provides a format-agnostic
framework for reading simulation trajectories, allowing rapid comparison of many
different computational methods. More information about the method can be found
in [Seyler2015]_.


Usage
=====

This tutorial demonstrates a straightforward application of PSA to a set of
transitions of the enzyme adenylate kinase (AdK) generated by a selection of methods
(for more background on this particular example see [Seyler2014]_). Two example scripts are
provided: a short version shows how to perform similarity analysis on a set of
trajectories that have been pre-processed for proper (frame-by-frame) structural
alignment; a full version additionally demonstrates, using the PSA framework,
how an alignment procedure would be performed prior to similarity analysis.

Analysis is performed with the ``psa_short.py`` or ``psa_full.py`` scripts,
which automatically read trajectories from the ``methods`` directory into a
PSA object, perform trajectory alignment (in the case of ``psa_full.py``),
generate discrete Fréchet distance matrices, and produce a heat map-dendrogram
plot representing the distance matrix after Ward hierarchical clustering.

The scripts can be run directly using

    python psa_short.py

or

    python psa_full.py

The user can also try adjusting settings at the top of each file to change the:

* path metric (default: discrete Fréchet [``discrete_frechet``])
* linkage algorithm for hierarchical clustering (default: ``Ward``)
* name and location of the plot (default: ``df_ward_psa-[short/full].pdf``)

These examples serve as a sufficient basis for understanding PSA's framework.
Some other techniques and analyses using PSA are described in [Seyler2015]_.


Help
====

If you have questions or problems using the package then ask on
the MDAnalysis user mailing list:
http://groups.google.com/group/mdnalysis-discussion


Implementation in MDAnalysis
============================

If you want to write your own code using PSA then use the ``MDAnalysis.analysis.psa``
module, which is part of MDAnalysis_ (since release 0.10.0) and have a look at
the `documentation of the PSA module`_.

.. _documentation of the PSA module: 
   http://devdocs.mdanalysis.org/documentation_pages/analysis/psa.html


References
==========

.. Links
.. -----

.. _MDAnalysis: http://www.mdanalysis.org

.. Articles
.. --------

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning,
   T. B. Woolf, and O. Beckstein. MDAnalysis: A toolkit for the
   analysis of molecular dynamics simulations. J Comp Chem,
   32:2319-2327, 2011. doi:`10.1002/jcc.21787`_. http://www.mdanalysis.org

.. _`10.1002/jcc.21787`: http://doi.org/10.1002/jcc.21787

.. [Seyler2014] S.L. Seyler and O. Beckstein, Sampling large conformational
   transitions: adenylate kinase as a testing ground. Mol Simul 40:855–877,
   2014. doi:`10.1080/08927022.2014.919497`_

.. _`10.1080/08927022.2014.919497`: http://dx.doi.org/10.1080/08927022.2014.919497

.. [Seyler2015] S.L. Seyler, A. Kumar, M.F. Thorpe, and O. Beckstein, Path
   Similarity Analysis: a Method for Quantifying Macromolecular Pathways.
   `arXiv:1505.04807`_ [q-bio.QM], 2015

.. _`arXiv:1505.04807`: http://arxiv.org/abs/1505.04807

.. IceBin documentation master file, created by
   sphinx-quickstart on Thu May  3 11:57:26 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview
========

*IceBin* is a C++ library, Python API and set of command-line utilites
 that provides the following:

#. Computation of regridding matrices involving *elevation classes*
#. Production of elevation-class enabled ModelE-specific input files.
#. Python API to work with elevation class data
#. Coupling of GCM with a two-way ice model using elevation classes.



ModelE Command-Line Utilities
=============================

These command-line utilities, included with IceBin, are used
ModelE-specific command line utilities included with IceBin.


.. toctree::
   :maxdepth: 2
   :caption: Theory and Infrastructure

   conservative_regridding
   sparse_matrices
   matrix_formats
   icebin_regridding
   mismatched_regridding


.. toctree::
   :maxdepth: 2
   :caption: ModelE Command-Line Utilities

   giss2nc
   etopo1_ice
   make_topoo
   global_ec
   combine_global_ec
   make_topoa



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. .. raw:: html
.. 
..    <object data="_static/ABCRegridding.svg" type="image/svg+xml"></object>

.. _giss2nc

giss2nc
=======

This program converts traditional GISS format Fortran files to NetCDF4.  For example:

.. code-block:: bash

   $ giss2nc Z1QX1N.BS1 Z1QX1N.BS1.nc

Some optional arguments can be used as needed; try ``giss2nc --help``
for more information:

* Use ``--endian=big`` if the Fortran-format
  file was written on a big-endian system; for example, older Sun
  systems.  Newer Intel-based systems write little-endian files by
  default.  See `here <https://en.wikipedia.org/wiki/Endianness>` for
  more information.

* By default, ``giss2nc`` reads single-precision ("float") and writes
  double-precision ("double") data.  It can be configured to
  read/write other data types as well, via the ``--type`` option.  The
  following types are current supported, and more can be added easily:

  * ``--type=float-double``: Read float, write double.
  * ``--type=int16``: Read and write 2-byte integers.

* Fortran-format files don't always provide a well-defined name for
  each variable.  The ``--names`` parameter can be used to specify
  these names, which will be used in the output NetCDF file.  For
  example, ``--names=A,B`` will store the first variable in the
  Fortran file as ``A`` and the second one as ``B``; if the input
  contains any further variables, they will be skipped.  The directive
  ``--names=A,_,C`` will process the first and third variable in the
  input file and skip the second.

Putting it together, the following is a typical command line:

.. code-block:: bash

   $ giss2nc ZETOPO1.NCEI ZETOPO1.NCEI.nc --endian=big \
     --names=FOCEAN,ZICTOP,ZSOLID --type=int16



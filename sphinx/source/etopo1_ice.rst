.. _etopo1_ice

etopo1_ice
==========

This program generates a global map of ice extent and elevation on a
1-minute grid.  It uses the following ModelE input files, found bin
searching the colon-separated path in the ``MODELE_FILE_PATH``
environment variable.

#. ``ZETOPO1.NCEI-SeparateGreenland.nc`` (1-minute): Extent of ice
   sheet, plus elevation, on global 1-minute grid.  This is a modified
   version of the standard ModelE input files ``ZETOPO1.NCEI``, in
   which ``FOCEAN`` (fractional ocean cover) is set to 1 for ocean
   grid cells, 0 for continent, an 2 for Greenland.  This convention
   allows ``etopo1_ice`` to remove Greenland from its output, if
   needed, for two-way coupled runs where a modeled Greenland will be
   used in its place.
#. ``ZNGDC1`` (1-degree): Ice sheet fractions in northern hemisphere.

``etopo_ice`` takes a singe argument, naming its output file (without
the ``.nc`` extension).  The ``-g`` flag may also be used to include
Greenland in the output; otherwise, the Greenland island and ice sheet
are removed, so a model-based Greenland may be used instead (for
two-way coupling).

For example, ``etopo1_ice -g etopo1_ice_g`` produces the files
``etopo1_ice_g1m.nc``, which contains global ice extent and elevation
at 1-minute resolution:

.. code-block:: none

   netcdf etopo1_ice_g1m {
   dimensions:
           jm1m = 10800 ;
           im1m = 21600 ;
   variables:
           short ZICETOP1m(jm1m, im1m) ;
                   ZICETOP1m:units = "m" ;
                   ZICETOP1m:source = "ETOPO1" ;
           short ZSOLG1m(jm1m, im1m) ;
                   ZSOLG1m:units = "m" ;
                   ZSOLG1m:source = "ETOPO1" ;
           short FOCEAN1m(jm1m, im1m) ;
                   FOCEAN1m:units = "m" ;
                   FOCEAN1m:source = "ETOPO1" ;
           short FGICE1m(jm1m, im1m) ;
                   FGICE1m:description = "Fractional ice cover (0 or 1)" ;
                   FGICE1m:units = "1" ;
                   FGICE1m:source = "etopo1_ice.cpp output" ;

   // global attributes:
                   :source = "etopo1_ice.cpp" ;
   }

Variables are as follows:

* ``ZICETOP``: Elevation of top of the ice
* ``ZSOLG``: Elevation of the ground.  Ice thickness is
  ``ZICETOP-ZSOLG``.  For areas not covered by ice,
  ``ZICETOP==ZSOLG``.
* ``FOCEAN``: Fraction of grid cell covered in ocean. Will be 0 or 1
  for 1-minute file.
* ``FGCIE``: Fraction of grid cell covered in ice. Will be 0 or 1 for
  1-minute file.



It also produces the file ``etopo1_ice_gh.nc``, which contains the
same data, averaged to 1/2-degree resolution.  This file is easier to
plot or otherwise visualize.



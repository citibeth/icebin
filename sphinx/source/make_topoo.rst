.. _make_topoo

make_topoo
==========

This program generates a ModelE TOPO file *on the ocean grid*; which
can later be regridded to the atmosphere grid using make_topoa_.  It
is based on the program ``Z1QX1N.BS1.F`` by Gary Russell.

.. note::

   A ModelE TOPO file describes land surface fractions for each ModelE
   land surface type).

Required inputs, which are found by searching the colon-separated
environment variable ``MODELE_FILE_PATH``, are:

#. ``etopo1_ice_g1m.nc``: Contains a high-resolution version of the
   variables ``FGICE1m`` (ground ice fraction), ``ZICETOP1m``
   (elevation, top of ice), ``ZSOLG1m`` (elevation, solid ground) and
   ``FOCEAN1m``.

#. ``Z10MX10M.nc``: Used for ``FLAKES``, that is lake extent.

``make_topoo`` takes a single command line argument, the name of the
output file.  For example:

.. code-block:: none

   make_topoo topoo.nc

The output file looks like this:

.. code-block:: none

   $ ncdump -h topoo.nc
   netcdf topoo {
   dimensions:
           jm = 180 ;
           im = 288 ;
   variables:
           int64 hspec ;
                   hspec:comment = "Parameters to instantiate a icebin::modele::HntrGrid description of this grid" ;
                   hspec:im = 288 ;
                   hspec:jm = 180 ;
                   hspec:offi = 0. ;
                   hspec:dlat = 60. ;
           double FOCEAN(jm, im) ;
                   FOCEAN:description = "0 or 1, Bering Strait 1 cell wide" ;
                   FOCEAN:units = "1" ;
                   FOCEAN:source = "GISS 1Qx1" ;
           double FLAKE(jm, im) ;
                   FLAKE:description = "Lake Surface Fraction" ;
                   FLAKE:units = "0:1" ;
                   FLAKE:sources = "GISS 1Qx1" ;
           double FGRND(jm, im) ;
                   FGRND:description = "Ground Surface Fraction" ;
                   FGRND:units = "0:1" ;
                   FGRND:sources = "GISS 1Qx1" ;
           double FGICE(jm, im) ;
                   FGICE:description = "Glacial Ice Surface Fraction" ;
                   FGICE:units = "0:1" ;
                   FGICE:sources = "GISS 1Qx1" ;
           double FGICE_greenland(jm, im) ;
                   FGICE_greenland:description = "Greenland Ice Surface Fraction" ;
                   FGICE_greenland:units = "0:1" ;
                   FGICE_greenland:sources = "GISS 1Qx1" ;
           double ZATMO(jm, im) ;
                   ZATMO:description = "Atmospheric Topography" ;
                   ZATMO:units = "m" ;
                   ZATMO:sources = "ETOPO2 1Qx1" ;
           double dZOCEN(jm, im) ;
                   dZOCEN:description = "Ocean Thickness" ;
                   dZOCEN:units = "m" ;
                   dZOCEN:sources = "ETOPO2 1Qx1" ;
           double dZLAKE(jm, im) ;
                   dZLAKE:description = "Lake Thickness" ;
                   dZLAKE:units = "m" ;
                   dZLAKE:sources = "ETOPO2 1Qx1" ;
           double dZGICE(jm, im) ;
                   dZGICE:description = "Glacial Ice Thickness" ;
                   dZGICE:units = "m" ;
                   dZGICE:sources = "Ekholm,Bamber" ;
           double ZSOLDG(jm, im) ;
                   ZSOLDG:description = "Solid Ground Topography" ;
                   ZSOLDG:units = "m" ;
                   ZSOLDG:sources = "ETOPO2 1Qx1" ;
           double ZICETOP(jm, im) ;
                   ZICETOP:description = "Atmospheric Topography (Ice-Covered Regions Only)" ;
                   ZICETOP:units = "m" ;
                   ZICETOP:sources = "ETOPO2 1Qx1" ;
           double ZSGLO(jm, im) ;
                   ZSGLO:description = "Lowest Solid Topography" ;
                   ZSGLO:units = "m" ;
                   ZSGLO:sources = "ETOPO2 1Qx1" ;
           double ZLAKE(jm, im) ;
                   ZLAKE:description = "Lake Surface Topography" ;
                   ZLAKE:units = "m" ;
                   ZLAKE:sources = "ETOPO2 1Qx1" ;
           double ZGRND(jm, im) ;
                   ZGRND:description = "Topography Break between Ground and GIce" ;
                   ZGRND:units = "m" ;
                   ZGRND:sources = "ETOPO2 1Qx1" ;
           double ZSGHI(jm, im) ;
                   ZSGHI:description = "Highest Solid Topography" ;
                   ZSGHI:units = "m" ;
                   ZSGHI:sources = "ETOPO2 1Qx1" ;
           double FOCEANF(jm, im) ;
                   FOCEANF:description = "Fractional ocean ocver" ;
                   FOCEANF:units = "1" ;
                   FOCEANF:sources = "GISS 1Qx1" ;
   }

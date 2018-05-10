.. _global_ec

global_ec
=========

*IceBin* was originally designed to define elevation classes for a
local ice sheet, eg. Greenland or Antarctica.  ``global_ec`` uses
*IceBin* to define elevation classes for global ice cover.

``global_ec`` works with four grids, all of them latitude/longitude
type grids:

#. ``A``: Atmosphere grid --- defined by the GCM
#. ``I``: Ice grid --- high-resolution
#. ``E``: Elevation grid --- defined by interactions between ``A`` and
   ``I``.

To define global elevation classes, ``global_ec`` requires the
following input:

#. ``<mask filename>``: A global map of ice cover and elevation (on the
   ice grid), with each grid cell either being entirely covered by ice
   or entirely ice-free.  It can be produced with the etopo1_ice_
   command.

``global_ec`` produces a single file containing the six regridding
matrices between ``A``, ``I`` and ``E``.

Command Line
------------

The ``global_ec`` command takes the following arguments:

.. code-block:: bash

   $ global_ec <atm. grid name> <ice grid name> <display ice grid name>
               <mask filename> <focean filename>
               -o <output filename>

For example:

.. code-block:: bash

   $ global_ec g1qx1 g1mx1m ghxh etopo1_ice_g1m.nc topoo.nc -o global_ec.nc

.. note::

   The ``<focean filename>`` argument is only used if *mismatched
   regridding* is turned on.  Otherwise, any filler is fine: ``foo``,
   ``_``, etc.


Grid Names
^^^^^^^^^^

The grid names are drawn from a list of pre-defined grids that
``global_ec`` is able to use (others may be added in the future).
They are:

* ``g1mx1m``: 1-minute (eg: ETOPO1)
* ``g2mx2m``: 2-minute (eg: ETOPO2)
* ``g10mx10m``: 10-minute
* ``ghxh``: :math:`\frac{1}{2} \times \frac{1}{2}` degree
* ``g1x1``: 1 degee
* ``g1qx1``: :math:`1 \frac{1}{4} \times 1` degree (eg: ModelE ocean)
* ``g2hx2``: :math:`2 \frac{1}{2} \times 2` degree (eg: ModelE atmosphere)
* ``g5x4`` : :math:`5 \times 4` degree (eg: Model3)

Elevations
^^^^^^^^^^

Elevation classes are defined at a discrete set of elevations; by
default, they go from :math:`-100 m` to :math:`5100 m`, every
:math:`200 m`.  This is specified on the command line with:

.. code-block:: none
   --elev-classes -100,5100,200

.. note::

   It is important that the highest elevation class is *higher* than
   the highest ice in your input file.

Radius of Earth
^^^^^^^^^^^^^^^

By default, the radius of the Earth is set to be the same as in ModelE
(:math:`6.371 \times 10^6 m`).  This can be changed with the ``-R`` flag.


NetCDF Variable Names
^^^^^^^^^^^^^^^^^^^^^

By default, ``global_ec`` looks for the variables ``FGICE1m`` and
``ZICETOP1m`` from the mask file.  ``global_ec`` can be directed to
look for these quantities under different variable names, using the
command line arguments ``--elev`` and ``--mask``.

Mismatched Regridding
---------------------

Due to properties of the ModelE ocean, the ice extent "seen" by ModelE
might not be the same as is provided by the high-resolution ice grid.
This is because ModelE must "round" ocean grid cells to be all ocean
or all non-ocean; and also because ice extent might chnage over time,
whereas the ModelE ocean cannot.

The ``--mismatched`` flag enables mismatched regridding for
``global_ec``.  If it is turned on, the behavior of ``global_ec`` changes:

#. Instead of providing the atmosphere grid on the command line, the
   user must provide the ocean grid, which will be 1/2 the resolution
   of the desired atmosphere grid.  All core calculations are done on
   the ocean grid, and then regridded to the atmosphere grid at the
   last step.

#. The ``<focean filename>`` parameter is required: A map of which
   ocean grid cells are used by the ModelE ocean (ocean grid), as
   produced by the make_topoo_ command.

Logistical Issues
-----------------

The computations involved in ``global_ec`` can require large amount of
RAM, depending on the size of the ice grid and the number of
ice-covered grid cells.  For example, the ETOPO1 grid contains 233
million grid cells, of which 26 million are covered in ice.
``global_ec`` is not able to run with this much ice, on a typical
computer with 16Gb RAM.

Therefore, ``global_ec`` works by splitting up the task into smaller
sub-matrices covering portions of the globe; and then combining the
matrices at the end.  The domain decomposition process is independent
of ice location or extent, and will work for exoplanets too.

When ``global_ec`` is run by the user, it examines the input files to
plan its subtasks.  It then writes the plan into a Makefile called
``<output>.nc.mk``.  For example, the following command line actually
produces the file ``global_ec.nc.mk``:

.. code-block:: bash

   $ global_ec g1qx1 g1mx1m ghxh etopo1_ice_g1m.nc topoo.nc -o global_ec.nc

The computation is then completed with the command:

.. code-block:: bash

   $ make -f global_ec.nc.mk

If you are using *make* to automate your workflow, the following
placed in your *makefile* will be sufficient:

.. code-block:: make

   # The six regridding matrices (compressed) for global ice
   global_ec.nc.mk : topoo.nc etopo1_ice_g1m.nc
       global_ec g1qx1 g1mx1m ghxh etopo1_ice_g1m.nc topoo.nc -o global_ec.nc

   global_ec.nc : global_ec.nc.mk
        make -f global_ec.nc.mk

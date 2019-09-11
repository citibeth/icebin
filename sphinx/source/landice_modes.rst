.. _landice_modes

Landice Modes
=============

ModelE (*e3/landice* branch) can run with three approaches to landice
(ice sheets): *Classic*, *Stieglitz* or *Coupled*:

* *Classic* (default): How ModelE has been running by default since
  the 1980's.  The top 3m of an ice sheet are modeled using a simple
  2-layer model built into ModelE [Alexander et al 2018].  If ice
  accumulates over 3m, it is "pushed down" into the *mdwnimp* and
  *edwnimp* variables; if it melts to less than 3m, it is "pulled up"
  from the same.

* *Stieglitz*: The 2-layer model was replaced with a Lynch-Stieglitz
  model, with 5 layers by default.  Otherwise, the "push down" and
  "pull up" behavior remains unchanged.

* *Coupled*: The Lynch-Stieglitz model is used for snow/firn, but it
  is placed on top of a dynamic ice sheet.  The ice sheet is modeled
  on a separate grid with a separate code (eg PISM), and *IceBin* is
  used to couple between the two.  In areas of the globe that are
  covered in ice but do not have a dynamic ice model, this mode is the
  same as *Stiegtliz*.

Instructions on using these different modes follow:

Classic Snow/Firn Model
-----------------------

The classic snow/firn model is enabled by default.  Nothing needs to
be done to the rundeck to enable this.  Note that Classic can be used
with or without elevation classes.

.. note::
   ``modele_control`` is required for:
   #. Code generation to dump constants

File Format: TOPO
`````````````````

The ``TOPO`` file provides land surface fraction and elevation class
information to ModelE; see `here
<http://twoway.readthedocs.io/en/latest/topo.html>` for more
information.  The number of elevation classes is picked up
automatically; if a classic ``TOPO`` file without elevation class
information is used, then ModelE will rever to running without
elevation classes.


File Format: GIC
````````````````

Variables ``snowli`` and ``tlandi`` in the ``GIC`` file may have an
extra ``nhc`` (elevation class) dimension, but they don't have to.  If
they do, the ``nhc`` from ``GIC`` must match that from ``TOPO``.
Otherwise, the values of ``snowli`` and ``tlandi`` are repeated for
every elevation class within each grid cell.



Stieglitz Snow/Firn Model
-------------------------

(Instructions not yet tested)

The Stieglitz snow/firn model may be enabled by making the following
changes to the ModelE rundeck:

#. Add ``lipluggable`` to the *Components* section of the rundeck.

#. Build with `modele_control
   <https://github.com/citibeth/modele-control>`.  If not, one must also:

#. The following rundeck parameters are now enabled.  See
   ``LISnowParams.F90`` for documentation:

   * ``lisnow_target_depth``
   * ``lisnow_rho_fresh_snow``
   * ``lisnow_rho_snow_firn_cutoff``
   * ``lisnow_max_fract_water``
   * ``lisnow_epsilon``
   * ``lisnow_min_snow_thickness``
   * ``lisnow_min_fract_cover``
   * ``lisnow_dump_forcing``
   * ``lisnow_percolate_nl``

File Format: TOPO
`````````````````

The ``TOPO`` file is the same for Stieglitz vs. Classic ice model.  As
with Classic, ModelE infers elevation classes for a run from the
``TOPO`` file.


File Format: GIC
````````````````

Stieglitz has a different internal state from the classic ice model.  Therefore, the variables it reads form ``GIC`` are different.  The variables ``snowli`` and ``tlandi`` are no longer needed.  They are replaced by (NetCDF index ordering):

* ``dz(nhc,j,i,nlice)`` [*m*]: Depth of each layer
* ``wsn(nhc,j,i,nlice)`` [*kg m-2*]: Mass of each layer
* ``hsn(nhc,j,i,nlice)`` [*J m-2*]: Enthalpy of each layer
* ``trsn(nhc,j,i,nlice,ntm)`` (OPTIONAL) Tracer state

.. note::

   #. Specific enthalpy [*J kg-1*] for ``hsn`` would have been more
      convenient than enthalpy because it corresponds directly to
      temperature.  That change might happen in the future.

   #. As with the Classic model, the ``nhc`` (elevation class)
      dimension is optional.  If it is missing, values for these
      variables provided by ``GIC`` are read into all elevation
      classes.


Coupled Mode
------------

(Instructions not yet tested)

Do everything needed for Stieglitz Model above, plus:

#. Build with ``modele-control``.

#. When setting up the Spack environment, make sure that the variants
   ``icebin+coupler+modele+pism`` and ``modele+icebin`` are used.

#. Add ``#define USE_ICEBIN`` to the *defines* section.  This enables
   coupled mode, but does not turn it on.

#. Set ``LIMODE=twoway`` in rundeck.

#. Spin up a PISM.  Note the name of its state file ``<pism-state>``

#. Createa ModelE run directory:

   .. code-block:: bash

      ectl setup <run-dir> ...

#. Somewhere do ``git clone https://github.com/citibeth/twoway.git``

#. In your ModelE run directory, do:

   .. code-block:: bash

      cd <modele-run-dir>
      # 20km grid
      python3 <twoway>/topo/modele_pism_inputs.py --out input --pism <pism-state>
      cp input/icebin.cdl config

   Try ``python3 <twoway>/topo/modele_pism_inputs.py --help`` for further
   details if you need something other than the SeaRise-style 20km
   grid.  (TODO: This program should really determine the grind from the PISM state file)

#. Edit ``config/icebin.cdl`` as appropriate.


(Where is this found?  How do I create one????)


netcdf OZ1QX1N.BS1 {
dimensions:
        lono = 288 ;
        lato = 180 ;
variables:
        float lono(lono) ;
                lono:units = "degrees_east" ;
        float lato(lato) ;
                lato:units = "degrees_north" ;
        float focean(lato, lono) ;
        float zatmo(lato, lono) ;
        float zocean(lato, lono) ;
}


Need to create TOPO_OC file for ModelE input from topoo.nc:
 1. Convert FOCEAN, ZATMO and dZOCEAN into focean, zatmo and zocean.
 2. 

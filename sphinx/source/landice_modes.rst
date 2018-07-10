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

Stieglitz Snow/Firn Model
-------------------------

(Instructions not yet tested)

The Stieglitz snow/firn model may be enabled by making the following
changes to the ModelE rundeck:

#. Add ``#define LIPLUGGABLE`` to the *defines* section.

#. Build with `modele_control
   <https://github.com/citibeth/modele-control>`.  This is required
   for:

   #. Adding the ``lipluggable`` code derectory when ``LIPLUGGABLE``
      is defined.

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

#. Spin up an PISM.

#. Use an appropriate *IceBin* configuration file that links to the
   *PISM* configuration file.


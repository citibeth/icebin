.. _icebin_regridding

IceBin Regridding
=================

*IceBin* uses the techniques of conservative_regridding_ to produce
regridding matrices for GCMs employing *elevation classes*.  *IceBin*
requires three grids to be defined:

* ``A`` = *atmosphere grid*, used by the GCM.
* ``I`` = *ice grid*, a high-resolution grid for an ice-covered region
  of the globe.  It may be defined in spherical or Cartesian geometry;
  if Cartesian, then overlap matrices will be computed with the
  atmosphere grid by projecting the atmosphere grid to the same
  Cartesian space.  Every grid cell in the ice grid has a single
  elevation; and must be either 100% ice-covered or 100% ice free.  If
  the GCM is coupled with a dynamic ice model, the ice model is run on
  this grid.
* ``E`` = *elevation grid*, an intermediate-resolution grid, derived
  from the atmosphere and ice grids.  The land/ice surface model is
  run on this grid.  See `Fischer and Nowicki, 2014
  <https://pubs.giss.nasa.gov/abs/fi03200y.html>` for more information
  on how this grid, and its basis functions, are derived and used.

Give these three grids, *IceBin* generates the six regridding matrices
between them, scaled or unscaled.


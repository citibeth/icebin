icebin22m:
    Puts a PISM ice sheet on the 2-minute lon/lat grid.  This was
    thought to be useful at one point, but is no longer needed.

giss2nc:
    Convert legacy GISS-format input files to NetCDF format.

make_topo_icebin:
    Adds ice sheet onto Greenland-free TOPO file.  This can be seen as
    a prototype of the GCMCoupler_ModelE::update_topo() that's now in
    IceBin.

make_topo:
    Modification of Gary's TOPO generator.

# Simplest possible ModelE data plotting demo

import netCDF4
import giss.basemap
import giss.modele
import sys

nc_name = sys.argv[1]
var_name = sys.argv[2]

basemap = giss.basemap.greenland_laea()

nc = netCDF4.Dataset(nc_name)
pp = giss.modele.plot_params(var_name, nc=nc)
giss.plot.plot_var(basemap=basemap, **pp)		# Plot, and show on screen

# Slightly more complex alternatives:
# Save figure:
# 	giss.plot.plot_var(fname='plottest1.png', **pp)
# Save figure and snow on screen
# 	giss.plot.plot_var(fname='plottest1.png', show=True, **pp)

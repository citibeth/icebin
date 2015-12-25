from icebin import ibgrid
from matplotlib import pyplot
import giss.basemap

# Loads an plots grids from an IceBin grid file
# See: https://github.com/citibob/icebin

figure = pyplot.figure(figsize=(11.,8.5))

# -------------- Plot ice grid
grid = ibgrid.Grid('/Users/rpfische/exp/151014-integrate/build/modele_ll_g2x2_5-searise_g20-40-PISM.nc', 'm.greenland.grid2')

ax = figure.add_subplot(121)
basemap = giss.basemap.greenland_laea()
grid.plot(basemap, ax=ax)

# -------------- Plot GCM grid
grid = ibgrid.Grid('/Users/rpfische/exp/151014-integrate/build/modele_ll_g2x2_5-searise_g20-40-PISM.nc', 'm.grid1')

ax = figure.add_subplot(122)
basemap = giss.basemap.greenland_laea()
grid.plot(basemap, ax=ax)

pyplot.show()


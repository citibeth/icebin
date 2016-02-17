import matplotlib
import netCDF4
import icebin
import sys

import giss.basemap

grid_fname = sys.argv[1]
vname = sys.argv[2]

# =============== Step 2: Plot grid outlines

# Set up the page
figure = matplotlib.pyplot.figure(figsize=(11,8.5))
# figure.set_size_inches(11., 8.5)# US Letter
figure.set_size_inches(8.267,11.692)# A4

# -------- First plot: grid1

# Read the grid
nc = netCDF4.Dataset(grid_fname + '.nc', 'r')
grid1 = icebin.Grid(nc, vname)
nc.close()

# Plot it!
ax = figure.add_subplot(111)
basemap = giss.basemap.greenland_laea(ax=ax)
grid1.plot(basemap)

# ------------ Render the page
figure.savefig(grid_fname + '.pdf', dpi=100, transparent=False)

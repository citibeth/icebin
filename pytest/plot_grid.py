import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import glint2
import sys

# Simple script to plot the polygons of a grid

fname = sys.argv[1]

nc = netCDF4.Dataset(fname)
grid = glint2.pyGrid(nc, 'grid')
nc.close()

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))
ax = figure.add_subplot(111)

mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
grid.plot(mymap)
mymap.drawcoastlines()

# Also show on screen
matplotlib.pyplot.show()

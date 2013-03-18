import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import glint2
import sys

# Usage: plot_grid <grid-file.nc> [var-name]
#    var-name default = 'grid'

# Simple script to plot the polygons of a grid
# m.greenland2.grid2

fname = sys.argv[1]

if len(sys.argv) > 2 :
	vname = sys.argv[2]
else :
	vname = 'grid'

nc = netCDF4.Dataset(fname)
grid = glint2.pyGrid(nc, vname)
nc.close()

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))
ax = figure.add_subplot(111)

mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
grid.plot(mymap, linewidth=.5)
mymap.drawcoastlines()

# Also show on screen
figure.savefig('x.pdf')
matplotlib.pyplot.show()

import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import glint2
import sys
import csv

# Simple script to plot the polygons of a grid
# m.greenland2.grid2

fname = sys.argv[1]


csvin = csv.reader(open(fname, 'r'))
xx = []
yy = []
for line in csvin :
	xx.append(float(line[1]))
	yy.append(float(line[2]))


# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))
ax = figure.add_subplot(111)

#mymap = giss.basemap.greenland_laea(ax)
##mymap = giss.basemap.north_laea(ax)
##grid.plot(mymap, linewidth=.5)
#mymap.drawcoastlines()

ax.plot(yy,xx, '.')

# Also show on screen
figure.savefig('x.ps')
matplotlib.pyplot.show()

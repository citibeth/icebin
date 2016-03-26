# GLINT2: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013 by Robert Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import glint2
import sys

# Simple script to plot the polygons of a grid
# m.greenland2.grid2

fname = sys.argv[1]

if len(sys.argv) > 2 :
    vname = sys.argv[2]
else :
    vname = 'grid'

nc = netCDF4.Dataset(fname)
grid = glint2.pyGrid(nc, vname)
grid1 = glint2.pyGrid(nc, 'm.grid1')
nc.close()

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))
ax = figure.add_subplot(111)

#mymap = giss.basemap.greenland_laea(ax)
mymap = giss.basemap.north_laea(ax)
grid.plot(mymap, linewidth=.5)
grid1.plot(mymap, linewidth=.5, color='r')
mymap.drawcoastlines()

# Also show on screen
figure.savefig('x.ps')
matplotlib.pyplot.show()

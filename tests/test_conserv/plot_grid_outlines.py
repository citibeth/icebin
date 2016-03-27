# IceBin: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

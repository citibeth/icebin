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


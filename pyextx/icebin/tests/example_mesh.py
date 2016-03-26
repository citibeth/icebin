# IceBin: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013-2016 by Elizabeth Fischer
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

from icebin import ibgrid, issm
from matplotlib import pyplot
import giss.basemap

# Create our own sample triangular mesh and other grids, and plot them
# around Greenland.

X0=-510000.
Y0=-2810000.
L = 1.e5
projection='+proj=stere +lon_0=-39 +lat_0=90 +lat_ts=71.0 +ellps=WGS84'
# ----------------------------------
# Atmopshere (GCM) grid

vertices = [ibgrid.Vertex(spec[0], X0+L*spec[1], Y0+L*spec[2]) for spec in [
    (0, -2., -2.),
    (1,  2., -2.),
    (2,  2.,  2.),
    (3, -2.,  2.)]]
vx = vertices
cells = [ibgrid.Cell(0, [vx[0], vx[1], vx[2], vx[3]], 0,0,0,0)]

gridA = ibgrid.Grid(vertices, cells,
    projection=projection,
    type='MESH',
    coordinates='XY')


# ----------------------------------
# Ice Mesh
vertices = [ibgrid.Vertex(spec[0], X0+L*spec[1], Y0+L*spec[2]) for spec in [
    (0, -1., -1.),
    (1,  1., -1.),
    (2,  1.,  1.),
    (3, -1.,  1.),
    (4,  0.,  0.)]]

vx = vertices
cells = [ibgrid.Cell(*spec) for spec in [
    (0, (vx[0], vx[1], vx[4]), 0,0,0,0),
    (1, (vx[1], vx[2], vx[4]), 0,0,0,0),
    (2, (vx[2], vx[3], vx[4]), 0,0,0,0),
    (3, (vx[3], vx[0], vx[4]), 0,0,0,0)]]

gridI = ibgrid.Grid(vertices, cells,
    projection=projection,
    type='L1',
    coordinates='XY')

# ----------------------------------
vx = vertices
cells = [ibgrid.Cell(*spec) for spec in [
    (0, (vx[0], vx[1], vx[4]), 0,0,0,0),
    (1, (vx[1], vx[2], vx[4]), 0,1,0,0),
    (2, (vx[2], vx[3], vx[4]), 0,2,0,0),
    (3, (vx[3], vx[0], vx[4]), 0,3,0,0)]]

gridX = ibgrid.Grid(vertices, cells,
    projection=projection,
    type='L1',
    coordinates='XY',
    grid1_ncells_full=gridA.cells_num_full,
    grid2_ncells_full=gridI.cells_num_full,
    grid1_nvertices_full=gridA.vertices_num_full,
    grid2_nvertices_full=gridI.vertices_num_full)



# ----------------------------------
figure = pyplot.figure(figsize=(8.5, 11))
ax = figure.add_subplot(111)
basemap = giss.basemap.greenland_laea()
gridA.plot(basemap=basemap, ax=ax, color='red')
gridX.plot(basemap=basemap, ax=ax, color='green')
gridI.plot(basemap=basemap, ax=ax, color='blue')

pyplot.show()

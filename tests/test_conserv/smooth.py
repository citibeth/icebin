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
from icebin import ibgrid, ibplotter
import sys
import rtree.index
import giss.basemap
import giss.modele
import pickle
import collections
import math
import matplotlib.pyplot
import scipy
import os
import numpy as np

ICEBIN_IN = 'icebin_in.nc'
grid_fname = 'sr_g20_pism'
ice_sheet = 'greenland'


# ========================================================


class Smoother(object):
    def __init__(self, inserts_iter, sigma):
        self.cells = collections.OrderedDict()
        self.rtree = rtree.index.Index()
        for tuple in inserts_iter:
            self.cells[tuple[0]] = tuple
            self.rtree.insert(*tuple)

        self.sigma = sigma


    def gaussian_weights(self, x0,y0,cell_mass0):
        """Gives the weights for a Gaussian centered at center=(x,y)"""
        two_sigma_squared_inv = 1./(2.*self.sigma*self.sigma)

#        index0,pos0,cell_mass0 = self.cells[index0]
#        x0 = pos0[0]
#        y0 = pos0[1]

        out_raw = list()
        mass_sum = 0.        # Sum of all the weights
        radius = 3.*self.sigma
        radius_squared = radius*radius
        for index in self.rtree.intersection((x0-radius, y0-radius, x0+radius, y0+radius)):
            index,pos,cell_mass = self.cells[index]

            x=pos[0]
            y=pos[1]

            # Theoretical weight for this cell's position
            # (don't worry about normalizing in theory, we will normalize ourselves)
            dx = (x-x0)
            dy = (y-y0)
            distance_squared = dx*dx + dy*dy
            if distance_squared < radius_squared:
                w0 = math.exp(-two_sigma_squared_inv * distance_squared)

                dmass = w0*cell_mass
                mass_sum += dmass
                out_raw.append((index, w0))

        # Normalize sum(w0*cell_mass) = cell_mass0
        factor = cell_mass0 / mass_sum
        out = list()
        for index,w in out_raw:
            out.append((index,w*factor))

        return out
# ==============================================================



# Get the regridding matrices
# (This tells us, among other thing, the mask and the weight of each gridcell)
mm = icebin.GCMRegridder(ICEBIN_IN)
rm = mm.regrid_matrices(ice_sheet)
with netCDF4.Dataset(ICEBIN_IN) as nc:
#    indexingA = ibgrid.Indexing(nc, 'm.gridA.indexing')
#    indexingHP = ibgrid.Indexing(nc, 'm.indexingHP')
    indexingI = ibgrid.Indexing(nc, 'm.{}.gridI.indexing'.format(ice_sheet))
    plotterI = ibplotter.read_nc(nc, 'm.{}.gridI'.format(ice_sheet))
#    plotterA = ibplotter.read_nc(nc, 'm.gridA')

IvE,wIvE = rm.regrid('IvE', scale=True)

wI = wIvE        # Numpy array
nI = len(wI)


# Read the grid
nc = netCDF4.Dataset(ICEBIN_IN, 'r')
pgridI = ibgrid.read_nc(nc, 'm.{}.gridI'.format(ice_sheet))
nc.close()


# Obtain a single point and weight for each basis functions.
# This is a simplification and approximate.  But it should work
# pretty well as long as grid cells are small.  And it will work
# for any kind of grid/mesh.
def inserts():
    n=0
    for cell in pgridI.cells.values():
        if wI[cell.index] > 0:
            centroid = cell.centroid()
            weight = wI[cell.index]
            yield  (cell.index, (centroid[0], centroid[1], centroid[0], centroid[1]), wI[cell.index])
            n += 1
#            if n == 1000:
#                break

# Construct sparse matrix...
M_data = list()
M_ii = list()
M_jj = list()

print('Making smoothing matrix...')
smoother = Smoother(inserts(), 60.*1000.)
n=0
for index0,pos0,mass0 in smoother.cells.values():
    gweights = smoother.gaussian_weights(pos0[0], pos0[1], mass0)
    for index1,weight in gweights:
        M_ii.append(index1)
        M_jj.append(index0)
        M_data.append(weight)
    n += 1
    if (n % 100) == 0:
        print('Finished {} cells'.format(n))

# This is our smoothing matrix!
M = scipy.sparse.coo_matrix((M_data, (M_ii, M_jj)), shape=(nI,nI))



# OK, now use it...
with netCDF4.Dataset(os.path.join(os.environ['HOME'], 'exp/151014-integrate/e4f40/pism_in.nc'), 'r') as nc:
    valI = nc.variables['glint2_massxfer_rate'][10,:,:] * 1000.

print('AA1')
valIs = M * valI.reshape(-1)
print('AA2')


print('Sums', np.sum(valI), np.sum(valIs))

figure = matplotlib.pyplot.figure(figsize=(11,8.5))

ax = figure.add_subplot(121)
basemap = giss.basemap.greenland_laea(ax=ax)
giss.plot.plot_var(ax=ax, basemap=basemap,
    **giss.modele.plot_params('I', val=valI, plotter=plotterI))

ax = figure.add_subplot(122)
basemap = giss.basemap.greenland_laea(ax=ax)
giss.plot.plot_var(ax=ax, basemap=basemap,
    **giss.modele.plot_params('I', val=valIs, plotter=plotterI))

print('Showing plot...')
matplotlib.pyplot.show()
# figure.savefig(fname, dpi=100, transparent=False)

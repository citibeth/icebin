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

import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import numpy as np

basemap = giss.basemap.greenland_laea()

fhc_nc = 'fhc-00.nc'

# Sum over all elevation classes
nc1 = netCDF4.Dataset(fhc_nc)
fhc1 = nc1.variables['fhc1h'][:]

# Should be zero everywhere
fhc1_sum = np.sum(fhc1, axis=(0)) - 1.0

pp = giss.modele.plot_params(val=fhc1_sum)
pp['title'] = 'fhc1_sum'
giss.plot.plot_var(basemap=basemap, fname='fhc1_sum.png', **pp)

# ====================================================

fall = nc1.variables['fgice1'][:] + nc1.variables['fgrnd1'][:] + nc1.variables['focean1'][:] + nc1.variables['flake1'][:] - 1
print 'F-all sum = %f' % np.sum(np.abs(fall))

pp = giss.modele.plot_params(val=fall)
pp['title'] = 'Total Surface Fraction-1 (should be 0)'
giss.plot.plot_var(basemap=basemap, fname='fall.png', **pp)


# ====================================================

for ihp in range(0,10) :

    nc0 = netCDF4.Dataset('GIC.144X90.DEC01.1.ext_hc.nc')
    fhc0 = nc0.variables['fhc'][ihp,:]

    nc1 = netCDF4.Dataset(fhc_nc)
    fhc1 = nc1.variables['fhc1h'][ihp,:]

    fhcdiff = fhc1 - fhc0
    # ----------------------------------------
    figure = matplotlib.pyplot.figure(figsize=(15,8.5))

    ax = figure.add_subplot(131)
    pp = giss.modele.plot_params('fhc', nc0, val=fhc0)
    giss.plot.plot_var(ax=ax, basemap=basemap, **pp)        # Plot, and show on screen


    ax = figure.add_subplot(132)
    pp = giss.modele.plot_params('fhc1h', nc1, val=fhc1)
    #pp = giss.modele.plot_params(val=fhcdiff)
    #pp['title'] = 'fhcdiff'
    giss.plot.plot_var(ax=ax, basemap=basemap, **pp)        # Plot, and show on screen

    ax = figure.add_subplot(133)
    pp = giss.modele.plot_params(val=fhcdiff * 1000)
    pp['title'] = 'fhcdiff (*1000)'
    giss.plot.plot_var(ax=ax, basemap=basemap, **pp)        # Plot, and show on screen

    figure.savefig('fhc-%d.png' % ihp)
#   matplotlib.pyplot.show()




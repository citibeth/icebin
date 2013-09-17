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

nc0 = netCDF4.Dataset('Z2HX2fromZ1QX1N_hc.nc')
fgice0 = nc0.variables['fgice'][:]

nc1 = netCDF4.Dataset('fhc.nc')
fgice1 = nc1.variables['fgice1'][:]

fgicediff = fgice1 - fgice0
fgicediff[fgicediff < -.004] = 0
# ----------------------------------------
basemap = giss.basemap.greenland_laea()
figure = matplotlib.pyplot.figure(figsize=(15,8.5))

ax = figure.add_subplot(131)
pp = giss.modele.plot_params('fgice', nc0)
giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen


ax = figure.add_subplot(132)
pp = giss.modele.plot_params('fgice1', nc1)
#pp = giss.modele.plot_params(val=fgicediff)
#pp['title'] = 'fgicediff'
giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen

ax = figure.add_subplot(133)
pp = giss.modele.plot_params(val=fgicediff * 1000)
pp['title'] = 'fgicediff (*1000)'
giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen


matplotlib.pyplot.show()

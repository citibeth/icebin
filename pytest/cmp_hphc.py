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
import numpy as np

basemap = giss.basemap.greenland_laea()

#fhc_nc = 'fhc-01.nc'

for rank in range(0,2) :
	fhc_nc = 'fhc-%02d.nc' % rank


	# Read a height-points field so we can play with it
	nc = netCDF4.Dataset('JUL1950.ijhchc1k225.nc')
	val1hp = nc.variables['impm_lndice'][:]		# C-style
	val1hp[np.abs(val1hp)>1e10] = 0
	nc.close()

	#val1hp[:] = 1.

	# Read the same thing on ice model
	nc = netCDF4.Dataset('dismal.nc')
	val2 = nc.variables['mass'][:]
	print nc.dimensions['nx']
	nx = len(nc.dimensions['nx'])
	ny = len(nc.dimensions['ny'])
	nhc = 40
	nc.close()

	# -----------------------------------------------------
	# Get the hp_to_hc matrix
	nc1 = netCDF4.Dataset(fhc_nc)
	rows_i = nc1.variables['rows_i'][:]
	rows_j = nc1.variables['rows_j'][:]
	rows_k = nc1.variables['rows_k'][:]
	cols_i = nc1.variables['cols_i'][:]
	cols_j = nc1.variables['cols_j'][:]
	cols_k = nc1.variables['cols_k'][:]
	vals = nc1.variables['vals'][:]
	nc1.close()

	val1hc = np.zeros(val1hp.shape)	# C-style

	# Multiply hp_to_hc * val1hp
	for i in range(0,len(rows_i)-1) :
	#	print i,vals[i],val1hp[cols_k[i]-1, cols_j[i]-1, cols_i[i]-1]
		val1hc[rows_k[i]-1, rows_j[i]-1, rows_i[i]-1] = \
			val1hc[rows_k[i]-1, rows_j[i]-1, rows_i[i]-1] + \
			vals[i] * val1hp[cols_k[i]-1, cols_j[i]-1, cols_i[i]-1]

	#print val1hp
	#print np.sum(np.abs(val1hp))

	#val1hc = val1hp

	# -------------------
	ncout = netCDF4.Dataset('val1hc-%02d.nc' % rank, 'w')
	ncout.createDimension('nhc', val1hc.shape[0])
	ncout.createDimension('jm', val1hc.shape[1])
	ncout.createDimension('im', val1hc.shape[2])
	val1hc_v = ncout.createVariable('val1hc', 'd', ('nhc', 'jm', 'im'))
	val1hc_v[:] = val1hc
	val1hp_v = ncout.createVariable('val1hp', 'd', ('nhc', 'jm', 'im'))
	val1hp_v[:] = val1hp
	ncout.close()
	# -------------------


	# -----------------------------------------------------
	# Compute height-classified answer we think we should get
	# (see unittests.py)
	# val1hp_0



# -----------------------------------------------------

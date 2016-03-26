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
import giss

nc = netCDF4.Dataset('hp2hc.nc')
hp2hc = giss.read_coo_matrix(nc, 'hp2hc')
nc.close()

#print hp2hc.__dict__

nc = netCDF4.Dataset('fhc-00.nc')
jm = len(nc.dimensions['jm'])
im = len(nc.dimensions['im'])
nhp = len(nc.dimensions['nhp'])
rows_i = nc.variables['rows_i'][:]
rows_j = nc.variables['rows_j'][:]
rows_k = nc.variables['rows_k'][:]
cols_i = nc.variables['cols_i'][:]
cols_j = nc.variables['cols_j'][:]
cols_k = nc.variables['cols_k'][:]
vals = nc.variables['vals'][:]
nc.close()

nele = len(rows_i)

for ix in range(0,nele) :
    vdiff = vals[ix] - hp2hc.data[ix]
    print '(%d - %d) --> ([%d %d %d] - [%d %d %d]) %f' % (hp2hc.row[ix], hp2hc.col[ix], rows_i[ix], rows_j[ix], rows_k[ix], cols_i[ix], cols_j[ix], cols_k[ix], vals[ix])

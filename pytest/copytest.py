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

import giss.ncutil
import netCDF4

nc0 = netCDF4.Dataset('mm.nc', 'r')
ncout = netCDF4.Dataset('mmx.nc', 'w')

cp = giss.ncutil.copy_nc(nc0, ncout)


# Define our own variables
ncout.variables['m.greenland.info'].ice_model = 'DISMAL'
#ncout.variables['m.greenland2.info'].ice_model = 'DISMAL'

cp.copy_data()

ncout.close()
nc0.close()

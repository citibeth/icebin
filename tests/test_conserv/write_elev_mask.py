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

import netCDF4
import sys
import giss.ncutil

# Equivalent to ncks -v

ifname = sys.argv[1]
ofname = sys.argv[2]

fin = netCDF4.Dataset(ifname)
fout = netCDF4.Dataset(ofname, 'w')

var_filter = lambda vname: vname if vname in {'topg', 'thk', 'mask'} else None

copy = giss.ncutil.copy_nc(fin, fout, var_filter=var_filter, attrib_filter=lambda x: True)
copy.define_vars()
copy.copy_data()

fout.close()
fin.close()

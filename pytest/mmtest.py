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

import glint2
import giss.searise
import numpy as np
import sys

hpdefs = np.array(range(0,40))*100.0 - 50.0
#hcmax = np.array(range(0,40))*100.0

#hpdefs = np.array([100,300,550, 850, 1150, 1450, 1800, 2250, 2750, 4000], dtype='d')
#hcmax = np.array([200,400,700,1000,1300,1600,2000,2500,3000,1000000], dtype='d')

mm = glint2.MatrixMaker()
mm.init('../greenland_2x2_5.nc', hpdefs) #, hcmax)  # no mask1


searise_fname = '../data/Greenland_5km_v1.1.nc'
(elev2, mask2) = giss.searise.read_elevation2_mask2(searise_fname)
greenland_id = mm.add_ice_sheet('../searise.nc', '../greenland_2x2_5-searise.nc',
    elev2, mask2=mask2, name='greenland')
#greenland2_id = mm.add_ice_sheet('../searise.nc', '../greenland_2x2_5-searise.nc',
#   elev2, mask2=mask2, name='greenland2')

print 'greenland_id = %d' % greenland_id
#print 'greenland2_id = %d' % greenland2_id

mm.realize()

print '***** Writing out MatrixMaker'
nc = glint2.NcFile('mm.nc', 'w')
mm.write(nc, 'm')
nc.close()


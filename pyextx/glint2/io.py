# pyGISS: GISS Python Library
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

import scipy.sparse

def read_overlap(nc, vname) :
	"""Reads the overlap matrix from the Exchange Grid file"""

	cells_ijk = nc.variables[vname + '.cells.ijk'][:]
	rows = cells_ijk[:,0]
	cols = cells_ijk[:,1]

	info_var = nc.variables[vname + '.info']
	cells_area = nc.variables[vname + '.cells.area'][:]

	shape = (info_var.__dict__['grid1.ncells_full'],
		info_var.__dict__['grid2.ncells_full'])
	return scipy.sparse.coo_matrix((cells_area, (rows, cols)), shape=shape)

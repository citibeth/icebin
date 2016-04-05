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

# Creates the base of a ICEBIN input file (a MatrixMaker output file)
# We will append stuff to this, based on the ice model used (PISM, etc).

import icebin
import ibmisc

import giss.pism
import numpy as np
import sys
import giss.ncutil
import netCDF4
import os
import os.path
import sys
import pickle

ICEBIN_IN = 'icebin_in.nc'
sheet_names = ['greenland']


mm = icebin.GCMRegridder(ICEBIN_IN)

# ========= Compute all regridding matrices for Python use
# Then store them in a Python Pickle-format file

matrices = dict()
for sheet_name in sheet_names:
    rm = mm.regrid_matrices(sheet_name)
    for mat_type in ('EvI', 'AvI', 'IvA', 'IvE', 'EvA', 'AvE'):
        key = (sheet_name, mat_type)
        print('-------- Computing', key)
        sys.stdout.flush()
        matrix,weights = rm.regrid(mat_type)
        matrices[key] = (matrix,weights)

with open('matrices.pik', 'wb') as fout:
    pickle.dump(matrices, fout)

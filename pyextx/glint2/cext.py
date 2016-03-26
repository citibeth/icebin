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

import _glint2
import scipy.sparse
import numpy as np
import giss.util

# Interface with C++ extension

# Pass-through class to glint2 module
NcFile = _glint2.NcFile
Grid = _glint2.Grid


# -----------------------------------------------------
@giss.util.inherit_docs
class MatrixMaker(_glint2.MatrixMaker) :
    """User-level python interface to MatrixMaker object."""

    # Used to load from existing GLINT2 config file
    def __init__(self, fname=None, vname='m', **kwds) :

        """Constructor, to be used either for loading existing Glint2
        config file, or for creating a new one from scratch.

        To create a new MatrixMaker from scratch:
            Use no arguments.

        To load existing config file, keyword arguments are:
            fname : str
                Name of Glint2 configuration file>
            vname (OPTIONAL):
                Name of variable to read in config file
            correct_area1 : int/bool (DEFAULT 1)
                if non-zero: account for projection error in transformations.
                For some reason, this is not stored in the Glint2 config file."""

        super(MatrixMaker, self).__init__(**kwds)
        self.fname = fname
        if fname is not None :
            # Used to load from GLINT2 config file
            super(MatrixMaker, self).load(fname, vname)
            super(MatrixMaker, self).realize()

    def init(self, *args, **kwargs) :
        if 'mask1' in kwargs :
            kwargs['mask1'] = giss.util.reshape_no_copy(kwargs['mask1'])
        super(MatrixMaker, self).init(*args, **kwargs)

    def add_ice_sheet(self, grid2_fname, exgrid_fname, elev2, **kwargs) :
        """When constructing a MatrixMaker via init() method, use to
        add an ice sheet to the system.

        grid2_fname : str
            Name of the grid file for the ice sheet.
        exgrid_fname : str
            Name of the grid file for the exchange grid and the GCM grid
            this MatrixMaker was instantiated with.
        elev2 : double[n2]
            Elevation of each ice grid cell
        mask2 : int[n2] (bool) (OPTIONAL)
            If non-zero, marks ice grid cells as unused
        name : str (OPTIONAL)
            Name used to identify this ice sheet later."""

        elev2 = giss.util.reshape_no_copy(elev2, -1)
        if 'mask2' in kwargs :
            kwargs['mask2'] = kwargs['mask2'].reshape(-1,)
        super(MatrixMaker, self).add_ice_sheet(grid2_fname, exgrid_fname, elev2, **kwargs)

    def hp_to_iceinterp(self, *args, **kwargs) :

        """Returns the regridding matrix to go from the elevation grid
        to an ice grid (E->I).

        sheetname : str
            The ice sheet to which to regrid.
        dest : glint2::IceInterp (str) (DEFAULT 'ICE')
            Specify whether the destination should be the ice grid or
            the interpolation grid (which might be the ice grid)
        returns : scipy.sparse.coo_matrix
        """
        if 'fill_masked' in kwargs:
            fill_masked = 1 if kwargs['fill_masked'] else 0
        else:
            fill_masked = 0
        tret = super(MatrixMaker, self).hp_to_iceinterp(*args, fill_masked=fill_masked)
        return _tuple_to_coo(tret)

    def hp_to_atm(self, *args) :
        """Returns the matrix E->A (elevation grid to atmosphere grid)
        Takes no arguments.  Returns scipy.sparse.coo_matrix""" 
        tret = super(MatrixMaker, self).hp_to_atm(*args)
        return _tuple_to_coo(tret)

    def iceinterp_to_atm(self, *args, **kwargs) :
        """Returns the matrix I->A, converting the ice to the atmosphere grid.

        Required arg:
        sheetname : str
            The ice sheet from which to regrid.

        Keyword args:
        src : glint2::IceInterp (OPTIONAL, DEFAULT 'ICE')
            {'ICE' | 'INTERP'}
            Specifies whether the source vector space is the ice grid or
            the interpolation grid (which might be the ice grid)
        returns : scipy.sparse.coo_matrix"""
        tret = super(MatrixMaker, self).iceinterp_to_atm(*args, **kwargs)
        return _tuple_to_coo(tret)

    def iceinterp_to_hp(self, f2s, *args, **kwargs):
        """Implements the I->E reverse transformation.

        f2s : dict(sheetname -> double[n2]
            Specifies the source field to regrid on each ice grid.
        initial3 : double[n3] (OPTIONAL)
            Initial guess at solution.
            (n3 = size of elevation grid)
        src : glint2::IceInterp (OPTIONAL, DEFAULT 'ICE')
            {'ICE' | 'INTERP'}
            Specifies whether the source vector space is the ice grid or
            the interpolation grid (which might be the ice grid)
        qp_algorithm : giss::QPAlgorithm (OPTIONAL, DEFAULT 'SINGLE_QP')
            The method used to create a QP program for this problem.
            'SINGLE_QP' :
                Create one single large QP program.
            'MULTI_QP' :
                Create a separate QP program for each atmosphere grid
                cell.  Only legal if the $RM$ matrix is local.
        """
        f2s_new = []
        for key, f2 in f2s.items() :
            f2s_new.append(  (key, giss.util.reshape_no_copy(f2, -1))  )

        nkwargs = kwargs.copy()
        if 'initial3' in nkwargs:
            nkwargs['initial3'] = nkwargs['initial3'].reshape(-1)

        nparray = super(MatrixMaker, self).iceinterp_to_hp(f2s_new, *args, **nkwargs)
        return nparray

    def ice_to_interp(self, sheetname, f2, *args, **kwargs):
        f2_1d = f2.reshape(-1)
        return super(MatrixMaker, self).ice_to_interp(sheetname, f2_1d, *args, **kwargs)

    def realize(self, *args) :
        super(MatrixMaker, self).realize(*args)

    def write(self, *args) :
        super(MatrixMaker, self).write(*args)
# -----------------------------------------------------

def _coo_to_tuple(coo) :
    return (coo._shape[0], coo._shape[1],
        coo.row, coo.col, coo.data)

def _tuple_to_coo(tuple) :
    nrow1 = tuple[0]
    ncol1 = tuple[1]
    rows1 = tuple[2]
    cols1 = tuple[3]
    data1 = tuple[4]
    return scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))

# -------------------------------------------------------
# Puts A*x into y, does not overwrite unused elements of y
# @param yy OUTPUT
def coo_matvec(coomat, xx, yy, ignore_nan=False) :
    yy = yy.reshape(-1)
    xx = xx.reshape(-1)

    _glint2.coo_matvec(_coo_to_tuple(coomat), xx, yy,
        ignore_nan=(1 if ignore_nan else 0))
    return

def coo_multiply(coomat, xx, fill=np.nan, ignore_nan=False) :
    xx = xx.reshape(-1)
    yy = np.zeros(coomat._shape[0])
    yy[:] = fill
    coo_matvec(coomat, xx, yy, ignore_nan)
    return yy


def multiply_bydiag(a1, a2) :
    if issubclass(type(a1), scipy.sparse.coo_matrix) :
        a1 = _coo_to_tuple(a1)
    else :
        a2 = _coo_to_tuple(a2)
    return _glint2.multiply_bydiag(a1, a2)

# --------------------------------------------------------

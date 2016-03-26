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

import giss.plot
import pyproj
import numpy as np
import giss.proj
import netCDF4
import giss.plot
import bisect

# Re-name variables for grid data files coming from ISSM
def standardize_names(ivars) :
    ovars = {}
    for vname,var in ivars.items() :
        print(vname)
        if vname != 'grid.info' and 'name' in var.__dict__ :
            ovars[var.name] = var
        else :
            ovars[vname] = var
    return ovars

# Local imports
import sys
#from overlap import *
from glint2.cext import *
glint2 = sys.modules[__name__]

class pyGrid(object) :
    """Utility class that reads a Grid file in pure python, rather
    than using C++.  Used mainly to plot grid outlines."""

    # Read Grid from a netCDF file
    def __init__(self, nc, vname) :
        """Read the Grid from a netCDF file.

            nc : netCDF4
                Open netCDF file handle.
            vname : str
                Name of variable from which to read the Grid."""
#       variables = standardize_names(nc.variables)
        variables = nc.variables

        info = variables[vname + '.info']

        # Read all attributes under .info:
        # name, type, cells.num_full, vertices.num_full
        self.__dict__.update(info.__dict__)
        self.cells_num_full = self.__dict__['cells.num_full']
        self.vertices_num_full = self.__dict__['vertices.num_full']

        self.vertices_index = variables[vname + '.vertices.index'][:]
        self.vertices_pos = dict()      # Position (by index) of each vertex in our arrays

        self.vertices_xy = variables[vname + '.vertices.xy'][:]
        self.cells_index = variables[vname + '.cells.index'][:]
        if vname + '.cells.ijk' in variables :
            self.cells_ijk = variables[vname + '.cells.ijk'][:]
        self.cells_area = variables[vname + '.cells.area'][:]
#       self.cells_proj_area = variables[vname + '.cells.proj_area'][:]
        self.cells_vertex_refs = variables[vname + '.cells.vertex_refs'][:]
        self.cells_vertex_refs_start = variables[vname + '.cells.vertex_refs_start'][:]

        # Compute vertices_pos
        for i in range(0, len(self.vertices_index)) :
            self.vertices_pos[self.vertices_index[i]] = i

        # Correct ISSM file
#       self.cells_vertex_refs = self.cells_vertex_refs - 1

        
        if self.coordinates == 'XY' :
            print('Projection = "' + self.projection + '"')
            print(type(self.projection))
            (self.llproj, self.xyproj) = giss.proj.make_projs(self.projection)
        else :
            self.xyproj = None

        if self.type == 'XY' :
            self.x_boundaries = variables[vname + '.x_boundaries'][:]
            self.y_boundaries = variables[vname + '.y_boundaries'][:]


    def plot(self, basemap, **kwargs) :
        """Plots the grid cell outlines
    
        Args:
            basemap: Map on which to plot
            **kwargs: Any options passed through to Matplotlib plotting
        """

        npoly = len(self.cells_vertex_refs_start)-1     # -1 for sentinel
        npoints = len(self.cells_vertex_refs)
#       xdata = np.zeros(npoints + npoly * 2)
#       ydata = np.zeros(npoints + npoly * 2)
        xdata = []
        ydata = []

        ipoint_dst = 0
        for ipoly in range(0,npoly) :
            iistart = self.cells_vertex_refs_start[ipoly]
            iinext = self.cells_vertex_refs_start[ipoly+1]
            npoints_this = iinext - iistart

#           print ipoint_dst, npoints_this, len(xdata)
#           print iistart, iinext, len(self.cells_vertex_refs)
#           print len(self.vertices_xy), self.cells_vertex_refs[iistart:iinext]

#           xdata[ipoint_dst:ipoint_dst + npoints_this] = \
#               self.vertices_xy[self.cells_vertex_refs[iistart:iinext], 0]
#           ydata[ipoint_dst:ipoint_dst + npoints_this] = \
#               self.vertices_xy[self.cells_vertex_refs[iistart:iinext], 1]

            refs = self.cells_vertex_refs[iistart:iinext]
            for i in range(0,len(refs)) :
                refs[i] = self.vertices_pos[refs[i]]
#           print refs
            xdata += list(self.vertices_xy[refs, 0])
            ydata += list(self.vertices_xy[refs, 1])

            ipoint_dst += npoints_this

            # Repeat the first point in the polygon
#           xdata[ipoint_dst] = \
#               self.vertices_xy[self.cells_vertex_refs[iistart], 0]
#           ydata[ipoint_dst] = \
#               self.vertices_xy[self.cells_vertex_refs[iistart], 1]

            xdata.append(self.vertices_xy[self.cells_vertex_refs[iistart], 0])
            ydata.append(self.vertices_xy[self.cells_vertex_refs[iistart], 1])


            ipoint_dst += 1

            # Add a NaN
#           xdata[ipoint_dst] = np.nan
#           ydata[ipoint_dst] = np.nan
            xdata.append(np.nan)
            ydata.append(np.nan)
            ipoint_dst += 1

        xdata = np.array(xdata)
        ydata = np.array(ydata)

        
        if self.xyproj is not None :    # translate xy->ll
            londata, latdata = pyproj.transform(self.xyproj, self.llproj, xdata, ydata)
            londata[np.isnan(xdata)] = np.nan
            latdata[np.isnan(ydata)] = np.nan

        else :      # Already in lon/lat coordinates
            londata = xdata
            latdata = ydata

        giss.basemap.plot_lines(basemap, londata, latdata, **kwargs)


    def plotter(self) :
        if self.type == 'XY' :
            return giss.plot.ProjXYPlotter(self.x_boundaries, self.y_boundaries, self.projection)
        return None


def _Grid_XY_read_plotter(nc, vname, transpose=False) :
    """Reads an plotter out of a netCDF file for a simple Cartesian grid"""

    # ======= Read our own grid2 info from the overlap file
    # Assumes an XY grid for grid2
    xb2 = nc.variables[vname + '.x_boundaries'][:]
    yb2 = nc.variables[vname + '.y_boundaries'][:]
    info_var = nc.variables[vname + '.info']
    sproj = info_var.projection
    return giss.plot.ProjXYPlotter(xb2, yb2, sproj, transpose=transpose)

def _Grid_LonLat_read_plotter(nc, vname, transpose=False) :
    lonb2 = nc.variables[vname + '.lon_boundaries'][:]
    latb2 = nc.variables[vname + '.lat_boundaries'][:]
    return giss.plot.LonLatPlotter(lonb2, latb2, True)

# -------------------------------
read_plotter_fn = {'XY' : _Grid_XY_read_plotter,
    'LONLAT' : _Grid_LonLat_read_plotter}

# Creates a plotter to plot data on an ice grid
# @param grid_nc Open netCDF file that has the ice grid
# @param vname Name of variable inside the netCDF file
# @param ice_sheet Name of ice sheet (works if variables follow std convention)
def Plotter2(nc=None, vname=None, fname=None, transpose=False) :
    if fname is not None :
        nc = netCDF4.Dataset(fname)
    stype = nc.variables[vname + '.info'].__dict__['type']
    read_fn = read_plotter_fn[stype]
    ret = read_fn(nc, vname, transpose=transpose)
    if fname is not None :
        nc.close()
    return ret

# ---------------------------------------------------
class Plotter1h(giss.plot.Plotter):
    # @param mmaker Instance of glint2.MatrixMaker
    # @param glint2_config Name of GLINT2 config file
    def __init__(self, glint2_config, ice_sheet, mmaker=None, **kwargs):
        self.ice_sheet = ice_sheet
        self.init_kwargs = kwargs

        if mmaker is None :
            self.glint2_config = glint2_config
            mmaker = glint2.MatrixMaker(glint2_config)
        else:
            # We'll try to reconstruct this mmaker on unpickle time
            self.glint2_config = mmaker.fname
        self.mat_1h_to_2 = mmaker.hp_to_iceinterp(ice_sheet, **kwargs)

        nc = netCDF4.Dataset(glint2_config)
        self.plotter2 = Plotter2(nc=nc, vname='m.' + ice_sheet + '.grid2')

        # Create a plotter1, to help determine coords when user clicks
        self.plotter1 = _Grid_LonLat_read_plotter(nc, 'm.grid1')

        # Used to determine nearest elevation point (for display)
        self.elev2 = nc.variables['m.'+ice_sheet+'.elev2'][:]
        self.elev2 = self.elev2.reshape(self.plotter2.ny2, self.plotter2.nx2)
        self.hpdefs = nc.variables['m.hpdefs'][:]

        # self.mask2 = nc.variables['m.' + ice_sheet + '.mask2'][:]
        nc.close()

    # For pickling
    def __getstate__(self):
        return ((self.glint2_config, self.ice_sheet), (self.init_kwargs))
    def __setstate__(self, state):
        self.__init__(*state[0], **state[1])


    def regrid2(self, val1h, mask=True):
        """mask: (default True)
            If set, use numpy masked_array for proper plotting of colorbar."""
        val1h = val1h.reshape(-1)
        val2 = glint2.coo_multiply(self.mat_1h_to_2, val1h, fill=np.nan, ignore_nan=False)  # Make np.nan
        if mask: mval2 = np.ma.masked_invalid(val2)
        return mval2


    def context(self, basemap, vals1h):
        context = giss.util.LazyDict()
        context['basemap'] = basemap
        context2 = self.plotter2.context(basemap, self.regrid2(vals1h))
        context['context2'] = context2
        return context      

    def plot(self, context, basemap_plot_fn, **plotargs):
        context2 = context['context2']
        return self.plotter2.plot(context2, basemap_plot_fn, **plotargs)

    def lookup(self, context, lon_d, lat_d):
        coords2,val2 = self.plotter2.lookup(context['context2'], lon_d, lat_d)

        # Find closest elevation point, based on our elevation on the ice grid.
        elev = self.elev2[coords2]
        ihp0 = bisect.bisect_left(self.hpdefs, elev)   # "rounds down"
        delta0 = abs(elev - self.hpdefs[ihp0])
        delta1 = abs(elev - self.hpdefs[ihp0+1])
        ihp = ihp0 if (delta0 <= delta1) else ihp0+1

        # Find enclosing grid cell on the GCM grid
        coords1 = self.plotter1.coords(lon_d, lat_d)

        return (ihp, coords1[0], coords1[1]), val2

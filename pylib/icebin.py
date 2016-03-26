# Try to import the Cython extension
try:
    from _icebin import *
except:
    pass

import giss.plot
import pyproj
import numpy as np
import giss.proj
import netCDF4
import bisect
import sys
import functools
import operator

class Indexing(object):
    def __init__(self, nc, vname):
        ncvar = nc.variables[vname]
        self.base = ncvar.base
        self.extent = ncvar.extent
        self.indices = ncvar.indices

        if not isinstance(self.base, np.ndarray):
            self.base = np.ndarray([self.base])
        if not isinstance(self.extent, np.ndarray):
            self.extent = np.ndarray([self.extent])
        if not isinstance(self.indices, np.ndarray):
            self.indices = np.ndarray([self.indices])

        self.size = 1
        for ex in self.extent:
            self.size *= ex

        # Shape of a row-major array in memory
        # Row-order ordered indices
        if (self.indices[0] == 0):
            self.shape = self.extent
        else:
            self.shape = (self.extent[1], self.extent[0])

    def __len__(self):
        return functools.reduce(operator.mul, self.extent)


class Grid(object):
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
        self.indexing = Indexing(nc, vname + '.indexing')

        # Read all attributes under .info:
        # name, type, cells.nfull, vertices.nfull
        self.__dict__.update(info.__dict__)
        for k in info.__dict__.keys():
            print(k)
        self.cells_nfull = self.__dict__['cells.nfull']
        self.vertices_nfull = self.__dict__['vertices.nfull']

        self.vertices_index = variables[vname + '.vertices.index'][:]
        self.vertices_pos = dict()      # Position (by index) of each vertex in our arrays

        self.vertices_xy = variables[vname + '.vertices.xy'][:]
        self.cells_index = variables[vname + '.cells.index'][:]
        if vname + '.cells.ijk' in variables :
            self.cells_ijk = variables[vname + '.cells.ijk'][:]
        self.cells_native_area = variables[vname + '.cells.native_area'][:]
#       self.cells_proj_area = variables[vname + '.cells.proj_area'][:]
        self.cells_vertex_refs = variables[vname + '.cells.vertex_refs'][:]
        self.cells_vertex_refs_start = variables[vname + '.cells.vertex_refs_start'][:]

        # Compute vertices_pos
        for i in range(0, len(self.vertices_index)) :
            self.vertices_pos[self.vertices_index[i]] = i
        
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
        xdata = []
        ydata = []

        ipoint_dst = 0
        for ipoly in range(0,npoly) :
            iistart = self.cells_vertex_refs_start[ipoly]
            iinext = self.cells_vertex_refs_start[ipoly+1]
            npoints_this = iinext - iistart

            refs = self.cells_vertex_refs[iistart:iinext]
            for i in range(0,len(refs)) :
                refs[i] = self.vertices_pos[refs[i]]
            xdata += list(self.vertices_xy[refs, 0])
            ydata += list(self.vertices_xy[refs, 1])

            ipoint_dst += npoints_this

            xdata.append(self.vertices_xy[self.cells_vertex_refs[iistart], 0])
            ydata.append(self.vertices_xy[self.cells_vertex_refs[iistart], 1])


            ipoint_dst += 1

            # Add a NaN
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


def _Grid_XY_read_plotter(nc, vname) :
    """Reads an plotter out of a netCDF file for a simple Cartesian grid"""

    sproj = nc.variables[vname + '.info'].projection
    xb2 = nc.variables[vname + '.x_boundaries'][:]
    yb2 = nc.variables[vname + '.y_boundaries'][:]
    indexing = Indexing(nc, vname + '.indexing')
    return giss.plot.ProjXYPlotter(xb2, yb2, sproj, indexing.indices[0] == 0)

def _Grid_LonLat_read_plotter(nc, vname) :
    lonb2 = nc.variables[vname + '.lon_boundaries'][:]
    latb2 = nc.variables[vname + '.lat_boundaries'][:]
    indexing = Indexing(nc, vname + '.indexing')
    return giss.plot.LonLatPlotter(lonb2, latb2, transpose=(indexing.indices[0] == 0), boundaries=True)

# -------------------------------
read_plotter_fn = {'XY' : _Grid_XY_read_plotter,
    'LONLAT' : _Grid_LonLat_read_plotter}

# Creates a plotter to plot data on an ice grid
# @param grid_nc Open netCDF file that has the ice grid
# @param vname Name of variable inside the netCDF file
# @param ice_sheet Name of ice sheet (works if variables follow std convention)
def read_plotter(nc=None, vname=None) :
    stype = nc.variables[vname + '.info'].__dict__['type']
    read_fn = read_plotter_fn[stype]
    ret = read_fn(nc, vname)
    return ret

# ---------------------------------------------------
class PlotterE(giss.plot.Plotter):
    # @param mmaker Instance of glint2.MatrixMaker
    # @param glint2_config Name of GLINT2 config file
    def __init__(self, icebin_config, ice_sheet, IvE=None, **kwargs):
        self.ice_sheet = ice_sheet
        self.init_kwargs = kwargs

        # Get the regridding matrix, if we don't already have it
        if not IvE:
            pass
        self.IvE = IvE

        with netCDF4.Dataset(glint2_config) as nc:
            self.plotterI = read_plotter(nc=nc, vname='m.' + ice_sheet + '.gridI')

            # Create a plotterA, to help determine coords when user clicks
            self.plotterA = _Grid_LonLat_read_plotter(nc, 'm.gridA')

            # Used to determine nearest elevation point (for display)
            self.elevI = nc.variables['m.'+ice_sheet+'.elevI'][:]
            self.elevI = self.elev2.reshape(self.plotterI.ny2, self.plotter2.nx2)
            self.hpdefs = nc.variables['m.hpdefs'][:]

            # self.mask2 = nc.variables['m.' + ice_sheet + '.mask2'][:]

    ## For pickling
    #def __getstate__(self):
    #   return ((self.glint2_config, self.ice_sheet), (self.init_kwargs))
    #def __setstate__(self, state):
    #   self.__init__(*state[0], **state[1])

    def regridI(self, valE, mask=True):
        """mask: (default True)
            If set, use numpy masked_array for proper plotting of colorbar."""
        valE = valE.reshape(-1)
        valI = np.zeros((self.IvE.shape[0],))
        valI[:] = np.nan
        valI = icebin.coo_multiply(self.IvE, valE, fill=np.nan, ignore_nan=False)   # Make np.nan
        if mask: mvalI = np.ma.masked_invalid(valI)
        return mvalI


    def context(self, basemap, valsE):
        context = giss.util.LazyDict()
        context['basemap'] = basemap
        contextI = self.plotterI.context(basemap, self.regridI(valsE))
        context['contextI'] = contextI
        return context      

    def plot(self, context, basemap_plot_fn, **plotargs):
        contextI = context['contextI']
        return self.plotterI.plot(contextI, basemap_plot_fn, **plotargs)

    def lookup(self, context, lon_d, lat_d):
        coordsI,valI = self.plotterI.lookup(context['contextI'], lon_d, lat_d)

        # Find closest elevation point, based on our elevation on the ice grid.
        elev = self.elevI[coordsI]
        ihp0 = bisect.bisect_left(self.hpdefs, elev)   # "rounds down"
        delta0 = abs(elev - self.hpdefs[ihp0])
        delta1 = abs(elev - self.hpdefs[ihp0+1])
        ihp = ihp0 if (delta0 <= delta1) else ihp0+1

        # Find enclosing grid cell on the GCM grid
        coords1 = self.plotterA.coords(lon_d, lat_d)

        return (ihp, coords1[0], coords1[1]), valI

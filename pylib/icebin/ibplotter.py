import giss.plot
import icebin
from icebin import ibgrid
import netCDF4
import numpy as np

"""Create plotters without reading the full grid."""


def _Grid_XY_read_plotter(nc, vname) :
    """Reads an plotter out of a netCDF file for a simple Cartesian grid"""

    sproj = nc.variables[vname + '.info'].sproj
    xb2 = nc.variables[vname + '.x_boundaries'][:]
    yb2 = nc.variables[vname + '.y_boundaries'][:]
    indexing = ibgrid.Indexing(nc, vname + '.indexing')
    return giss.plot.ProjXYPlotter(xb2, yb2, sproj, indexing.indices[0] == 0)

def _Grid_LonLat_read_plotter(nc, vname) :
    lonb2 = nc.variables[vname + '.lon_boundaries'][:]
    latb2 = nc.variables[vname + '.lat_boundaries'][:]
    indexing = ibgrid.Indexing(nc, vname + '.indexing')
    return giss.plot.LonLatPlotter(lonb2, latb2, transpose=(indexing.indices[0] == 0), boundaries=True)

# ---------------------------------------------------
read_plotter_fn = {'XY' : _Grid_XY_read_plotter,
    'LONLAT' : _Grid_LonLat_read_plotter}

def read_nc(nc, vname):
    """Creates a plotter to plot data on an ice grid
    @param grid_nc Open netCDF file that has the ice grid
    @param vname Name of variable inside the netCDF file
    @param ice_sheet Name of ice sheet (works if variables follow std convention)"""
    stype = nc.variables[vname + '.info'].__dict__['type']
    read_fn = read_plotter_fn[stype]
    ret = read_fn(nc, vname)
    return ret
# ---------------------------------------------------
class PlotterE(giss.plot.Plotter):
    def __init__(self, icebin_config, ice_sheet, elevmaskI=None, IvE=None, correctA=True):
        """icebin_config: str
            Name of IceBin configuration file
        ice_sheet: str
            Name of ice sheet to create plotter for
        elevmaskI:
            Elevation on ice grid; NaN if no ice
            OPTIONAL: Only needed for interactive plot (lookup())
        IvE: sparse matrix
            Matrix to convert [Ice Grid <-- Elevation Grid]"""
        self.ice_sheet = ice_sheet

        # Get the regridding matrix, if we don't already have it
        if IvE is None:
            mm = icebin.GCMRegridder(icebin_config)
            rm = mm.regrid_matrices(ice_sheet)
            IvE,wIvE = rm.regrid('IvE', scale=True)#, correctA=correctA)
        self.IvE = IvE

        with netCDF4.Dataset(icebin_config) as nc:
            self.plotterI = read_nc(nc=nc, vname='m.' + ice_sheet + '.agridI')

            # Create a plotterA, to help determine coords when user clicks
            self.plotterA = _Grid_LonLat_read_plotter(nc, 'm.agridA')

            # Used to determine nearest elevation point (for display)
            # (convert sparse to dense vector)
            #elevI = nc.variables['m.'+ice_sheet+'.elevI'][:]

            self.elevmaskI = elevmaskI
            if elevmaskI is not None:
                self.elevmaskI = elevmaskI.reshape(self.plotterI.ny2, self.plotterI.nx2)
            self.hcdefs = nc.variables['m.hcdefs'][:]

            # self.mask2 = nc.variables['m.' + ice_sheet + '.mask2'][:]

    ## For pickling
    #def __getstate__(self):
    #   return ((self.icebin_config, self.ice_sheet), (self.init_kwargs))
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
        context = giss.giutil.LazyDict()
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
        elev = self.elevmaskI[coordsI] if self.elevmaskI is not None else 0.
        ihp0 = bisect.bisect_left(self.hcdefs, elev)   # "rounds down"
        delta0 = abs(elev - self.hcdefs[ihp0])
        delta1 = abs(elev - self.hcdefs[ihp0+1])
        ihp = ihp0 if (delta0 <= delta1) else ihp0+1

        # Find enclosing grid cell on the GCM grid
        coords1 = self.plotterA.coords(lon_d, lat_d)

        return (ihp, coords1[0], coords1[1]), valI


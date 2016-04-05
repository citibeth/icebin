import icebin
from icebin import ibgrid, element_l1, ibplotter
from matplotlib import pyplot
import giss.basemap
import unittest
import numpy as np
import netCDF4
import sys
import scipy.sparse
import matplotlib
import giss.modele

SEARISE_GRID = 'sr_g20_pism.nc'
OVERLAP_FILE = 'sr_g20_pism-ISSM_mesh.nc'

class TestISSMRegrid(unittest.TestCase):

    def setUp(self):
        with netCDF4.Dataset(SEARISE_GRID, 'r') as nc:
            self.plotterA = ibplotter.read_nc(nc, 'grid')
            gridA = ibgrid.read_nc(nc, 'grid')
            self.indexingA = ibgrid.Indexing(nc, 'grid.indexing')

        with netCDF4.Dataset(OVERLAP_FILE, 'r') as nc:
            gridI = ibgrid.read_nc(nc, 'gridI')
            exgrid = ibgrid.read_nc(nc, 'exgrid')


        nA = gridA.cells_nfull
        nI = gridI.vertices_nfull
        AvI,weightsA,weightsI = element_l1.compute_AvI(exgrid, nA, gridI)

        # Diagonal scale matrices based on regrid weights
        scaleA = scipy.sparse.dia_matrix( ([1. / weightsA], [0]), shape=(nA,nA))
        scaleI = scipy.sparse.dia_matrix( ([1. / weightsI], [0]), shape=(nI,nI))

        self.AvI = (scaleA * AvI).tocoo()
        self.IvA = (scaleI * AvI.transpose()).tocoo()


    def diagonalA(self):
        """Sample function"""
        shapeA = self.indexingA.shape
        valA = np.zeros(shapeA)
        for i in range(0,shapeA[0]):
            for j in range(0,shapeA[1]):
                valA[i,j] = i + j
        return valA.reshape(-1)

    def constantA(self):
        """Sample function"""
        valA = np.zeros(self.indexingA.shape) + 1
        return valA.reshape(-1)


    def test_regrids(self):
        """Does a sample regrid A <- I <- A, and plots in A-space.
        A = SeaRISE Grid
        I = ISSM Grid
        I don't know how to plot ISSM meshes directly."""

        for fname,valA_fn in [('test_diagonal.png', self.diagonalA), ('test_constant.png', self.constantA)]:

            valA = valA_fn()

            # A <--> I
            valIA = icebin.coo_multiply(self.IvA, valA, fill=np.nan)
            valAIA = icebin.coo_multiply(self.AvI, valIA, fill=np.nan)

            # Plot it...
            figure = matplotlib.pyplot.figure(figsize=(17, 11))

            ax = figure.add_subplot(121)
            basemap = giss.basemap.greenland_laea(ax=ax)
            giss.plot.plot_var(ax=ax, basemap=basemap,
                **giss.modele.plot_params('A', val=valA, plotter=self.plotterA))

            ax = figure.add_subplot(122)
            basemap = giss.basemap.greenland_laea(ax=ax)
            giss.plot.plot_var(ax=ax, basemap=basemap,
                **giss.modele.plot_params('A', val=valAIA, plotter=self.plotterA))

            figure.savefig(fname, dpi=100, transparent=False)

            print('----------- Range Results -----------')
            print('valA range = [{}, {}]'.format(np.nanmin(valA), np.nanmax(valA)))
            print('valIA range = [{}, {}]'.format(np.nanmin(valIA), np.nanmax(valIA)))
            print('valAIA range = [{}, {}]'.format(np.nanmin(valAIA), np.nanmax(valAIA)))



if __name__ == '__main__':
    unittest.main()


from icebin import ibgrid, element_l1
from matplotlib import pyplot
import giss.basemap
import unittest
import numpy as np

projection='+proj=stere +lon_0=-39 +lat_0=90 +lat_ts=71.0 +ellps=WGS84'
class TestMesh(unittest.TestCase):

    def setUp(self):
        # -------------- Ice Grid
        vertices = [ibgrid.Vertex(spec[0], spec[1], spec[2]) for spec in [
            (0, -1., -1.),
            (1,  1., -1.),
            (2,  1.,  1.),
            (3, -1.,  1.),
            (4,  0.,  0.)]]

        vx = vertices
        cells = [ibgrid.Cell(*spec) for spec in [
            (0, (vx[0], vx[1], vx[4]), 0,0,0,0),
            (1, (vx[1], vx[2], vx[4]), 0,0,0,0),
            (2, (vx[2], vx[3], vx[4]), 0,0,0,0),
            (3, (vx[3], vx[0], vx[4]), 0,0,0,0)]]

        self.gridI = ibgrid.Grid(vertices, cells,
            projection=projection,
            type='MESH',
            coordinates='XY')

    def test_A1(self):
        """A has a single grid cell covering everything."""

        # ------- Ice Grid
        gridI = self.gridI

        # ------- Atmosphere Grid
        vertices = [ibgrid.Vertex(spec[0], spec[1], spec[2]) for spec in [
            (0, -2., -2.),
            (1,  2., -2.),
            (2,  2.,  2.),
            (3, -2.,  2.)]]
        vx = vertices
        cells = [ibgrid.Cell(0, [vx[0], vx[1], vx[2], vx[3]], 0,0,0,0)]

        gridA = ibgrid.Grid(vertices, cells,
            projection=projection,
            type='MESH',
            coordinates='XY')

        # ------ Exchange Grid
        vertices = dict(gridI.vertices)
        vx = vertices
        cells = [ibgrid.Cell(*spec) for spec in [
            (0, (vx[0], vx[1], vx[4]), 0,0,0,0),
            (1, (vx[1], vx[2], vx[4]), 0,1,0,0),
            (2, (vx[2], vx[3], vx[4]), 0,2,0,0),
            (3, (vx[3], vx[0], vx[4]), 0,3,0,0)]]
                
        gridX = ibgrid.Grid(vx, cells,
            projection=projection,
            type='MESH',
            coordinates='XY',
            grid1_ncells_full=gridA.cells_num_full,
            grid2_ncells_full=gridI.cells_num_full,
            grid1_nvertices_full=gridA.vertices_num_full,
            grid2_nvertices_full=gridI.vertices_num_full)

        AvI,weightsA,weightsI = element_l1.compute_AvI(gridX, gridI)
        print(AvI.todense())
        print(weightsA)
        print(weightsI)

        self.assertAlmostEqual(2./3., weightsI[0])
        self.assertAlmostEqual(2./3., weightsI[1])
        self.assertAlmostEqual(2./3., weightsI[2])
        self.assertAlmostEqual(2./3., weightsI[3])
        self.assertAlmostEqual(1.+1./3., weightsI[4])

        M = np.asarray(AvI.todense())
        self.assertAlmostEqual(M[0,0], weightsI[0])
        self.assertAlmostEqual(M[0,1], weightsI[1])
        self.assertAlmostEqual(M[0,2], weightsI[2])
        self.assertAlmostEqual(M[0,3], weightsI[3])
        self.assertAlmostEqual(M[0,4], weightsI[4])

        self.assertAlmostEqual(sum(weightsA), sum(weightsI))


    def test_A2(self):
        """A has a two grid cells."""

        # ------- Ice Grid
        gridI = self.gridI

        # ------- Atmosphere Grid
        vertices = [ibgrid.Vertex(spec[0], spec[1], spec[2]) for spec in [
            (0, -2., -2.),
            (1,  0., -2.),
            (2,  2., -2.),
            (3,  2.,  2.),
            (4,  0.,  2.),
            (5, -2.,  2.),]]
        vx = vertices
        cells = [ibgrid.Cell(*spec) for spec in [
            (0, (vx[0], vx[1], vx[4], vx[5]), 0,0,0,0),
            (1, (vx[1], vx[2], vx[3], vx[4]), 0,0,0,0)]]

        gridA = ibgrid.Grid(vertices, cells,
            projection=projection,
            type='MESH',
            coordinates='XY')

        # ------ Exchange Grid
        vertices = [ibgrid.Vertex(spec[0], spec[1], spec[2]) for spec in [
            (0, -1., -1.),      # Original points from gridI
            (1,  1., -1.),
            (2,  1.,  1.),
            (3, -1.,  1.),
            (4,  0.,  0.),      # Middle point
            (5,  0., -1.),      # Split points created by overlap
            (6,  0.,  1.)]]

        vx = vertices
        cells = [ibgrid.Cell(*spec) for spec in [
            (0, (vx[0], vx[5], vx[4]), 0,0,0,0),
            (1, (vx[5], vx[1], vx[4]), 1,0,0,0),
            (2, (vx[1], vx[2], vx[4]), 1,1,0,0),
            (3, (vx[2], vx[6], vx[4]), 1,2,0,0),
            (4, (vx[6], vx[3], vx[4]), 0,2,0,0),
            (5, (vx[3], vx[0], vx[4]), 0,3,0,0)]]

        gridX = ibgrid.Grid(vx, cells,
            projection=projection,
            type='MESH',
            coordinates='XY',
            grid1_ncells_full=gridA.cells_num_full,
            grid2_ncells_full=gridI.cells_num_full,
            grid1_nvertices_full=gridA.vertices_num_full,
            grid2_nvertices_full=gridI.vertices_num_full)

        AvI,weightsA,weightsI = element_l1.compute_AvI(gridX, gridI)
        #print(AvI.todense())
        #print(weightsA)
        #print(weightsI)

        self.assertAlmostEqual(sum(weightsA), sum(weightsI))


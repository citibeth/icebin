import unittest
import scipy.sparse
import itertools

import icebin
import matplotlib
#import basemap
import numpy as np
import netCDF4
import pickle
import sys
import giss.basemap
import giss.modele
import giss.plot
import giss.pism

ice_sheet = 'greenland'
ICEBIN_IN = 'icebin_in.nc'

regrid_names = ('IvA', 'AvI', 'IvE', 'EvI', 'EvA', 'AvE')

class RegridTests(unittest.TestCase):

	def setUp(self):
		self.mm = icebin.GCMRegridder(ICEBIN_IN)
		self.rm = self.mm.regrid_matrices('greenland')
		with netCDF4.Dataset(ICEBIN_IN) as nc:
			self.indexingA = icebin.Indexing(nc, 'm.gridA.indexing')
			self.indexingHP = icebin.Indexing(nc, 'm.indexingHP')
			self.indexingI = icebin.Indexing(nc, 'm.{}.gridI.indexing'.format(ice_sheet))
			self.plotterI = icebin.read_plotter(nc, 'm.{}.gridI'.format(ice_sheet))
			self.plotterA = icebin.read_plotter(nc, 'm.gridA')

		self.elevI, self.maskI = giss.pism.read_elevI_maskI('elev_mask.nc')

		self.AvI,self.wAvI = self.rm.regrid('AvI', scale=True)
		self.IvA,self.wIvA = self.rm.regrid('IvA', scale=True)
		self.EvI,self.wEvI = self.rm.regrid('EvI', scale=True)
		self.IvE,self.wIvE = self.rm.regrid('IvE', scale=True)
		self.EvA,self.wEvA = self.rm.regrid('EvA', scale=True)
		self.AvE,self.wAvE = self.rm.regrid('AvE', scale=True)



	# -----------------------------------------------------
	def assert_equal_sparse(self, rname, _A, _B):
		print('Regrid {}:'.format(rname))
		A = scipy.sparse.coo_matrix(_A)
		B = scipy.sparse.coo_matrix(_B)

		print('    A.shape = ', A.shape, A.nnz)
		print('    B.shape = ', B.shape, B.nnz)

		At = sorted(zip(A.row, A.col, A.data))
		Bt = sorted(zip(B.row, B.col, B.data))

		for a,b in zip(At,Bt):
			self.assertEqual(a[0:2],b[0:2])
			self.assertAlmostEqual(1.,a[2]/b[2])

	def assert_equal_np(self, A, B):
		for a,b in zip(A,B):
			ratio = a/b
			if not np.isnan(ratio):
				try:
					self.assertAlmostEqual(1.,ratio)
				except:
					print('A={}, B={}'.format(a,b))
					raise

	# -----------------------------------------------------
	def run_test_matrix(self, rname):
		BvA_s, w0 = self.rm.regrid(rname, scale=True)
		BvA,   w1 = self.rm.regrid(rname, scale=False)

		# Weights un-scale the matrix
		w0_2 = w0.reshape(1,len(w0))

		dia_w0 = scipy.sparse.dia_matrix((w0_2, [0]), shape=(len(w0), len(w0)))
		BvA2 = dia_w0 * BvA_s

		self.assert_equal_sparse(rname, BvA, BvA2)

	def test_matrices(self):
		"""Check that the weighted and unweighted versions of the
		matrices are consistent with the weight vectors.
		     We should have:   BvA = wBvA * BvA(SCALED)"""

		self.run_test_matrix('IvA')
		self.run_test_matrix('AvI')
		self.run_test_matrix('IvE')
		self.run_test_matrix('EvI')
		self.run_test_matrix('EvA')
		self.run_test_matrix('AvE')
	# -----------------------------------------------------
	# Sample functions

	def diagonalI(self):
		shapeI = self.indexingI.shape
		valI = np.zeros(shapeI)
		for i in range(0,shapeI[0]):
			for j in range(0,shapeI[1]):
				valI[i,j] = i + j
		return valI.reshape(-1)

	def constantI(self):
		valI = np.zeros(self.indexingI.shape) + 1
		return valI.reshape(-1)

	def elevationI(self):
		return self.elevI.reshape(-1)

	# -----------------------------------------------------

	def test_constant_regrid(self):
		"""Regrid the constant function, check that answer is always
		what we started with."""

		valI = self.constantI()

		# A <-- I
		valAI = icebin.coo_multiply(self.AvI, valI, fill=np.nan)		# Most recent space is on the left
		valIAI = icebin.coo_multiply(self.IvA, valAI, fill=np.nan)		# Read history right-to-left.
		self.assert_equal_np(valI, valIAI)

		# Make sure A<--I DOES have projection scaling
		# Look at valAI, make sure it's not all the same
		x = valAI - 1
		xsum = np.nansum(x*x)
		# This will be 0 if valAI is 1 everywhere, non-zero if valAI varies.
		# valAI should vary because the projection used (Stereographic) is non-area-preserving.
		self.assertNotEqual(0., xsum)
			

		# E <-- I
		valEI = icebin.coo_multiply(self.EvI, valI, fill=np.nan)		# Most recent space is on the left
		valIEI = icebin.coo_multiply(self.IvE, valEI, fill=np.nan)		# Read history right-to-left.
		self.assert_equal_np(valI, valIEI)

		# A <-- E
		valE = np.zeros(len(self.indexingHP)) + 1
		valE = valE.reshape(-1)

		valAE = icebin.coo_multiply(self.AvE, valE, fill=np.nan)		# Most recent space is on the left
		valEAE = icebin.coo_multiply(self.EvA, valAE, fill=np.nan)		# Read history right-to-left.
		self.assert_equal_np(valE, valEAE)

		# Make sure A<-E and E<-A have no net projection scaling
		valA = np.zeros(len(self.indexingA)) + 1
		valA = valA.reshape(-1)
		valEA = icebin.coo_multiply(self.EvA, valA, fill=np.nan)
		valAEA = icebin.coo_multiply(self.AvE, valEA, fill=np.nan)
		self.assert_equal_np(valA, valAEA)
		self.assert_equal_np(valA, valAE)
		self.assert_equal_np(valEA, valE)

	def assert_eq_weighted(self, A, wA, B, wB):
		Asum = np.nansum(np.multiply(A, wA))
		Bsum = np.nansum(np.multiply(B, wB))
		self.assertAlmostEqual(1., Asum/Bsum)


	def test_conserv(self):
		"""Tests conservation of mass.  Checks that the amount of
		stuff (== value * weight) remains constant."""

		for fname,valI_fn in (('test_diagonal.png', self.diagonalI), ('test_elevation.png', self.elevationI)):

			valI = valI_fn()

			# A <--> I
			valAI = icebin.coo_multiply(self.AvI, valI, fill=np.nan)
			self.assert_eq_weighted(valI, self.wIvA, valAI, self.wAvI)
			valIAI = icebin.coo_multiply(self.IvA, valAI, fill=np.nan)
			self.assert_eq_weighted(valIAI, self.wIvA, valAI, self.wAvI)

			# E <--> I
			valEI = icebin.coo_multiply(self.EvI, valI, fill=np.nan)
			self.assert_eq_weighted(valI, self.wIvE, valEI, self.wEvI)
			valIEI = icebin.coo_multiply(self.IvE, valEI, fill=np.nan)
			self.assert_eq_weighted(valIEI, self.wIvE, valEI, self.wEvI)

			# A <--> E
			valE = valEI
			valAE = icebin.coo_multiply(self.AvE, valE, fill=np.nan)
			self.assert_eq_weighted(valE, self.wEvA, valAE, self.wAvE)
			valEAE = icebin.coo_multiply(self.EvA, valAE, fill=np.nan)
			self.assert_eq_weighted(valEAE, self.wEvA, valAE, self.wAvE)

			# Plot it...
			figure = matplotlib.pyplot.figure(figsize=(22,17))

			ax = figure.add_subplot(231)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('I', val=valI, plotter=self.plotterI))

			ax = figure.add_subplot(232)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('AI', val=valAI, plotter=self.plotterA))

			ax = figure.add_subplot(233)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('IAI', val=valIAI, plotter=self.plotterI))

			ax = figure.add_subplot(234)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('IEI', val=valIEI, plotter=self.plotterI))

			valAEI = icebin.coo_multiply(self.AvE, valEI, fill=np.nan)
			ax = figure.add_subplot(235)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('AEI', val=valAEI, plotter=self.plotterA))

			valEAI = icebin.coo_multiply(self.EvA, valAI, fill=np.nan)
			valIEAI = icebin.coo_multiply(self.IvE, valEAI, fill=np.nan)
			ax = figure.add_subplot(236)
			basemap = giss.basemap.greenland_laea(ax=ax)
			giss.plot.plot_var(ax=ax, basemap=basemap,
				**giss.modele.plot_params('IEAI', val=valIEAI, plotter=self.plotterI))

			figure.savefig(fname, dpi=100, transparent=False)



if __name__ == '__main__':
    unittest.main()

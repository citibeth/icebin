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

ice_sheet = 'greenland'
ICEBIN_IN = sys.argv[1]

mm = icebin.GCMRegridder(ICEBIN_IN)
rm = mm.regrid_matrices('greenland')
with netCDF4.Dataset(ICEBIN_IN) as nc:
	indexingA = icebin.Indexing(nc, 'm.gridA.indexing')
	indexingHP = icebin.Indexing(nc, 'm.indexingHP')
	indexingI = icebin.Indexing(nc, 'm.{}.gridI.indexing'.format(ice_sheet))
	plotterI = icebin.read_plotter(nc, 'm.{}.gridI'.format(ice_sheet))
	plotterA = icebin.read_plotter(nc, 'm.gridA')

# -----------------------------------

def exampleI():
	shapeI = indexingI.shape
	print('shapeI', shapeI)
	valI = np.zeros(shapeI)
	for i in range(0,shapeI[0]):
		for j in range(0,shapeI[1]):
			valI[i,j] = i + j
	return valI.reshape(-1)

def plot_IG():
	valI = exampleI()
	GvI = rm.regrid('GvI')
	sGvI = rm.scale('sGvI')
	valG = np.multiply(sGvI, icebin.coo_multiply(GvI, valI))

	IvG = GvI.transpose()		# This is same as rm.regrid('IvG')
	sIvG = rm.scale('sIvG')
	valI2 = np.multiply(sIvG, icebin.coo_multiply(IvG, valG))

	print('valI numnan = {}'.format(np.count_nonzero(np.isnan(valI))))
	print('valG numnan = {}'.format(np.count_nonzero(np.isnan(valG))))
	print('valI2 numnan = {}'.format(np.count_nonzero(np.isnan(valI2))))

	# Plot our value
	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	ax = figure.add_subplot(121)
	basemap = giss.basemap.greenland_laea(ax=ax)
	pp = giss.modele.plot_params('valI', val=valI, plotter=plotterI)
	print(pp)
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)

	# Plot our value
	ax = figure.add_subplot(122)
	basemap = giss.basemap.greenland_laea(ax=ax)
	pp = giss.modele.plot_params('valI', val=valI2, plotter=plotterI)
	print(pp)
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)

	figure.savefig('example_IG.png', dpi=100, transparent=False)




def conserv_IvA():
	valI = exampleI()

	wI = rm.scale('wIvG')
	wA = rm.scale('wAvG')
	wE = rm.scale('wEvG')
#	wA = rm.scale('wA')
	print('wA', wA.shape, wA, wA[10698])
	print('wI', wI.shape, wI)
	IvA = rm.regrid('IvA(PARTIAL_CELL)')

	AvI = rm.regrid('AvI(PARTIAL_CELL)')
	print('AvI', AvI.shape, AvI.nnz)
	valA = icebin.coo_multiply(AvI, valI, fill=np.nan)
	#valA = AvI * valI
	print('valA', valA.shape, valA)

	# Plot our value
	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	ax = figure.add_subplot(121)
	basemap = giss.basemap.greenland_laea(ax=ax)
	pp = giss.modele.plot_params('valA', val=valA, plotter=plotterA)
	print(pp)
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)

	valI2 = icebin.coo_multiply(IvA, valA, fill=np.nan)
	valI[np.isnan(valI2)] = np.nan

	print('sum', np.nansum(valI), np.nansum(valI2))

	# Plot our value
	ax = figure.add_subplot(122)
	basemap = giss.basemap.greenland_laea(ax=ax)
	pp = giss.modele.plot_params('valI', val=valI2, plotter=plotterI)
	print(pp)
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)

	figure.savefig('AxI.png', dpi=100, transparent=False)

plot_IG()
conserv_IvA()

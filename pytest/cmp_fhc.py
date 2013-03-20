import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot

basemap = giss.basemap.greenland_laea()

# Sum over all elevation classes
nc1 = netCDF4.Dataset('fhc.nc')
fhc1 = nc1.variables['fhc1h'][:]

fhc1_sum = np.sum(fhc1, axis=(0))

pp = giss.modele.plot_params(val=fhc1_sum)
pp['title'] = 'fhc1_sum'
giss.plot.plot_var(basemap=basemap, **pp)		# Plot, and show on screen


for ihc in range(0,10) :

	nc0 = netCDF4.Dataset('GIC.144X90.DEC01.1.ext_hc.nc')
	fhc0 = nc0.variables['fhc'][ihc,:]

	nc1 = netCDF4.Dataset('fhc.nc')
	fhc1 = nc1.variables['fhc1h'][ihc,:]

	fhcdiff = fhc1 - fhc0
	# ----------------------------------------
	figure = matplotlib.pyplot.figure(figsize=(15,8.5))

	ax = figure.add_subplot(131)
	pp = giss.modele.plot_params('fhc', nc0, val=fhc0)
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen


	ax = figure.add_subplot(132)
	pp = giss.modele.plot_params('fhc1h', nc1, val=fhc1)
	#pp = giss.modele.plot_params(val=fhcdiff)
	#pp['title'] = 'fhcdiff'
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen

	ax = figure.add_subplot(133)
	pp = giss.modele.plot_params(val=fhcdiff * 1000)
	pp['title'] = 'fhcdiff (*1000)'
	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen

	figure.savefig('fhc-%d.png' % ihc)
#	matplotlib.pyplot.show()




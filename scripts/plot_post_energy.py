import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import mpl_toolkits.basemap
import glint2
import numpy as np
import sys


def pp_posneg(name, val):
	global plotno
	plotno += 1
	ax = figure.add_subplot(1,2,plotno)
	pp = giss.modele.plot_params(name, val=val, nc=nc, plotter=plotter2)
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=True).cmap
	pp['plot_args']['norm'] = giss.plot.AsymmetricNormalize()

	dmin = np.nanmin(pp['val'])
	dmax = np.nanmax(pp['val'])
	print 'dminmax: %s' % name,dmin,dmax
	ticks = [dmin,0,dmax]
	pp['cb_args']['ticks'] = ticks
	pp['cb_args']['format']='%d'

	pp['ax'] = ax
	pp['basemap'] = basemap
	return pp

def pp_pos(name, val):
	global plotno
	plotno += 1
	ax = figure.add_subplot(1,2,plotno)
	pp = giss.modele.plot_params(name, val=val, nc=nc, plotter=plotter2)
#	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=True).cmap
#	pp['plot_args']['norm'] = giss.plot.AsymmetricNormalize()

	dmin = np.nanmin(pp['val'])
	dmax = np.nanmax(pp['val'])
	print 'dminmax: %s' % name,dmin,dmax
	ticks = [dmin,dmax]
	pp['cb_args']['ticks'] = ticks
	pp['cb_args']['format']='%d'

	pp['ax'] = ax
	pp['basemap'] = basemap
	return pp

ice_density = 916.6
def compute_enth2(z, Enth3, ice_thickness):
	print z.shape
	print Enth3.shape
	print ice_thickness.shape
	nx = Enth3.shape[0]
	ny = Enth3.shape[1]
	nz = Enth3.shape[2]

	enth2 = np.zeros((nx, ny))
	for i in xrange(0, nx):
		for j in xrange(0, ny):
			for k in xrange(0, nz):
				z0 = z[k]
				z1 = z[k+1]
				if z1 > ice_thickness[i,j]:
					enth2[i,j] += (ice_thickness[i,j] - z0) * Enth3[i,j,k]
					break
				else:
					pass
#					enth2[i,j] += (z1 - z0) * Enth3[i,j,k]
			enth2[i,j] *= ice_density

	return enth2

grid_fname = '/Users/rpfische/git/glint2-attachments/example2/modele_ll_g2x2_5-searise_g20-40-PISM.nc'
data_fname = 'post_energy.nc'

nc = netCDF4.Dataset(grid_fname, 'r')
plotter2 = glint2.Plotter2(nc=nc, vname='m.greenland.grid2', transpose=True)
nc.close()


# Use a custom basemap
basemap = giss.basemap.greenland_laea()

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(14,11))


nc = netCDF4.Dataset(data_fname)

plotno = 0

# -------------------------------------------
ix0 = 0
ix1 = 1
t0 = nc.variables['time'][ix0]
t1 = nc.variables['time'][ix1]
z = nc.variables['z'][:]
by_dt = 1.0/(t1 - t0)

if False:
	# -------------------------------------------
	vname='total.mass'
	val0 = nc.variables[vname][ix0,:]
	val1 = nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=(val1-val0)*by_dt)
	pp['cb_args']['format']='%1.2g'
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='surface_mass_balance.mass'
	val=nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='internal_advection.mass'
	val=nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	giss.plot.plot_var(**pp)
	# ------------------------------------------
	vname='epsilon.mass'
	val=nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	giss.plot.plot_var(**pp)
	# ------------------------------------------
		
# ======================================
# *.enth


ice_thickness_0 = nc.variables['thk'][ix0,:]
Enth3_0 = nc.variables['enthalpy'][ix0,:]
enth2c_0 = compute_enth2(z, Enth3_0, ice_thickness_0)
enth2_0 = nc.variables['total.enth'][ix0,:]

ice_thickness_1 = nc.variables['thk'][ix1,:]
Enth3_1 = nc.variables['enthalpy'][ix1,:]
enth2c_1 = compute_enth2(z, Enth3_1, ice_thickness_1)
enth2_1 = nc.variables['total.enth'][ix1,:]


# -------------------------------------------
pp = pp_posneg('thk_diff', val=ice_thickness_1 - ice_thickness_0)
pp['cb_args']['format']='%1.2g'
pp['title'] = pp['title']
pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
giss.plot.plot_var(**pp)

# -------------------------------------------
pp = pp_posneg('enth_diff', val=(enth2c_1 - enth2c_0) * by_dt)
pp['cb_args']['format']='%1.2g'
pp['title'] = pp['title']
pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
giss.plot.plot_var(**pp)

# -------------------------------------------

if False:
	## -------------------------------------------
	#vname='total.enth'
	#val0 = nc.variables[vname][ix0,:]
	#val1 = nc.variables[vname][ix1,:]
	##val=(val1-val0) * by_dt
	#val=val1
	#pp = pp_posneg(vname, val=val)
	#pp['cb_args']['format']='%1.2f'
	#pp['title'] = pp['title']
	#pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	#giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='total.enth'
	val0 = nc.variables[vname][ix0,:]
	val1 = nc.variables[vname][ix1,:]
	#val=(val1-val0) * by_dt
	val=val0
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2f'
	pp['title'] = pp['title']
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	#vname='enth2c'
	#pp = pp_posneg(vname, val=enth2c)
	#pp['cb_args']['format']='%1.2f'
	#pp['title'] = pp['title']
	#pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	#giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='surface_mass_balance.enth'
	val=nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	pp['title'] = pp['title']
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='internal_advection.enth'
	val=nc.variables[vname][ix1,:]
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	pp['title'] = pp['title']
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname='epsilon.enth'
	val=nc.variables[vname][ix1,:]
	val[val == float('inf')] = 0
	pp = pp_posneg(vname, val=val)
	pp['cb_args']['format']='%1.2g'
	pp['title'] = pp['title']
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
	giss.plot.plot_var(**pp)
	# -------------------------------------------

if False:
	vname = 'strain_heating'
	val = nc.variables[vname][ix1,:]
	val[abs(val) > 1e100] = 0
	pp = pp_pos(vname, val)
	pp['cb_args']['format']='%1.2g'
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname = 'geothermal_flux'
	pp = pp_pos(vname, nc.variables[vname][ix1,:])
	pp['cb_args']['format']='%f'
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	vname = 'upward_geothermal_flux'
	pp = pp_pos(vname, nc.variables[vname][ix1,:])
	giss.plot.plot_var(**pp)
	# -------------------------------------------
	#vname = 'basal_frictional_heating'
	#pp = pp_pos(vname, nc.variables[vname][ix1,:])
	#giss.plot.plot_var(**pp)
	# -------------------------------------------
		





# Save to a file as png
figure.savefig('fig.png', dpi=300, transparent=True)

# Also show on screen
matplotlib.pyplot.show()

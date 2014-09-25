import netCDF4
import giss.basemap
import giss.modele
import matplotlib.pyplot
import mpl_toolkits.basemap
import glint2
import numpy as np



def pp_posneg(name, val):
	global plotno
	plotno += 1
	ax = figure.add_subplot(3,4,plotno)
	pp = giss.modele.plot_params(name, val=val, nc=nc, plotter=plotter2)
	pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=True).cmap
	pp['plot_args']['norm'] = giss.plot.AsymmetricNormalize()

	dmin = np.nanmin(pp['val'])
	dmax = np.nanmax(pp['val'])
#	print dmin,dmax
	ticks = [dmin,0,dmax]
	pp['cb_args']['ticks'] = ticks
	pp['cb_args']['format']='%d'

	pp['ax'] = ax
	pp['basemap'] = basemap
	return pp

def pp_pos(name, val):
	global plotno
	plotno += 1
	ax = figure.add_subplot(3,4,plotno)
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
ix0 = 1
ix1 = 2
t0 = nc.variables['time'][ix0]
t1 = nc.variables['time'][ix1]
by_dt = 1.0/(t1 - t0)
# -------------------------------------------
vname='total.mass'
val0 = nc.variables[vname][ix0,:]
val2 = nc.variables[vname][ix1,:]
pp = pp_posneg(vname, val=(val2-val0)*by_dt)
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
# -------------------------------------------
vname='total.enth'
val0 = nc.variables[vname][ix0,:]
val2 = nc.variables[vname][ix1,:]
val=(val2-val0) * by_dt
pp = pp_posneg(vname, val=val)
pp['cb_args']['format']='%1.2f'
pp['title'] = pp['title']
pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
giss.plot.plot_var(**pp)
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
pp = pp_posneg(vname, val=val)
pp['cb_args']['format']='%1.2g'
pp['title'] = pp['title']
pp['plot_args']['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=False).cmap
giss.plot.plot_var(**pp)
# -------------------------------------------
#vname = 'basal_frictional_heating'
#pp = pp_posneg(vname, nc.variables[vname][ix1,:])
#giss.plot.plot_var(**pp)
# -------------------------------------------
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






# Save to a file as png
figure.savefig('fig.png', dpi=300, transparent=True)

# Also show on screen
#matplotlib.pyplot.show()

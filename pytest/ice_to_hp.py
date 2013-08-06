# See what kind of fill-in we get when we compute A^T A,
# where A is hp_to_ice().

import glint2
import netCDF4
import numpy as np
import giss.modele
import sys
import os.path
import giss.plot
import matplotlib.pyplot as plt

#mm = glint2.MatrixMaker(correct_area1=False)
glint2_config_file = '/Users/rpfische/cmrun/rpfische/modele_ll_g2x2_5-searise_g20-40-DISMAL.nc'
mm = glint2.MatrixMaker(glint2_config_file, 'm', correct_area1=True)

interp_grid = 'EXCH'
mm.set_interp_grid(interp_grid)

RUNDIR1='/Users/rpfische/savedisk/e4f40-hc40'
RUNDIR2='/Users/rpfische/savedisk/e4f40-hc40'

#nc = netCDF4.Dataset('../data/JUL1950.ijhchc1k225.nc')
#nc = netCDF4.Dataset('../data/JUL1956.ijhce4f40-hc40.nc')
nc = netCDF4.Dataset(os.path.join(RUNDIR1, 'JUN1950.ijhce4f40-hc40.nc'))
impm3 = giss.modele.read_ncvar(nc, 'impm_lndice')
impm3 = impm3[1:,:,:]			# Remove legacy height point
frac3 = giss.modele.read_ncvar(nc, 'frac')
frac3 = frac3[1:,:,:]			# Remove legacy height point
print 'frac3.shape = ', frac3.shape
nc.close()

RM = mm.hp_to_atm()
XM = mm.hp_to_iceinterp('greenland', dest='ICE')
M = mm.hp_to_iceinterp('greenland', dest='INTERP')

impms = []
impms.append(RM.dot(impm3.reshape(-1)))
impm2s = []
impm2s.append(XM.dot(impm3.reshape(-1)))

# The first one is OFFICIAL!!!
impms.append(RM.dot(impm3.reshape(-1)))


if interp_grid == 'EXCH' :
	impm3b = mm.iceinterp_to_hp(
		{'greenland' : M.dot(impm3.reshape(-1))},
		initial3=impm3.reshape(-1),
		src='INTERP',
		qp_algorithm='MULTI_QP')
	impms.append(RM.dot(impm3b))
	impm2s.append(XM.dot(impm3b))

impm4 = M.dot(impm3.reshape(-1))
impm3b = mm.iceinterp_to_hp(
	{'greenland' : impm4},
	initial3=impm3.reshape(-1),
	src='INTERP',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))
impm2s.append(XM.dot(impm3b))

impm2 = XM.dot(impm3.reshape(-1))
impm4b = mm.ice_to_interp('greenland', impm2)
impm3b = mm.iceinterp_to_hp(
	{'greenland' : impm4b},
	initial3=impm3.reshape(-1),
	src='INTERP',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))
impm2s.append(XM.dot(impm3b))

impm3b = mm.iceinterp_to_hp(
	{'greenland' : impm2},
	initial3=impm3.reshape(-1),
	src='ICE',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))
impm2s.append(XM.dot(impm3b))

print 'impm2.shape = ',impm2.shape
S = mm.iceinterp_to_atm('greenland', src='ICE')
impms.append(S.dot(impm2))

R = mm.iceinterp_to_atm('greenland', src='INTERP')
impms.append(R.dot(impm4b))
impms.append(R.dot(impm4))

# ----------------------------------------
# Test atm_to_hp

impm1 = RM.dot(impm3.reshape(-1))
impm3b = mm.atm_to_hp(impm1)
impms.append(RM.dot(impm3b))
impm2s.append(XM.dot(impm3b))




M = mm.hp_to_iceinterp('greenland', dest='INTERP')

area1 = mm.area1('greenland')
area1xy = area1.reshape(impm3.shape[1:])
print RM.shape
for (j,i) in {(75,55)} :
	for k in range(0,len(impms)) :
		impm1 = impms[k]
		impm1xy = impm1.reshape(impm3.shape[1:])
		total = np.dot(impm1, area1)
#		impm1 = RM.dot(impms[k].reshape(-1)).reshape(impm3.shape[1:])
		print '%d: conserv[%d,%d] = %g / %g (total=%1.15f)' % (k, j, i, impm1xy[j,i], impm1xy[j,i] * area1xy[j,i], total)



plotter2 = glint2.Plotter2(fname=glint2_config_file, vname='m.greenland.grid2')
for k in range(0,len(impm2s)) :
	fig = plt.figure(figsize=(8.5,11))
	ax = fig.add_subplot(111)

	basemap = giss.basemap.greenland_laea()

	var_name = 'Try %d' % k
	pp = giss.modele.plot_params(var_name, val=impm2s[k], plotter=plotter2)
	print pp['plot_args']
	#pp = {}
	# pp['var_name'] = 'Try %d' % k
	#pp['val'] = nc.variables[var_name][:]
	#pp['title'] = "%s:%s" % (nc_name, var_name)
	#pp['plot_boundaries'] = _default_plot_boundaries
	pp['plotter'] = plotter2

	plot_args = {}
	pp['plot_args'] = plot_args
	pp['plot_args']['vmax'] = .00033
	#plot_args['norm'] = giss.plot.AsymmetricNormalize()
	#_reverse_scale = {'mass'}
	#reverse = (var_name in _reverse_scale)
	#plot_args['cmap'] = giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=reverse).cmap

	giss.plot.plot_var(ax=ax, basemap=basemap, **pp)		# Plot, and show on screen

	# Save to a file as png
	fig.savefig('fig%d.png' % k, dpi=300, transparent=True)

	# Also show on screen
#	plt.show()


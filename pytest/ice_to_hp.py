# See what kind of fill-in we get when we compute A^T A,
# where A is hp_to_ice().

import glint2
import netCDF4
import numpy as np
import giss.modele
import sys
import os.path

#mm = glint2.MatrixMaker(correct_area1=False)
mm = glint2.MatrixMaker('/Users/rpfische/cmrun/rpfische/modele_ll_g2x2_5-searise_g20-40-DISMAL.nc', 'm', correct_area1=False)


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

impms = []
impms.append(RM.dot(impm3.reshape(-1)))

M = mm.hp_to_ice('greenland', dest='EXCH')
impm3b = mm.ice_to_hp(
	{'greenland' : M.dot(impm3.reshape(-1))},
	initial3=impm3.reshape(-1),
	src='EXCH',
	qp_algorithm='MULTI_QP')
impms.append(RM.dot(impm3b))

M = mm.hp_to_ice('greenland', dest='EXCH')
impm4 = M.dot(impm3.reshape(-1))
impm3b = mm.ice_to_hp(
	{'greenland' : impm4},
	initial3=impm3.reshape(-1),
	src='EXCH',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))

XM = mm.hp_to_ice('greenland', dest='ICE')
impm2 = XM.dot(impm3.reshape(-1))
impm4b = mm.ice_to_exch('greenland', impm2)
impm3b = mm.ice_to_hp(
	{'greenland' : impm4b},
	initial3=impm3.reshape(-1),
	src='EXCH',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))

impm3b = mm.ice_to_hp(
	{'greenland' : impm2},
	initial3=impm3.reshape(-1),
	src='ICE',
	qp_algorithm='SINGLE_QP')
impms.append(RM.dot(impm3b))


S = mm.ice_to_atm('greenland', src='ICE')
impms.append(S.dot(impm2))

R = mm.ice_to_atm('greenland', src='EXCH')
impms.append(R.dot(impm4b))
impms.append(R.dot(impm4))

area1 = mm.area1('greenland')
area1xy = area1.reshape(impm3.shape[1:])
print RM.shape
for (j,i) in {(75,55)} :
	for k in range(0,len(impms)) :
		impm1 = impms[k]
		impm1xy = impm1.reshape(impm3.shape[1:])
		total = np.dot(impm1, area1)
#		impm1 = RM.dot(impms[k].reshape(-1)).reshape(impm3.shape[1:])
		print '%d: conserv[%d,%d] = %g (%1.15f)' % (k, j, i, impm1xy[j,i] * area1xy[j,i], total)

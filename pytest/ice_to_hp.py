# See what kind of fill-in we get when we compute A^T A,
# where A is hp_to_ice().

import glint2
import netCDF4
import numpy as np
import giss.modele
import sys
import os.path

#mm = glint2.MatrixMaker(correct_area1=False)
mm = glint2.MatrixMaker('/Users/rpfische/cmrun/rpfische/modele_ll_g2x2_5-searise_g20-40-DISMAL.nc', 'm', correct_area1=True)

XM = mm.hp_to_ice('greenland')
print 'XM', XM.shape, XM.nnz

XMt = XM.transpose()
print 'XMt', XMt.shape, XMt.nnz

PP = XMt * XM
print 'PP', PP.shape, PP.nnz

# This is the wrong way, and it gets lots of fill-in to boot
#PP = XM * XMt
#print 'PP', PP.shape, PP.nnz


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


impm2 = glint2.coo_multiply(XM, impm3)
print 'impm2', impm2.shape
f2s = dict()
f2s['greenland'] = impm2

print 'impm2.shape',impm2.shape
print 'impm3.shape',impm3.shape

impm3r = impm3.reshape(-1)
impm2b = np.zeros(impm2.shape)
for row,col,val in zip(XM.row, XM.col, XM.data) :
	impm2b[row] += impm3r[col] * val
#	print 'XM-A,row,col,val,impm3r[col]'

#for i in range(0,impm2.shape[0]) :
#	print impm2[i], impm2b[i],impm2[i] - impm2b[i]

#sys.exit(0)


#initial = (impm3 + np.random.normal(scale=1e-6,size=impm3.shape)).reshape(-1)
initial = impm3.reshape(-1)
#initial[:] = 0
impm3b = mm.ice_to_hp(f2s, initial).reshape(impm3.shape)
#impm3b = impm3


# ------------------------------------------------------
RM = mm.hp_to_atm()

conserv1 = glint2.coo_multiply(RM, impm3)
conserv1b = glint2.coo_multiply(RM, impm3b)

# ------------------------------------------------------
nc = netCDF4.Dataset('./ice_to_hp.nc', 'w')
print impm3.shape
nc.createDimension('nhc', impm3.shape[0])
nc.createDimension('jm', impm3.shape[1])
nc.createDimension('im', impm3.shape[2])

dims = ('nhc', 'jm', 'im')
v = nc.createVariable('impm3', 'd', dims)
v[:] = impm3[:]
v.missing_value = -1e30

v = nc.createVariable('impm3b', 'd', dims)
v[:] = impm3b[:]
v.missing_value = -1e30

diff = impm3b - impm3
v = nc.createVariable('diff', 'd', dims)
v[:] = diff[:]
v.missing_value = -1e30


dims = ('jm', 'im')
v = nc.createVariable('conserv1', 'd', dims)
v[:] = conserv1[:]
v.missing_value = -1e30

v = nc.createVariable('conserv1b', 'd', dims)
v[:] = conserv1b[:]
v.missing_value = -1e30

conserv_diff = conserv1b - conserv1
v = nc.createVariable('conserv_diff', 'd', dims)
v[:] = conserv_diff[:]
v.missing_value = -1e30




#v = nc.createVariable('nan3', 'd', dims)
#v[:] = nan3[:]
#v.missing_value = -1e30

nc.close()

#print impm3b - impm3

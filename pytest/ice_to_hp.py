# See what kind of fill-in we get when we compute A^T A,
# where A is hp_to_ice().

import glint2
import netCDF4
import numpy as np

mm = glint2.MatrixMaker()
print '***** Reading new MatrixMaker'
mm.load('mm.nc', 'm')

XM = mm.hp_to_ice('greenland')
print 'XM', XM.shape, XM.nnz

XMt = XM.transpose()
print 'XMt', XMt.shape, XMt.nnz

PP = XMt * XM
print 'PP', PP.shape, PP.nnz

# This is the wrong way, and it gets lots of fill-in to boot
#PP = XM * XMt
#print 'PP', PP.shape, PP.nnz



nc = netCDF4.Dataset('../data/JUL1950.ijhchc1k225.nc')
impm3 = np.zeros(nc.variables['impm_lndice'].shape, dtype='d')
impm3[:] = nc.variables['impm_lndice'][:]
nc.close()

impm2 = glint2.coo_multiply(XM, impm3)
print 'impm2', impm2.shape

f2s = dict()
f2s['greenland'] = impm2
impm3b = mm.ice_to_hp(f2s).reshape(impm3.shape)

nc = netCDF4.Dataset('./ice_to_hp.nc', 'w')
print impm3.shape
nc.createDimension('nhc', impm3.shape[0])
nc.createDimension('jm', impm3.shape[1])
nc.createDimension('im', impm3.shape[2])

dims = ('nhc', 'jm', 'im')
nc.createVariable('impm3', 'd', dims)[:] = impm3[:]
nc.createVariable('impm3b', 'd', dims)[:] = impm3b[:]

nc.close()

#print impm3b - impm3

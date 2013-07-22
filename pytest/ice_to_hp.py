# See what kind of fill-in we get when we compute A^T A,
# where A is hp_to_ice().

import glint2
import netCDF4
import numpy as np
import giss.modele
import sys

mm = glint2.MatrixMaker()
print '***** Reading new MatrixMaker'
mm.load('../data/mmx.nc', 'm')

XM = mm.hp_to_ice('greenland')
print 'XM', XM.shape, XM.nnz

XMt = XM.transpose()
print 'XMt', XMt.shape, XMt.nnz

PP = XMt * XM
print 'PP', PP.shape, PP.nnz

# This is the wrong way, and it gets lots of fill-in to boot
#PP = XM * XMt
#print 'PP', PP.shape, PP.nnz



#nc = netCDF4.Dataset('../data/JUL1950.ijhchc1k225.nc')
nc = netCDF4.Dataset('../data/JUL1956.ijhce4f40-hc40.nc')
impm3 = giss.modele.read_ncvar(nc, 'impm_lndice')
impm3 = impm3[1:,:]			# Remove legacy height point
frac3 = giss.modele.read_ncvar(nc, 'frac')
frac3 = frac3[1:,:]			# Remove legacy height point
print 'frac3.shape = ', frac3.shape
nc.close()

#nc = netCDF4.Dataset('../data/mmx.nc')
#nc.close()


impm2 = glint2.coo_multiply(XM, impm3)
print 'impm2', impm2.shape

f2s = dict()
f2s['greenland'] = impm2
initial = impm3.reshape(-1)
print initial[initial != 0]
#impm3b = mm.ice_to_hp(f2s, initial).reshape(impm3.shape)
impm3b = impm3

# ------------------------------------------------------
# Plot our list of nan
idx3 = [36775, 49735, 50171, 62695, 63131, 75655, 76091, 88615, 88905, 89051, 101575, 101865, 102011, 114535, 114971 ]
nan3 = np.zeros(initial.shape, dtype='d')
nan3[idx3] = 1
nan3 = nan3.reshape(impm3.shape)

nhc = impm3.shape[0]
jm = impm3.shape[1]
im = impm3.shape[2]
n1 = im*jm
for i3 in idx3 :
	ihc = i3 / n1
	i2 = i3 - ihc * n1
	jj = i2 / im
	ii = i2 - jj * im

#	x = i3 / im
#	ii = i3 - x*im
#	ihc = x / jm
#	jj = x - ihc*jm
#	ii = 59
#	jj = 81
	print '%d -> (%d, %d=(%d, %d))' % (i3, ihc,i2,jj,ii)
#	sys.stdout.write('    ')
#	for f in frac3[:,jj,ii] :
#		sys.stdout.write('%e ' % f)
#	print
	print '    ',frac3[:,jj,ii]

# ------------------------------------------------------
nc = netCDF4.Dataset('./ice_to_hp.nc', 'w')
print impm3.shape
nc.createDimension('nhc', impm3.shape[0])
nc.createDimension('jm', impm3.shape[1])
nc.createDimension('im', impm3.shape[2])

dims = ('nhc', 'jm', 'im')
nc.createVariable('impm3', 'd', dims)[:] = impm3[:]
nc.createVariable('impm3b', 'd', dims)[:] = impm3b[:]
nc.createVariable('nan3', 'd', dims)[:] = nan3[:]

nc.close()

#print impm3b - impm3

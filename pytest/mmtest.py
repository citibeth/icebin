import glint2
import giss.searise
import numpy as np
import sys

hpdefs = np.array(range(0,40))*100.0 - 50.0
hcmax = np.array(range(0,40))*100.0

#hpdefs = np.array([100,300,550, 850, 1150, 1450, 1800, 2250, 2750, 4000], dtype='d')
#hcmax = np.array([200,400,700,1000,1300,1600,2000,2500,3000,1000000], dtype='d')

mm = glint2.MatrixMaker()
mm.init('../greenland_2x2_5.nc', hpdefs, hcmax)	# no mask1


searise_fname = '../data/Greenland_5km_v1.1.nc'
(elev2, mask2) = giss.searise.read_elevation2_mask2(searise_fname)
greenland_id = mm.add_ice_sheet('../searise.nc', '../greenland_2x2_5-searise.nc',
	elev2, mask2=mask2, name='greenland')
#greenland2_id = mm.add_ice_sheet('../searise.nc', '../greenland_2x2_5-searise.nc',
#	elev2, mask2=mask2, name='greenland2')

print 'greenland_id = %d' % greenland_id
#print 'greenland2_id = %d' % greenland2_id

mm.realize()

print '***** Writing out MatrixMaker'
nc = glint2.NcFile('mm.nc', 'w')
mm.write(nc, 'm')
nc.close()


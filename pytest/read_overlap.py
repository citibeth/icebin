import glint2
import netCDF4
import numpy as np

np.set_printoptions(threshold='nan')

fname = '../greenland_4x5-searise.nc'

nc = netCDF4.Dataset(fname)

overlap = glint2.read_overlap(nc, 'grid')

print overlap.__dict__

# Looks good!

sum_by_row = np.array(overlap.sum(axis=1)).reshape(-1)
sum_by_col = np.array(overlap.sum(axis=0)).reshape(-1)

print len(sum_by_row)
print sum_by_row

print len(sum_by_col)
print sum_by_col

nc.close()

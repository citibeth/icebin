import sys
import shutil
import netCDF4
import numpy as np

# Constructs a Greenland-only TOPO file from a Whole-Earth TOPO file
# and an Earth-Minus-Greenland TOPO file

rootname = sys.argv[1]

topo_withg_fname = rootname+'-withgr.nc'
topo_nog_fname = rootname+'-nogr.nc'
ofname = rootname+'-onlygr.nc'

shutil.copyfile(topo_withg_fname, ofname)

ncg = netCDF4.Dataset(ofname, 'a')
nc = netCDF4.Dataset(topo_nog_fname, 'r')

varg = {}
vnames = list(ncg.variables.keys())
for vname in vnames:
    varg[vname] = ncg.variables[vname][:]

# Take out just the continental cells, mask the rest
diff = varg['FOCEAN'] - nc.variables['FOCEAN'][:]
mask = (diff != 0.)
notmask = np.logical_not(mask)

varg['FOCEAN'][mask] = 0.0

for vname in vnames:
    varg[vname][notmask] = np.nan
    ncg.variables[vname][:] = varg[vname]

ncg.close()




## See where we chaned the ocean mask
#
#
#print('MASK', np.sum(np.sum(mask)), np.sum(np.sum(notmask)))
#
#varg['FOCEAN'][mask] = 0.
#
#
## By definition, we only affect land cells
#varg['dZOCEN'][:] = 0
#
#val = ncg.variables['ZSOLDG'][:]
#varg['ZSOLDG'][mask] = val[mask]
#
#val = ncg.variables['ZSGLO'][:]
#varg['ZSGLO'][mask] = val[mask]
#
#
#

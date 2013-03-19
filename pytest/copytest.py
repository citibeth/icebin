import giss.ncutil
import netCDF4

nc0 = netCDF4.Dataset('mm.nc', 'r')
ncout = netCDF4.Dataset('mmx.nc', 'w')

cp = giss.ncutil.copy_nc(nc0, ncout)


# Define our own variables
ncout.variables['m.greenland.info'].ice_model = 'DISMAL'
ncout.variables['m.greenland2.info'].ice_model = 'DISMAL'

cp.copy_data()

ncout.close()
nc0.close()

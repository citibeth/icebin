import netCDF4
import giss

nc = netCDF4.Dataset('hp2hc.nc')
hp2hc = giss.read_coo_matrix(nc, 'hp2hc')
nc.close()

#print hp2hc.__dict__

nc = netCDF4.Dataset('fhc-00.nc')
jm = len(nc.dimensions['jm'])
im = len(nc.dimensions['im'])
nhc = len(nc.dimensions['nhc'])
rows_i = nc.variables['rows_i'][:]
rows_j = nc.variables['rows_j'][:]
rows_k = nc.variables['rows_k'][:]
cols_i = nc.variables['cols_i'][:]
cols_j = nc.variables['cols_j'][:]
cols_k = nc.variables['cols_k'][:]
vals = nc.variables['vals'][:]
nc.close()

nele = len(rows_i)

for ix in range(0,nele) :
	vdiff = vals[ix] - hp2hc.data[ix]
	print '(%d - %d) --> ([%d %d %d] - [%d %d %d]) %f' % (hp2hc.row[ix], hp2hc.col[ix], rows_i[ix], rows_j[ix], rows_k[ix], cols_i[ix], cols_j[ix], cols_k[ix], vals[ix])

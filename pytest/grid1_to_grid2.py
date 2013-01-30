import sys
import glint2
import netCDF4
import numpy as np
import matplotlib.pyplot
import giss.basemap
import numpy.ma as ma
import giss.util
import giss.modele
import sys

np.set_printoptions(threshold='nan')

#grid1_fname = '../greenland_4x5.nc'
#grid2_fname = '../searise50.nc'
#exgrid_fname = '../greenland_4x5-searise50.nc'
#data_fname = '../data/JUL1954.ijhcncar45.nc'

grid1_fname = '../greenland_2x2_5.nc'
grid2_fname = '../searise.nc'
exgrid_fname = '../greenland_2x2_5-searise.nc'
data_fname = '../data/JUL1954.ijhchctest225.nc'

# ---------- Read from Exchange Grid
nc = netCDF4.Dataset(exgrid_fname)
overlap = glint2.read_overlap(nc, 'grid')
sproj = nc.variables['grid.info'].projection
nc.close()

# ---------- Read from grid1
grid1 = glint2.Grid(grid1_fname, 'grid')
n1 = grid1.ncells_full

proj_area1 = grid1.get_proj_area(sproj)
native_area1 = grid1.get_native_area()

nc = netCDF4.Dataset(grid1_fname)
info_var = nc.variables['grid.info']
n1 = info_var.__dict__['cells.num_full']
nlon = info_var.nlon
nlat = info_var.nlat
#vname = 'grid'
#cells_index = nc.variables[vname + '.cells.index'][:]
#cells_area = nc.variables[vname + '.cells.area'][:]
nc.close()

#print '****** Area Ratios'
#good = np.logical_not(np.logical_or(np.isnan(proj_area1), np.isnan(native_area1)))
#print native_area1[good] / proj_area1[good]

# ========= Read some data
nc = netCDF4.Dataset(data_fname)
val1h = giss.modele.read_ncvar(nc, 'impm_lndice') * 86400.
nc.close()
nhc = val1h.shape[0]
val1 = val1h[0,:].astype(np.double)

# Re-make our data for testing purposes
if True :
	for j in range(0,val1.shape[0]) :
		for i in range(0, val1.shape[1]) :
			val1[j,i] = i+j
			if  j == 41 : val1[j,i] = np.nan

# =========== Mask out overlap matrix
# (to get correct conservation numbers)
mask1 = (np.isnan(val1)).astype(np.int32)
mask2 = None
overlap_m = glint2.mask_out(overlap, mask1, mask2)

# ========== Now try it out

area1p = giss.util.sum_by_rows(overlap_m)
area2 = giss.util.sum_by_cols(overlap_m)

# Matrix to go from PROJECTED grid1 to grid2
mat = glint2.grid1_to_grid2(overlap_m)

# Diagonal matrix converts from a vector of values in the native
# grid to a vector of values in the projected grid.
diag = glint2.proj_native_area_correct(grid1, sproj, 'n2p')
glint2.multiply_bydiag(mat, diag)

# Regrid
val1 = val1.reshape(-1)
val2 = glint2.coo_multiply(mat, val1, fill=np.nan)

# =============================================

total1n = np.nansum(area1p * (native_area1 / proj_area1) * val1)
total1p = np.nansum(area1p * val1)
total2 = np.nansum(area2 * val2)

# These should be the same for conservatie regridding
print 'total1n = %f' % total1n
print 'total1p = %f' % total1p
print 'total2  = %f' % total2

# =============================================

#sys.exit(0)

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))

# ---------- Plot the GCM Grid
ax = figure.add_subplot(121)
plotter = glint2.Grid_read_plotter(grid1_fname, 'grid')
mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
im = plotter.pcolormesh(mymap, val1)
mymap.colorbar(im, "right", size='5%')
mymap.drawcoastlines()



# ---------- Plot the ice grid
ax = figure.add_subplot(122)
plotter = glint2.Grid_read_plotter(grid2_fname, 'grid')
mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
im = plotter.pcolormesh(mymap, val2)
mymap.colorbar(im, "right", size='5%')
mymap.drawcoastlines()



# Also show on screen
matplotlib.pyplot.show()


# Looks good!



sys.exit(0)


sum_by_row = np.array(overlap_m.sum(axis=1)).reshape(-1)
sum_by_col = np.array(overlap_m.sum(axis=0)).reshape(-1)

print len(sum_by_row)
print sum_by_row

print len(sum_by_col)
print sum_by_col


import glint2
import netCDF4
import numpy as np
import matplotlib.pyplot
import giss.basemap
import numpy.ma as ma
import giss.util
import giss.modele
import sys

#np.set_printoptions(threshold='nan')

grid1_fname = '../greenland_4x5.nc'
grid2_fname = '../searise.nc'
exgrid_fname = '../greenland_4x5-searise.nc'
data_fname = '../data/JUL1954.ijhcncar45.nc'

nc = netCDF4.Dataset(exgrid_fname)
overlap = glint2.read_overlap(nc, 'grid')
nc.close()


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

area1 = giss.util.sum_by_rows(overlap_m)
area2 = giss.util.sum_by_cols(overlap_m)

mat = glint2.grid1_to_grid2(overlap_m)

# Regrid
val1 = val1.reshape(-1)
val2 = glint2.coo_multiply(mat, val1, fill=np.nan)

# =============================================

area1 = giss.util.sum_by_rows(overlap_m)
area2 = giss.util.sum_by_cols(overlap_m)

total1 = np.nansum(area1 * val1)
total2 = np.nansum(area2 * val2)

# These should be the same for conservatie regridding
print 'total1 = %f' % total1
print 'total2 = %f' % total2

# =============================================

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


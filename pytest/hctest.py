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
import giss.searise

np.set_printoptions(threshold='nan')

#grid1_fname = '../greenland_4x5.nc'
#data_fname = '../data/JUL1954.ijhcncar45.nc'
#grid2_fname = '../searise.nc'
#searise_fname = '../data/Greenland_5km_v1.1.nc'
#exgrid_fname = '../greenland_4x5-searise.nc'

grid1_fname = '../greenland_2x2_5.nc'
grid2_fname = '../searise.nc'
exgrid_fname = '../greenland_2x2_5-searise.nc'
data_fname = '../data/JUL1950.ijhchc1k225.nc'
searise_fname = '../data/Greenland_5km_v1.1.nc'

#hcmax = np.array([200,400,700,1000,1300,1600,2000,2500,3000,10000], dtype='d')
#hpdefs = np.array([100,300,550,850,1150,1450,1800,2250,2750,4000], dtype='d')

hcmax = [25.]
hcmax.extend(np.array(range(1,40))*100.0)
hcmax = np.array(hcmax)
print hcmax

hpdefs = np.array([
           0.,    50.,   150.,   250.,   350.,   450.,   550.,   650.,
         750.,   850.,   950.,  1050.,  1150.,  1250.,  1350.,  1450.,
        1550.,  1650.,  1750.,  1850.,  1950.,  2050.,  2150.,  2250.,
        2350.,  2450.,  2550.,  2650.,  2750.,  2850.,  2950.,  3050.,
        3150.,  3250.,  3350.,  3450.,  3550.,  3650.,  3750.,  3850.])


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

# ----------- Read from Searise
(elev2, mask2) = giss.searise.read_elevation2_mask2(searise_fname)


# ========= Read some data
nc = netCDF4.Dataset(data_fname)
val1h = giss.modele.read_ncvar(nc, 'tsurf_lndice').astype('d') # * 86400.
nc.close()
nhc = val1h.shape[0]

# ========== Height-classify overlap matrix
print len(hcmax), nhc
overlaph = glint2.height_classify(overlap, elev2, hcmax)

# =========== Mask out overlap matrix
# (to get correct conservation numbers)
mask1h = (np.isnan(val1h)).astype(np.int32)
#mask2 = mask2
overlaph_m = glint2.mask_out(overlaph, mask1h, mask2)

# ========== Now try it out

area1hc = giss.util.sum_by_rows(overlaph_m)
area2 = giss.util.sum_by_cols(overlaph_m)

# Matrix to go from PROJECTED grid1 to grid2
mat_1h_to_2 = glint2.grid1_to_grid2(overlaph_m)

# ...and back
mat_2_to_1h = glint2.grid2_to_grid1(overlaph_m)


# Diagonal matrix converts from a vector of values in the native
# grid to a vector of values in the projected grid.
if False :
	diag = glint2.proj_native_area_correct(grid1, sproj, 'n2p')
	diag = np.tile(diag, (1,nhc)).reshape(-1)
	glint2.multiply_bydiag(mat_1h_to_2, diag)

# Regrid
val1h = val1h.reshape(-1)
val2_hc = glint2.coo_multiply(mat_1h_to_2, val1h, fill=np.nan)

# ========== Interpolate to ice grid, using height point smoothing
# Collapse height-classified mask1 into non-height-classified version
mask1 = np.array(mask1h.sum(axis=0)).reshape(mask1h.shape[1:]).astype(np.int32)
# Mask overlap matrix using this (for correct conservation checking)
overlap_m = glint2.mask_out(overlap, mask1, mask2)
mat_1p_to_2 = glint2.hp_interp(overlap_m, elev2, hpdefs)
val2_hp = glint2.coo_multiply(mat_1p_to_2, val1h, fill=np.nan)

val1hx = glint2.coo_multiply(mat_2_to_1h, val2_hp, fill=np.nan)
val1hx_plt = glint2.coo_multiply(mat_1h_to_2, val1hx, fill=np.nan)

#print 'Computing AlmostI'
#AlmostI = mat_2_to_1h * mat_1p_to_2
#print 'Done Computing AlmostI'
#print type(AlmostI)
#print AlmostI.shape


# =============================================

# A bit more complication extending native_area1 to full hc scheme
factor1h = np.tile(native_area1 / proj_area1, (nhc,1)).reshape(-1)

total1n = np.nansum(area1hc * factor1h * val1h)

# These should be the same for conservatie regridding
print 'total1_native = %f' % total1n
print 'total1_proj   = %f' % np.nansum(area1hc * val1h)
print 'total2        = %f' % np.nansum(area2 * val2_hc)
print
print 'total2_hp     = %f' % np.nansum(area2 * val2_hp)
print 'total1hx      = %f' % np.nansum(area1hc * val1hx)

# =============================================

# Plot multiple plots on one page
figure = matplotlib.pyplot.figure(figsize=(8.5,11))

# ---------- Plot field with height classes
ax = figure.add_subplot(121)
plotter = glint2.Grid_read_plotter(grid2_fname, 'grid')
mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
im = plotter.pcolormesh(mymap, val2_hc)
mymap.colorbar(im, "right", size='5%')
mymap.drawcoastlines()

# ---------- Plot with height classes, after going through ice grid
ax = figure.add_subplot(122)
plotter = glint2.Grid_read_plotter(grid2_fname, 'grid')
mymap = giss.basemap.greenland_laea(ax)
#mymap = giss.basemap.north_laea(ax)
im = plotter.pcolormesh(mymap, val1hx_plt)
mymap.colorbar(im, "right", size='5%')
mymap.drawcoastlines()

## ---------- Plot the ice grid
#ax = figure.add_subplot(122)
#plotter = glint2.Grid_read_plotter(grid2_fname, 'grid')
#mymap = giss.basemap.greenland_laea(ax)
##mymap = giss.basemap.north_laea(ax)
#im = plotter.pcolormesh(mymap, val2_hp)
#mymap.colorbar(im, "right", size='5%')
#mymap.drawcoastlines()



# Also show on screen
matplotlib.pyplot.show()


# Looks good!



sys.exit(0)


sum_by_row = np.array(overlaph_m.sum(axis=1)).reshape(-1)
sum_by_col = np.array(overlaph_m.sum(axis=0)).reshape(-1)

print len(sum_by_row)
print sum_by_row

print len(sum_by_col)
print sum_by_col


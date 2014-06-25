import glint2
import netCDF4
import giss.basemap
import matplotlib
import subprocess
import sys

steps = set([int(x) for x in sys.argv[1:]])

if 1 in steps :
	# ================ Step 1: Create some grids
	subprocess.call(['../../sbin/greenland_2x2_5'])
	subprocess.call(['../../sbin/searise_g', '20'])
	subprocess.call(['../../sbin/searise_g', '5'])
	subprocess.call(['../../sbin/overlap', 'greenland_2x2_5.nc', 'searise_g5.nc'])
	subprocess.call(['../../sbin/overlap', 'greenland_2x2_5.nc', 'searise_g20.nc'])


if 2 in steps :
	# =============== Step 2: Plot grid outlines

	# Set up the page
	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	# figure.set_size_inches(11., 8.5)		# US Letter
	figure.set_size_inches(8.267,11.692)	# A4

	# -------- First plot: grid1

	# Read the grid
	nc = netCDF4.Dataset('greenland_2x2_5.nc', 'r')
	grid1 = glint2.pyGrid(nc, 'grid')
	nc.close()

	# Plot it!
	ax = figure.add_subplot(131)
	basemap = giss.basemap.greenland_laea(ax=ax)
	grid1.plot(basemap)

	# ---------- Second plot: grid2
	# Read the grid
	nc = netCDF4.Dataset('searise_g20.nc', 'r')
	grid2 = glint2.pyGrid(nc, 'grid')
	nc.close()

	ax = figure.add_subplot(132)
	basemap = giss.basemap.greenland_laea(ax=ax)
	grid2.plot(basemap, linewidth=.5, color='gray')


	# ---------- Third plot: exchange grid
	nc = netCDF4.Dataset('greenland_2x2_5-searise_g20.nc', 'r')
	grid4 = glint2.pyGrid(nc, 'grid')
	nc.close()

	ax = figure.add_subplot(133)
	basemap = giss.basemap.greenland_laea(ax=ax)
	grid4.plot(basemap, linewidth=.5)


	# ------------ Render the page
	figure.savefig('grids.pdf', dpi=100, transparent=False)

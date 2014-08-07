import sys
import subprocess
import urllib
import netCDF4
import numpy as np
import matplotlib
import glint2
import giss.basemap
import giss.searise
import giss.ncutil
import giss.modele

steps = set([int(x) for x in sys.argv[1:]])

def regrid_5_to_20(ff) :
	"""Regrids from the 5km SeaRISE grid to the 20km SeaRISE grid
	Just sub-sample, don't do anything fancy"""

	print (ff.shape[0]/4+1, ff.shape[1]/4+1)
	ret = np.zeros((ff.shape[0]/4+1, ff.shape[1]/4+1), dtype=ff.dtype)
	for j in range(0, ret.shape[0]) :
		for i in range(0, ret.shape[1]) :
			ret[j,i] = ff[j*4,i*4]
	print 'regrid: %s to %s' % (ff.shape, ret.shape)
	return ret
# --------------------------------------------------------------------
class Glint2Files(object) :
	"""Sets up names of Glint2 files"""
	def __init__(self, grid1, grid2, ice_dx, nhc, ice_model, searise_fname) :
		self.grid1 = grid1
		self.grid2 = grid2
		self.ice_dx = ice_dx
		self.nhc = nhc
		self.searise_fname = searise_fname

		self.grid1_fname = grid1 + '.nc'
		self.grid2_fname = grid2 + '.nc'
		self.grid4_fname = '%s-%s.nc' % (grid1, grid2)	# Exchange grid file
		self.glint2_in_fname = '%s-%s-%d-%s.nc' % (grid1, grid2, nhc, ice_model)

		self.gcm_ijhc_fname = '130516-Regrid_examples_e4f40-hc40-g5_JUL1956.ijhce4f40-hc40.nc'
		self.gcm_aij_fname = '130516-Regrid_examples_e4f40-hc40-g5_JUL1956.aije4f40-hc40-prec.nc'


# --------------------------------------------------------------------
def make_glint2_in(gf) :
	"""
	gfiles : Glint2Files
	"""

	# DISMAL = "Demonstration Ice Sheet Model and Landice"
	# DISMAL does nothing, it just writes its inputs out to disk every time it is called.
	ice_model = 'DISMAL'

	# Output: glint2 config file
	tmp_ofname = 'tmp.nc'


	# ========== Set up grid1 and height points
	grid1_name = 'modele_ll_g2x2_5'
	hpdefs = np.array(range(0,gf.nhc))*(4000. / gf.nhc) - (4000. / (gf.nhc * 2))
	mm = glint2.MatrixMaker()
	print hpdefs
	mm.init(gf.grid1_fname, 'MODELE', hpdefs) # no mask1



	# ========= Add each grid2

	# --- Greenland
	(elev2_5, mask2_5) = giss.searise.read_elevation2_mask2(searise_fname)

	# Sub-sample the 5km data if we really wanted 20km
	if gf.ice_dx == 20 :
		elev2 = regrid_5_to_20(elev2_5)
		mask2 = regrid_5_to_20(mask2_5)
	else :
		elev2 = elev2_5
		mask2 = mask2_5

	print 'gf = ', gf.__dict__
	greenland_id = mm.add_ice_sheet(gf.grid2_fname, gf.grid4_fname,
	        elev2, mask2=mask2, name='greenland')

	# ========== Finish up and write out

	mm.realize()

	print '***** Writing MatrixMaker'
	nc = glint2.NcFile(tmp_ofname, 'w')
	mm.write(nc, 'm')
	nc.close()

	# ==========================================================
	# Add to the MatrixMaker netCDF file to obtain a full GLINT2 configuration

	nc_write = []

	ncout = netCDF4.Dataset(gf.glint2_in_fname, 'w')


	# ----------- Copy our existing netCDF file
	nc0 = netCDF4.Dataset(tmp_ofname, 'r')

	cp = giss.ncutil.copy_nc(nc0, ncout,
		attrib_filter = lambda x : not x.endswith('comment'))	# Remove comments from final file for better readibility
	cp.define_vars()

	def _write_icemodel() :
		cp.copy_data()
		nc0.close()
	nc_write.append(_write_icemodel)

	# ------------ Add the name of the ice model we wish to use
	# Define our own variables
	greenland_info = ncout.variables['m.greenland.info']
	greenland_info.ice_model = ice_model
	greenland_info.interp_grid = 'EXCH'
	greenland_info.interp_style = 'Z_INTERP'
	# TODO: This should be set to 'true'.  ModelE should be re-run, using the PISM spun-up ice model.
	greenland_info.update_elevation = 'false';
	ncout.variables['m.info'].hc_index_type = 'MODELE'
	ncout.variables['m.info'].gcm_out_file = 'modele_out.nc'

	# ------------ Add DISMAL parameters
	if True :
		dismal = ncout.createVariable('m.greenland.dismal', 'd')
		dismal.output_dir = 'dismal_out2';

	# ------------ Add PISM parameters
	if True:
		pism = ncout.createVariable('m.greenland.pism', 'd')
		pism.i = 'std-greenland/g20km_10ka.nc'
		pism.skip = ''
		pism.skip_max = '10'
		pism.surface = 'given'
		pism.surface_given_file = 'std-greenland/pism_Greenland_5km_v1.1.nc'
		pism.calving = 'ocean_kill'
		pism.ocean_kill_file = 'std-greenland/pism_Greenland_5km_v1.1.nc'
		pism.sia_e = '3.0'
		pism.ts_file = 'ts_g20km_10ka_run2.nc'
		pism.ts_times = '0:1:1000'
		pism.extra_file = 'ex_g20km_10ka_run2.nc'
		pism.extra_times = '0:.1:1000'
		pism.extra_vars = 'climatic_mass_balance,ice_surface_temp,diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,csurf,mask,thk,topg,usurf,discharge_flux_cumulative'
		pism.o = 'g20km_10ka_run2.nc'

	if True :
		modele = ncout.createVariable('m.greenland.modele', 'd')
		# modele.coupling_type = 'NEUMANN_BC'		# Heat flux over timestep specified
		modele.coupling_type = 'DIRICHLET_BC'		# Constant T over timestep

	# ----------------- Write data and finish up
	for fn in nc_write : fn()
	ncout.close()
# ------------------------------------------------------------------



# ========== Set-up
grid1_name = 'greenland_2x2_5'
grid2_name = 'searise_g5'
grid4_name = grid1_name + '-' + grid2_name
searise_fname = 'Greenland_5km_dev1.2.nc'

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
	nc = netCDF4.Dataset(grid1_name + '.nc', 'r')
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

if 3 in steps :
	# =========== Step 3: Download SeaRISE file
	url = 'http://websrv.cs.umt.edu/isis/images/e/e9/%s' % searise_fname
	urllib.urlretrieve(url, searise_fname)


if 4 in steps :
	# ============ Step 4: Read and plot elevations
	nc = netCDF4.Dataset(searise_fname, 'r')
	thk = nc.variables['thk'][0,:]
	topg = nc.variables['topg'][0,:]
	nc.close()
	elev2 = thk + topg

	# Get a plotter for the ice grid
	plotter2 = glint2.Plotter2(fname=grid2_name + '.nc', vname='grid')

	# Plot it!
	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	ax = figure.add_subplot(111)
	basemap = giss.basemap.greenland_laea(ax=ax)

	giss.plot.plot_var(basemap=basemap,
		plotter=plotter2, val=elev2, title='Ice Elevation',
		cb_args = {})

	figure.savefig('elev2.png', transparent=True, dpi=300)


if 5 in steps :
	# ============ Step 5: Create Glint2 Configuration File
	gf = Glint2Files(grid1_name, grid2_name, 5, 40, 'DISMAL', searise_fname)
	make_glint2_in(gf)

if 6 in steps:
	# ============ Step 6: Run the GCM to create output file
	gf = Glint2Files(grid1_name, grid2_name, 5, 40, 'DISMAL', searise_fname)
	print \
		"""We will not run the climate model here, just
		use output from a previous run.
		See the file: %s""" % gf.gcm_ijhc_fname

if 7 in steps:
	# ============ Step 7: Load GCM output, regrid and plot
	gf = Glint2Files(grid1_name, grid2_name, 5, 40, 'DISMAL', searise_fname)
	mm = glint2.MatrixMaker(gf.glint2_in_fname, 'm', correct_area1=False)
	nc = netCDF4.Dataset(gf.glint2_in_fname)
	mask2 = (nc.variables['m.greenland.mask2'][:] != 0)
	nc.close()

	nc = netCDF4.Dataset(gf.gcm_ijhc_fname, 'r')
	# Convert to non-SI units,kg/(day m^2)
	impm3 = giss.modele.read_ncvar(nc, 'impm_lndice') * 86400.
	impm3 = impm3[1:,:,:]			# Remove legacy height point
	frac3 = giss.modele.read_ncvar(nc, 'frac')
	frac3 = frac3[1:,:,:]			# Remove legacy height point
	nc.close()

	nc = netCDF4.Dataset(gf.gcm_aij_fname, 'r')
	prec1 = giss.modele.read_ncvar(nc, 'prec')	# precipitation
	mask1 = (np.sum(frac3,0) == 0)
	prec1[mask1] = 0	# Ignore Antarctica, the tropics, etc
	nc.close()

	# These are already set in the config file, but they can be
	# re-set on the fly.
	# mm.set_interp_grid('EXCH')	# Interpolate to the exchange grid, not directly to ice
	# mm.set_interp_style('Z_INTERP')

	# Convert impm3 to ice grid (2)
	XM = mm.hp_to_iceinterp('greenland', dest='ICE')
	impm32 = XM.dot(impm3.reshape(-1))

	# ...and to atmosphere grid
	RM = mm.hp_to_atm()
	impm31 = RM.dot(impm3.reshape(-1))
	RXp = mm.iceinterp_to_atm('greenland', src='ICE')
	impm321 = RXp.dot(impm32)

	# Now plot the two
	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	ax = figure.add_subplot(131)
	basemap = giss.basemap.greenland_laea(ax=ax)

	cb_args = {}
	plot_args = {
		'cmap' : giss.plot.cpt('giss-cpt/BlRe.cpt', reverse=True).cmap,
		'vmin' : -15.,
		'vmax' : 5.,
		'norm' : giss.plot.AsymmetricNormalize()
	}

	plotter2 = glint2.Plotter2(fname=gf.glint2_in_fname, vname='m.greenland.grid2')
	val = impm32
	val[mask2] = np.nan
	ret = giss.plot.plot_var(basemap=basemap,
		plotter=plotter2, val=val, title='SMB (mm/d): E->I',
		plot_args = plot_args, cb_args = cb_args)

	ax = figure.add_subplot(132)
	basemap = giss.basemap.greenland_laea(ax=ax)

	plotter1 = giss.modele.plotters.get_byname('2x2.5')
	val = impm31
	#val[mask1.reshape(-1)] = np.nan
	giss.plot.plot_var(basemap=basemap,
		plotter=plotter1, val=val, title='SMB (mm/d): E->A',
		plot_args = plot_args, cb_args = cb_args)

	ax = figure.add_subplot(133)
	basemap = giss.basemap.greenland_laea(ax=ax)
	val = impm321
	#val[mask1.reshape(-1)] = np.nan
	giss.plot.plot_var(basemap=basemap,
		plotter=plotter1, val=val, title='SMB (mm/d): E->I->A',
		plot_args = plot_args, cb_args = cb_args)


	figure.savefig('impm_a.png', transparent=True, dpi=300)

	# ------------------------------------------
	# Regrid back to elevation grid

	impm323 = mm.iceinterp_to_hp(
		{'greenland' : impm32},
		src='ICE',
		qp_algorithm='SINGLE_QP')
	impm313 = mm.atm_to_hp(impm31)

	# Plot things in "plain" elevaition space
	mm.set_interp_style('ELEV_CLASS_INTERP')
	XM_ec = mm.hp_to_iceinterp('greenland', dest='ICE')

	figure = matplotlib.pyplot.figure(figsize=(11,8.5))
	ax = figure.add_subplot(121)
	basemap = giss.basemap.greenland_laea(ax=ax)

	val = XM_ec.dot(impm313)
	val[mask2] = np.nan
	ret = giss.plot.plot_var(basemap=basemap,
		plotter=plotter2, val=val, title='SMB (mm/d): E->A->E',
		plot_args = plot_args, cb_args = cb_args)

	ax = figure.add_subplot(122)
	basemap = giss.basemap.greenland_laea(ax=ax)

	val = XM_ec.dot(impm323)
	giss.plot.plot_var(basemap=basemap,
		plotter=plotter2, val=val, title='SMB (mm/d): E->I->E',
		plot_args = plot_args, cb_args = cb_args)

	figure.savefig('impm_b.png', transparent=True, dpi=300)




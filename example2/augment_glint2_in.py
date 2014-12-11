import glint2
import giss.searise
import numpy as np
import sys
import giss.ncutil
import netCDF4
import os
import os.path

ice_dx=20		# 20km ice grid
#ice_dx=5		# 5km ice grid
ice_model = 'PISM'
#ice_model = 'DISMAL'

ifname = os.path.join('/Users/rpfische/exp/130613-pism', 'modele_ll_g2x2_5-searise_g%d-40-base.nc' % ice_dx)
ofname = 'modele_ll_g2x2_5-searise_g%d-40-%s.nc' % (ice_dx, ice_model)



nc_write = []

ncout = netCDF4.Dataset(ofname, 'w')

# ----------- Copy our existing netCDF file
nc0 = netCDF4.Dataset(ifname, 'r')

cp = giss.ncutil.copy_nc(nc0, ncout)
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
ncout.variables['m.info'].gcm_in_file = 'modele_in.nc'

# ------------ Add PISM parameters
#if ice_model == 'DISMAL' :
if True :
	dismal = ncout.createVariable('m.greenland.dismal', 'd')
	dismal.output_dir = 'dismal_out2';

# PISM dev branch
#if True and ice_model == 'PISM' :
if True:
	pism = ncout.createVariable('m.greenland.pism', 'd')
	pism.i = 'std-greenland/g20km_10ka.nc'
	pism.skip = ''
	pism.skip_max = '10'
#	pism.ys = '0'
#	pism.ye = '100'
	pism.surface = 'given'
	pism.surface_given_file = 'std-greenland/pism_Greenland_5km_v1.1.nc'
	pism.calving = 'ocean_kill'
	pism.ocean_kill_file = 'std-greenland/pism_Greenland_5km_v1.1.nc'
	pism.sia_e = '3.0'
	pism.ts_file = 'ts_g20km_10ka_run2.nc'
	pism.ts_times = '0:1:1000'
	pism.extra_file = 'ex_g20km_10ka_run2.nc'
	pism.extra_times = '0:.1:1000'
	pism.extra_vars = 'climatic_mass_balance,ice_surface_temp,diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,csurf,mask,thk,topg,usurf'
	pism.o = 'g20km_10ka_run2.nc'

if True :
	modele = ncout.createVariable('m.greenland.modele', 'd')
#	modele.coupling_type = 'NEUMANN_BC'			# Heat flux over timestep specified
	modele.coupling_type = 'DIRICHLET_BC'		# Constant T over timestep


# PISM stable0.5 branch
if False and ice_model == 'PISM' :
	pism = ncout.createVariable('m.greenland.pism', 'd')
	pism.bed_def = 'lc'
	pism.config_override = 'searise-greenland/searise_config.nc'
	pism.ssa_sliding = ''
	pism.thk_eff = ''
	pism.topg_to_phi = '5.0,20.0,-300.0,700.0'
	pism.ocean_kill = ''
	pism.acab_cumulative = ''
	pism.skip = '200'
	pism.i = 'searise-greenland/g20km_0_ftt.nc'
	pism.ocean = 'constant'
	pism.atmosphere = 'searise_greenland'
	pism.surface = 'given'
	pism.surface_given_file = 'pism_in.nc'		# This is what we don't need...
	pism.o = 'pism_out2.nc'
#	pism.ys = '0'
#	pism.ye = '100'
#	pism.reference_date = '1951-01-01'
	pism.extra_file = 'pism_ex.nc'
	pism.extra_times = '0:.1:1000'
	pism.extra_vars = 'climatic_mass_balance,ice_surface_temp,usurf,topg,thk,bmelt,bwat,bwp,mask,velsurf,wvelsurf,velbase,wvelbase,tempsurf,tempbase,diffusivity,acab_cumulative,cbase,csurf,tempicethk_basal,tauc,temppabase'
	pism.ts_file = 'pism_ts.nc'
	pism.ts_times = '0:1:1000'
	pism.ts_vars = 'ivol,iareag,iareaf'


# ----------------- Write data and finish up
for fn in nc_write : fn()
ncout.close()

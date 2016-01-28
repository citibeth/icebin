# Creates the base of a ICEBIN input file (a MatrixMaker output file)
# We will append stuff to this, based on the ice model used (PISM, etc).

import icebin
import ibmisc

import giss.pism
import numpy as np
import sys
import giss.ncutil
import netCDF4
import os
import os.path


def make_icebin_in_base(grid_dir, gridA_name, grid2_name, pism_spinup_fname, ofname):
	# gridA_name = 'modele_ll_g2x2_5'
	# grid2_name = 'searise_g%d' % ice_dx
	# DATA_PATH = os.environ['DATA_PATH']
	# 	pism_spinup_fname = os.path.join(DATA_PATH, 'searise/Greenland_5km_v1.1.nc')

	ice_dx=20		# 20km ice grid
	#ice_dx=5		# 5km ice grid


#	ofname = 'modele_ll_g2x2_5-searise_g%d-40-base.nc' % ice_dx


	# Regrids data from the 5km SeaRISE grid to the 20km SeaRISE grid
	# Just sub-sample, don't do anything fancy
	def regrid_5_to_20(ff) :
		print (ff.shape[0]/4+1, ff.shape[1]/4+1)
		ret = np.zeros((ff.shape[0]/4+1, ff.shape[1]/4+1), dtype=ff.dtype)
		for j in range(0, ret.shape[0]) :
			for i in range(0, ret.shape[1]) :
				ret[j,i] = ff[j*4,i*4]
		print('regrid: %s to %s' % (ff.shape, ret.shape))
		return ret

	# ========== Set up gridA and height points
	gridA_fname = os.path.join(grid_dir, gridA_name + '.nc')
	hpdefs = np.array(range(0,40))*100.0 - 50.0
	mm = icebin.GCMRegridder(gridA_fname, 'grid', hpdefs, True)


	# ========= Add each grid2
	if False:
		# --- Greenland
		# pism_spinup_fname = os.path.join(DATA_PATH, 'searise/Greenland_5km_v1.1.nc')
		print('PISM spinup file: {}'.format(pism_spinup_fname))
		(elev2, mask2) = giss.pism.read_elevation2_mask2(pism_spinup_fname)
	
		grid2_fname = os.path.join(grid_dir, '%s.nc' % grid2_name)
		overlap_fname = os.path.join(grid_dir, '%s-%s.nc' % (gridA_name, grid2_name))
	
		print('mask2',mask2.shape)
		greenland_id = mm.add_ice_sheet(grid2_fname, overlap_fname,
		        elev2, mask2=mask2, name='greenland')

	# ========== Finish up and write out
	print('Writing: {}'.format(ofname))
	ncio = ibmisc.NcIO(ofname, 'replace')
	mm.ncio(ncio, 'm')
	ncio.close()

make_icebin_in_base(
	'/Users/rpfische/git/icebin/build/grids',
	'modele_ll_g2x2_5',
	'sr_g20_pism',
	'pism_spinup.nc',
	'x.nc')

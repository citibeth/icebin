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
import sys
import pickle

def make_icebin_in_base(grid_dir, gridA_name, gridI_name, pism_spinup_fname, ofname):
    # gridA_name = 'modele_ll_g2x2_5'
    # gridI_name = 'searise_g%d' % ice_dx
    # DATA_PATH = os.environ['DATA_PATH']
    #   pism_spinup_fname = os.path.join(DATA_PATH, 'searise/Greenland_5km_v1.1.nc')

    ice_dx=20       # 20km ice grid
    #ice_dx=5       # 5km ice grid

    # ========== Set up gridA and height points
    gridA_fname = os.path.join(grid_dir, gridA_name + '.nc')
    hpdefs = np.array(range(0,40))*100.0 - 50.0
    mm = icebin.GCMRegridder(gridA_fname, 'grid', hpdefs, True)


    # ========= Add each Ice Sheet

    # --- Greenland
    print('PISM spinup file: {}'.format(pism_spinup_fname))
    (elevI, maskI) = giss.pism.read_elevI_maskI(pism_spinup_fname)
    print('len(elevI) = ', elevI.shape)

    gridI_fname = os.path.join(grid_dir, '%s.nc' % gridI_name)
    overlap_fname = os.path.join(grid_dir, '%s-%s.nc' % (gridA_name, gridI_name))
    
    print('maskI',maskI.shape)
    mm.add_sheet('greenland',
        gridI_fname, 'grid',
        overlap_fname, 'exgrid',
        'Z_INTERP',
        elevI, maskI)

    # ========== Finish up and write out
    print('Writing: {}'.format(ofname))
    ncio = ibmisc.NcIO(ofname, 'replace')
    mm.ncio(ncio, 'm')
    ncio.close()

make_icebin_in_base( \
    '.',
    'modele_ll_g2x2_5',
    'sr_g20_pism',
    'elev_mask.nc',
    './icebin_in.nc')
#   './modele_ll_g2x2_5-sr_g20_pism-base.nc')

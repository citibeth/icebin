
# Set to your Python 3 command.  Should not need to be set if an
# appropriate Python3 environment module is already loaded.
PYTHON = python3

# Set this to the directory where your PISM spinup file resides.
# (NOT necessary to run the tests in this directory)
PISM_SPINUP_DIR = /Users/rpfische/exp/151014-integrate/build/std-greenland

# --------------------------------------------------------------

# Run the main test.
all: test_diagonal.png test_elevation.png sr_g20_pism.pdf modele_ll_g2x2_5.pdf modele_ll_g2x2_5-sr_g20_pism.pdf matrices.pik

clean:
	rm *.nc *.pdf *.png

# =========== Generate stuff we check into git
# This requires external software (PISM), and should not normally need
# to be run.

# Files that will be checked into git
git_files : elev_mask.cdl

.INTERMEDIATE : _elev_mask.nc

# Extract topg, thk and mass from a PISM state file.
_elev_mask.nc :
	$(PYTHON) write_elev_mask.py $(PISM_SPINUP_DIR)/g20km_10ka.nc _elev_mask.nc

# Convert netCDF to text-based format
elev_mask.cdl : _elev_mask.nc
	ncdump -l 10 _elev_mask.nc >elev_mask.cdl

# ================ Generate input files for the example, from text antecedents

# Reconstitute topg, thk and mass
elev_mask.nc : elev_mask.cdl
	echo '*****************************************************************'
	echo '[makefile] Reconstituting elevation and land mask data (from PISM)'
	ncgen -o elev_mask.nc elev_mask.cdl

# Generate grid: SeaRISE-Greenland-20km
sr_g20_pism.nc :
	echo '*****************************************************************'
	echo '[makefile] Generating sr_g20_pism.nc: SeaRISE Greenland 20km grid (PISM-style indexing)'
	searise

# Generate grid: ModelE - LonLat - Greenland - 2x2.5 degrees
modele_ll_g2x2_5.nc :
	echo '*****************************************************************'
	echo '[makefile] Generating modele_ll_g2x2_5.nc: ModelE Atmosphere Grid'
	modele_ll --zone g

# Overlap two grids
modele_ll_g2x2_5-sr_g20_pism.nc : modele_ll_g2x2_5.nc sr_g20_pism.nc
	echo '*****************************************************************'
	echo '[makefile] Overlapping ModelE and SeaRISE Grids'
	overlap modele_ll_g2x2_5.nc sr_g20_pism.nc

icebin_in.nc : modele_ll_g2x2_5-sr_g20_pism.nc elev_mask.nc
	echo '*****************************************************************'
	echo '[makefile] Assembling IceBin Input File from grids.'
	$(PYTHON) write_icebin_in_base.py

matrices.pik : icebin_in.nc
	echo '*****************************************************************'
	echo '[makefile] Writing regridding matrices to matrices.pik'
	$(PYTHON) write_matrices.py icebin_in.nc

test_diagonal.png test_elevation.png : icebin_in.nc
	echo '*****************************************************************'
	echo '[makefile] Testing that regridding matrices conserve; see test_diagonal.png and test_elevation.png'
	$(PYTHON) test_conserv.py

# =============== The actual example

# Grid outlines
sr_g20_pism.pdf : sr_g20_pism.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for SeaRISE Grid'
	$(PYTHON) plot_grid_outlines.py sr_g20_pism grid

modele_ll_g2x2_5.pdf : modele_ll_g2x2_5.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for ModelE Grid'
	$(PYTHON) plot_grid_outlines.py modele_ll_g2x2_5 grid

modele_ll_g2x2_5-sr_g20_pism.pdf : modele_ll_g2x2_5-sr_g20_pism.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for Exchange Grid'
	$(PYTHON) plot_grid_outlines.py modele_ll_g2x2_5-sr_g20_pism exgrid



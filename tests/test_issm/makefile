PYTHON = python3

all: issm3.nc sr_g20_pism.nc sr_g20_pism-ISSM_mesh.nc sr_g20_pism.pdf issm3.pdf sr_g20_pism-ISSM_mesh.pdf test_diagonal.png test_constant.png

clean:
	rm *.nc *.pdf *.png


# =========== Generate stuff we check into git
# This requires external software (PISM), and should not normally need
# to be run.

# Files that will be checked into git
git_files : issm3.cdl

.INTERMEDIATE : _issm3.nc

# Convert netCDF to text-based format
issm3.cdl : _issm3.nc
	ncdump -l 10 _issm3.nc >issm3.cdl


# ================ Generate input files for the example, from text antecedents

# Reconstitute topg, thk and mass
issm3.nc : issm3.cdl
	echo '*****************************************************************'
	echo '[makefile] Reconstituting elevation and land mask data (from PISM)'
	ncgen -o issm3.nc issm3.cdl


# Generate grid: SeaRISE-Greenland-20km
sr_g20_pism.nc :
	echo '*****************************************************************'
	echo '[makefile] Generating sr_g20_pism.nc: SeaRISE Greenland 20km grid (PISM-style indexing)'
	searise

sr_g20_pism-ISSM_mesh.nc : sr_g20_pism.nc issm3.nc
	echo '*****************************************************************'
	echo '[makefile] Overlapping the grids'
	overlap sr_g20_pism.nc issm3.nc 

# =============== The actual example

# Grid outlines
sr_g20_pism.pdf : sr_g20_pism.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for SeaRISE Grid'
	$(PYTHON) plot_grid_outlines.py sr_g20_pism grid

issm3.pdf : issm3.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for ISSM Grid'
	$(PYTHON) plot_grid_outlines.py issm3 grid

sr_g20_pism-ISSM_mesh.pdf : sr_g20_pism-ISSM_mesh.nc
	echo '*****************************************************************'
	echo '[makefile] Plotting grid outlines for Exchange Grid'
	$(PYTHON) plot_grid_outlines.py sr_g20_pism-ISSM_mesh exgrid

# The regridding test
test_diagonal.png test_constant.png : sr_g20_pism.nc sr_g20_pism-ISSM_mesh.nc
	$(PYTHON) test_regrid_issm.py

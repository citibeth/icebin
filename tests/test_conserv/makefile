PYTHON = python3

# Run the main test.
all:
	$(PYTHON) test_conserv.py

clean:
	rm *.nc *.pdf *.png

# =========== Generate stuff we check into git
# This requires external software (PISM), and should not normally need
# to be run.

# Files that will be checked into git
git_files : elev_mask.cdl

PISM_SPINUP_DIR = /Users/rpfische/exp/151014-integrate/build/std-greenland

.INTERMEDIATE : _elev_mask.nc

# Extract topg, thk and mass from a PISM state file.
_elev_mask.nc :
	$(PYTHON) write_elev_mask.py $(PISM_SPINUP_DIR)/g20km_10ka.nc _elev_mask.nc

# Convert netCDF to text-based format
elev_mask.cdl : _elev_mask.nc
	ncdump -l 10 _elev_mask.nc >elev_mask.cdl

# ================ To generate input files for the example, from text antecedents

# Reconstitute topg, thk and mass
elev_mask.nc : elev_mask.cdl
	ncgen -o elev_mask.nc elev_mask.cdl

# Generate grid: SeaRISE-Greenland-20km
sr_g20_pism.nc :
	searise

# Generate grid: ModelE - LonLat - Greenland - 2x2.5 degrees
modele_ll_g2x2_5.nc :
	modele_ll --zone g

# Overlap two grids
modele_ll_g2x2_5-sr_g20_pism.nc : modele_ll_g2x2_5.nc sr_g20_pism.nc
	$(GRIDS_BUILD)/overlap modele_ll_g2x2_5.nc sr_g20_pism.nc

icebin_in.nc : modele_ll_g2x2_5-sr_g20_pism.nc elev_mask.nc
	$(PYTHON) write_icebin_in_base.py

#matrices.pik : modele_ll_g2x2_5-sr_g20_pism-base.nc
#	$(PYTHON) write_matrices.py modele_ll_g2x2_5-sr_g20_pism-base.nc

# =============== The actual example

# Grid outlines
sr_g20_pism.png : sr_g20_pism.nc
	$(PYTHON) plot_grid_outlines.py sr_g20_pism grid

modele_ll_g2x2_5.png : modele_ll_g2x2_5.nc
	$(PYTHON) plot_grid_outlines.py modele_ll_g2x2_5 grid

modele_ll_g2x2_5-sr_g20_pism.png : modele_ll_g2x2_5-sr_g20_pism.nc
	$(PYTHON) plot_grid_outlines.py modele_ll_g2x2_5-sr_g20_pism exgrid



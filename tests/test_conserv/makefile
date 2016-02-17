GRIDS_BUILD = ../../build/grids
SLIB_BUILD = ../../build/slib
PYLIB_BUILD = ../../build/pylib
ENV = $(PYLIB_BUILD)/env
PYTHON = $(ENV) python

# Run the main test.
all:
	$(PYTHON) test_conserv.py

clean:
	rm *.nc *.pdf *.png

# =========== To generate stuff we check into git
git_files : elev_mask.cdl

PISM_SPINUP_DIR = /Users/rpfische/exp/151014-integrate/build/std-greenland

.INTERMEDIATE : _elev_mask.nc

_elev_mask.nc :
	python write_elev_mask.py $(PISM_SPINUP_DIR)/g20km_10ka.nc _elev_mask.nc

elev_mask.cdl : _elev_mask.nc
	ncdump -l 10 _elev_mask.nc >elev_mask.cdl

# ================ To generate input files for the example, from text antecedents
elev_mask.nc : elev_mask.cdl
	ncgen -o elev_mask.nc elev_mask.cdl

sr_g20_pism.nc :
	$(GRIDS_BUILD)/searise

modele_ll_g2x2_5.nc :
	$(GRIDS_BUILD)/modele_ll --zone g

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



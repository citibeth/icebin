from icebin import ibgrid, element_l1

ISSM_IN = 'modele_ll_g2x2_5-ISSM_mesh.nc'

with netCDF4.Dataset(ISSM_IN) as nc:
    gridA = ibgrid.read_nc(nc, 'gridA')    # ISMIP6 Grid
    gridI = ibgrid.read_nc(nc, 'gridI')    # ISSM Grid
    exgrid = ibgrid.read_nc(nc, 'exgrid')    # Exchange Grid

AvI,weightsA,weightsI = element_l1.compute_AvI(gridX, gridI)

valueI = ...
valueA = (1. / weightsA) * AvI * valueI

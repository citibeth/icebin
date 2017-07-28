

# Generate grid: SeaRISE-Greenland-20km
sr_g20_searise.nc :
	echo '*****************************************************************'
	echo '[makefile] Generating sr_g20_searise.nc: SeaRISE Greenland 20km grid (PISM-style indexing)'
	searise --icemodel searise

# Generate the 2 minute grid
modele_ll_g2mx2m.nc:
	modele_ll --zone g --grid 2mx2m --pole-caps=0

# Overlap two grids
modele_ll_g2mx2m-sr_g20_searise.nc: modele_ll_g2mx2m.nc sr_g20_searise.nc
	overlap modele_ll_g2mx2m.nc sr_g20_searise.nc


# ------------------------------------------------------------------

# Convert TOPO-generation input files to NetCDF, for easy viewing
Z2MX2M_NGDC = $(shell ./searchenv Z2MX2M.NGDC MODELE_FILE_PATH)
Z2MX2M.NGDC.nc: $(Z2MX2M_NGDC)
	giss2nc $(Z2MX2M_NGDC) Z2MX2M.NGDC.nc --endian=big --names=FOCEN2,ZETOP2


Z10MX10M = $(shell ./searchenv Z10MX10M MODELE_FILE_PATH)
Z10MX10M.nc: $(Z10MX10M)
	giss2nc $(Z10MX10M) Z10MX10M.nc --endian=big --names=_,FLAKES

ZICEHXH = $(shell ./searchenv ZICEHXH MODELE_FILE_PATH)
ZICEHXH.nc: $(ZICEHXH)
	giss2nc $(ZICEHXH) ZICEHXH.nc --endian=big --names=dZGICH,FGICEH,ZSOLDH

ZNGDC1 = $(shell ./searchenv ZNGDC1 MODELE_FILE_PATH)
ZNGDC1.nc: $(ZNGDC1)
	giss2nc $(ZNGDC1) ZNGDC1.nc --endian=big --names=_,_,_,FCONT1,_,FGICE1



#giss2nc $DIR/Z10MX10M Z10MX10M.nc --endian=big \
#	--names=_,FLAKES
#
#giss2nc $DIR/ZICEHXH ZICEHXH.nc --endian=big \
#	--names=dZGICH,FGICEH,ZSOLDH
#
#giss2nc $DIR/ZNGDC1 ZNGDC1.nc --endian=big \
#	--names=_,_,_,FCONT1,FGICE1

# Make sure you've done configure and make in the above directory
# before trying to run the tests.


all : greenland_2x2_5-searise.nc

data :
	./fetch_data

greenland_2x2_5-searise.nc : greenland_2x2_5.nc searise.nc
	../sbin/overlap greenland_2x2_5.nc searise.nc

greenland_2x2_5.nc :
	../sbin/greenland_2x2_5

searise.nc :
	../sbin/searise



import giss.proj
import pyproj

sproj = "+proj=stere +lon_0=-39 +lat_0=90 +lat_ts=71.0 +ellps=WGS84"

xyproj = pyproj.Proj(sproj)
llproj = xyproj.to_latlong()

x=15.0
y=15.0

lon, lat = pyproj.transform(xyproj, llproj, x, y)
#lon, lat = pyproj.transform(llproj, xyproj, x, y)

print lon
print lat

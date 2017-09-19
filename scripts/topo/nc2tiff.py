import numpy as np
import PIL
import netCDF4

from PIL import Image

#fname = 'ZNGDC1.nc'
#vname = 'FCONT1'

fname = 'Z2MX2M.NGDC.nc'
vname = 'FOCEN2'

with netCDF4.Dataset(fname) as nc:
    val_d = nc.variables[vname][:]

val_i = np.zeros(val_d.shape, dtype='i')
val_i[:] = val_d[:]
#val_i = val_i.transpose()

#val_i[:] = 1
#val_i[0:20,:] = 0


# http://www.effbot.org/imagingbook/image.htm
#img = PIL.Image.fromarray(val_i, mode='1')

img = PIL.Image.new('P', val_i.shape)
for i in range(0,val_i.shape[0]):
    print('i=%d of %d' % (i, val_i.shape[0]))
    for j in range(0,val_i.shape[1]):
        img.putpixel((i,j), (val_i[i,j],) )

img = img.transpose(PIL.Image.ROTATE_90)

# Create initial grayscale palette
palette = np.zeros((256,3), dtype='i')
for i in range(0,256):
    palette[i,:] = 255-i

# Replace certain colors
palette[0,:] = (0,0,255)        # FCONT=0 = Ocean
palette[1,:] = (0,255,0)        # 1 = Land
palette[2,:] = (255,255,255)    # 2 = Greenland


pl = list(palette.reshape(-1))
print(pl)
img.putpalette(pl)





img.save('img.tiff')

print(img.getbands())

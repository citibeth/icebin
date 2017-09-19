import numpy as np
import PIL
import netCDF4

from PIL import Image

img = PIL.Image.open('img2.tiff')
img = img.transpose(PIL.Image.ROTATE_270)


val_d = np.zeros(img.size)
for i in range(0,val_d.shape[0]):
    print('%d of %d' % (i, val_d.shape[0]))
    for j in range(0,val_d.shape[1]):
        val_d[i,j] = img.getpixel((i,j))


print(val_d.shape)

with netCDF4.Dataset('x.nc', 'a') as nc:
    nc.variables['FOCEN2'][:] = val_d[:]

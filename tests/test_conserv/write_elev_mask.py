import netCDF4
import sys
import giss.ncutil

# Equivalent to ncks -v

ifname = sys.argv[1]
ofname = sys.argv[2]

fin = netCDF4.Dataset(ifname)
fout = netCDF4.Dataset(ofname, 'w')

var_filter = lambda vname: vname if vname in {'topg', 'thk', 'mask'} else None

copy = giss.ncutil.copy_nc(fin, fout, var_filter=var_filter, attrib_filter=lambda x: True)
copy.define_vars()
copy.copy_data()

fout.close()
fin.close()

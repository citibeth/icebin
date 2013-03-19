import glint2

mm2 = glint2.MatrixMaker()
print '***** Reading new MatrixMaker'
mm2.load('mm.nc', 'm')

nc = glint2.NcFile('mm2.nc', 'w')
mm2.write(nc, 'm')
nc.close()

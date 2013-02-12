import glint2

mm2 = glint2.MatrixMaker()
print '***** Reading new MatrixMaker'
mm2.load('mm.nc', 'm')
mm2.save('mm2.nc', 'm')

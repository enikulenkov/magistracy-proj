#!/usr/bin/python2

import numpy as np
from fortran_file import FortranFile

f = FortranFile('FCC.dat')
reals = f.readReals(prec='d')
matr = reals.reshape((6369/3, 3))
np.savetxt('out.mat', matr)


#>>> f = FortranFile('tetrahedron.dat','@','i','w')
#>>> f.writeReals([1,-1,-1,1,1,-1,1,-1,1,1,-1,-1],'d')
#>>> f.close()

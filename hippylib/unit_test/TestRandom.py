import dolfin as dl
import numpy as np
import sys
sys.path.append( "../../" )
from hippylib import *

nx = 10
ny = 10
mesh = dl.UnitSquareMesh(nx, ny)
    
rank = dl.MPI.rank(mesh.mpi_comm())
nproc = dl.MPI.size(mesh.mpi_comm())

#Random.split(rank, nproc, 10, 1)

#np.random.seed(1)
pp=np.random.RandomState(1)
z = np.zeros(20)
for i in np.arange(z.shape[0]):
    z[i] = Random.normal(0., 1.)
    #z[i] = np.random.randn(1)
    #z[i] = pp.standard_normal(1)
    
print "PID {0}:".format(rank), z, z.mean(), z.std()
    

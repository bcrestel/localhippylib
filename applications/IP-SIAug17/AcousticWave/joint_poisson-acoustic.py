import dolfin as dl

from exple_poisson import model_poisson
from exple_acousticinversion import model_acoustic

from fenicstools.mpicomm import create_communicators

from hippylib import ZeroPrior,\
STATE, ADJOINT, PARAMETER

dl.set_log_active(False)





mpicomm_local, mpicomm_global = create_communicators()
mpiworldrank = dl.MPI.rank(dl.mpi_comm_world())
#PRINT = (mpiworldrank == 0)
PRINT = True

nx, ny = 20, 20
mesh = dl.UnitSquareMesh(mpicomm_local, nx, ny)
Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
Vh = [None, None, None]
Vh[PARAMETER] = Vh1
Vh[STATE] = Vh2
Vh[ADJOINT] = Vh2

modelpoisson = model_poisson(Vh, ZeroPrior(Vh[PARAMETER]), PRINT)
modelacoustic = model_acoustic(mpicomm_local, mpicomm_global, Vh,\
ZeroPrior(Vh[PARAMETER]), PRINT)

#TODO: continue here

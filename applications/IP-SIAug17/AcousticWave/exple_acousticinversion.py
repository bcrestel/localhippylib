"""
Solve acoustic wave inverse problem in terms of parameter a only
"""
import dolfin as dl

from hippylib import *

from fenicstools.acousticwave import AcousticWave
from fenicstools.sourceterms import PointSources, RickerWavelet
from fenicstools.observationoperator import TimeObsPtwise
from fenicstools.regularization import TVPD, TV
from fenicstools.mpicomm import create_communicators, partition_work
from fenicstools.examples.acousticwave.mediumparameters import \
targetmediumparameters, initmediumparameters, loadparameters

# Create local and global communicators
mpicomm_local, mpicomm_global = create_communicators()
mpiworldrank = dl.MPI.rank(dl.mpi_comm_world())
PRINT = (mpiworldrank == 0)
mpicommbarrier = dl.mpi_comm_world()

Nxy, Dt, fpeak, t0, t1, t2, tf = loadparameters(False)
h = 1./Nxy
if PRINT:
    print 'Nxy={} (h={}), Dt={}, fpeak={}, t0,t1,t2,tf={}'.format(\
    Nxy, h, Dt, fpeak, [t0,t1,t2,tf])

X, Y = 1.0, 1.0
mesh = dl.UnitSquareMesh(mpicomm_local, Nxy, Nxy)
Vl = dl.FunctionSpace(mesh, 'Lagrange', 1)

# Source term:
Ricker = RickerWavelet(fpeak, 1e-6)
r = 2   # polynomial degree for state and adj
V = dl.FunctionSpace(mesh, 'Lagrange', r)
y_src = 0.1 # 1.0->reflection, 0.1->transmission
#Pt = PointSources(V, [[0.1*ii*X-0.05, y_src] for ii in range(1,11)])
Pt = PointSources(V, [[0.1,y_src], [0.5,y_src], [0.9,y_src]])
#Pt = PointSources(V, [[0.5, y_src]])
srcv = dl.Function(V).vector()

# Boundary conditions:
class ABCdom(dl.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[1] < Y)

Wave = AcousticWave({'V':V, 'Vm':Vl}, 
{'print':False, 'lumpM':True, 'timestepper':'backward'})
Wave.set_abc(mesh, ABCdom(), lumpD=False)


at, bt,_,_,_ = targetmediumparameters(Vl, X)
a0, _,_,_,_ = initmediumparameters(Vl, X)
b0 = bt
Wave.update({'b':bt, 'a':at, 't0':0.0, 'tf':tf, 'Dt':Dt,\
'u0init':dl.Function(V), 'utinit':dl.Function(V)})
if PRINT:
    print 'nb of src={}, nb of timesteps={}'.format(len(Pt.src_loc), Wave.Nt)

sources, timesteps = partition_work(mpicomm_local, mpicomm_global, \
len(Pt.src_loc), Wave.Nt)

mpilocalrank = dl.MPI.rank(mpicomm_local)
mpiglobalrank = dl.MPI.rank(mpicomm_global)
mpiworldsize = dl.MPI.size(dl.mpi_comm_world())
print 'mpiworldrank={}, mpiglobalrank={}, mpilocalrank={}, sources={}, timestep=[{},{}]'.format(\
mpiworldrank, mpiglobalrank, mpilocalrank, sources,\
timesteps[0], timesteps[-1])

obspts = [[ii*float(X)/float(Nxy), Y] for ii in range(1,Nxy)]
#obspts = [[0.0, ii/10.] for ii in range(1,10)] + \
#[[1.0, ii/10.] for ii in range(1,10)] + \
#[[ii/10., 0.0] for ii in range(1,10)] + \
#[[ii/10., 1.0] for ii in range(1,10)]
tfilterpts = [t0, t1, t2, tf]
obsop = TimeObsPtwise({'V':V, 'Points':obspts}, tfilterpts)

reg = TVPD({'Vm':Vl, 'eps':1.0, 'k':1e-6, 'print':PRINT})

model = ModelAcoustic(mpicomm_global, Wave, [Ricker, Pt, srcv], sources,
timesteps, obsop, at, reg)

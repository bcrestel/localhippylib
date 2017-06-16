"""
Solve acoustic wave inverse problem in terms of parameter a only
"""
import dolfin as dl

from hippylib.model_acousticinversiona import ModelAcoustic
from hippylib import ReducedSpaceNewtonCG
from hippylib import STATE, PARAMETER, ADJOINT

from fenicstools.acousticwave import AcousticWave
from fenicstools.sourceterms import PointSources, RickerWavelet
from fenicstools.observationoperator import TimeObsPtwise
from fenicstools.regularization import TVPD, TV
from fenicstools.mpicomm import create_communicators, partition_work
from fenicstools.examples.acousticwave.mediumparameters import \
targetmediumparameters, initmediumparameters, loadparameters

dl.set_log_active(False)


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
tfilterpts = [t0, t1, t2, tf]
obsop = TimeObsPtwise({'V':V, 'Points':obspts}, tfilterpts)

reg = TVPD({'Vm':Vl, 'eps':1.0, 'k':1e-6, 'print':PRINT})

if PRINT:   print 'Create model'
model = ModelAcoustic(mpicomm_global, Wave, [Ricker, Pt, srcv], sources,
timesteps, obsop, reg, at, bt)

if PRINT:   print 'Generate synthetic data'
model.generate_synthetic_obs(20.0)

out = model.generate_vector(PARAMETER)
x = model.generate_vector("ALL")
x[PARAMETER] = a0
model.solveFwd(out, x)
_, costreg0, costmisf0 = model.cost(x)
x[PARAMETER] = at
model.solveFwd(out, x)
_, costregt, costmisft = model.cost(x)
if PRINT:   
    print 'misfit at target={:.4e}, at initial state={:.4e}'.format(\
    costmisft, costmisf0)
    print 'Regularization at target={:.2e}, at initial state={:.2e}'.format(\
    costregt, costreg0)

if PRINT:   print 'Solve inverse problem'
solver = ReducedSpaceNewtonCG(model)
solver.parameters["rel_tolerance"] = 1e-10
solver.parameters["abs_tolerance"] = 1e-12
solver.parameters["inner_rel_tolerance"] = 1e-15
solver.parameters["gda_tolerance"] = 1e-24
solver.parameters["c_armijo"] = 5e-5
solver.parameters["max_backtracking_iter"] = 20
solver.parameters["GN_iter"] = 20
solver.parameters["max_iter"] = 500
solver.parameters["print_level"] = 0
if not PRINT:   solver.parameters["print_level"] = -1

x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[1e-4, 1.0])

minat = at.vector().min()
maxat = at.vector().max()
minbt = bt.vector().min()
maxbt = bt.vector().max()
mina0 = a0.vector().min()
maxa0 = a0.vector().max()
minb0 = b0.vector().min()
maxb0 = b0.vector().max()
a = model.objacoustic.PDE.a
b = model.objacoustic.PDE.b
mina = a.vector().min()
maxa = a.vector().max()
minb = b.vector().min()
maxb = b.vector().max()
if PRINT:
    print '\ntarget: min(a)={}, max(a)={}'.format(minat, maxat)
    print 'init: min(a)={}, max(a)={}'.format(mina0, maxa0)
    print 'MAP: min(a)={}, max(a)={}'.format(mina, maxa)

    print '\ntarget: min(b)={}, max(b)={}'.format(minbt, maxbt)
    print 'init: min(b)={}, max(b)={}'.format(minb0, maxb0)
    print 'MAP: min(b)={}, max(b)={}'.format(minb, maxb)

mdmis, mdmisperc = model.mediummisfit(x[PARAMETER])
if PRINT:
    print 'medmisft={:e} ({:.1f}%)'.format(mdmis, mdmisperc)
    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost

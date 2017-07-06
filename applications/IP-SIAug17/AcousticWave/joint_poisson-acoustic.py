import sys
import dolfin as dl

from exple_poisson import model_poisson
from exple_acousticinversion import model_acoustic
from targetmedium import targetmediumparameters, initmediumparameters

from fenicstools.mpicomm import create_communicators
from fenicstools.regularization import TVPD
from fenicstools.jointregularization import \
SumRegularization, V_TVPD, NuclearNormSVD2D
from fenicstools.plotfenics import PlotFenics

from hippylib import ZeroPrior, ReducedSpaceNewtonCG, BFGS,\
STATE, ADJOINT, PARAMETER
from hippylib.jointmodel import JointModel

dl.set_log_active(False)


# Command-line argument
try:
    k = float(sys.argv[1])
except:
    k = 2e-8
#######################


mpicomm_local, mpicomm_global = create_communicators()
mpiworldrank = dl.MPI.rank(dl.mpi_comm_world())
PRINT = (mpiworldrank == 0)

nx, ny = 20, 20
if PRINT:   'nx={}, ny={}'.format(nx, ny)
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

# regularization
eps = 1e-3
#regpoisson = TVPD({'Vm':Vh[PARAMETER], 'k':2e-8, 'eps':eps, 'print':PRINT})
#regacoustic = TVPD({'Vm':Vh[PARAMETER], 'k':2e-8, 'eps':eps, 'print':PRINT})
#jointregul = SumRegularization(regpoisson, regacoustic, \
#coeff_cg=k,\
#coeff_ncg=0.0, parameters_ncg={'eps':1e-5},\
#coeff_vtv=0.0, isprint=PRINT)
#prefix = 'SumRegul'

#jointregul = V_TVPD(Vh[PARAMETER], {'k':k, 'eps':eps,\
#'rescaledradiusdual':1.0, 'print':PRINT})
#prefix = 'VTV'

jointregul = NuclearNormSVD2D(mesh, {'eps':eps, 'k':k}, isprint=PRINT)
prefix = 'NN'


jointmodel = JointModel(modelpoisson, modelacoustic, jointregul,\
parameters={'print':PRINT, 'splitassign':True})

if PRINT:   print '\nSolve joint inverse problem'
if prefix == 'NN':
    solver = BFGS(jointmodel)
    solver.parameters["H0inv"] = 'd0'
else:
    solver = ReducedSpaceNewtonCG(jointmodel)
solver.mpicomm_global = mpicomm_global
solver.parameters["rel_tolerance"] = 1e-10
solver.parameters["abs_tolerance"] = 1e-12
solver.parameters["inner_rel_tolerance"] = 1e-15
solver.parameters["gda_tolerance"] = 1e-24
solver.parameters["c_armijo"] = 5e-5
solver.parameters["max_backtracking_iter"] = 20
solver.parameters["GN_iter"] = 20
solver.parameters["check_param"] = 10
solver.parameters["max_iter"] = 5000
solver.parameters["print_level"] = 0
if not PRINT:   solver.parameters["print_level"] = -1

a0poisson,_ = initmediumparameters(Vh[PARAMETER], 0.0)
a0acoustic,_ = initmediumparameters(Vh[PARAMETER], 1.0)
a0 = dl.interpolate(dl.Constant(("0.0","0.0")), jointmodel.Vh[PARAMETER])
dl.assign(a0.sub(0), a0poisson)
dl.assign(a0.sub(1), a0acoustic)
if prefix == 'NN':
    x = solver.solve(a0.vector(), bounds_xPARAM=[1e-8, 10.0])
else:
    x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[1e-8, 10.0])

xP, xAW = jointmodel.splitvector(x)
minxP = xP[PARAMETER].min()
maxxP = xP[PARAMETER].max()
mmfP, mmfPperc = modelpoisson.mediummisfit(xP[PARAMETER])
minxAW = xAW[PARAMETER].min()
maxxAW = xAW[PARAMETER].max()
mmfAW, mmfAWperc = modelacoustic.mediummisfit(xAW[PARAMETER])
atP, btP = targetmediumparameters(Vh[PARAMETER], 0.0)
minatP = atP.vector().min()
maxatP = atP.vector().max()
atA, btA = targetmediumparameters(Vh[PARAMETER], 1.0)
minatA = atA.vector().min()
maxatA = atA.vector().max()
if PRINT:
    print '\ntarget: min(a)={}, max(a)={}'.format(minatP, maxatP)
    print 'MAP-poisson: min(a)={}, max(a)={}'.format(minxP, maxxP)
    print 'medmisfit={:e} ({:.1f}%)'.format(mmfP, mmfPperc)
    print '\ntarget: min(a)={}, max(a)={}'.format(minatA, maxatA)
    print 'MAP-acoustic: min(a)={}, max(a)={}'.format(minxAW, maxxAW)
    print 'medmisfit={:e} ({:.1f}%)'.format(mmfAW, mmfAWperc)

    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost

    myplot = PlotFenics(comm = mesh.mpi_comm(),\
    Outputfolder='joint_poisson-acoustic/plots')
    myplot.set_varname(prefix+'-poisson-MAP_k' + str(k) + '_e' + str(eps))
    myplot.plot_vtk(vector2Function(xP[PARAMETER], Vh[PARAMETER]))
    myplot.set_varname(prefix+'-acoustic-MAP_k' + str(k) + '_e' + str(eps))
    myplot.plot_vtk(vector2Function(xAW[PARAMETER], Vh[PARAMETER]))

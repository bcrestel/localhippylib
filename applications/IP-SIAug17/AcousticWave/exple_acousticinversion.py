"""
Solve acoustic wave inverse problem in terms of parameter a only
"""
import sys
import dolfin as dl

from hippylib.model_acousticinversiona import ModelAcoustic
from hippylib import ReducedSpaceNewtonCG
from hippylib import STATE, PARAMETER, ADJOINT

from fenicstools.acousticwave import AcousticWave
from fenicstools.sourceterms import PointSources, RickerWavelet
from fenicstools.observationoperator import TimeObsPtwise
from fenicstools.regularization import TVPD
from fenicstools.plotfenics import PlotFenics
from fenicstools.mpicomm import create_communicators, partition_work
from targetmedium import targetmediumparameters, initmediumparameters, loadparameters

dl.set_log_active(False)

LARGELOADPARAMETERS=True


def model_acoustic(mpicomm_local, mpicomm_global, Vh, reg, PRINT=False):
    X, Y = 1.0, 1.0
    V = Vh[STATE]
    Vl = Vh[PARAMETER]

    # source locations:
    y_src = 0.1 # 1.0->reflection, 0.1->transmission
    #Pt = PointSources(V, [[0.1,y_src], [0.25,y_src], [0.4,y_src],\
    #[0.6,y_src], [0.75,y_src], [0.9,y_src]])
    #Pt = PointSources(V, [[0.2,y_src], [0.5,y_src], [0.8,y_src]])
    Pt = PointSources(V, [[0.5,y_src]])

    # Absorbing Boundary Conditions on left, bott, & right:
    class ABCdom(dl.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (x[1] < Y)

    Wave = AcousticWave({'V':V, 'Vm':Vl}, 
    {'print':False, 'lumpM':True, 'timestepper':'backward'})
    Wave.set_abc(Vl.mesh(), ABCdom(), lumpD=False)

    _, Dt, fpeak, t0, t1, t2, tf = loadparameters(LARGELOADPARAMETERS)
    at, bt = targetmediumparameters(Vl, 1.0)
    Wave.update({'b':bt, 'a':at, 't0':t0, 'tf':tf, 'Dt':Dt,\
    'u0init':dl.Function(V), 'utinit':dl.Function(V)})
    if PRINT:
        print '[acoustic] nb of src={}, nb of timesteps={}'.format(len(Pt.src_loc), Wave.Nt)

    sources, timesteps = partition_work(mpicomm_local, mpicomm_global, \
    len(Pt.src_loc), Wave.Nt)

    mpilocalrank = dl.MPI.rank(mpicomm_local)
    mpiglobalrank = dl.MPI.rank(mpicomm_global)
    mpiworldrank = dl.MPI.rank(dl.mpi_comm_world())
    print '[acoustic] mpiworldrank={}, mpiglobalrank={}, mpilocalrank={}, sources={}, timestep=[{},{}]'.format(\
    mpiworldrank, mpiglobalrank, mpilocalrank, sources,\
    timesteps[0], timesteps[-1])

    obspts = [[ii*float(X)/20.0, Y] for ii in range(1,20)]
    tfilterpts = [t0, t1, t2, tf]
    obsop = TimeObsPtwise({'V':V, 'Points':obspts}, tfilterpts)

    # Source term:
    Ricker = RickerWavelet(fpeak, 1e-6)
    srcv = dl.Function(V).vector()

    model = ModelAcoustic(mpicomm_global, Wave, [Ricker, Pt, srcv], sources,
    timesteps, obsop, reg, at, bt)

    model.generate_synthetic_obs(20.0)

    out = model.generate_vector(PARAMETER)
    x = model.generate_vector("ALL")
    a0,_ = initmediumparameters(Vl, 1.0)
    x[PARAMETER] = a0
    model.solveFwd(out, x)
    _, costreg0, costmisf0 = model.cost(x)
    x[PARAMETER] = at
    model.solveFwd(out, x)
    _, costregt, costmisft = model.cost(x)
    if PRINT:   
        print '[acoustic] misfit at target={:.4e}, at initial state={:.4e}'.format(\
        costmisft, costmisf0)
        print '[acoustic] Regularization at target={:.2e}, at initial state={:.2e}'.format(\
        costregt, costreg0)

    return model




if __name__ == "__main__":
    # Command-line argument
    try:
        k = float(sys.argv[1])
        eps = float(sys.argv[2])
    except:
        k = 2e-7
        eps = 1e-3
    #######################

    PLOT = True

    # Create local and global communicators
    mpicomm_local, mpicomm_global = create_communicators()
    mpiworldrank = dl.MPI.rank(dl.mpi_comm_world())
    PRINT = (mpiworldrank == 0)
    mpicommbarrier = dl.mpi_comm_world()

    Nxy, Dt, fpeak, t0, t1, t2, tf = loadparameters(LARGELOADPARAMETERS)
    h = 1./Nxy
    if PRINT:
        print 'Nxy={} (h={}), Dt={}, fpeak={}, t0,t1,t2,tf={}'.format(\
        Nxy, h, Dt, fpeak, [t0,t1,t2,tf])

    mesh = dl.UnitSquareMesh(mpicomm_local, Nxy, Nxy)
    Vl = dl.FunctionSpace(mesh, 'Lagrange', 1)
    r = 2   # polynomial degree for state and adj
    V = dl.FunctionSpace(mesh, 'Lagrange', r)
    Vh = [None, None, None]
    Vh[PARAMETER] = Vl
    Vh[STATE] = V
    Vh[ADJOINT] = V

    reg = TVPD({'Vm':Vl, 'eps':eps, 'k':k, 'print':PRINT})

    model = model_acoustic(mpicomm_local, mpicomm_global, Vh, reg, PRINT)

    if PRINT:   print 'Solve inverse problem'
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-10
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 20
    solver.parameters["GN_iter"] = 20
    solver.parameters["max_iter"] = 5000
    solver.parameters["print_level"] = 0
    if not PRINT:   solver.parameters["print_level"] = -1

    a0, b0 = initmediumparameters(Vl, 1.0)
    x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[1e-4, 1.0])

    at = model.atrue
    bt = model.btrue
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

        if PLOT:
            myplot = PlotFenics(comm = mesh.mpi_comm(),\
            Outputfolder='exple_acousticinversion/plots')
            if LARGELOADPARAMETERS: suff = '4Hz'
            else:   suff = '2Hz'
            model.objacoustic._plotab(myplot, \
            'acoustic' + suff + '-MAP_k' + str(k) + '_e' + str(eps))

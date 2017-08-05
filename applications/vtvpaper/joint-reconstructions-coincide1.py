import sys, os
import numpy as np
import dolfin as dl
from dolfin import Expression

from hippylib import ZeroPrior, ReducedSpaceNewtonCG, amg_method, Random, Model,\
STATE, ADJOINT, PARAMETER, ReducedHessian, CGSolverSteihaug, vector2Function
from hippylib.jointmodel import JointModel
from hippylib.bfgs import H0invdefault

from fenicstools.jointregularization import V_TVPD, V_TV
from fenicstools.plotfenics import PlotFenics

from targetmedium_coincide1 import targetmedium
from definePDE import pdes
from definemisfit import defmisfit

PLOT = False
TESTPC = False
EIG = False

if __name__ == "__main__":
    dl.set_log_active(False)
    nx = 64
    ny = 64

    ######################################
    try:
        k = float(sys.argv[1])
        eps = float(sys.argv[2])
    except:
        k = 3e-7
        eps = 1e-3
    ######################################
    
    mesh = dl.UnitSquareMesh(nx, ny)
    
    mpicomm = mesh.mpi_comm()
    mpirank = dl.MPI.rank(mpicomm)
    mpisize = dl.MPI.size(mpicomm)

    PRINT = (mpirank == 0)
    
    if mpisize > 1:
        Random.split(mpirank, mpisize, 1000000, 1)
    
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    # Target medium parameters
    a1true, a2true = targetmedium(Vh, PARAMETER)
    if PLOT:
        PltFen = PlotFenics(mpicomm, os.path.splitext(sys.argv[0])[0] + '/plots')
        PltFen.set_varname('target1')
        PltFen.plot_vtk(a1true)
        PltFen.set_varname('target2')
        PltFen.plot_vtk(a2true)

    # Define PDE
    amg_pde = amg_method()
    pde1, pde2 = pdes(Vh, STATE, amg_pde)
 
    # Define misfit functions
    misfit1, misfit2, targets1, targets2 = defmisfit(Vh, STATE)

    # Generate synthetic observations
    rel_noise_level = 0.02
    utrue1 = pde1.generate_state()
    utrue2 = pde2.generate_state()
    for misfit, atrue, targets, utrue, pde in zip([misfit1, misfit2], \
    [a1true, a2true], [targets1, targets2], [utrue1, utrue2], [pde1, pde2]):
        x = [utrue, atrue.vector(), None]
        minatrue = dl.MPI.min(mesh.mpi_comm(), np.amin(atrue.vector().array()))
        maxatrue = dl.MPI.max(mesh.mpi_comm(), np.amax(atrue.vector().array()))
        if PRINT:
            print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
        pde.solveFwd(x[STATE], x, 1e-9)
        noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
        Random.normal(x[STATE], noise_level, False)
        misfit.B.mult(x[STATE], misfit.d)
        misfit.noise_variance = np.sqrt(targets.shape[0])   # hack to compare both models

    # Define models
    model1 = Model(pde1, ZeroPrior(Vh[PARAMETER]), misfit1, a1true.vector())
    model2 = Model(pde2, ZeroPrior(Vh[PARAMETER]), misfit2, a2true.vector())
    x[STATE].zero()
    c1, r1, m1 = model1.cost(x)
    c2, r2, m2 = model2.cost(x)
    if PRINT:
        print 'Cost @ MAP for m1: cost={}, misfit={}, reg={}'.format(c1, m1, r1)
        print 'Cost @ MAP for m2: cost={}, misfit={}, reg={}'.format(c2, m2, r2)

    ############ Regularization #############
    #jointregul = V_TV(Vh[PARAMETER], {'k':k, 'eps':eps, 'amg':'petsc_amg',\
    #'print':PRINT, 'GNhessian':True})
    jointregul = V_TVPD(Vh[PARAMETER], {'k':k, 'eps':eps, 'amg':'petsc_amg',\
    'rescaledradiusdual':1.0, 'print':PRINT, 'PCGN':False})
    filename = 'coincide1-k' + str(k) + '_e' + str(eps)
    filename += '-VTVPD_petsc_amg'
    #########################################

    jointmodel = JointModel(model1, model2, jointregul,
    parameters={'print':PRINT, 'splitassign':True})

    SOLVER = 'Newton'
    #########################################
    if SOLVER == 'Newton':
        solver = ReducedSpaceNewtonCG(jointmodel)
        solver.parameters['PC'] = 'BFGS'
        solver.parameters['memory_limit'] = 1000
        solver.parameters['H0inv'] = 'Rinv'
        filename += '-Newton_PCBFGSRinv'
    elif SOLVER == 'BFGS':
        solver = BFGS(jointmodel)
        solver.parameters['H0inv'] = 'default'
        solver.parameters['memory_limit'] = 5000
        filename += '-BFGS'
    #########################################

    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 25 
    solver.parameters["cg_coarse_tolerance"] = 0.5
    solver.parameters["GN_iter"] = 10
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if not PRINT:
        solver.parameters["print_level"] = -1

    pde1.PDEcounts = 0
    pde2.PDEcounts = 0
    
    a0 = dl.interpolate(dl.Constant(("0.0","0.0")),jointmodel.Vh[PARAMETER])
    if SOLVER == 'Newton':
        x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[-10., 20.])
    elif SOLVER == 'BFGS':
        x = solver.solve(a0.vector(), bounds_xPARAM=[-10., 20.])

    x1, x2 = jointmodel.splitvector(x)
    minaf1 = x1[PARAMETER].min()
    maxaf1 = x1[PARAMETER].max()
    md1mis, md1misperc = model1.mediummisfit(x1[PARAMETER])
    minaf2 = x2[PARAMETER].min()
    maxaf2 = x2[PARAMETER].max()
    md2mis, md2misperc = model2.mediummisfit(x2[PARAMETER])
    if PRINT:
        print 'min(af1)={}, max(af1)={}, med1misft={:e} ({:.1f}%)'.format(\
        minaf1, maxaf1, md1mis, md1misperc)
        print 'min(af2)={}, max(af2)={}, med2misft={:e} ({:.1f}%)'.format(\
        minaf2, maxaf2, md2mis, md2misperc)
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"

        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost

    
    if PLOT:
        PltFen.set_varname('jointmodel1-' + filename)
        PltFen.plot_vtk(vector2Function(x1[PARAMETER], Vh[PARAMETER]))
        PltFen.set_varname('jointmodel2-' + filename)
        PltFen.plot_vtk(vector2Function(x2[PARAMETER], Vh[PARAMETER]))


    if TESTPC:
        if PRINT:
            print '\n\nCompare different preconditioners at MAP point'

        dl.MPI.barrier(dl.mpi_comm_world())
        # BFGS-Rinv
        if PRINT:
            print '\nPreconditioner: BFGS_H0=Rinv'
        ahat = jointmodel.generate_vector(PARAMETER)
        mg = jointmodel.generate_vector(PARAMETER)
        jointmodel.solveFwd(x[STATE], x, solver.parameters["inner_rel_tolerance"])
        jointmodel.solveAdj(x[ADJOINT], x, solver.parameters["inner_rel_tolerance"])
        jointmodel.setPointForHessianEvaluations(x)
        gradnorm = jointmodel.evalGradientParameter(x, mg)
        bfgsPC = solver.bfgsPC
        bfgsPC.set_H0inv(jointmodel.Prior.getprecond())
        H = ReducedHessian(jointmodel, solver.parameters["inner_rel_tolerance"], 
        gauss_newton_approx=False, misfit_only=False)
        solvercg = CGSolverSteihaug()
        solvercg.set_operator(H)
        solvercg.set_preconditioner(bfgsPC)
        solvercg.parameters["rel_tolerance"] = 1e-24
        solvercg.parameters["zero_initial_guess"] = True
        solvercg.parameters["max_iter"] = 5000
        solvercg.parameters["print_level"] = 1
        if not PRINT:
            solvercg.parameters["print_level"] = -1
        ahat.zero()
        solvercg.solve(ahat, -mg)
        c, r, m = jointmodel.cost(x)
        normahat = ahat.norm('l2')
        normmg = mg.norm('l2')
        if PRINT:
            print 'Nb of Hessian-vect={}'.format(H.ncalls)
            print 'check--|ahat|={}, |mg|={}, gradnorm={}'.format(normahat, normmg, gradnorm)
            print 'check--c={}, r={}, m={}'.format(c, r, m)


        dl.MPI.barrier(dl.mpi_comm_world())
        # Prior-preconditioned
        if PRINT:
            print '\nPreconditioner: prior'
        ahat = jointmodel.generate_vector(PARAMETER)
        mg = jointmodel.generate_vector(PARAMETER)
        jointmodel.solveFwd(x[STATE], x, solver.parameters["inner_rel_tolerance"])
        jointmodel.solveAdj(x[ADJOINT], x, solver.parameters["inner_rel_tolerance"])
        jointmodel.setPointForHessianEvaluations(x)
        gradnorm = jointmodel.evalGradientParameter(x, mg)
        H = ReducedHessian(jointmodel, solver.parameters["inner_rel_tolerance"], 
        gauss_newton_approx=False, misfit_only=False)
        solvercg = CGSolverSteihaug()
        solvercg.set_operator(H)
        solvercg.set_preconditioner(jointmodel.Rsolver())
        solvercg.parameters["rel_tolerance"] = 1e-24
        solvercg.parameters["zero_initial_guess"] = True
        solvercg.parameters["max_iter"] = 5000
        solvercg.parameters["print_level"] = 1
        if not PRINT:
            solvercg.parameters["print_level"] = -1
        ahat.zero()
        solvercg.solve(ahat, -mg)
        c, r, m = jointmodel.cost(x)
        normahat = ahat.norm('l2')
        normmg = mg.norm('l2')
        if PRINT:
            print 'Nb of Hessian-vect={}'.format(H.ncalls)
            print 'check--|ahat|={}, |mg|={}, gradnorm={}'.format(normahat, normmg, gradnorm)
            print 'check--c={}, r={}, m={}'.format(c, r, m)


#        dl.MPI.barrier(dl.mpi_comm_world())
#        # BFGS-d0
#        if PRINT:
#            print '\nPreconditioner: BFGS_H0=d0.I'
#        ahat = jointmodel.generate_vector(PARAMETER)
#        mg = jointmodel.generate_vector(PARAMETER)
#        jointmodel.solveFwd(x[STATE], x, solver.parameters["inner_rel_tolerance"])
#        jointmodel.solveAdj(x[ADJOINT], x, solver.parameters["inner_rel_tolerance"])
#        jointmodel.setPointForHessianEvaluations(x)
#        gradnorm = jointmodel.evalGradientParameter(x, mg)
#        bfgsPC = solver.bfgsPC
#        bfgsPC.isH0invdefault = True
#        bfgsPC.H0inv = H0invdefault()
#        if PRINT:
#            print 'd0={}'.format(bfgsPC.H0inv.d0)
#        bfgsPC.updated0()
#        bfgsPC.isupdated0 = False
#        if PRINT:
#            print 'd0={}'.format(bfgsPC.H0inv.d0)
#        H = ReducedHessian(jointmodel, solver.parameters["inner_rel_tolerance"], 
#        gauss_newton_approx=False, misfit_only=False)
#        solvercg = CGSolverSteihaug()
#        solvercg.set_operator(H)
#        solvercg.set_preconditioner(bfgsPC)
#        solvercg.parameters["rel_tolerance"] = 1e-24
#        solvercg.parameters["zero_initial_guess"] = True
#        solvercg.parameters["max_iter"] = 5000
#        solvercg.parameters["print_level"] = 1
#        if not PRINT:
#            solvercg.parameters["print_level"] = -1
#        ahat.zero()
#        solvercg.solve(ahat, -mg)
#        c, r, m = jointmodel.cost(x)
#        normahat = ahat.norm('l2')
#        normmg = mg.norm('l2')
#        if PRINT:
#            print 'Nb of Hessian-vect={}'.format(H.ncalls)
#            print 'check--|ahat|={}, |mg|={}, gradnorm={}'.format(normahat, normmg, gradnorm)
#            print 'check--c={}, r={}, m={}'.format(c, r, m)



    if EIG:
        jointmodel.setPointForHessianEvaluations(x)
        H = ReducedHessian(jointmodel, solver.parameters["inner_rel_tolerance"], 
        gauss_newton_approx=False, misfit_only=False)
        k = x[PARAMETER].size()
        p = 20
        if rank == 0:
            print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)

        Omega = MultiVector(x[PARAMETER], k+p)
        for i in range(k+p):
            Random.normal(Omega[i], 1., True)

        d, U = doublePassG(H, jointmodel.M[PARAMETER], jointmodel.Msolver, Omega, k, s=1, check=False)
        if rank == 0:
            np.savetxt('eigenvaluesatMAP' + filename + '.txt', d)

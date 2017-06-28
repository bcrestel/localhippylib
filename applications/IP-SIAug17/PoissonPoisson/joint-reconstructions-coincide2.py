import sys, os
import numpy as np
import dolfin as dl
from dolfin import Constant

from hippylib import Random, ZeroPrior, Model, ReducedSpaceNewtonCG, BFGS,\
STATE, ADJOINT, PARAMETER, amg_method, vector2Function
from hippylib.jointmodel import JointModel

from fenicstools.regularization import TV, TVPD
from fenicstools.jointregularization import \
SumRegularization, V_TVPD, NuclearNormSVD2D
from fenicstools.plotfenics import PlotFenics

from targetmedium_coincide2 import targetmedium
from definePDE import pdes
from definemisfit import defmisfit

PLOT = True
SOLVER = 'Newton'
suffix = '-c2'

if __name__ == "__main__":
    dl.set_log_active(False)

    try:
        k = float(sys.argv[1])
    except:
        k = 1e-7


    nx = 64
    ny = 64
    mesh = dl.UnitSquareMesh(nx, ny)
    
    rank = dl.MPI.rank(mesh.mpi_comm())
    nproc = dl.MPI.size(mesh.mpi_comm())
    
    if nproc > 1:
        Random.split(rank, nproc, 1000000, 1)
    
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    # Target medium parameters
    a1true, a2true = targetmedium(Vh, PARAMETER)
    if PLOT:
        #PltFen = PlotFenics(Outputfolder=os.path.splitext(sys.argv[0])[0] + '/Plots', 
        PltFen = PlotFenics(Outputfolder='output/plots', 
        comm=mesh.mpi_comm())
        PltFen.set_varname('model1' + suffix)
        PltFen.plot_vtk(a1true)
        PltFen.set_varname('model2' + suffix)
        PltFen.plot_vtk(a2true)

    # Define PDE
    pde1, pde2 = pdes(Vh, STATE, amg_method)
 
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
        if rank == 0:
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
    if rank == 0:
        print 'Cost @ MAP for m1: cost={}, misfit={}, reg={}'.format(c1, m1, r1)
        print 'Cost @ MAP for m2: cost={}, misfit={}, reg={}'.format(c2, m2, r2)

    ############ Regularization #############
    if SOLVER == 'Newton':
        reg1 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':3e-7, 
        'rescaledradiusdual':1.0, 'print':(not rank)})
        reg2 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':4e-7, 
        'rescaledradiusdual':1.0, 'print':(not rank)})
    elif SOLVER == 'BFGS':
        reg1 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':3e-7, 
        'GNhessian':True, 'print':(not rank)})
        reg2 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':4e-7, 
        'GNhessian':True, 'print':(not rank)})
    jointregul = SumRegularization(reg1, reg2, 
    coeff_cg=0.0,
    coeff_ncg=k, parameters_ncg={'eps':1e-5},
    coeff_vtv=0.0, parameters_vtv={'eps':1e-3, 'k':5e-9, 'rescaledradiusdual':1.0},
    isprint=(not rank))
    suffix += '-TV-e' + str(reg1.parameters['eps']) \
    + '-CG' + str(jointregul.coeff_cg) \
    + '-NCG' + str(jointregul.coeff_ncg) 

    """
    jointregul = V_TVPD(Vh[PARAMETER], {'k':3e-7, 'eps':1e-3, \
    'rescaledradiusdual':1.0, 'print':not rank})
    suffix += '-VTV-k' + str(jointregul.parameters['k']) \
    + '-e' + str(jointregul.parameters['eps'])
    """

    """
    jointregul = NuclearNormSVD2D(mesh, {'eps':1e-3, 'k':3e-7}, isprint=(not rank))
    suffix += '-NN-k' + str(jointregul.parameters['k']) \
    + '-e' + str(jointregul.parameters['eps'])
    """
    #########################################


    jointmodel = JointModel(model1, model2, jointregul,
    parameters={'print':(not rank), 'splitassign':True})


    #########################################
    if SOLVER == 'Newton':
        solver = ReducedSpaceNewtonCG(jointmodel)
        suffix += '-Newton'
    elif SOLVER == 'BFGS':
        solver = BFGS(jointmodel)
        #solver.apply_H0 = solver.apply_Rinv #note: for TV+CG only
        suffix += '-BFGS'
    #########################################

    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 25 # !!! very large
    solver.parameters["GN_iter"] = 10
    solver.parameters["max_iter"] = 10000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    a0 = dl.interpolate(dl.Constant(("0.0","0.0")),jointmodel.Vh[PARAMETER])

    if SOLVER == 'Newton':
        x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[-10., 10.])
    elif SOLVER == 'BFGS':
        x = solver.solve(a0.vector(), bounds_xPARAM=[-10., 10.])

    x1, x2 = jointmodel.splitvector(x)
    minaf1 = dl.MPI.min(mesh.mpi_comm(), np.amin(x1[PARAMETER].array()))
    maxaf1 = dl.MPI.max(mesh.mpi_comm(), np.amax(x1[PARAMETER].array()))
    md1mis, md1misperc = model1.mediummisfit(x1[PARAMETER])
    minaf2 = dl.MPI.min(mesh.mpi_comm(), np.amin(x2[PARAMETER].array()))
    maxaf2 = dl.MPI.max(mesh.mpi_comm(), np.amax(x2[PARAMETER].array()))
    md2mis, md2misperc = model2.mediummisfit(x2[PARAMETER])
    if rank == 0:
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
        PltFen.set_varname('joint-model1' + suffix)
        PltFen.plot_vtk(vector2Function(x1[PARAMETER], Vh[PARAMETER]))
        PltFen.set_varname('joint-model2' + suffix)
        PltFen.plot_vtk(vector2Function(x2[PARAMETER], Vh[PARAMETER]))



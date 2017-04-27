import sys
import numpy as np
import matplotlib.pyplot as plt
import dolfin as dl
from dolfin import Expression

from hippylib import *

from fenicstools.prior import LaplacianPrior
from fenicstools.regularization import TV, TVPD
from fenicstools.jointregularization import \
SumRegularization, VTV, V_TV, V_TVPD, NuclearNormformula, NuclearNormSVD2D
from fenicstools.plotfenics import PlotFenics


PLOT = True

def u_boundary(x, on_boundary):
    return on_boundary

if __name__ == "__main__":
    dl.set_log_active(False)
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
    # coincide:
    a1true = dl.interpolate(dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '4*(x[0]<=0.5) + 8*(x[0]>0.5) ))'), Vh[PARAMETER])
#    # coincide2:
#    a1true = dl.interpolate(dl.Expression('log(10 - ' + \
#    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * 8 )'), Vh[PARAMETER])
    a2true = dl.interpolate(dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '8*(x[0]<=0.5) + 4*(x[0]>0.5) ))'), Vh[PARAMETER])
    if PLOT:
        PltFen = PlotFenics(comm=mesh.mpi_comm())
        PltFen.set_varname('jointa1')
        PltFen.plot_vtk(a1true)
        PltFen.set_varname('jointa2')
        PltFen.plot_vtk(a2true)

    # Define PDE
    f = dl.Constant(1.0)
    u_bdr = dl.Constant(0.0)
    u_bdr0 = dl.Constant(0.0)
    bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)
    
    def pde_varf(u,a,p):
        return dl.exp(a)*dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - f*p*dl.dx
    
    pde1 = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    pde2 = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    for pde in [pde1, pde2]:
        pde.solver = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver.parameters["relative_tolerance"] = 1e-15
        #pde.solver.parameters["absolute_tolerance"] = 1e-20
        pde.solver_fwd_inc = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver_fwd_inc.parameters = pde.solver.parameters
        pde.solver_adj_inc = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver_adj_inc.parameters = pde.solver.parameters
 
    # Define misfit functions
    nbobsperdir=50
    targets1 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range((nbobsperdir+2)/2, nbobsperdir+1) \
    for j in range((nbobsperdir+2)/2, nbobsperdir+1)])
    targets2 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])
    misfit1 = PointwiseStateObservation(Vh[STATE], targets1)
    misfit2 = PointwiseStateObservation(Vh[STATE], targets2)

    #TODO: apply same noise in serial and parallell
    # could save noise vector, along with corresponding coordinates to a file,
    # then read noise in serial/parallel for each observation point
    # 1. generate noise in serial for each obs points
    # 2. save to file: x|y|noise1|noise2    
    # 3. read file from each processor
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
    #reg1 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #reg2 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})

    #reg1 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})
    #reg2 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})

    #reg1 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':3e-7, 
    #'rescaledradiusdual':1.0, 'print':(not rank)})
    #reg2 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':4e-7, 
    #'rescaledradiusdual':1.0, 'print':(not rank)})

    #jointregul = SumRegularization(reg1, reg2, mesh.mpi_comm(), 
    #coeff_cg=0.0,
    #coeff_ncg=1e-4, parameters_ncg={'eps':1e-4},
    #coeff_vtv=0.0, parameters_vtv={'eps':1e-3, 'k':5e-9, 'rescaledradiusdual':1.0},
    #isprint=(not rank))
    #jointregul = Tikhonovab({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #jointregul = VTV(Vh[PARAMETER], {'k':4e-7, 'eps':1e-2})
    #jointregul = V_TV(Vh[PARAMETER], {'k':4e-7, 'eps':1e-3})
    #jointregul = V_TVPD(Vh[PARAMETER], {'k':4e-7, 'eps':1e-3, \
    #'rescaledradiusdual':1.0, 'print':not rank})
    #jointregul = NuclearNormformula(mesh, {'eps':1000., 'k':1e-7}, 
    #isprint=(not rank))
    jointregul = NuclearNormSVD2D(mesh, {'eps':1e-8, 'k':1e-7}, isprint=(not rank))
    #########################################
    try:
        if jointregul.coeff_cg > 0.0:
            plot_suffix = 'TVPD+' + str(jointregul.coeff_cg) + 'CG'
        elif jointregul.coeff_ncg > 0.0:
            plot_suffix = 'TVPD+' + str(jointregul.coeff_ncg) + 'NCG'
    except:
        #plot_suffix = 'NN'
        #plot_suffix = 'TVPD'
        plot_suffix = 'NNsvd'
        plot_suffix += '-e' + str(jointregul.parameters['eps']) \
        + '-k' + str(jointregul.parameters['k'])
    #######################

    jointmodel = JointModel(model1, model2, jointregul,
    parameters={'print':(not rank), 'splitassign':True})

    #solver = ReducedSpaceNewtonCG(jointmodel)
    solver = BFGS(jointmodel)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
#    solver.parameters["rel_tolerance"] = 1e-10
#    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 20 # !!! very large
    solver.parameters["GN_iter"] = 10
    solver.parameters["max_iter"] = 10000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    a0 = dl.interpolate(dl.Expression(("0.0","0.0")),jointmodel.Vh[PARAMETER])
    #x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[-10., 25.])
    x = solver.solve(a0.vector(), bounds_xPARAM=[-10., 25.])

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
        PltFen.set_varname('jointsolution1-' + plot_suffix)
        PltFen.plot_vtk(vector2Function(x1[PARAMETER], Vh[PARAMETER]))
        PltFen.set_varname('jointsolution2-' + plot_suffix)
        PltFen.plot_vtk(vector2Function(x2[PARAMETER], Vh[PARAMETER]))

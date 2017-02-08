# Copyright (c) 2016, The University of Texas at Austin & University of
# California, Merced.
#
# All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the hIPPYlib library. For more information and source code
# availability see https://hippylib.github.io.
#
# hIPPYlib is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License (as published by the Free
# Software Foundation) version 3.0 dated June 2007.

import dolfin as dl
import math
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append( "../../" )
from hippylib import *

from fenicstools.prior import LaplacianPrior
from fenicstools.plotfenics import PlotFenics

def u_boundary(x, on_boundary):
    return on_boundary

            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    ndim = 2
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
    
    ndofs = [Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim()]
    if rank == 0:
        print sep, "Set up the mesh and finite element spaces", sep
        print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(*ndofs)
    
    # Target medium paramters
    a1true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '4*(x[0]<=0.5) + 8*(x[0]>0.5) ))')
    a2true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '8*(x[0]<=0.5) + 4*(x[0]>0.5) ))')
        
    # Initialize Expressions
    f = dl.Constant(1.0)
    u_bdr = dl.Constant(0.0)
    u_bdr0 = dl.Constant(0.0)
    bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)
    
    def pde_varf(u,a,p):
        return dl.exp(a)*dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - f*p*dl.dx
    
    pde = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    pde.solver = dl.PETScKrylovSolver("cg", amg_method())
    pde.solver.parameters["relative_tolerance"] = 1e-15
    #pde.solver.parameters["absolute_tolerance"] = 1e-20    #TODO: uncomment line
    pde.solver_fwd_inc = dl.PETScKrylovSolver("cg", amg_method())
    pde.solver_fwd_inc.parameters = pde.solver.parameters
    pde.solver_adj_inc = dl.PETScKrylovSolver("cg", amg_method())
    pde.solver_adj_inc.parameters = pde.solver.parameters
 
    nbobsperdir=50
    targets1 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range((nbobsperdir+2)/2, nbobsperdir+1) \
    for j in range((nbobsperdir+2)/2, nbobsperdir+1)])
    targets2 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])

    #TODO: add medium misfit to model, and compare old and new models
    SELECTMODEL = 2 ###### CHANGE THIS VALUE ########
    PltFen = PlotFenics()
    if SELECTMODEL == 1:
        atrue = dl.interpolate(a1true, Vh[PARAMETER])
        targets = targets1
        PltFen.set_varname('a1')
        PltFen.plot_vtk(atrue)
        PltFen.set_varname('solutionptwise1')
    else:
        atrue = dl.interpolate(a2true, Vh[PARAMETER])
        targets = targets2
        PltFen.set_varname('a2')
        PltFen.plot_vtk(atrue)
        PltFen.set_varname('solutionptwise2')

    misfit = PointwiseStateObservation(Vh[STATE], targets)
    
    prior = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':5e-8, 'beta':5e-8})
    
    #Generate synthetic observations
    utrue = pde.generate_state()
    x = [utrue, atrue.vector(), None]
    minatrue = dl.MPI.min(mesh.mpi_comm(), np.amin(atrue.vector().array()))
    maxatrue = dl.MPI.max(mesh.mpi_comm(), np.amax(atrue.vector().array()))
    if rank == 0:
        print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
    pde.solveFwd(x[STATE], x, 1e-9)
    rel_noise_level = 0.02
    noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
    Random.normal(x[STATE], noise_level, False)
    misfit.B.mult(x[STATE], misfit.d)
    misfit.noise_variance = np.sqrt(targets.shape[0])   # hack to compare both models
    
    model = Model(pde, prior, misfit)
    
#    if rank == 0:
#        print sep, "Test the gradient and the Hessian of the model", sep
#    
#    a0 = dl.interpolate(dl.Expression("sin(x[0])", element=Vh[PARAMETER].ufl_element() ), Vh[PARAMETER])
#    modelVerify(model, a0.vector(), 1e-12, is_quadratic = False, verbose = (rank == 0) )

    if rank == 0:
        print sep, "Find the MAP point", sep
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 12
    solver.parameters["GN_iter"] = 0
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    InexactCG = 1
    GN = True
    a0 = dl.interpolate(dl.Expression("0.0"),Vh[PARAMETER])
    x = solver.solve(a0.vector(), InexactCG, GN)
    
    if rank == 0:
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"

        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost

    # Plot reconstruction
    PltFen.plot_vtk(vector2Function(x[PARAMETER], Vh[PARAMETER]))
        
    if nproc == 1:
        xx = [vector2Function(x[i], Vh[i]) for i in range(len(Vh))]
        dl.plot(xx[STATE], title = "State")
        dl.plot(dl.exp(xx[PARAMETER]), title = "exp(Parameter)")
        dl.interactive()

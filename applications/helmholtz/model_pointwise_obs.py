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

#from fenicstools.prior import LaplacianPrior
from fenicstools.regularization import TV, TVPD
from fenicstools.plotfenics import PlotFenics
from fenicstools.sourceterms import PointSources

def u_boundary(x, on_boundary):
    return on_boundary and not x[1] > 1.0 - 1e-16

            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    ndim = 2
    nx = 50
    ny = 50
    
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
    
    # Target medium parameters
    freq = 5.0
    w = 2*np.pi*freq
    wsq = dl.Constant(w*w)
    c = 3.0
    k = w/c
    l = 2*np.pi/k
    print 'f={}, w={}, c={},\nk={}, l={}, h_max={}'.format(freq, w, c, k, l, 2*np.pi/(10.*k))
    atrue = dl.interpolate(dl.Expression('log(m)', m=1./(c*c)), Vh[PARAMETER])
    dl.plot(atrue, title="atrue")
        
    # Define PDE
    #f = dl.Expression("10*(pow( pow(x[0]-0.5,2) + pow(x[1]-0.5,2), 0.5 ) < 0.1)")
    f = dl.Expression("sin(x[0]*pi)*sin(x[1]*pi)")
    #test = dl.TestFunction(Vh[STATE])
    #f = dl.assemble(dl.Constant(0.0)*test*dl.dx)
    #ptsrc = dl.PointSource(Vh[STATE], dl.Point(2, np.array([0.5, 0.5])))
    #ptsrc.apply(f)
    u_bdr = dl.Constant(0.0)
    u_bdr0 = dl.Constant(0.0)
    bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)
    
    def pde_varf(u,a,p):
        return dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - dl.inner(wsq*dl.exp(a)*u, p)*dl.dx - f*p*dl.dx

    # Helmholtz operator is indefinite -> change solver (LU or gmres)
    pde = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    pde.solver = dl.PETScLUSolver("petsc")
    pde.solver.parameters['symmetric'] = True
#    pde.solver = dl.PETScKrylovSolver("cg", amg_method())
#    pde.solver.parameters["relative_tolerance"] = 1e-15
#    pde.solver.parameters["absolute_tolerance"] = 1e-20
#    pde.solver_fwd_inc = dl.PETScKrylovSolver("cg", amg_method())
    pde.solver_fwd_inc = dl.PETScLUSolver("petsc")
    pde.solver_fwd_inc.parameters = pde.solver.parameters
#    pde.solver_adj_inc = dl.PETScKrylovSolver("cg", amg_method())
    pde.solver_adj_inc = dl.PETScLUSolver("petsc")
    pde.solver_adj_inc.parameters = pde.solver.parameters

    # Define misfit functions
    nbobsperdir=20
    targets = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])
    misfit = PointwiseStateObservation(Vh[STATE], targets)

    # Generate synthetic observations
    rel_noise_level = 0.1
    utrue = pde.generate_state()
    x = [utrue, atrue.vector(), None]
    minatrue = dl.MPI.min(mesh.mpi_comm(), np.amin(atrue.vector().array()))
    maxatrue = dl.MPI.max(mesh.mpi_comm(), np.amax(atrue.vector().array()))
    if rank == 0:
        print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
    pde.solveFwd(x[STATE], x, 1e-9)
    noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
    dl.plot(vector2Function(x[STATE], Vh[STATE]), title="u")
    Random.normal(x[STATE], noise_level, False)
    dl.plot(vector2Function(x[STATE], Vh[STATE]), title="u_noisy")
    dl.interactive()
    misfit.B.mult(x[STATE], misfit.d)
    misfit.noise_variance = 1.0
    
    # Regularization
    prior = LaplacianPrior(Vh[PARAMETER], 1e-2, 1e-2, atrue.vector())
    #prior = TVPD({'Vm':Vh[PARAMETER], 'k':5e-7, 'eps':1e-3})
    
    model = Model(pde, prior, misfit, atrue.vector())
    x[STATE].zero()
    c, r, m = model.cost(x)
    if rank == 0:
        print 'Cost @ MAP: cost={}, misfit={}, reg={}'.format(c, m, r)
    
    if rank == 0:
        print sep, "Test the gradient and the Hessian of the model", sep
    
    a0 = dl.interpolate(dl.Expression("sin(x[0])", element=Vh[PARAMETER].ufl_element() ), Vh[PARAMETER])
    modelVerify(model, a0.vector(), 1e-12, is_quadratic = False, verbose = (rank == 0) )

    if rank == 0:
        print sep, "Find the MAP point", sep
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 12
    solver.parameters["GN_iter"] = 5
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    #a0 = dl.interpolate(dl.Expression("0.0"),Vh[PARAMETER])
    x = solver.solve(atrue.vector(), InexactCG=1, GN=False)
    
    minaf = dl.MPI.min(mesh.mpi_comm(), np.amin(x[PARAMETER].array()))
    maxaf = dl.MPI.max(mesh.mpi_comm(), np.amax(x[PARAMETER].array()))
    mdmis, mdmisperc = model.mediummisfit(x[PARAMETER])
    if rank == 0:
        print 'min(af)={}, max(af)={}, medmisft={:e} ({:.1f}%)'.format(\
        minaf, maxaf, mdmis, mdmisperc)
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"

        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost

    # Plot reconstruction
    dl.plot(vector2Function(x[PARAMETER], Vh[PARAMETER]), "asol")

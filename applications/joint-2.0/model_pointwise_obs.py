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
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append( "../../" )
from hippylib import *

from fenicstools.prior import LaplacianPrior
from fenicstools.regularization import TV, TVPD
from fenicstools.plotfenics import PlotFenics

def u_boundary(x, on_boundary):
    return on_boundary

class Poisson:
    def __init__(self, mesh, Vh, atrue, targets, prior, noiselevel, alphareg=1.0):
        """
        Construct a model by proving
        - the mesh
        - the finite element spaces for the STATE/ADJOINT variable and the PARAMETER variable
        - the Prior information
        """
        self.mesh = mesh
        self.Vh = Vh
        
        # Initialize Expressions
        self.atrue = atrue
        self.f = dl.Constant(1.0)
        self.u_o = dl.Vector()
        
        self.u_bdr = dl.Constant(0.0)
        self.u_bdr0 = dl.Constant(0.0)
        self.bc = dl.DirichletBC(self.Vh[STATE], self.u_bdr, u_boundary)
        self.bc0 = dl.DirichletBC(self.Vh[STATE], self.u_bdr0, u_boundary)
                
        # Assemble constant matrices      
        self.Prior = prior
        self.B = assemblePointwiseObservation(Vh[STATE], targets)

        self.alphareg = alphareg
                
        _ = self.computeObservation(self.u_o, noiselevel, mesh.mpi_comm())
                
        self.A = []
        self.At = []
        self.C = []
        self.Raa = []
        self.Wau = []
        
    def generate_vector(self, component="ALL"):
        """
        Return the list x=[u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the parameter.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component is STATE, PARAMETER, or ADJOINT return x[component]
        """
        if component == "ALL":
            x = [dl.Vector(), dl.Vector(), dl.Vector()]
            self.B.init_vector(x[STATE], 1)
            self.Prior.init_vector(x[PARAMETER], 0)
            self.B.init_vector(x[ADJOINT], 1)
        elif component == STATE:
            x = dl.Vector()
            self.B.init_vector(x, 1)
        elif component == PARAMETER:
            x = dl.Vector()
            self.Prior.init_vector(x, 0)
        elif component == ADJOINT:
            x = dl.Vector()
            self.B.init_vector(x, 1)
            
        return x
    
    def init_parameter(self, a):
        """
        Reshape a so that it is compatible with the parameter variable
        """
        self.Prior.init_vector(a, 0)
        
    def assembleA(self, x, assemble_adjoint = False, assemble_rhs = False):
        """
        Assemble the matrices and rhs for the forward/adjoint problems
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test = dl.TestFunction(self.Vh[STATE])
        c = vector2Function(x[PARAMETER], self.Vh[PARAMETER])
        Avarf = dl.inner(dl.exp(c)*dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        if not assemble_adjoint:
            bform = dl.inner(self.f, test)*dl.dx
            Matrix, rhs = dl.assemble_system(Avarf, bform, self.bc)
        else:
            # Assemble the adjoint of A (i.e. the transpose of A)
            s = vector2Function(x[STATE], self.Vh[STATE])
            bform = dl.inner(dl.Constant(0.), test)*dl.dx
            Matrix, _ = dl.assemble_system(dl.adjoint(Avarf), bform, self.bc0)
            Bu = -(self.B*x[STATE])
            Bu += self.u_o
            rhs = dl.Vector()
            self.B.init_vector(rhs, 1)
            self.B.transpmult(Bu,rhs)
            #rhs *= 1.0/self.noise_variance
            
        if assemble_rhs:
            return Matrix, rhs
        else:
            return Matrix
    
    def assembleC(self, x):
        """
        Assemble the derivative of the forward problem with respect to the parameter
        """
        trial = dl.TrialFunction(self.Vh[PARAMETER])
        test = dl.TestFunction(self.Vh[STATE])
        s = vector2Function(x[STATE], Vh[STATE])
        c = vector2Function(x[PARAMETER], Vh[PARAMETER])
        Cvarf = dl.inner(dl.exp(c) * trial * dl.nabla_grad(s), dl.nabla_grad(test)) * dl.dx
        C = dl.assemble(Cvarf)
        self.bc0.zero(C)
        return C
                
    def assembleWau(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the state
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        a = vector2Function(x[ADJOINT], Vh[ADJOINT])
        c = vector2Function(x[PARAMETER], Vh[PARAMETER])
        varf = dl.inner(dl.exp(c)*dl.nabla_grad(trial),dl.nabla_grad(a))*test*dl.dx
        Wau = dl.assemble(varf)
        Wau_t = Transpose(Wau)
        self.bc0.zero(Wau_t)
        Wau = Transpose(Wau_t)
        return Wau
    
    def assembleRaa(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the parameter (Newton method)
        """
        trial = dl.TrialFunction(self.Vh[PARAMETER])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        s = vector2Function(x[STATE], Vh[STATE])
        c = vector2Function(x[PARAMETER], Vh[PARAMETER])
        a = vector2Function(x[ADJOINT], Vh[ADJOINT])
        varf = dl.inner(dl.nabla_grad(a),dl.exp(c)*dl.nabla_grad(s))*trial*test*dl.dx
        return dl.assemble(varf)

        
    def computeObservation(self, u_o, rel_noise_level, mpi_comm):
        """
        Compute the syntetic observation
        """
        self.at = dl.interpolate(self.atrue, self.Vh[PARAMETER])
        rank = dl.MPI.rank(mpi_comm)
        minatrue = dl.MPI.min(mpi_comm, np.amin(self.at.vector().array()))
        maxatrue = dl.MPI.max(mpi_comm, np.amax(self.at.vector().array()))
        if rank == 0:
            print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
        x = [self.generate_vector(STATE), self.at.vector(), None]
        self.solveFwd(x[STATE], x, tol=1e-9)
        
        # Create noisy data, ud
        noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(self.Vh[PARAMETER].dim())
        Random.normal(x[STATE], noise_level, False)
        
        self.B.init_vector(u_o,0)
        self.B.mult(x[STATE], u_o)

        x[STATE].zero()
        c, r, m = self.cost(x)
        if rank == 0:
            print 'Cost @ MAP: cost={}, misfit={}, reg={}'.format(c, m, r)
        
        return noise_level*noise_level
        
    def mediummisfit(self, m):
        """
        Compute medium misfit
        """
        diff = m - self.at.vector()
        nd = dl.norm(diff)
        return nd, 100.*nd/dl.norm(self.at.vector())
    
    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, parameter, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        assert x[STATE] != None
                       
        diff = self.B*x[STATE]
        diff -= self.u_o
        misfit = .5 * diff.inner(diff)
        #misfit = (.5/self.noise_variance) * diff.inner(diff)
        
        reg = self.Prior.costvect(x[PARAMETER])
        
        c = misfit + reg
        
        return c, reg, misfit
    
    def solveFwd(self, out, x, tol=1e-9):
        """
        Solve the forward problem.
        """
        A, b = self.assembleA(x, assemble_rhs = True)
        A.init_vector(out, 1)
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.parameters["relative_tolerance"] = tol
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["nonzero_initial_guess"] = False
        solver.set_operator(A)
        nit = solver.solve(out,b)
        
    
    def solveAdj(self, out, x, tol=1e-9):
        """
        Solve the adjoint problem.
        """
        At, badj = self.assembleA(x, assemble_adjoint = True, assemble_rhs = True)
        At.init_vector(out, 1)
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.parameters["relative_tolerance"] = tol
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["nonzero_initial_guess"] = False
        solver.set_operator(At)
        nit = solver.solve(out,badj)
        
    
    def evalGradientParameter(self,x, mg):
        """
        Evaluate the gradient for the variation parameter equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the parameter space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        C = self.assembleC(x)

        self.Prior.init_vector(mg,0)
        C.transpmult(x[ADJOINT], mg)

        mg.axpy(self.alphareg, self.Prior.gradvect(x[PARAMETER]))
        
        g = dl.Vector()
        self.Prior.init_vector(g,1)
        
        self.Prior.Msolver.solve(g, mg)
        g_norm = dl.sqrt( g.inner(mg) )
        
        return g_norm
        
    
    def setPointForHessianEvaluations(self, x):  
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        """      
        self.A  = self.assembleA(x)
        self.At = self.assembleA(x, assemble_adjoint=True )
        self.C  = self.assembleC(x)
        self.Wau = self.assembleWau(x)
        self.Raa = self.assembleRaa(x)
        self.Prior.assemble_hessian(x[PARAMETER])

        
    def solveFwdIncremental(self, sol, rhs, tol):
        """
        Solve the incremental forward problem for a given rhs
        """
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.set_operator(self.A)
        solver.parameters["relative_tolerance"] = tol
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["nonzero_initial_guess"] = False
        self.A.init_vector(sol,1)
        nit = solver.solve(sol,rhs)
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        """
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.set_operator(self.At)
        solver.parameters["relative_tolerance"] = tol
        solver.parameters["error_on_nonconvergence"] = True
        solver.parameters["nonzero_initial_guess"] = False
        self.At.init_vector(sol,1)
        nit = solver.solve(sol, rhs)
    
    def applyC(self, da, out):
        self.C.mult(da,out)
    
    def applyCt(self, dp, out):
        self.C.transpmult(dp,out)
    
    def applyWuu(self, du, out, gn_approx=False):
        help = dl.Vector()
        self.B.init_vector(help, 0)
        self.B.mult(du, help)
        self.B.transpmult(help, out)
        #out *= 1./self.noise_variance
    
    def applyWua(self, da, out):
        self.Wau.transpmult(da,out)
    
    def applyWau(self, du, out):
        self.Wau.mult(du, out)
    
    def applyR(self, da, out):
        out.zero()
        out.axpy(self.alphareg, self.Prior.hessian(da))
        
    def Rsolver(self):        
        return self.Prior.getprecond()
    
    def applyRaa(self, da, out):
        self.Raa.mult(da, out)
            
if __name__ == "__main__":
    dl.set_log_active(False)
    nx, ny = 64, 64 
    mesh = dl.UnitSquareMesh(nx, ny)
    
    rank = dl.MPI.rank(mesh.mpi_comm())
    nproc = dl.MPI.size(mesh.mpi_comm())
    
    if nproc > 1:
        Random.split(rank, nproc, 1000000, 1)
    
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    #Prior = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #Prior = TV({'Vm':Vh[PARAMETER], 'k':1e-8, 'eps':1e-3, 'GNhessian':False})
    Prior = TVPD({'Vm':Vh[PARAMETER], 'k':1e-9, 'eps':1e-3})

    a1true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '4*(x[0]<=0.5) + 8*(x[0]>0.5) ))')
    a2true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '8*(x[0]<=0.5) + 4*(x[0]>0.5) ))')

    nbobsperdir=50
    targets = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])

    model1 = Poisson(mesh, Vh, a1true, targets, Prior, noiselevel=0.02, alphareg=1.0)
    model2 = Poisson(mesh, Vh, a2true, targets, Prior, noiselevel=0.02, alphareg=1.0)
    PltFen = PlotFenics()
    PltFen.set_varname('a1')
    PltFen.plot_vtk(model1.at)
    PltFen.set_varname('a2')
    PltFen.plot_vtk(model2.at)

    # modify here! #######
    model = model2
    PltFen.set_varname('solutionptwise2-k1e-9')
    ######################
        
    if rank == 0 and Prior.isTV():
        print 'TV parameters: k={}, eps={}, alphareg={}'.format(\
        Prior.parameters['k'], Prior.parameters['eps'], model.alphareg)

    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-10
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 12
    solver.parameters["GN_iter"] = 0
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    InexactCG = 0
    GN = True
    a0 = dl.interpolate(dl.Expression("0.0"),Vh[PARAMETER])
    x = solver.solve(a0.vector(), InexactCG, GN)

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
    
    PltFen.plot_vtk(vector2Function(x[PARAMETER], Vh[PARAMETER]))

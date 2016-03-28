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
import sys
sys.path.append( "../" )
from hippylib import *
import numpy as np

def u_boundary(x, on_boundary):
    return on_boundary

class Poisson:
    def __init__(self, mesh, Vh, atrue, targets, prior, rel_noise_level):
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
        self.f = dl.Expression("1.0")
        self.u_o = dl.Vector()
        
        self.u_bdr = dl.Expression("0.0")
        self.u_bdr0 = dl.Expression("0.0")
        self.bc = dl.DirichletBC(self.Vh[STATE], self.u_bdr, u_boundary)
        self.bc0 = dl.DirichletBC(self.Vh[STATE], self.u_bdr0, u_boundary)
                
        # Assemble constant matrices      
        self.prior = prior
        self.B = assemblePointwiseObservation(self.Vh[STATE],targets)
                
        self.noise_variance = self.computeObservation(self.u_o, rel_noise_level)
        print "Noise variance:", self.noise_variance
                
        self.A = []
        self.At = []
        self.C = []
        self.Raa = []
        self.Wau = []
        
    def generate_vector(self, component="ALL"):
        """
        Return the list x=[u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the parameter variable.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component is STATE, PARAMETER, or ADJOINT return x[component]
        """
        if component == "ALL":
            x = [dl.Vector(), dl.Vector(), dl.Vector()]
            self.B.init_vector(x[STATE],1)
            self.prior.init_vector(x[PARAMETER],0)
            self.B.init_vector(x[ADJOINT], 1)
        elif component == STATE:
            x = dl.Vector()
            self.B.init_vector(x,1)
        elif component == PARAMETER:
            x = dl.Vector()
            self.prior.init_vector(x,0)
        elif component == ADJOINT:
            x = dl.Vector()
            self.B.init_vector(x,1)
            
        return x
    
    def init_parameter(self, a):
        """
        Reshape a so that it is compatible with the parameter variable
        """
        self.prior.init_vector(a,0)
        
    def assembleA(self,x, assemble_adjoint = False, assemble_rhs = False):
        """
        Assemble the matrices and rhs for the forward/adjoint problems
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test = dl.TestFunction(self.Vh[STATE])
        c = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        Avarf = dl.inner(dl.exp(c)*dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        if not assemble_adjoint:
            bform = dl.inner(self.f, test)*dl.dx
            Matrix, rhs = dl.assemble_system(Avarf, bform, self.bc)
        else:
            # Assemble the adjoint of A (i.e. the transpose of A)
            s = dl.Function(self.Vh[STATE], x[STATE])
            bform = dl.inner(dl.Constant(0.), test)*dl.dx
            Matrix, _ = dl.assemble_system(dl.adjoint(Avarf), bform, self.bc0)
            Bu = -(self.B*x[STATE])
            Bu += self.u_o
            rhs = dl.Vector()
            self.B.init_vector(rhs, 1)
            self.B.transpmult(Bu,rhs)
            rhs *= 1.0/self.noise_variance
            
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
        s = dl.Function(self.Vh[STATE], x[STATE])
        c = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        Cvarf = dl.inner(dl.exp(c) * trial * dl.nabla_grad(s), dl.nabla_grad(test)) * dl.dx
        C = dl.assemble(Cvarf)
#        print "||c||", x[PARAMETER].norm("l2"), "||s||", x[STATE].norm("l2"), "||C||", C.norm("linf")
        self.bc0.zero(C)
        return C
        
    def assembleWau(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the state
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        a = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        c = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        varf = dl.inner(dl.exp(c)*dl.nabla_grad(trial),dl.nabla_grad(a))*test*dl.dx
        Wau = dl.assemble(varf)
        dummy = dl.Vector()
        Wau.init_vector(dummy,0)
        self.bc0.zero_columns(Wau, dummy)
        return Wau
    
    def assembleRaa(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the parameter (Newton method)
        """
        trial = dl.TrialFunction(self.Vh[PARAMETER])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        s = dl.Function(self.Vh[STATE], x[STATE])
        c = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        a = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        varf = dl.inner(dl.nabla_grad(a),dl.exp(c)*dl.nabla_grad(s))*trial*test*dl.dx
        return dl.assemble(varf)

        
    def computeObservation(self, u_o, rel_noise_level):
        """
        Compute the syntetic observation
        """
        x = [self.generate_vector(STATE), self.atrue, None]
        A, b = self.assembleA(x, assemble_rhs = True)
        
        dl.solve(A, x[STATE], b)
        
        # Create noisy data, ud
        MAX = x[STATE].norm("linf")
        noise_level = rel_noise_level * MAX
        noise = noise_level * np.random.normal(0, 1, len(x[STATE].array()))
        x[STATE].set_local(x[STATE].array() + noise)
        
        self.B.init_vector(u_o,0)
        self.B.mult(x[STATE], u_o)
        
        return noise_level*noise_level
        
    
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
        misfit = (.5/self.noise_variance) * diff.inner(diff)
        
        Rdiff_x = dl.Vector()
        self.prior.init_vector(Rdiff_x,0)
        diff_x = x[PARAMETER] - self.prior.mean
        self.prior.R.mult(diff_x, Rdiff_x)
        reg = .5 * diff_x.inner(Rdiff_x)
        
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
        solver.set_operator(A)
        nit = solver.solve(out,b)
        
#        print "FWD", (self.A*out - b).norm("l2")/b.norm("l2"), nit

    
    def solveAdj(self, out, x, tol=1e-9):
        """
        Solve the adjoint problem.
        """
        At, badj = self.assembleA(x, assemble_adjoint = True,assemble_rhs = True)
        At.init_vector(out, 1)
        
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.parameters["relative_tolerance"] = tol
        solver.set_operator(At)
        nit = solver.solve(out,badj)
        
#        print "ADJ", (self.At*out - badj).norm("l2")/badj.norm("l2"), nit
    
    def evalGradientParameter(self,x, mg):
        """
        Evaluate the gradient for the variational parameter equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the parameter space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        C = self.assembleC(x)

        self.prior.init_vector(mg,0)
        C.transpmult(x[ADJOINT], mg)
        Rdx = dl.Vector()
        self.prior.init_vector(Rdx,0)
        dx = x[PARAMETER] - self.prior.mean
        self.prior.R.mult(dx, Rdx)   
        mg.axpy(1., Rdx)
        
        g = dl.Vector()
        self.prior.init_vector(g,1)
        
        self.prior.Msolver.solve(g, mg)
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

        
    def solveFwdIncremental(self, sol, rhs, tol):
        """
        Solve the incremental forward problem for a given rhs
        """
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.set_operator(self.A)
        solver.parameters["relative_tolerance"] = tol
        self.A.init_vector(sol,1)
        nit = solver.solve(sol,rhs)
#        print "FwdInc", (self.A*sol-rhs).norm("l2")/rhs.norm("l2"), nit
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        """
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.set_operator(self.At)
        solver.parameters["relative_tolerance"] = tol
        self.At.init_vector(sol,1)
        nit = solver.solve(sol, rhs)
#        print "AdjInc", (self.At*sol-rhs).norm("l2")/rhs.norm("l2"), nit
    
    def applyC(self, da, out):
        self.C.mult(da,out)
    
    def applyCt(self, dp, out):
        self.C.transpmult(dp,out)
    
    def applyWuu(self, du, out, gn_approx=False):
        help = dl.Vector()
        self.B.init_vector(help, 0)
        self.B.mult(du, help)
        self.B.transpmult(help, out)
        out *= 1./self.noise_variance
    
    def applyWua(self, da, out):
        self.Wau.transpmult(da,out)

    
    def applyWau(self, du, out):
        self.Wau.mult(du, out)
    
    def applyR(self, da, out):
        self.prior.R.mult(da, out)
        
    def Rsolver(self):        
        return self.prior.Rsolver
    
    def applyRaa(self, da, out):
        self.Raa.mult(da, out)


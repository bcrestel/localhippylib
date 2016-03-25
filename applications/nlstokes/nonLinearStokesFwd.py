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
# Software Foundation) version 2.1 dated February 1999.

import dolfin as dl
import math
import sys
sys.path.append( "../../" )
from hippylib import *

class IceProblem(PDEProblem):
    
    def __init__(self, Vh, ds_robin, RateFactor, GlenExp, density, gravity, eps, pen=1e3):
        """ Construct an IceProblem with the following paramenters:
        - Vh: Finite Element Space for the state variable
        - ds_robin: Boundary integrate for the Robin Boundary Condition
        - RateFactor: The (temperature dependent) flow rate factor 
        - GlenExp: The Glen's flow law exponent
        - density: The ice density
        - gravity: The gravity field.
        - eps: a regularization parameter for the Glen Flow effective viscosity. 
        """
        
        self.Vh = Vh
        self.ds_robin = ds_robin
        self.RateFactor = RateFactor
        self.GlenExp = GlenExp
        self.density = density
        self.gravity = gravity
        self.eps = eps
        self.pen = pen
        
        up = dl.TrialFunction(self.Vh[STATE])
        vq = dl.TestFunction(self.Vh[STATE])
        self.M = dl.assemble( dl.inner(up,vq)*dl.dx )
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        
        self.newton_parameters = {}
        self.newton_parameters["rel_tolerance"] = 1e-6
        self.newton_parameters["abs_tolerance"] = 1e-9
        self.newton_parameters["gdu_tolerance"] = 1e-20
        self.newton_parameters["lin_tolerance"] = 1e-12
        self.newton_parameters["max_iter"]      = 100
        self.newton_parameters["c_armijo"]      = 1e-5
        self.newton_parameters["max_backtracking_iter"] = 10
        
        self.use_Nitsche = False
        
        self.KKT = {}
        
    def generate_state(self):
        return dl.Function(self.Vh[STATE]).vector()
    
    def normF(self, F):
        f = dl.Vector()
        self.M.init_vector(f,1)
        self.Msolver.solve(f,F)
        return math.sqrt( f.inner(F) )
    
    def FNewtonean(self,up,a,vq, eta_noA):
        n = dl.FacetNormal(self.Vh[STATE].mesh())
        h = dl.CellSize(self.Vh[STATE].mesh())
        
        uh,ph = dl.split(up)
        vh,qh = dl.split(vq)
        
        myA = self.RateFactor**(-1./self.GlenExp)
        eta = myA*eta_noA
                
        def t(u,p,n): return dl.dot(eta*dl.sym(dl.nabla_grad(u)),n) - p*n
        
        a11 = eta*dl.inner(dl.sym( dl.nabla_grad(uh) ), dl.sym( dl.nabla_grad(vh) ))*dl.dx
        a12 = -dl.nabla_div(vh)*ph*dl.dx
        a21 = -dl.nabla_div(uh)*qh*dl.dx
        robinvarf  = dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), vh) * self.ds_robin
        if self.use_Nitsche:
            weak_bc = -dl.dot(n, t(uh, ph, n) )*dl.dot(vh,n)*self.ds_robin \
                      - dl.dot(n, t(vh, qh, n) )*dl.dot(uh,n)*self.ds_robin \
                      + self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        else:
            weak_bc = self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        
        f1 = dl.inner(vh, self.gravity*self.density)*dl.dx
        f2 = dl.Constant(0.)*qh*dl.dx

        if type(up) == dl.Function:        
            return a11+a12+a21+robinvarf+weak_bc - f1-f2
        elif type(up) == dl.Argument:
            return a11+a12+a21+robinvarf+weak_bc
        else:
            raise Exception()

    def Energy(self,up, a):
        
        if self.use_Nitsche:
            return None
        else:
            n = dl.FacetNormal(self.Vh[STATE].mesh())
            h = dl.CellSize(self.Vh[STATE].mesh())
            
            uh,ph = dl.split(up)
            
            myA = self.RateFactor**(-1./self.GlenExp)
            myExp = dl.Constant(.5)*( dl.Constant(1.) + self.GlenExp)/self.GlenExp
            
            Phi = dl.Constant(1.)/myExp*myA*( dl.Constant(.5)*dl.inner( dl.sym(dl.nabla_grad(uh)), dl.sym(dl.nabla_grad(uh)) )  + self.eps )**myExp
            
            a11 = Phi*dl.dx
            robinvarf  = dl.Constant(.5)*dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), uh) * self.ds_robin
            weak_bc = dl.Constant(.5)*self.pen/h*dl.dot(uh,n)*dl.dot(uh, n)*self.ds_robin
            
            f1 = dl.inner(uh, self.gravity*self.density)*dl.dx
            
            return a11+robinvarf+weak_bc - f1
        

                
    def F(self,up,a,vq):
        n = dl.FacetNormal(self.Vh[STATE].mesh())
        h = dl.CellSize(self.Vh[STATE].mesh())
        
        uh,ph = dl.split(up)
        vh,qh = dl.split(vq)
        
        myA = self.RateFactor**(-1./self.GlenExp)
        myExp = dl.Constant(.5)*( dl.Constant(1.) - self.GlenExp)/self.GlenExp
        
        eta = myA*( dl.Constant(.5)*dl.inner( dl.sym(dl.nabla_grad(uh)), dl.sym(dl.nabla_grad(uh)) )  + self.eps )**myExp
        
        def t(u,p,n): return dl.dot(eta*dl.sym(dl.nabla_grad(u)),n) - p*n
        
        a11 = eta*dl.inner(dl.sym( dl.nabla_grad(uh) ), dl.sym( dl.nabla_grad(vh) ))*dl.dx
        a12 = -dl.nabla_div(vh)*ph*dl.dx
        a21 = -dl.nabla_div(uh)*qh*dl.dx
        robinvarf  = dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), vh) * self.ds_robin
        if self.use_Nitsche:
            weak_bc = -dl.dot(n, t(uh, ph, n) )*dl.dot(vh,n)*self.ds_robin \
                      - dl.dot(n, t(vh, qh, n) )*dl.dot(uh,n)*self.ds_robin \
                      + self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        else:
            weak_bc = self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        
        f1 = dl.inner(vh, self.gravity*self.density)*dl.dx
        f2 = dl.Constant(0.)*qh*dl.dx

        if type(up) == dl.Function:        
            return a11+a12+a21+robinvarf+weak_bc - f1-f2
        elif type(up) == dl.Argument:
            return a11+a12+a21+robinvarf+weak_bc
        else:
            raise Exception()
        
    def solveFwd(self, state, x, mytol):
        # Assume u0 is div free
        up = dl.Function(self.Vh[STATE], state)
        a  = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        vq = dl.TestFunction(self.Vh[STATE])
        
        E0 = dl.assemble( self.Energy(up, a) )
        F0 = dl.assemble( self.F(up, a, vq) )
        norm_F0 = self.normF(F0)
        
        if mytol is None:
            rtol = self.newton_parameters["rel_tolerance"]
        else:
            rtol = mytol
        atol = self.newton_parameters["abs_tolerance"]
        gdu_tol  = self.newton_parameters["gdu_tolerance"]
        lin_tol = self.newton_parameters["lin_tolerance"]
        maxiter = self.newton_parameters["max_iter"]
        c_armijo = self.newton_parameters["c_armijo"]
        maxiter_backtracking = self.newton_parameters["max_backtracking_iter"]
        converged = False
        
        tol = max( norm_F0*rtol, atol)
        
        E_prev = E0
        up_prev = dl.Function(self.Vh[STATE])
        du = self.generate_state()
        for it in range(maxiter):
            f_form = self.F(up,a,vq)
            J_form = dl.derivative(f_form, up)
            J, F = dl.assemble_system(J_form, f_form)
            
            norm_F = self.normF(F)
            
            if norm_F < tol:
                print "Converged ||F||_L2 = {0:5e}.".format(norm_F)
                converged = True
                break
            
            Jsolver = dl.PETScLUSolver()
            Jsolver.set_operator(J)
            Jsolver.solve(du, F)
            
            duF = F.inner(du)
            
            if duF < gdu_tol:
                print "Converged (du,F) = {0:5e}.".format(duF)
                converged = True
                break
            
            bt_converged = False
            alpha = 1.
            up_prev.assign(up)
            
            for j in range(maxiter_backtracking):
                up.vector().axpy(-alpha, du)
                E_cur = dl.assemble( self.Energy(up, a) )
                if E_cur < E_prev - alpha*c_armijo*duF:
                    bt_converged = True
                    E_prev = E_cur
                    break
                else:
                    up.assign(up_prev)
                    alpha *= .5
            
            print "{0:3d} {1:15e} {2:15e} {3:15e} {4:15e}".format(it, E_cur, norm_F, duF, alpha)       
            if not bt_converged:
                print "Maximum number of backtracking reached."
                break
                        
        #assert converged
        
    def solveAdj(self, adj, x, adj_rhs, tol):
        up = dl.Function(self.Vh[STATE], x[STATE], name="up")
        a  = dl.Function(self.Vh[PARAMETER], x[PARAMETER], name="a")
        vq = dl.Function(self.Vh[STATE], name="vq")  
        f_form = self.F(up,a,vq)
        adj_form = dl.derivative( dl.derivative(f_form, up), vq )
        Aadj = dl.assemble(adj_form) 
        Aadjsolver = dl.PETScLUSolver()
        Aadjsolver.set_operator(Aadj)
        Aadjsolver.solve(adj, adj_rhs) 
        
    def eval_da(self, x, out):
        up = dl.Function(self.Vh[STATE], x[STATE], name="up")
        a  = dl.Function(self.Vh[PARAMETER], x[PARAMETER], name="a")
        vq = dl.Function(self.Vh[STATE], x[ADJOINT], name="vq")
        f_form = self.F(up,a,vq)
        da_form = dl.derivative(f_form, a)
        dl.assemble(da_form, tensor=out)
                    
    def setLinearizationPoint(self,x):
        up = dl.Function(self.Vh[STATE], x[STATE], name="up")
        a  = dl.Function(self.Vh[PARAMETER], x[PARAMETER], name="a")
        vq = dl.Function(self.Vh[STATE], x[ADJOINT], name="vq")
        
        x_fun = [up,a,vq]
        
        f_form = self.F(up,a,vq)
        
        g_form = [None,None,None]
        for i in range(3):
            g_form[i] = dl.derivative(f_form, x_fun[i])

        kkt_form = {}
        for i in range(3):
            for j in range(3):
                kkt_form[i,j] = dl.derivative(g_form[i], x_fun[j])
                self.KKT[i,j] = dl.assemble(kkt_form[i,j])
                
    def apply_ij(self,i,j, dir, out):
        self.KKT[i,j].mult(dir, out)
            
    def solveIncremental(self, out, rhs, is_adj, mytol):
        solver = dl.PETScLUSolver()
        if is_adj:
            solver.set_operator(self.KKT[STATE,ADJOINT])
        else:
            solver.set_operator(self.KKT[ADJOINT,STATE])
            
        solver.solve(out, rhs)

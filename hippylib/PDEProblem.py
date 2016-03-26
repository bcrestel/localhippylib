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
from variables import STATE, PARAMETER, ADJOINT

class PDEProblem:
    """ Consider the PDE Problem:
        Given a, find u s.t. 
        F(u,a,p) = ( f(u,a), p) = 0 for all p.
        Here F is linear in p, but it may be non-linear in u and a.
    """
        
    def generate_state(self):
        """ return a vector in the shape of the state """
        
    def generate_parameter(self):
        """ return a vector in the shape of the parameter """
        
    def init_parameter(self, a):
        """ initialize the parameter """
    
    def solveFwd(self, state, x, tol):
        """ Solve the possibly nonlinear Fwd Problem:
        Given a, find u such that
        \delta_p F(u,a,p;\hat_p) = 0 \for all \hat_p"""
        
    def solveAdj(self, state, x, adj_rhs, tol):
        """ Solve the linear Adj Problem: 
            Given a, u; find p such that
            \delta_u F(u,a,p;\hat_u) = 0 \for all \hat_u
        """
     
    def eval_da(self, x, out):
        """Given u,a,p; eval \delta_a F(u,a,p; \hat_a) \for all \hat_a """
         
    def setLinearizationPoint(self,x):
        """ Set the values of the state and parameter
            for the incremental Fwd and Adj solvers """
        
    def solveIncremental(self, out, rhs, is_adj, mytol):
        """ If is_adj = False:
            Solve the forward incremental system:
            Given u, a, find \tilde_u s.t.:
            \delta_{pu} F(u,a,p; \hat_p, \tilde_u) = rhs for all \hat_p.
            
            If is_adj = True:
            Solve the adj incremental system:
            Given u, a, find \tilde_p s.t.:
            \delta_{up} F(u,a,p; \hat_u, \tilde_p) = rhs for all \delta_u.
        """
    
    def apply_ij(self,i,j, dir, out):   
        """
            Given u, a, p; compute 
            \delta_{ij} F(u,a,p; \hat_i, \tilde_j) in the direction \tilde_j = dir for all \hat_i
        """

class PDEVariationalProblem(PDEProblem):
    def __init__(self, Vh, varf_handler, bc, bc0):
        self.Vh = Vh
        self.varf_handler = varf_handler
        self.bc = bc
        self.bc0 = bc0
        
        self.A  = []
        self.At = []
        self.C = []
        self.Wau = []
        self.Waa = []
        self.Wuu = []
        
    def generate_state(self):
        """ return a vector in the shape of the state """
        return dl.Function(self.Vh[STATE]).vector()
    
    def generate_parameter(self):
        """ return a vector in the shape of the parameter """
        return dl.Function(self.Vh[PARAMETER]).vector()
    
    def init_parameter(self, a):
        """ initialize the parameter """
        dummy = self.generate_parameter()
        a.init( dummy.mpi_comm(), dummy.local_range() )
    
    def solveFwd(self, state, x, tol):
        """ Solve the possibly nonlinear Fwd Problem:
        Given a, find u such that
        \delta_p F(u,a,p;\hat_p) = 0 \for all \hat_p"""
        state.set_local(x[STATE].array())
        u = dl.Function(self.Vh[STATE], state)
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.TestFunction(self.Vh[ADJOINT])
        res_form = self.varf_handler(u,a,p)
        dl.solve(res_form == 0, u, self.bc)
        
    def solveAdj(self, adj, x, adj_rhs, tol):
        """ Solve the linear Adj Problem: 
            Given a, u; find p such that
            \delta_u F(u,a,p;\hat_u) = 0 \for all \hat_u
        """
        u = dl.Function(self.Vh[STATE], x[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.Function(self.Vh[ADJOINT])
        du = dl.TestFunction(self.Vh[STATE])
        dp = dl.TrialFunction(self.Vh[ADJOINT])
        varf = self.varf_handler(u,a,p)
        adj_form = dl.derivative( dl.derivative(varf, u, du), p, dp )
        Aadj, dummy = dl.assemble_system(adj_form, dl.Constant(0.)*du*dl.dx, self.bc0)
        solver = dl.PETScLUSolver()
        solver.set_operator(Aadj)
        solver.solve(adj, adj_rhs)
     
    def eval_da(self, x, out):
        """Given u,a,p; eval \delta_a F(u,a,p; \hat_a) \for all \hat_a """
        u = dl.Function(self.Vh[STATE], x[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        da = dl.TestFunction(self.Vh[PARAMETER])
        res_form = self.varf_handler(u,a,p)
        out.zero()
        dl.assemble( dl.derivative(res_form, a, da), tensor=out)
         
    def setLinearizationPoint(self,x):
        """ Set the values of the state and parameter
            for the incremental Fwd and Adj solvers """
        u = dl.Function(self.Vh[STATE], x[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        x_fun = [u,a,p]
        
        f_form = self.varf_handler(u,a,p)
        
        g_form = [None,None,None]
        for i in range(3):
            g_form[i] = dl.derivative(f_form, x_fun[i])
            
        self.A, dummy = dl.assemble_system(dl.derivative(g_form[ADJOINT],u), g_form[ADJOINT], self.bc0)
        self.At, dummy = dl.assemble_system(dl.derivative(g_form[STATE],p),  g_form[STATE], self.bc0)
        self.C = dl.assemble(dl.derivative(g_form[ADJOINT],a))
        self.bc0.zero(self.C)
        self.Wau = dl.assemble(dl.derivative(g_form[PARAMETER],u))
        self.bc0.zero_columns(self.Wau, dummy)
        
        self.Wuu = dl.assemble(dl.derivative(g_form[STATE],u))
        self.bc0.zero(self.Wuu)
        self.bc0.zero_columns(self.Wuu, dummy)
        
        self.Waa = dl.assemble(dl.derivative(g_form[PARAMETER],a))
        
                
    def solveIncremental(self, out, rhs, is_adj, mytol):
        """ If is_adj = False:
            Solve the forward incremental system:
            Given u, a, find \tilde_u s.t.:
            \delta_{pu} F(u,a,p; \hat_p, \tilde_u) = rhs for all \hat_p.
            
            If is_adj = True:
            Solve the adj incremental system:
            Given u, a, find \tilde_p s.t.:
            \delta_{up} F(u,a,p; \hat_u, \tilde_p) = rhs for all \delta_u.
        """
        solver = dl.PETScLUSolver()
        if is_adj:
            solver.set_operator(self.At)
        else:
            solver.set_operator(self.A)
            
        solver.solve(out, rhs)
    
    def apply_ij(self,i,j, dir, out):   
        """
            Given u, a, p; compute 
            \delta_{ij} F(u,a,p; \hat_i, \tilde_j) in the direction \tilde_j = dir for all \hat_i
        """
        KKT = {}
        KKT[STATE,STATE] = self.Wuu
        KKT[PARAMETER, STATE] = self.Wau
        KKT[PARAMETER, PARAMETER] = self.Waa
        KKT[ADJOINT, STATE] = self.A
        KKT[ADJOINT, PARAMETER] = self.C
        
        if i >= j:
            KKT[i,j].mult(dir, out)
        else:
            KKT[j,i].transpmult(dir, out)

    
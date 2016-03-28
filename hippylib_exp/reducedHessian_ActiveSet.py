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

from dolfin import Vector, PETScKrylovSolver
import sys
sys.path.append( "../")
from hippylib import *

class ReducedHessianActiveSet:
    """
    Experimental implementation of ActiveSet constraints
    """
    def __init__(self, model, innerTol, gauss_newton_approx=False):
        self.H = ReducedHessian(model, innerTol, gauss_newton_approx)
        self.I_set = model.getIdentityMatrix(PARAMETER)
        self.A_set = model.getIdentityMatrix(PARAMETER)
        self.A_set.zero()
        
        self.index_active_set = None
        self.x_inactive = Vector()
        self.y_inactive = Vector()
        self.I_diag = Vector()
        
        self.I_set.init_vector(self.x_inactive, 0)
        self.I_set.init_vector(self.y_inactive, 0)
        self.I_set.init_vector(self.I_diag, 0)
        self.active_set = False
        
    def init_vector(self, x, dim):
        self.H.init_vector(x, dim)
        
    def setActiveSet(self, index_active_set):
        self.index_active_set = index_active_set
        
        self.I_diag.set_local(abs( index_active_set ))
                       
        self.A_set.set_diagonal(self.I_diag)
        self.I_set.set_diagonal(1. - self.I_diag )
        
        self.active_set = True
        
    def applyBounds(self, x, lb, ub):
        # Only works for lb = 0, ub = None 
        if self.index_active_set is None:
            return
        
        if lb is not None:
            index = self.index_active_set == 1.
            x.array()[index] = lb[index] 
        if ub is not None:
            index = self.index_active_set == -1.
            x.array()[index] = ub[index] 
            
    def getPreconditioner(self):
        
        if self.active_set == False:
            return self.H.model.RPreconditioner()
        
        R_inactive = MatPtAP(self.H.model.Prior.R, self.I_set)
        R_inactive.ident_zeros()
#        R_inactive.axpy(1., self.A_set, False)
        
        Rprec = PETScKrylovSolver("richardson", "amg")
        Rprec.parameters["maximum_iterations"] = 1
        Rprec.parameters["error_on_nonconvergence"] = False
        Rprec.parameters["nonzero_initial_guess"] = False
        Rprec.set_operator(R_inactive)
        
        return Rprec
              
    def mult(self,x,y):
        if self.active_set:
            # y(active_set) = x(active_set)
            self.A_set.mult(x,y)
            # y(inactive_set)= H(inactive_set,inactive_set)*x(inactive_set)
            self.I_set.mult(x, self.x_inactive)
            self.H.mult(self.x_inactive, self.y_inactive)
            self.I_set.mult(self.y_inactive, self.x_inactive)
            y.axpy(1., self.x_inactive)
        else:                  
            self.H.mult(x, y)   
            

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

import math
import numpy as np
from variables import PARAMETER
from cgsolverSteihaug import CGSolverSteihaug
from reducedHessian import ReducedHessianActiveSet


class PDActiveSet:
    """
    Primal-dual Active Set Method for box constraints.
    This function is sperimental.
    """
    termination_reasons = [
                           "Maximum Number of Iterations Reached",
                           "Norm of the gradient less then tolerance"
                           ]
    
    def __init__(self, model):
        """
        Initialize the Primal-Dual Active Set with the following parameters.
        abs_tolerance         --> we achieve convernge when sqrt(g,g) <= abs_tolerance
        max_iter              --> maximum number of iterations
        inner_rel_tolerance   --> relative tolerance used for the solution of the
                                  forward, adjoint, and incremental (fwd,adj) problems
        print_level           --> Print info on screen
        GN_iter               --> Number of Gauss Newton iterations before switching to Newton
        cg_coarse_tolerance   --> Coarsest tolerance for the 
        """
        self.model = model
        
        self.parameters = {}
        self.parameters["abs_tolerance"]         = 1e-6
        self.parameters["max_iter"]              = 20
        self.parameters["inner_rel_tolerance"]   = 1e-9
        self.parameters["print_level"]           = 0
        self.parameters["GN_iter"]               = 5
        self.parameters["cg_coarse_tolerance"]   = .5
        self.parameters["cg_tolerance"]   = 1e-8
        
        self.c = 100.
        
        self.it = 0
        self.converged = False
        self.total_cg_iter = 0
        self.ncalls = 0
        self.reason = 0
        self.final_grad_norm = 0
        
    def solve_linear(self,a0, lb = None, ub = None):
        """
        Solve the constrained optimization problem with initial guess a0.
        Return the solution [u,a,p] 
        """
        tol = self.parameters["abs_tolerance"]
        max_iter = self.parameters["max_iter"]
        innerTol = self.parameters["inner_rel_tolerance"]
        print_level = self.parameters["print_level"]
        
        [u,a,p] = self.model.generate_vector()
        self.model.solveFwd(u, [u, a0, p], innerTol)
        self.model.solveAdj(p, [u,a0,p], innerTol)
        
        mg = self.model.generate_vector(PARAMETER)    
        self.model.setPointForHessianEvaluations([u,a0,p])
        self.model.evalGradientParameter([u,a0,p], mg)
        
        self.it = 0
        self.converged = False
        self.ncalls += 1
                
        HessApply = ReducedHessianActiveSet(self.model, innerTol)
        solver = CGSolverSteihaug()
        solver.set_operator(HessApply)
        solver.set_preconditioner(HessApply.getPreconditioner())
        solver.parameters["zero_initial_guess"] = True
        solver.parameters["print_level"] = print_level-1
        solver.parameters["rel_tolerance"] = self.parameters["cg_tolerance"]
            
        solver.solve(a, -mg)
        
        l = self.model.generate_vector(PARAMETER)
        l.zero()
        res = self.model.generate_vector(PARAMETER)
        
        while self.it < max_iter:
            self.it += 1
            
            solver.parameters["zero_initial_guess"] = False
            index_active = np.zeros(l.array().shape, dtype = l.array().dtype)
            self.update_active_set(l, a, lb, ub, index_active)
            HessApply.setActiveSet(index_active)
            x = self.model.generate_vector(PARAMETER)
            x.axpy(-1., mg)
            
            HessApply.applyBounds(x, lb, ub)
            solver.set_operator(HessApply)
            solver.set_preconditioner(HessApply.getPreconditioner())
            solver.solve(a, x)
            HessApply.active_set = False
            HessApply.mult(a, res)
            res *= -1.
            res.axpy(-1., mg)
            l.zero()
            l[index_active != 0] = res[index_active != 0]
            
       
        return [u,a,p,l]
    
    def update_active_set(self, l, a, lb, ub, index_active):
        if lb is not None:
            tmp = l + self.c*(a - lb)         
            index_active[ tmp.array()[:] < 0.] = 1
            
        if ub is not None:
            tmp = l + self.c*(a - ub)
            index_active[ tmp.array()[:] > 0.] = -1

    

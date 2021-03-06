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

import sys
import dolfin as dl
import numpy as np
import math

from variables import PARAMETER
from cgsolverSteihaug import CGSolverSteihaug
from reducedHessian import ReducedHessian
from linalg import vector2Function, MPIAllReduceVector
from bfgs import BFGS_operator

from fenicstools.plotfenics import PlotFenics
from fenicstools.linalg.miscroutines import compute_eigfenics

class ReducedSpaceNewtonCG:
    
    """
    Inexact Newton-CG method to solve constrained optimization problems in the reduced parameter space.
    The Newton system is solved inexactly by early termination of CG iterations via Eisenstat-Walker
    (to prevent oversolving) and Steihaug (to avoid negative curvature) criteria.
    Globalization is performed using the armijo sufficient reduction condition (backtracking).
    The stopping criterion is based on a control on the norm of the gradient and a control of the
    inner product between the gradient and the Newton direction.
       
    The user must provide a model that describes the forward problem, cost functionals, and all the
    derivatives for the gradient and the Hessian.
    
    More specifically the model object should implement following methods:
       - generate_vector() -> generate the object containing state, parameter, adjoint
       - cost(x) -> evaluate the cost functional, report regularization part and misfit separately
       - solveFwd(out, x,tol) -> solve the possibly non linear Fwd Problem up a tolerance tol
       - solveAdj(out, x,tol) -> solve the linear adj problem
       - evalGradientParameter(x, out) -> evaluate the gradient of the parameter and compute its norm
       - setPointForHessianEvaluations(x) -> set the state to perform hessian evaluations
       - solveFwdIncremental(out, rhs, tol) -> solve the linearized forward problem for a given rhs
       - solveAdjIncremental(out, rhs, tol) -> solve the linear adjoint problem for a given rhs
       - applyC(da, out)    --> Compute out = C_x * da
       - applyCt(dp, out)   --> Compute out = C_x' * dp
       - applyWuu(du,out)   --> Compute out = Wuu_x * du
       - applyWua(da, out)  --> Compute out = Wua_x * da
       - applyWau(du, out)  --> Compute out = Wau * du
       - applyR(da, out)    --> Compute out = R * da
       - applyRaa(da,out)   --> Compute out = Raa * out
       - Rsolver()          --> A solver for the regularization term
       
    Type help(ModelTemplate) for additional information
    """
    termination_reasons = [
                           "Maximum number of Iteration reached",      #0
                           "Norm of the gradient less than tolerance", #1
                           "Maximum number of backtracking reached",   #2
                           "Norm of (g, da) less than tolerance"       #3
                           ]
    
    def __init__(self, model):
        """
        Initialize the ReducedSpaceNewtonCG with the following parameters.
        rel_tolerance         --> we converge when sqrt(g,g)/sqrt(g_0,g_0) <= rel_tolerance
        abs_tolerance         --> we converge when sqrt(g,g) <= abs_tolerance
        gda_tolerance         --> we converge when (g,da) <= gda_tolerance
        max_iter              --> maximum number of iterations
        inner_rel_tolerance   --> relative tolerance used for the solution of the
                                  forward, adjoint, and incremental (fwd,adj) problems
        c_armijo              --> Armijo constant for sufficient reduction
        max_backtracking_iter --> Maximum number of backtracking iterations
        print_level           --> Print info on screen
        GN_iter               --> Number of Gauss Newton iterations before switching to Newton
        cg_coarse_tolerance   --> Coarsest tolerance for the CG method (Eisenstat-Walker)
        PC                    --> type of preconditioner for inexact CG ('prior' or 'BFGS')
        parameters for BFGS (optional):
        memory_limit
        H0inv
        """
        self.model = model
        
        self.parameters = {}
        self.parameters["rel_tolerance"]         = 1e-6
        self.parameters["abs_tolerance"]         = 1e-12
        self.parameters["gda_tolerance"]         = 1e-18
        self.parameters["max_iter"]              = 20
        self.parameters["inner_rel_tolerance"]   = 1e-9
        self.parameters["c_armijo"]              = 1e-4
        self.parameters["max_backtracking_iter"] = 10
        self.parameters["print_level"]           = 0
        self.parameters["GN_iter"]               = 5
        self.parameters["check_param"]           = 10
        self.parameters["cg_coarse_tolerance"]   = .5
        self.parameters['PC']                    = 'prior' 
        # BFGS parameters:
        self.parameters['memory_limit']          = 50
        self.parameters['H0inv']                 = 'Rinv'
        
        self.it = 0
        self.converged = False
        self.total_cg_iter = 0
        self.ncalls = 0
        self.reason = 0
        self.final_grad_norm = 0
        self.mpicomm_global = dl.mpi_comm_self()
        self.bfgsPC = None
        
    # add callback function for ReducedHessian (allow user-defined)
    def solve(self, a0, InexactCG=0, GN=False, bounds_xPARAM=None):
        """
        Solve the constrained optimization problem with initial guess a0.
        InexactCG: use Inexact CG method to solve Newton system; tol_cg = sqrt(||grad||)
        GN: use GN Hessian
        bounds_xPARAM: set bounds for parameter a in line search to avoid
        potential instabilities
        Return the solution [u,a,p] 
        """
        rel_tol = self.parameters["rel_tolerance"]
        abs_tol = self.parameters["abs_tolerance"]
        max_iter = self.parameters["max_iter"]
        innerTol = self.parameters["inner_rel_tolerance"]
        c_armijo = self.parameters["c_armijo"]
        max_backtracking_iter = self.parameters["max_backtracking_iter"]
        print_level = self.parameters["print_level"]
        GN_iter = self.parameters["GN_iter"]
        check_param = self.parameters["check_param"]
        cg_coarse_tolerance = self.parameters["cg_coarse_tolerance"]
        PC = self.parameters['PC']

        mpirank = dl.MPI.rank(a0.mpi_comm())

        try:
            self.model.mediummisfit(a0)
            self.mm = True
        except:
            self.mm = False

        if PC == 'BFGS':
            self.bfgsPC = BFGS_operator(self.parameters)
            H0inv = self.bfgsPC.parameters['H0inv']

        [u,a,p] = self.model.generate_vector()
        self.model.solveFwd(u, [u, a0, p], innerTol)
        
        self.it = 0
        self.converged = False
        self.ncalls += 1
        
        ahat = self.model.generate_vector(PARAMETER)    
        mg = self.model.generate_vector(PARAMETER)
        
        cost_old, reg_old, misfit_old = self.model.cost([u,a0,p])
        c, r, m = cost_old, reg_old, misfit_old
        alpha = 0.0

        if self.mm:
            medmisf, perc = self.model.mediummisfit(a0)
        else:
            medmisf, perc = -99, -99
        if(print_level >= 0):
            print "\n {:3} {:5} {:5} {:15} {:15} {:15} {:15} {:14} {:14} {:14} {:14}".format(
            "It", "cg_it", "nbPDE", "cost", "misfit", "reg", "(g,da)", "||g||L2", "alpha", "tolcg", "medmisf")
            #print "{:3d} {:5} {:3d} {:15e} {:15e} {:15e} {:15} {:14} {:14} {:14} {:14e} ({:3.1f}%)".format(
            #self.it, "", self.model.getPDEcounts(), cost_old, misfit_old, reg_old, "", "", "", "", medmisf, perc)
        
        while (self.it < max_iter) and (self.converged == False):
            self.model.solveAdj(p, [u,a0,p], innerTol)
            
            self.model.setPointForHessianEvaluations([u,a0,p])
            mg_old = mg.copy()
            gradnorm = self.model.evalGradientParameter([u,a0,p], mg)
            
            if PC == 'BFGS':
                # update s and y
                if self.it > 0:
                    s = ahat * alpha
                    y = mg - mg_old
                    theta = self.bfgsPC.update(s, y)
                else:
                    theta = 1.0
                # update H0
                if H0inv == 'Rinv':
                    self.bfgsPC.set_H0inv(self.model.Prior.getprecond())
                elif H0inv == 'Minv':
                    self.bfgsPC.set_H0inv(self.model.Prior.Msolver)

            if self.it == 0:
                gradnorm_ini = gradnorm
                tol = max(abs_tol, gradnorm_ini*rel_tol)

            # check if solution is reached
            if (gradnorm < tol) and (self.it > 0):
                self.converged = True
                self.reason = 1
                if self.mm:
                    medmisf, perc = self.model.mediummisfit(a0)
                else:
                    medmisf, perc = -99, -99
                if print_level >= 0:
                    print "{:3d} {:5} {:3d} {:15e} {:15e} {:15e} {:15} {:14e} {:14e} {:14} {:14e} ({:3.1f}%)".format(
                    self.it, "", self.model.getPDEcounts(), 
                    c, m, r, "", gradnorm, alpha, "", medmisf, perc)
                break
            
            if InexactCG==1:
                tolcg = min(cg_coarse_tolerance, math.sqrt(gradnorm/gradnorm_ini))
            elif InexactCG==2:
                tolcg = min(cg_coarse_tolerance, gradnorm/gradnorm_ini)
            else:
                tolcg = 1e-12
            
            if GN:
                HessApply = ReducedHessian(self.model, innerTol, True)
            else:
                HessApply = ReducedHessian(self.model, innerTol, self.it < GN_iter)
            solver = CGSolverSteihaug()
            solver.set_operator(HessApply)
            if PC == 'prior':
                solver.set_preconditioner(self.model.Rsolver())
            elif PC == 'BFGS':
                solver.set_preconditioner(self.bfgsPC)
            solver.parameters["rel_tolerance"] = tolcg
            solver.parameters["zero_initial_guess"] = True
            solver.parameters["print_level"] = print_level-1
            
            try:
                solver.solve(ahat, -mg)
            except RuntimeError as err:
                if mpirank == 0:    
                    print str(err)
                    print 'Could not solve Newton system'
                    print 'plot a1, a2, mg1, mg2'
                a0fun = vector2Function(a0, self.model.Vh[PARAMETER])
                a1, a2 = a0fun.split(deepcopy=True)
                mgfun = vector2Function(mg, self.model.Vh[PARAMETER])
                mg1, mg2 = mgfun.split(deepcopy=True)
                plt = PlotFenics(a0.mpi_comm(), 'Output-failure-NewtonCG')
                plt.set_varname('a1')
                plt.plot_vtk(a1)
                plt.set_varname('a2')
                plt.plot_vtk(a2)
                plt.set_varname('mg1')
                plt.plot_vtk(mg1)
                plt.set_varname('mg2')
                plt.plot_vtk(mg2)
                if mpirank == 0:    print 'compute and print eigenvalues preconditioner'
                compute_eigfenics(self.model.Prior.precond,
                'Output-failure-NewtonCG/eigvaluesPrecond.txt')
                sys.exit(1)
            self.total_cg_iter += HessApply.ncalls
            mg_ahat = mg.inner(ahat)

            if self.mm:
                medmisf, perc = self.model.mediummisfit(a0)
            else:
                medmisf, perc = -99, -99
            if print_level >= 0:
                print "{:3d} {:5d} {:3d} {:15e} {:15e} {:15e} {:15e} {:14e} {:14e} {:14e} {:14e} ({:3.1f}%)".format(
                self.it, HessApply.ncalls, self.model.getPDEcounts(), 
                c, m, r, mg_ahat, gradnorm, alpha, tolcg, medmisf, perc)
            
            alpha = 1.0
            descent = 0
            n_backtrack = 0
            
            while descent == 0 and n_backtrack < max_backtracking_iter:
                # update a and u
                a.zero()
                a.axpy(1., a0)
                a.axpy(alpha, ahat)
                if bounds_xPARAM is not None:
                    amin = a.min()
                    amax = a.max()
                    if amin < bounds_xPARAM[0] or amax > bounds_xPARAM[1]:
                        u.zero()
                        n_backtrack += 1
                        alpha *= 0.5
                        continue
                try:
                    self.model.solveFwd(u, [u, a, p], innerTol)
                except RuntimeError as err:
                    if mpirank == 0:    
                        print 'Forward solve failed during line search'
                        print 'plot a1, a2'
                        print str(err)
                    try:
                        afun = vector2Function(a, self.model.Vh[PARAMETER])
                        a1, a2 = afun.split(deepcopy=True)
                        a1min = a1.vector().min()
                        a1max = a1.vector().max()
                        a2min = a2.vector().min()
                        a2max = a2.vector().max()
                        if mpirank == 0:
                            print 'min(a1)={}, max(a1)={}, min(a2)={}, max(a2)={}'.format(\
                            a1min, a1max, a2min, a2max)
                        plt = PlotFenics('Output-failure-NewtonCG')
                        plt.set_varname('a1')
                        plt.plot_vtk(a1)
                        plt.set_varname('a2')
                        plt.plot_vtk(a2)
                    except:
                        amin = a.min()
                        amax = a.max()
                        if mpirank == 0:
                            print 'min(a)={}, max(a)={}'.format(amin, amax)
                        plt = PlotFenics('Output-failure-NewtonCG')
                        plt.set_varname('a')
                        plt.plot_vtk(vector2Function(a, self.model.Vh[PARAMETER]))
                    sys.exit(1)

                cost_new, reg_new, misfit_new = self.model.cost([u,a,p])
                c, r, m = cost_new, reg_new, misfit_new
                
                # Check if armijo conditions are satisfied
                if (cost_new < cost_old + alpha * c_armijo * mg_ahat) or (-mg_ahat <= self.parameters["gda_tolerance"]):
                    cost_old = cost_new
                    descent = 1
                    a0.zero()
                    a0.axpy(1., a)
                else:
                    n_backtrack += 1
                    alpha *= 0.5

            if self.it % check_param == 0:
                self.compareparam(a0)

            self.it += 1

            # for primal-dual Newton method only:
            if self.model.Prior.isPD():
                self.model.Prior.update_w(ahat, alpha)

            if n_backtrack == max_backtracking_iter:
                self.converged = False
                self.reason = 2
                break
            
            if -mg_ahat <= self.parameters["gda_tolerance"]:
                self.converged = True
                self.reason = 3
                break
                            
        self.final_grad_norm = gradnorm
        self.final_cost      = cost_new
        return [u,a0,p]


    def compareparam(self, a0):
        a0_recv = a0.copy()
        a0_recv.zero()
        na0 = np.linalg.norm(a0.array())
        MPIAllReduceVector(a0, a0_recv, self.mpicomm_global)
        a0_recv /= dl.MPI.size(self.mpicomm_global)
        diff = a0_recv - a0
        reldiff = np.linalg.norm(diff.array())/na0
        assert reldiff < 2e-16, 'Diff in a0 across proc: {:.2e}'.format(reldiff)

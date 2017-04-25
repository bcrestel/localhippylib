
import sys
import dolfin as dl

import math
from variables import PARAMETER
from linalg import vector2Function

from fenicstools.plotfenics import PlotFenics
from fenicstools.linalg.miscroutines import compute_eigfenics


#TODO: convert into L-BFGS
#TODO: create mode to compute BFGS approx to data misfit Hessian only (???)
#TODO: count PDE solves
class BFGS:
    """
    Implement BFGS technique with backtracking inexact line search and damped updating
    See Nocedal & Wright (06), $6.2, $7.3, $18.3

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
        Initialize the BFGS model with the following parameters.
        rel_tolerance         --> we converge when sqrt(g,g)/sqrt(g_0,g_0) <= rel_tolerance
        abs_tolerance         --> we converge when sqrt(g,g) <= abs_tolerance
        gda_tolerance         --> we converge when (g,da) <= gda_tolerance
        max_iter              --> maximum number of iterations
        inner_rel_tolerance   --> relative tolerance used for the solution of the
                                  forward, adjoint, and incremental (fwd,adj) problems
        c_armijo              --> Armijo constant for sufficient reduction
        max_backtracking_iter --> Maximum number of backtracking iterations
        print_level           --> Print info on screen
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
        self.parameters["BFGS_damping"]          = 0.2
        
        self.it = 0
        self.converged = False
        self.ncalls = 0
        self.reason = 0
        self.final_grad_norm = 0

        self.S, self.Y, self.R = [], [], []
        self.apply_H0 = self.apply_H0_default
        self.d0 = 1.0   # H0 = d0 x Identity_matrix



    def apply_H0_default(self, x_in, Hx_out):
        """
        Default application of component H0
        """
        xcopy = x_in.copy()
        Hx_out.zero()
        Hx_out.axpy(self.d0, xcopy)

    def apply_Rinv(self, x_in, Hx_out):
        xcopy = x_in.copy()
        Rinvsolver = self.model.Prior.getprecond()
        Rinvsolver.solve(Hx_out, xcopy)

    def apply_Minv(self, x_in, Hx_out):
        xcopy = x_in.copy()
        Minvsolver = self.model.Prior.Msolver
        Minvsolver.solve(Hx_out, xcopy)



    def apply_Hk(self, x_in, Hx_out):
        """
        Compute matvec between current Hk matrix (approx to inverse of Hessian)
        and a vector x_in:  
                            Hx_out = Hk * x_in
        The Hk matrix is the BFGS approximation to the inverse of the Hessian
        Computation done via double-loop algorithm
        Input:
            x_in = vector (Fenics)  [in]
            Hx_out = vector(Fenics) [out]
        """
        A = []
        Hx_out.zero()
        Hx_out.axpy(1.0, x_in)

        for s, y, r in zip(reversed(self.S), reversed(self.Y), reversed(self.R)):
            a = r * s.inner(Hx_out)
            A.append(a)
            Hx_out.axpy(-a, y)

        self.apply_H0(Hx_out, Hx_out)

        for s, y, r, a in zip(self.S, self.Y, self.R, reversed(A)):
            b = r * y.inner(Hx_out)
            Hx_out.axpy(a - b, s)



    def solve(self, a0, bounds_xPARAM=None):
        """
        Solve the constrained optimization problem with initial guess a0.
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
        damp = self.parameters["BFGS_damping"]
        mpirank = dl.MPI.rank(a0.mpi_comm())

        try:
            self.model.mediummisfit(a0)
            self.mm = True
        except:
            self.mm = False

        [u,a,p] = self.model.generate_vector()
        self.model.solveFwd(u, [u, a0, p], innerTol)
        
        self.it = 0
        self.converged = False
        self.ncalls += 1
        
        ahat = self.model.generate_vector(PARAMETER)    
        mg = self.model.generate_vector(PARAMETER)
        Hy = self.model.generate_vector(PARAMETER)
        
        cost_old, reg_old, misfit_old = self.model.cost([u,a0,p])

        if self.mm:
            medmisf, perc = self.model.mediummisfit(a0)
        else:
            medmisf, perc = -99, -99
        if(print_level >= 0):
            print "\n{:3} {:15} {:15} {:15} {:15} {:14} {:14} {:14} {:14}".format(
            "It", "cost", "misfit", "reg", "(g,da)", "||g||L2", "alpha", "theta", "medmisf")
            print "{:3d} {:15e} {:15e} {:15e} {:15} {:14} {:14} {:14} {:14e} ({:3.1f}%)".format(
            self.it, cost_old, misfit_old, reg_old, "", "", "", "", medmisf, perc)
        
        while (self.it < max_iter) and (self.converged == False):
            self.model.solveAdj(p, [u,a0,p], innerTol)
            
            self.model.setPointForHessianEvaluations([u,a0,p])

            mg_old = mg.copy()
            gradnorm = self.model.evalGradientParameter([u,a0,p], mg)
            # Update BFGS with damped update
            theta = 1.0 # default value
            if self.it > 0:
                s = ahat * alpha
                y = mg - mg_old
                sy = s.inner(y)
                self.apply_Hk(y, Hy)
                yHy = y.inner(Hy)
                if mpirank == 0:
                    print 'sy={}, yHy={}, damp*yHy={}'.format(sy, yHy, damp*yHy)
                if sy < damp*yHy:
                    theta = (1.0-damp)*yHy/(yHy-sy)
                    s *= theta
                    s.axpy(1-theta, Hy)
                    sy = s.inner(y)
                assert(sy > 0.)
                rho = 1./sy
                self.S.append(s.copy())
                self.Y.append(y.copy())
                self.R.append(rho)
                # re-scale H0 based on first update
                if self.it == 1:
                    self.d0 = sy/y.inner(y)
            
            if self.it == 0:
                gradnorm_ini = gradnorm
                tol = max(abs_tol, gradnorm_ini*rel_tol)
                
            # check if solution is reached
            if (gradnorm < tol) and (self.it > 0):
                self.converged = True
                self.reason = 1
                break
            
            self.it += 1

            # compute search direction with BFGS:
            self.apply_Hk(-mg, ahat)
            
            # backtracking line-search
            alpha = 1.0
            descent = 0
            n_backtrack = 0
            mg_ahat = mg.inner(ahat)
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
                    print str(err)
                    sys.exit(1)

                cost_new, reg_new, misfit_new = self.model.cost([u,a,p])
                
                # Check if armijo conditions are satisfied
                if (cost_new < cost_old + alpha * c_armijo * mg_ahat) or (-mg_ahat <= self.parameters["gda_tolerance"]):
                    cost_old = cost_new
                    descent = 1
                    a0.zero()
                    a0.axpy(1., a)
                else:
                    n_backtrack += 1
                    alpha *= 0.5

            if self.mm:
                medmisf, perc = self.model.mediummisfit(a)
            else:
                medmisf, perc = -99, -99
            if print_level >= 0:
                print "{:3d} {:15e} {:15e} {:15e} {:15e} {:14e} {:14e} {:14e} {:14e} ({:3.1f}%)".format(
                self.it, cost_new, misfit_new, reg_new, mg_ahat, gradnorm, alpha, theta, medmisf, perc)
                
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
